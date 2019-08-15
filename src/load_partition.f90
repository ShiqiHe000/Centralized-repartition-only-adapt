

MODULE LOAD_PARTITION

USE MESH, ONLY: NEL_TOTAL, ELE_PROC_NUM, NEL_LOCAL, X_LOCAL, Y_LOCAL
USE ADAPT, ONLY: LPOLY_ORDER, H_DEPTH
USE SIZE
USE MPI
USE FIELDS

IMPLICIT NONE

!-----------------------------------------------------------------------
INTEGER, ALLOCATABLE, DIMENSION(:) :: W ! COMPUTATIONAL LOAD OF EACH ELEMENT
!INTEGER :: PREFIX_SUM_OFFSET    ! GLOBAL PREFIX_SUM OF WEIGHTS
INTEGER, ALLOCATABLE, DIMENSION(:) :: PREFIX_SUM   ! PREFIX SUM 
INTEGER, ALLOCATABLE, DIMENSION(:) :: LOAD_LOCAL ! LOCAL LOAD IN EACH PROCESSOR
INTEGER :: B    !BOTTLENECK

REAL(KIND=8) :: B_OPT   ! OPTIMAL BOTTLENECK
REAL(KIND=8) :: LOAD_BALANCING  ! LOAD BALANCING RATIO
!-----------------------------------------------------------------------

CONTAINS

SUBROUTINE GET_WEIGHT(NEL)

    INTEGER :: IEL
    INTEGER, INTENT(IN) :: NEL
    
    ! ALLOCATE----------------------------------------------------------
    ALLOCATE(W(NEL))
    !-------------------------------------------------------------------
    
    ! INITIALIZE--------------------------------------------------------
    W=0
    !-------------------------------------------------------------------
    
    IF ( LP_SERIAL ) THEN   ! NOT SERIAL LOAD BALANCING
        DO IEL=1, NEL
            W(IEL) = (GPOLY_ORDER(IEL))**4
        ENDDO
    ELSE    ! SERIAL LOAD BALANCING
        DO IEL=1, NEL
            W(IEL) = (LPOLY_ORDER(IEL))**4 
        ENDDO
    ENDIF

END SUBROUTINE GET_WEIGHT

SUBROUTINE H1_H2_PARTITION_SERIAL
    !-------------------------------------------------------------------
    INTEGER :: IEL, CUT_I, I
    INTEGER :: DIS1, DIS2   ! THE DISTANCE BETWEEN THE TWO ELEMENTS WEIGHT AND THE B_OPT*CUT_I
    INTEGER, ALLOCATABLE, DIMENSION(:) :: CUTTING_POINT
    INTEGER, ALLOCATABLE, DIMENSION(:) :: SENDCOUNTS
    INTEGER, ALLOCATABLE, DIMENSION(:) :: DISPLS1, DISPLS2
    
    !-------------------------------------------------------------------

    ! SET LP_SERIAL TO BE TRUE------------------------------------------
    LP_SERIAL= .TRUE.
    !-------------------------------------------------------------------
    
    IF (RANK == 0) THEN ! ONLY PARTITION IN ROOT -----------------------
        
        ! ALLOCATE------------------------------------------------------
        ALLOCATE(PREFIX_SUM(0:NEL_TOTAL))
        ALLOCATE(CUTTING_POINT(0:NUM_PROC))
        ALLOCATE(LOAD_LOCAL(NUM_PROC))
        ALLOCATE(SENDCOUNTS(NUM_PROC))
        ALLOCATE(DISPLS1(NUM_PROC), DISPLS2(NUM_PROC))
        !---------------------------------------------------------------
        
        ! INITIALIZATION------------------------------------------------
        PREFIX_SUM=0
        CUTTING_POINT=0
        CUT_I=1
        LOAD_LOCAL=0
        SENDCOUNTS=0
        DISPLS1=0; DISPLS2=0
        !---------------------------------------------------------------
    
        ! GET WEIGHT----------------------------------------------------
        CALL GET_WEIGHT(NEL_TOTAL)
        !---------------------------------------------------------------
        
        ! CALCULATE THE PREFIX_SUM--------------------------------------
        DO IEL=1, NEL_TOTAL
            PREFIX_SUM(IEL)=PREFIX_SUM(IEL-1)+W(IEL)
        ENDDO
        !---------------------------------------------------------------
        
        ! OPTIMAL BOTTLENECK--------------------------------------------
        B_OPT=REAL(PREFIX_SUM(NEL_TOTAL), KIND=8)/REAL(NUM_PROC, KIND=8)
        !---------------------------------------------------------------
        
!        PRINT *, "B_OPT", B_OPT
!        PRINT *, "PREFIX_SUM", PREFIX_SUM

        ! HI PARTITION--------------------------------------------------
        DO IEL=1, NEL_TOTAL
            
            IF(PREFIX_SUM(IEL) > (B_OPT*CUT_I)) THEN
                CUTTING_POINT(CUT_I) = IEL-1
                CUT_I=CUT_I+1
            ENDIF
        
            IF(CUT_I > (NUM_PROC-1)) EXIT
        ENDDO
        
        CUTTING_POINT(NUM_PROC) = NEL_TOTAL
        
        IF(CUTTING_POINT(1) == 0) THEN
            CUTTING_POINT(1)=1
        ENDIF
        !---------------------------------------------------------------
        
!        PRINT *, CUTTING_POINT

        ! H2 PARITITON--------------------------------------------------
        DO CUT_I=1, NUM_PROC-1
            DIS1=B_OPT*CUT_I - PREFIX_SUM(CUTTING_POINT(CUT_I))
            DIS2=PREFIX_SUM(CUTTING_POINT(CUT_I+1)) - B_OPT*CUT_I
            
            IF(DIS2 < DIS1) THEN
                CUTTING_POINT(CUT_I)=CUTTING_POINT(CUT_I)+1
            ENDIF
        ENDDO
        !---------------------------------------------------------------
        
         ! GET LOCAL LOAD ON EACH PROCESSOR AND ELEMENT PROC NUMBER-----
        DO I=1, NUM_PROC
            LOAD_LOCAL(I) = PREFIX_SUM(CUTTING_POINT(I))- &
                            PREFIX_SUM(CUTTING_POINT(I-1))
            GELE_PROC_NUM(CUTTING_POINT(I-1)+1: CUTTING_POINT(I)) = I
        ENDDO
        !---------------------------------------------------------------
    
        ! GET BOTTLENECK------------------------------------------------
        B = MAXVAL(LOAD_LOCAL)
        !---------------------------------------------------------------
        
!        PRINT *, "B", B
        
        ! COMPUTE LOAD BALANCING RATIO----------------------------------
        LOAD_BALANCING=B_OPT/REAL(B, KIND=8)*100D0
        OPEN(UNIT=3, FILE='load_balancing_serial.txt', ACCESS='APPEND')
        WRITE(3, 20) B_OPT, B, LOAD_BALANCING
20 FORMAT(F11.3, 2X, I7, 2X, F10.5)
        CLOSE(UNIT=3)
!        PRINT *, "H2_LOAD_BALANCING", LOAD_BALANCING
        !---------------------------------------------------------------
        
        ! PREPARE FOR MPI_SCATTERV--------------------------------------
        DO I=1, NUM_PROC
            SENDCOUNTS(I)=CUTTING_POINT(I)-CUTTING_POINT(I-1)
            DISPLS1(I) = CUTTING_POINT(I-1)*4
            DISPLS2(I) = CUTTING_POINT(I-1)
        ENDDO
        !---------------------------------------------------------------
!        PRINT *, "SENDCOUNTS", SENDCOUNTS
!        PRINT *, "NEL_TOTAL", NEL_TOTAL
!        PRINT *, "CUTTING_POINT", CUTTING_POINT

    ENDIF
    !--------------------------------------------------------------------
    
    ! DEALLOCATE--------------------------------------------------------
    DEALLOCATE(X_LOCAL, Y_LOCAL)
    DEALLOCATE(ELE_PROC_NUM)
    DEALLOCATE(H_DEPTH)
    DEALLOCATE(LPOLY_ORDER)
    !-------------------------------------------------------------------
    
    ! SEND ELEMENT NUMBER TO EACH PROC----------------------------------
    CALL MPI_SCATTER(SENDCOUNTS, 1, MPI_INT, NEL_LOCAL, 1, &
                    MPI_INT, 0, MPI_COMM_WORLD, IERROR)
    !-------------------------------------------------------------------
    
    ! ALLOCATE----------------------------------------------------------
    ALLOCATE(X_LOCAL(4, NEL_LOCAL), Y_LOCAL(4, NEL_LOCAL))
    ALLOCATE(ELE_PROC_NUM(NEL_LOCAL))
    ALLOCATE(H_DEPTH(NEL_LOCAL))
    ALLOCATE(LPOLY_ORDER(NEL_LOCAL))
    !-------------------------------------------------------------------
    
!    PRINT *, "RANK", RANK, "NEL", NEL_LOCAL
    
    ! TRANSFORM THE ELEMENT TO EACH PROC--------------------------------
    CALL MPI_SCATTERV(X_GLOBAL, SENDCOUNTS*4, DISPLS1, MPI_DOUBLE_PRECISION, &
                       X_LOCAL, NEL_LOCAL*4, MPI_DOUBLE_PRECISION, &
                        0, MPI_COMM_WORLD, IERROR )
                        
    CALL MPI_SCATTERV(Y_GLOBAL, SENDCOUNTS*4, DISPLS1, MPI_DOUBLE_PRECISION, &
                        Y_LOCAL, NEL_LOCAL*4, MPI_DOUBLE_PRECISION, &
                        0, MPI_COMM_WORLD, IERROR )
                        
    CALL MPI_SCATTERV(GELE_PROC_NUM, SENDCOUNTS, DISPLS2, MPI_INT, &
                        ELE_PROC_NUM, NEL_LOCAL, MPI_INT, &
                        0, MPI_COMM_WORLD, IERROR)
                        
    CALL MPI_SCATTERV(GH_DEPTH, SENDCOUNTS, DISPLS2, MPI_INT, &
                       H_DEPTH, NEL_LOCAL, MPI_INT, &
                        0, MPI_COMM_WORLD, IERROR)
                    
    CALL MPI_SCATTERV(GPOLY_ORDER, SENDCOUNTS, DISPLS2, MPI_INT, &
                       LPOLY_ORDER, NEL_LOCAL, MPI_INT, &
                        0, MPI_COMM_WORLD, IERROR)
    !-------------------------------------------------------------------
    
!    IF(RANK==0) THEN
!        PRINT *, X_LOCAL(:, 1)
    
!    ENDIF
    
    ! DEALLOCATE--------------------------------------------------------
    DEALLOCATE(X_GLOBAL, Y_GLOBAL)
    DEALLOCATE(GELE_PROC_NUM)
    DEALLOCATE(GH_DEPTH)
    DEALLOCATE(GPOLY_ORDER)

    IF (RANK == 0) THEN
        DEALLOCATE(DISPLS1, DISPLS2)
        DEALLOCATE(PREFIX_SUM)
        DEALLOCATE(CUTTING_POINT)
        DEALLOCATE(LOAD_LOCAL)
        DEALLOCATE(SENDCOUNTS)
       DEALLOCATE(W)
    ENDIF
    !-------------------------------------------------------------------
    
    ! SYNCHRONIZE-------------------------------------------------------
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !-------------------------------------------------------------------

END SUBROUTINE H1_H2_PARTITION_SERIAL

SUBROUTINE H1_PARTITION_LOCAL
    !-------------------------------------------------------------------
    INTEGER :: IEL, I
    INTEGER :: PREFIX_SUM_OFFSET    ! GLOBAL PREFIX_SUM OF WEIGHTS
    INTEGER, ALLOCATABLE, DIMENSION(:) :: P_MAPPING    ! PROCESSOR MAPPING, LOCAL 
    INTEGER, ALLOCATABLE, DIMENSION(:) :: SDISPLS1, RDISPLS1, SDISPLS2, RDISPLS2
    INTEGER, ALLOCATABLE, DIMENSION(:) :: EXCHANGE_SEND ! ELEMENT EXCHANGE INFORMATION
    INTEGER, ALLOCATABLE, DIMENSION(:) :: EXCHANGE_RECV ! ELEMENT RECV INFORMATION
    INTEGER, ALLOCATABLE, DIMENSION(:) :: H_DEPTH_NEW  
    INTEGER, ALLOCATABLE, DIMENSION(:) :: LPOLY_ORDER_NEW  
    INTEGER, ALLOCATABLE, DIMENSION(:) :: W_NEW

    REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :) :: XX, YY    ! X, Y COORDINATE AFTER EXCHANGE    
    
    !-------------------------------------------------------------------
    
    ! ALLOCATE----------------------------------------------------------
    ALLOCATE(PREFIX_SUM(0:NEL_LOCAL))
    ALLOCATE(P_MAPPING(NEL_LOCAL))
    ALLOCATE(LOAD_LOCAL(1))
    ALLOCATE(SDISPLS1(NUM_PROC), RDISPLS1(NUM_PROC))
    ALLOCATE(SDISPLS2(NUM_PROC), RDISPLS2(NUM_PROC))
    ALLOCATE(EXCHANGE_SEND(NUM_PROC))
    ALLOCATE(EXCHANGE_RECV(NUM_PROC))
    !-------------------------------------------------------------------
    
    ! INITIALIZATION----------------------------------------------------
    PREFIX_SUM_OFFSET=0
    LOAD_LOCAL=0
    PREFIX_SUM=0
    P_MAPPING=0
    SDISPLS1=0; RDISPLS1=0
    SDISPLS2=0; RDISPLS2=0
    EXCHANGE_SEND=0
    EXCHANGE_RECV=0
    !-------------------------------------------------------------------
    
    ! GET WEIGHT--------------------------------------------------------
    CALL GET_WEIGHT(NEL_LOCAL)
    !-------------------------------------------------------------------

    ! CALCULATE THE LOCAL PREFIX SUM -----------------------------------
    DO IEL=1, NEL_LOCAL
        PREFIX_SUM(IEL)=PREFIX_SUM(IEL-1)+W(IEL)
    ENDDO
    !-------------------------------------------------------------------
    
    
    ! CACULATE THE GLOBAL PREFIX_SUM OFFSET (EXCLUSIVE)-----------------
    CALL MPI_EXSCAN(PREFIX_SUM(NEL_LOCAL), PREFIX_SUM_OFFSET, 1, &
                    MPI_INT, MPI_SUM, MPI_COMM_WORLD, IERROR)
    !-------------------------------------------------------------------
    
    ! LOCAL LOAD + LOAD OFFSET------------------------------------------
    PREFIX_SUM(1:NEL_LOCAL)=PREFIX_SUM(1:NEL_LOCAL)+PREFIX_SUM_OFFSET
    !-------------------------------------------------------------------
    
    ! COMPUTE THE OPTIMAL BOTTLENECK------------------------------------
    IF (RANK == NUM_PROC-1) THEN
        B_OPT = REAL(PREFIX_SUM(NEL_LOCAL), KIND=8)/REAL(NUM_PROC, KIND=8)
!        IF (MAXVAL(W) > B_OPT) THEN
!                PRINT *, " "
!                PRINT *, "PERFECT BALANCE CANNOT BE ACHIEVED"
!                PRINT *, " "
!        ENDIF
!    PRINT *, "B_OPT=", B_OPT
    ENDIF
    !-------------------------------------------------------------------
    
    ! BROADCAST B_OPT---------------------------------------------------
    CALL MPI_BCAST(B_OPT, 1, MPI_DOUBLE_PRECISION, NUM_PROC-1, MPI_COMM_WORLD, IERROR)
    !-------------------------------------------------------------------
    
    ! GET PROCESSOR MAPPING---------------------------------------------
    DO IEL=1, NEL_LOCAL
        P_MAPPING(IEL) = (PREFIX_SUM(IEL)-1)/B_OPT
    ENDDO
    !-------------------------------------------------------------------
    
    ! GET ELEMNT MAPPING INDICE------------------------------------------
    DO IEL=1, NEL_LOCAL
        IF(P_MAPPING(IEL) < RANK) THEN  ! THIS ELEMENT SHOULD BE MAPPED TO FORMER PROC
            EXCHANGE_SEND(RANK)=EXCHANGE_SEND(RANK)+1
        ELSEIF(P_MAPPING(IEL) == RANK) THEN ! THIS ELEMENT SHOULD STAY IN THE SAME PROC
            EXCHANGE_SEND(RANK+1)=EXCHANGE_SEND(RANK+1)+1
        ELSE    ! THIS ELEMENT SHOULD BE MAPPED TO NEXT PROC
            EXCHANGE_SEND(RANK+2)=EXCHANGE_SEND(RANK+2)+1
        ENDIF
    ENDDO
    !-------------------------------------------------------------------
    
    ! BROADCAST THE ELEMENT EXCHANGE INFORMATION TO EACH PROC-----------
    CALL MPI_ALLTOALL(EXCHANGE_SEND, 1, MPI_INT, EXCHANGE_RECV, &
                        1, MPI_INT, MPI_COMM_WORLD, IERROR)
    !-------------------------------------------------------------------

    ! UPDATE NEL_LOCAL--------------------------------------------------
    NEL_LOCAL=SUM(EXCHANGE_RECV)
    !-------------------------------------------------------------------
    
    ! ALLOCATE ---------------------------------------------------------
    ALLOCATE(XX(4, NEL_LOCAL), YY(4, NEL_LOCAL))
    ALLOCATE(H_DEPTH_NEW(NEL_LOCAL))
    ALLOCATE(LPOLY_ORDER_NEW(NEL_LOCAL))
    ALLOCATE(W_NEW(NEL_LOCAL))
    !-------------------------------------------------------------------
    
    ! INITILAIZE--------------------------------------------------------
    XX=0.0D0
    YY=0.0D0
    H_DEPTH_NEW=0
    LPOLY_ORDER_NEW=0
    W_NEW=0
    !-------------------------------------------------------------------
    
    ! ELEMENT TRANSFER--------------------------------------------------
    DO I=2, NUM_PROC
        SDISPLS1(I) = SDISPLS1(I-1)+EXCHANGE_SEND(I-1)*4
        RDISPLS1(I) = RDISPLS1(I-1)+EXCHANGE_RECV(I-1)*4
        
        SDISPLS2(I) = SDISPLS2(I-1)+EXCHANGE_SEND(I-1)
        RDISPLS2(I) = RDISPLS2(I-1)+EXCHANGE_RECV(I-1)
    
    ENDDO
    
    CALL MPI_ALLTOALLV(X_LOCAL, EXCHANGE_SEND*4, SDISPLS1, MPI_DOUBLE_PRECISION, &
                        XX, EXCHANGE_RECV*4, RDISPLS1, MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD, IERROR)
                        
    CALL MPI_ALLTOALLV(Y_LOCAL, EXCHANGE_SEND*4, SDISPLS1, MPI_DOUBLE_PRECISION, &
                        YY, EXCHANGE_RECV*4, RDISPLS1, MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD, IERROR)
                        
    CALL MPI_ALLTOALLV(H_DEPTH, EXCHANGE_SEND, SDISPLS2, MPI_INT, &
                        H_DEPTH_NEW, EXCHANGE_RECV, RDISPLS2, MPI_INT, &
                        MPI_COMM_WORLD, IERROR)
                        
    CALL MPI_ALLTOALLV(LPOLY_ORDER, EXCHANGE_SEND, SDISPLS2, MPI_INT, &
                        LPOLY_ORDER_NEW, EXCHANGE_RECV, RDISPLS2, MPI_INT, &
                        MPI_COMM_WORLD, IERROR)
                        
    CALL MPI_ALLTOALLV(W, EXCHANGE_SEND, SDISPLS2, MPI_INT, &
                        W_NEW, EXCHANGE_RECV, RDISPLS2, MPI_INT, &
                        MPI_COMM_WORLD, IERROR)
    !-------------------------------------------------------------------
    
    ! DEALLOCATE--------------------------------------------------------
    DEALLOCATE(X_LOCAL, Y_LOCAL)
    DEALLOCATE(H_DEPTH)
    DEALLOCATE(LPOLY_ORDER)
    DEALLOCATE(ELE_PROC_NUM)
    DEALLOCATE(W)
    !-------------------------------------------------------------------
    
    ! ALLOCATE----------------------------------------------------------
    ALLOCATE(X_LOCAL(4, NEL_LOCAL), Y_LOCAL(4, NEL_LOCAL))
    ALLOCATE(ELE_PROC_NUM(NEL_LOCAL))
    ALLOCATE(H_DEPTH(NEL_LOCAL))
    ALLOCATE(LPOLY_ORDER(NEL_LOCAL))
    ALLOCATE(W(NEL_LOCAL))
    !-------------------------------------------------------------------
    
    ! UPDATE LOCAL INFORMATION------------------------------------------
    X_LOCAL=XX
    Y_LOCAL=YY
    LPOLY_ORDER=LPOLY_ORDER_NEW
    H_DEPTH=H_DEPTH_NEW
    W=W_NEW
    ELE_PROC_NUM=RANK+1
    !-------------------------------------------------------------------
    
    ! DEALLOCATE--------------------------------------------------------
    DEALLOCATE(XX, YY)
    DEALLOCATE(LPOLY_ORDER_NEW)
    DEALLOCATE(H_DEPTH_NEW)
    DEALLOCATE(W_NEW)
    !-------------------------------------------------------------------

    ! GET LOCAL LOAD ON EACH PROCESSOR----------------------------------
    DO IEL=1, NEL_LOCAL
        LOAD_LOCAL=LOAD_LOCAL+W(IEL)
    ENDDO
    !-------------------------------------------------------------------
    
    ! GET BOTTLENECK----------------------------------------------------
    CALL MPI_REDUCE(LOAD_LOCAL, B, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD, IERROR)
    !-------------------------------------------------------------------
    
    ! COMPUTE LOAD BALANCING RATIO--------------------------------------
    IF (RANK == 0) THEN
        LOAD_BALANCING=B_OPT/REAL(B, KIND=8)*100D0
        PRINT *, "-----------------------------------------------------"
        PRINT *, "H1_LOAD_BALANCING", LOAD_BALANCING
            
        ! WRITE LOAD BALANCING TO FILE----------------------------------
        OPEN(UNIT=3, FILE="load_balancing_parallel.txt", ACCESS="APPEND")
        WRITE(3, 30) B_OPT, B, LOAD_BALANCING
30 FORMAT(F11.3, 2X, I7, 2X, F10.5)
        CLOSE(UNIT=3)
    ENDIF
    !-------------------------------------------------------------------
    
!    PRINT *, "B_OPT", B_OPT
!    PRINT *, "rank", RANK, "B", B, "LOAD_LOCAL", LOAD_LOCAL
!    PRINT *, "LOAD_LOCAL", LOAD_LOCAL

    ! DEALLOCATE--------------------------------------------------------
    DEALLOCATE(P_MAPPING)
    DEALLOCATE(SDISPLS1, SDISPLS2)
    DEALLOCATE(RDISPLS1, RDISPLS2)
    DEALLOCATE(LOAD_LOCAL)
    DEALLOCATE(EXCHANGE_RECV, EXCHANGE_SEND)
    DEALLOCATE(PREFIX_SUM)
    DEALLOCATE(W)
    !-------------------------------------------------------------------

    
!    ! SYNCHRONIZE-------------------------------------------------------
!    CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
!    !-------------------------------------------------------------------
    
END SUBROUTINE H1_PARTITION_LOCAL




END MODULE LOAD_PARTITION
