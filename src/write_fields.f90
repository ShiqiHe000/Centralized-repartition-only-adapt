
MODULE FIELDS

USE MPI
USE SIZE

IMPLICIT NONE

INTEGER :: FRAME=1
INTEGER, ALLOCATABLE, DIMENSION(:) :: GPOLY_ORDER   ! GLOBAL POINT POLYNOMIAL ORDER
INTEGER, ALLOCATABLE, DIMENSION(:) :: GELE_PROC_NUM ! GLOBAL ELEMENT PROCESSOR NUMBER
INTEGER, ALLOCATABLE, DIMENSION(:) :: GH_DEPTH      ! GLOBAL ELEMENT H-REFINEMENT DEPTH

REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :) :: X_GLOBAL, Y_GLOBAL


CONTAINS

SUBROUTINE WRITE_FIELDS
    USE MESH, ONLY: ELE_PROC_NUM, X_LOCAL, Y_LOCAL, NEL_LOCAL, NEL_TOTAL
    USE ADAPT, ONLY: H_DEPTH, LPOLY_ORDER
    
    ! VARIABLE----------------------------------------------------------
    CHARACTER(LEN=16) :: FILENAME
    
    INTEGER :: I, IEL
    INTEGER :: ELEM
    INTEGER, ALLOCATABLE, DIMENSION(:) :: RECVCOUNT     ! RANK 0 RECV X_LOCAL ARRAIES LENGTH
    INTEGER, ALLOCATABLE, DIMENSION(:) :: DISPLS1        ! THE DISPLACEMANT RELATIVE TO RECVBUF AT WHICH TO PLACE THE INCOMING DATA
    INTEGER, ALLOCATABLE, DIMENSION(:) :: DISPLS2        ! THE DISPLACEMANT RELATIVE TO RECVBUF AT WHICH TO PLACE THE INCOMING DATA
    !-------------------------------------------------------------------
    
    ! INITIALIZE -------------------------------------------------------
    ELEM=1
    !-------------------------------------------------------------------
    
    ! UPDATE NEL_TOTAL--------------------------------------------------
    CALL MPI_ALLREDUCE(NEL_LOCAL, NEL_TOTAL, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, IERROR)
    !-------------------------------------------------------------------
    
    ! ALLOCATE ---------------------------------------------------------
    ALLOCATE(X_GLOBAL(4, NEL_TOTAL), Y_GLOBAL(4, NEL_TOTAL))
    ALLOCATE(GPOLY_ORDER(NEL_TOTAL))
    ALLOCATE(GELE_PROC_NUM(NEL_TOTAL))
    ALLOCATE(GH_DEPTH(NEL_TOTAL))
    ALLOCATE(RECVCOUNT(NUM_PROC))
    ALLOCATE(DISPLS1(NUM_PROC))
    ALLOCATE(DISPLS2(NUM_PROC))
    !-------------------------------------------------------------------
    
    ! GATHRE NEL_LOCAL ON PROCESS 0-------------------------------------
    CALL MPI_GATHER(NEL_LOCAL, 1, MPI_INT, RECVCOUNT, 1, MPI_INT, 0, &
                    MPI_COMM_WORLD, IERROR)
    !-------------------------------------------------------------------
    
    ! THE POSITION TO PLACE DATA FROM DIFFERENT PROCESSES---------------
    IF (RANK == 0) THEN
        DISPLS1=0    ! INITIALIZE
        DISPLS2=0    ! INITIALIZE
        DO I=2, NUM_PROC
            DISPLS1(I)=DISPLS1(I-1)+RECVCOUNT(I-1)*4
            DISPLS2(I)=DISPLS2(I-1)+RECVCOUNT(I-1)
        ENDDO
        
        RECVCOUNT=RECVCOUNT*4   ! 4 NODES
    ENDIF
    !-------------------------------------------------------------------

    
    ! GATHER LOCAL INFORMATION-------------------------------------------
    CALL MPI_GATHERV(X_LOCAL, 4*NEL_LOCAL, MPI_DOUBLE_PRECISION, X_GLOBAL, &
            RECVCOUNT, DISPLS1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERROR)
            
    CALL MPI_GATHERV(Y_LOCAL, 4*NEL_LOCAL, MPI_DOUBLE_PRECISION, Y_GLOBAL, &
            RECVCOUNT, DISPLS1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERROR)

    CALL MPI_GATHERV(LPOLY_ORDER, NEL_LOCAL, MPI_INT, GPOLY_ORDER, &
            RECVCOUNT/4, DISPLS2, MPI_INT, 0, MPI_COMM_WORLD, IERROR)
            
    CALL MPI_GATHERV(ELE_PROC_NUM, NEL_LOCAL, MPI_INT, GELE_PROC_NUM, &
            RECVCOUNT/4, DISPLS2, MPI_INT, 0, MPI_COMM_WORLD, IERROR)
            
    CALL MPI_GATHERV(H_DEPTH, NEL_LOCAL, MPI_INT, GH_DEPTH, &
            RECVCOUNT/4, DISPLS2, MPI_INT, 0, MPI_COMM_WORLD, IERROR)
    !-------------------------------------------------------------------
    
    
    ! WRITE TO FILE ----------------------------------------------------
    IF (RANK == 0) THEN ! ONLY PROCESSOR ONE WRITE
    
        WRITE(FILENAME,FMT='(''aoutput'',I5.5,''.dat'')') FRAME

        OPEN(2,FILE=FILENAME)
        
        WRITE(2,FMT='(''variables="X","Y",'')',advance='no')
        WRITE(2,FMT='(''"GELE_PROC_NUM"'','','')',advance='no')
        WRITE(2,FMT='(''"GH_DEPTH"'','','')',advance='no')
        WRITE(2,FMT='(''"GPOLY_ODER"'','','')')
        

        DO IEL=1, NEL_TOTAL
        
            WRITE(2,FMT='(''ZONE T="IEL'',I6,''", i='',I2,'', j='',I2, &
                    '' F=POINT'')') &
                        ELEM, 2, 2
            ELEM=ELEM+1
            
!            DO I=1, 4
!                WRITE(2, 10) X_GLOBAL(I, IEL), Y_GLOBAL(I, IEL), GELE_PROC_NUM(IEL), &
!                                GH_DEPTH(IEL), GPOLY_ORDER(IEL)
!!                WRITE(*, 10) X_GLOBAL(I, IEL), Y_GLOBAL(I, IEL), GELE_PROC_NUM(IEL), &
!!                                GH_DEPTH(IEL), GPOLY_ORDER(IEL)
!            ENDDO
             WRITE(2, 10) X_GLOBAL(1, IEL), Y_GLOBAL(1, IEL), GELE_PROC_NUM(IEL), &
                                GH_DEPTH(IEL), GPOLY_ORDER(IEL)
             WRITE(2, 10) X_GLOBAL(2, IEL), Y_GLOBAL(2, IEL), GELE_PROC_NUM(IEL), &
                                GH_DEPTH(IEL), GPOLY_ORDER(IEL)
             WRITE(2, 10) X_GLOBAL(4, IEL), Y_GLOBAL(4, IEL), GELE_PROC_NUM(IEL), &
                                GH_DEPTH(IEL), GPOLY_ORDER(IEL)
             WRITE(2, 10) X_GLOBAL(3, IEL), Y_GLOBAL(3, IEL), GELE_PROC_NUM(IEL), &
                                GH_DEPTH(IEL), GPOLY_ORDER(IEL)
        ENDDO
10 FORMAT(F10.5, 2X, F10.5, 2X, I2, 2X, I2, 2X, I2)

    ENDIF

    !-------------------------------------------------------------------
    
    ! SYNCHRONIZE-------------------------------------------------------
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !-------------------------------------------------------------------
    
    ! UPDATE FRAME------------------------------------------------------
    FRAME=FRAME+1
    !-------------------------------------------------------------------
    
    ! DEALLOCATE--------------------------------------------------------
    DEALLOCATE(DISPLS1, DISPLS2)
    DEALLOCATE(RECVCOUNT)
    !-------------------------------------------------------------------
    
    ! IF NO SERIAL LOAD PARITION THEN DEALLOCATE GLOBAL ARRAY-----------
!    IF (.NOT. LP_SERIAL) THEN
!        CALL DEALLOCATE_GLOBAL
!    ENDIF
    !-------------------------------------------------------------------
    
    
END SUBROUTINE WRITE_FIELDS

SUBROUTINE DEALLOCATE_GLOBAL

    ! DEALLOCATE--------------------------------------------------------
    DEALLOCATE(X_GLOBAL, Y_GLOBAL)
    DEALLOCATE(GPOLY_ORDER)
    DEALLOCATE(GELE_PROC_NUM)
    DEALLOCATE(GH_DEPTH)
    !-------------------------------------------------------------------

END SUBROUTINE DEALLOCATE_GLOBAL





END MODULE FIELDS
