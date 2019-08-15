

MODULE ADAPT

USE MPI
USE SIZE

IMPLICIT NONE

INTEGER, ALLOCATABLE, DIMENSION(:) :: LPOLY_ORDER   ! LOCAL POLYNOMIAL ORDER
INTEGER, ALLOCATABLE, DIMENSION(:) :: H_DEPTH   ! ELEMENT H-REFINEMENT DEPTH


CONTAINS

SUBROUTINE PADAPT
    !-------------------------------------------------------------------
    ! PADAPTATION
    ! POLYNOMIAL ORDER= 4, 6, 8, 10, 12
    ! EACH ELEMENT HAS 30% OF CHANCE TO INCREASING POLYNOMIAL ORDER 
    !-------------------------------------------------------------------
    USE MESH, ONLY: NEL_LOCAL
    
    !--------------------------------------------------------------------
    INTEGER :: IEL
    
    REAL(KIND=8) :: RAND_NUM
    !-------------------------------------------------------------------
    
    DO IEL=1, NEL_LOCAL ! ----------------------------------------------
        CALL RANDOM_NUMBER(RAND_NUM)
        
        IF ((RAND_NUM <= 0.3) .AND. (LPOLY_ORDER(IEL) < 12) &
            .AND. (H_DEPTH(IEL) >= 2)) THEN ! P-REFINEMENT
            
            LPOLY_ORDER(IEL) = LPOLY_ORDER(IEL)+2
        
        ENDIF
        !---------------------------------------------------------------
    ENDDO
    !-------------------------------------------------------------------
    
    ! SET LP_SERIAL-----------------------------------------------------
!    LP_SERIAL = .TRUE.
    !-------------------------------------------------------------------
    
    
!    CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)

END SUBROUTINE PADAPT


SUBROUTINE HADAPT !-----------------------------------------------------
    !-------------------------------------------------------------------
    ! H-REFINEMENT
    ! 30% REFINE
    ! REFINEMENT DEPTH <= 4 (EACH PROCESSOR <= 256 ELEMENT)
    !-------------------------------------------------------------------
    USE MESH, ONLY: NEL_LOCAL, NELMAX, X_LOCAL, Y_LOCAL, ELE_PROC_NUM
    
    !-------------------------------------------------------------------
    INTEGER :: IEL, I
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ELE_PROC_NUM_NEW
    INTEGER, ALLOCATABLE, DIMENSION(:) :: H_DEPTH_NEW   ! ELEMENT H-REFINEMENT DEPTH
    INTEGER, ALLOCATABLE, DIMENSION(:) :: LPOLY_ORDER_NEW   ! ELEMENT POLYMOIAL ORDER, LOCAL
    INTEGER :: NEL_LOCAL_NEW    ! LOCAL ELEMENT NUMBER AFTER ADAPTATION
    
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :) :: X_NEW ! RECORD ELEMENT COORDINATE AFTER H-REFINEMENT
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :) :: Y_NEW
    REAL(KIND=8), PARAMETER :: THRESHOLD=0.7D0   ! THRESHOLD FOR ADAPTATION
    REAL(KIND=8) :: RAND_NUM  ! RANDOM NUMBER 
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: XX, YY
    !-------------------------------------------------------------------
    
    ! INITIALIZATION----------------------------------------------------
    NEL_LOCAL_NEW=1
    !-------------------------------------------------------------------
    
    ! ALLOCATE ---------------------------------------------------------
    ALLOCATE(X_NEW(4, NELMAX))
    ALLOCATE(Y_NEW(4, NELMAX))
    ALLOCATE(ELE_PROC_NUM_NEW(NELMAX))
    ALLOCATE(H_DEPTH_NEW(NELMAX))
    ALLOCATE(LPOLY_ORDER_NEW(NELMAX))
    ALLOCATE(XX(3), YY(3))
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    DO IEL=1, NEL_LOCAL ! LOOP THROUGH ALL ELEMENTS
        
        ! GENERATE A RANDOM NUMBER
        CALL RANDOM_NUMBER(RAND_NUM)
        
        ! 30% ELEMENT H REFINE
        IF (RAND_NUM <= 0.3 .AND. H_DEPTH(IEL) < 4) THEN ! H-REFINEMENT------------------------
            XX(1)=X_LOCAL(1, IEL)
            XX(2)=(X_LOCAL(1, IEL)+X_LOCAL(2, IEL))/2.0D0
            XX(3)=X_LOCAL(2, IEL)
            
            YY(1)=Y_LOCAL(1, IEL)
            YY(2)=(Y_LOCAL(2, IEL)+Y_LOCAL(3, IEL))/2.0D0
            YY(3)=Y_LOCAL(3, IEL)
            
            ! P REFINEMENT GEENERATES 4 NEW ELEMENTS--------------------
            DO I=1,4
                IF ((I == 1) .OR. (I == 3)) THEN
                    X_NEW(1, NEL_LOCAL_NEW)=XX(1)
                    X_NEW(4, NEL_LOCAL_NEW)=XX(1)
                    X_NEW(2, NEL_LOCAL_NEW)=XX(2)
                    X_NEW(3, NEL_LOCAL_NEW)=XX(2)
                ELSE
                    X_NEW(1, NEL_LOCAL_NEW)=XX(2)
                    X_NEW(4, NEL_LOCAL_NEW)=XX(2)
                    X_NEW(2, NEL_LOCAL_NEW)=XX(3)
                    X_NEW(3, NEL_LOCAL_NEW)=XX(3)
                    
                ENDIF
                
                IF ((I == 1) .OR. (I == 2)) THEN
                    Y_NEW(1, NEL_LOCAL_NEW)=YY(1)
                    Y_NEW(2, NEL_LOCAL_NEW)=YY(1)
                    Y_NEW(3, NEL_LOCAL_NEW)=YY(2)
                    Y_NEW(4, NEL_LOCAL_NEW)=YY(2)
                ELSE
                    Y_NEW(1, NEL_LOCAL_NEW)=YY(2)
                    Y_NEW(2, NEL_LOCAL_NEW)=YY(2)
                    Y_NEW(3, NEL_LOCAL_NEW)=YY(3)
                    Y_NEW(4, NEL_LOCAL_NEW)=YY(3)
                ENDIF
                
                ! NEW GLOBAL ELEMENT PROCESSOR NUMEBR
                ELE_PROC_NUM_NEW(NEL_LOCAL_NEW) = ELE_PROC_NUM(IEL)
                
                H_DEPTH_NEW(NEL_LOCAL_NEW) = H_DEPTH(IEL) + 1
                
                LPOLY_ORDER_NEW(NEL_LOCAL_NEW) = LPOLY_ORDER(IEL) 
                
                NEL_LOCAL_NEW=NEL_LOCAL_NEW+1
            
            ENDDO
            !-----------------------------------------------------------
        ELSE !------NO REFINEMENT, STORE THE ELEMENT COORDINATE IN THE NEW GLOBAL ARRAY
        
            X_NEW(:, NEL_LOCAL_NEW)=X_LOCAL(:, IEL)
            Y_NEW(:, NEL_LOCAL_NEW)=Y_LOCAL(:, IEL)
            
            ELE_PROC_NUM_NEW(NEL_LOCAL_NEW) = ELE_PROC_NUM(IEL)
            
            H_DEPTH_NEW(NEL_LOCAL_NEW) = H_DEPTH(IEL) 
            
            LPOLY_ORDER_NEW(NEL_LOCAL_NEW) = LPOLY_ORDER(IEL) 
            
            NEL_LOCAL_NEW=NEL_LOCAL_NEW+1
        
        ENDIF
        !---------------------------------------------------------------
        
    ENDDO
    !--------------------------------------------------------------------
    
    ! NEW TOTAL NUMBER OF ELEMENT---------------------------------------
    NEL_LOCAL_NEW=NEL_LOCAL_NEW-1
    NEL_LOCAL=NEL_LOCAL_NEW
    !-------------------------------------------------------------------
    
    ! UPDATE GLOBAL ELEMENT COORDINATES---------------------------------
    DEALLOCATE(X_LOCAL, Y_LOCAL)
    ALLOCATE(X_LOCAL(4, NEL_LOCAL), Y_LOCAL(4, NEL_LOCAL))
    
    X_LOCAL(:, :) = X_NEW(:, 1:NEL_LOCAL)
    Y_LOCAL(:, :) = Y_NEW(:, 1:NEL_LOCAL)
    !-------------------------------------------------------------------
    
    ! UPDATE GLOBAL PROCESSOR NUMBER------------------------------------
    DEALLOCATE(ELE_PROC_NUM)
    ALLOCATE(ELE_PROC_NUM(NEL_LOCAL))
    ELE_PROC_NUM(:) = ELE_PROC_NUM_NEW(1: NEL_LOCAL)
    !-------------------------------------------------------------------
    
    ! UPDATE GLOBAL ELEMENT H-REFINEMENT DEPTH -------------------------
    DEALLOCATE(H_DEPTH)
    ALLOCATE(H_DEPTH(NEL_LOCAL))
    H_DEPTH(:) = H_DEPTH_NEW(1 : NEL_LOCAL)
    !-------------------------------------------------------------------
    
    !UPDATE GLOBAL ELEMENT POLYNOMIAL DEGREE----------------------------
    DEALLOCATE(LPOLY_ORDER)
    ALLOCATE(LPOLY_ORDER(NEL_LOCAL))
    LPOLY_ORDER(:) = LPOLY_ORDER_NEW(1:NEL_LOCAL)
    !-------------------------------------------------------------------
    
    ! DEALLOCATE--------------------------------------------------------
    DEALLOCATE(X_NEW, Y_NEW)
    DEALLOCATE(LPOLY_ORDER_NEW)
    DEALLOCATE(ELE_PROC_NUM_NEW)
    DEALLOCATE(H_DEPTH_NEW)
    DEALLOCATE(XX, YY)
    !-------------------------------------------------------------------
    

END SUBROUTINE HADAPT !-------------------------------------------------

SUBROUTINE HPADAPT_INIT !------------------------------------------------
    !-------------------------------------------------------------------
    ! INITIALIZE POLYNOMIAL ORDER
    ! POLYNOMIAL ORDER =4, 6, 8, 10, 12, START WITH 4
    ! H-REFINEMENT DEPTH =1, 2, 3, 4, START WITH 1
    !-------------------------------------------------------------------
    USE MESH, ONLY: NEL_LOCAL
    
    ALLOCATE(LPOLY_ORDER(NEL_LOCAL))
    ALLOCATE(H_DEPTH(NEL_LOCAL))
    
    ! POLYNOMIAL ORDER START WITH 4
    LPOLY_ORDER(:) = 4
    H_DEPTH(:) = 1
    

END SUBROUTINE HPADAPT_INIT !--------------------------------------------



END MODULE ADAPT
