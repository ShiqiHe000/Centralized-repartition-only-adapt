
MODULE MESH

USE MPI
USE SIZE

IMPLICIT NONE

INTEGER, ALLOCATABLE, DIMENSION(:) :: ELE_PROC_NUM   ! ELEMENT PROCCESSOR NUMBER
INTEGER, PARAMETER :: NELMAX=256       ! MAX NUMBER OF ELEMENT EACH PROCESSOR
INTEGER :: NEL_LOCAL ! NUMBER OF ELEMENT IN EACH PROC
INTEGER :: NEL_TOTAL                ! TOTAL NUMBER OF ELEMENTS
 
REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :) :: X_LOCAL, Y_LOCAL ! NODE COORDINATES 

CONTAINS

SUBROUTINE READ_MESH(MESHFILE)
    IMPLICIT NONE
    
    ! VARIABLES---------------------------------------------------------
    INTEGER :: NDIM                     ! DIMENSION
    INTEGER :: NGROUPS                  ! TOTAL GROUP NUMBER
    INTEGER :: IEL, IC     
    INTEGER :: GROUP  
    INTEGER, ALLOCATABLE, DIMENSION(:) :: NEL_GLOBAL ! ELEMENT NUMBER IN EACH PROC 
    
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :) :: X_GLOBAL, Y_GLOBAL ! NODE COORDINATES 
    
    CHARACTER(LEN=*), INTENT(IN) :: MESHFILE
    !-------------------------------------------------------------------

    ! ONLY PROC 0 READ--------------------------------------------------
    IF (RANK ==0) THEN
    
        ! READ MESH FILE OPENING--------------------------------------------
        WRITE(*,*)'Opening .rea file: ', MESHFILE
        OPEN (UNIT=1,FILE=MESHFILE,STATUS='OLD')
        
        WRITE(*,'(32A)')'Reading .rea file...'
        
        READ(1,*)! READ dummy line

        READ(1,*)  NEL_TOTAL,NDIM,NGROUPS !Read mesh size and dimension
        
        
        ! CHECK NGROUP <= NPROC --------------------------------------------
        IF (NGROUPS.GT.NUM_PROC)THEN
            PRINT '("ERROR: Mesh file has ",I0," groups and requires at least ",&
                    & I0," processors. Only ", I0," launched")', NGROUPS,NGROUPS,NUM_PROC
            CALL MPI_FINALIZE(IERROR)
            STOP 
        ENDIF
        !--------------------------------------------------------------------
        
        ! ALLOCATE ----------------------------------------------------------
        ALLOCATE(NEL_GLOBAL(NUM_PROC))         ! ELEMENT INDEX IN EACH PROC
        ALLOCATE(X_GLOBAL(4,NEL_TOTAL))
        ALLOCATE(Y_GLOBAL(4,NEL_TOTAL))
        !--------------------------------------------------------------------
        
        ! INITIALIZATION-----------------------------------------------------
        NEL_GLOBAL=0.0D0
        X_GLOBAL=0.0D0
        Y_GLOBAL=0.0D0
        !-------------------------------------------------------------------
        
        
        ! CHECK-------------------------------------------------------------
        IF (NEL_TOTAL .GT. NELMAX*NUM_PROC) THEN
            PRINT *, "TOTAL NUMBER OF ELEMENT EXCEEDS THE PROCESSOR ABILITY."
            STOP
        ENDIF
        !--------------------------------------------------------------------
     
        ! READ THE ELEMENT COORDINATES--------------------------------------
        DO IEL=1,NEL_TOTAL
            READ(1,30) GROUP    !Read group number of current element
            
            ! ELEMENT NUMBER IN EACH PROC
            NEL_GLOBAL(GROUP) = NEL_GLOBAL(GROUP)+1 

!            ! GLOBAL PROCESSOR NUMBER
!            PROC_GLOBAL(IEL)= GROUP
            
            READ(1,*) (X_GLOBAL(IC, IEL),IC=1,4) !Read x coordinates

            READ(1,*) (Y_GLOBAL(IC, IEL),IC=1,4) !Read y coordinates
            
    30  FORMAT(43X,I5)

        ENDDO
        !-------------------------------------------------------------------
        
        ! CHECK --------------------------------------------------------------
        DO IC=1,NUM_PROC
            IF (NEL_GLOBAL(IC).GT.NELMAX) THEN
                PRINT '("Number of elements in .rea file is too large (Elem/proc)")'
                CALL MPI_FINALIZE(IERROR)
                STOP 
            ENDIF
        ENDDO
        !--------------------------------------------------------------------
    
    
    ENDIF
    !----------------------------------------------------------------------
    
    ! BROADCAST INFORMATION---------------------------------------------
    CALL MPI_BCAST(NDIM, 1, MPI_INT, 0, MPI_COMM_WORLD, IERROR)
    CALL MPI_BCAST(NEL_TOTAL, 1, MPI_INT, 0, MPI_COMM_WORLD, IERROR)
    CALL MPI_SCATTER(NEL_GLOBAL, 1, MPI_INT, NEL_LOCAL, 1, MPI_INT, 0, MPI_COMM_WORLD, IERROR)
    !-------------------------------------------------------------------
    
    ! SCATTER COORDINATE------------------------------------------------
    ALLOCATE(X_LOCAL(4, NEL_TOTAL/NUM_PROC))
    ALLOCATE(Y_LOCAL(4, NEL_TOTAL/NUM_PROC))
    CALL MPI_SCATTER(X_GLOBAL, 4, MPI_DOUBLE_PRECISION, X_LOCAL, 4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERROR)
    CALL MPI_SCATTER(Y_GLOBAL, 4, MPI_DOUBLE_PRECISION, Y_LOCAL, 4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERROR)
    !-------------------------------------------------------------------
        
    ! RECORD ELEMENT PROCESS NUMBER-------------------------------------
    CALL GET_PROC_NUMBER
    !-------------------------------------------------------------------
        

    ! DEALLOCATE--------------------------------------------------------
    IF (RANK == 0) THEN
        DEALLOCATE(X_GLOBAL, Y_GLOBAL)
    ENDIF
    !-------------------------------------------------------------------

    ! SYNCHRONIZE-------------------------------------------------------
!    CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !-------------------------------------------------------------------
    
END SUBROUTINE READ_MESH

SUBROUTINE GET_PROC_NUMBER

    INTEGER :: I
    
    ! GET ELEMENT PROCESSOR NUMBER--------------------------------------

    ALLOCATE(ELE_PROC_NUM(NEL_LOCAL))
    ELE_PROC_NUM=RANK+1
    !-------------------------------------------------------------------

END SUBROUTINE GET_PROC_NUMBER


END MODULE MESH
