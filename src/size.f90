
MODULE SIZE
!-----------------------------------------------------------------------
!PARAMETERS WHEN USING MPI
!-----------------------------------------------------------------------
USE MPI
IMPLICIT NONE

INTEGER :: RANK     ! PROCESSOR RANK
INTEGER :: NUM_PROC ! NUMBER OF PROCESSOR
INTEGER :: IERROR 

LOGICAL :: LP_SERIAL = .FALSE.  ! USE SERIAL LP

CONTAINS 

SUBROUTINE START_MPI

    ! MPI INITIALIZE 
    CALL MPI_INIT(IERROR)
    
    ! START UP MPI_ENVIRONEMENT
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, RANK, IERROR)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NUM_PROC, IERROR)

END SUBROUTINE START_MPI



END MODULE SIZE
