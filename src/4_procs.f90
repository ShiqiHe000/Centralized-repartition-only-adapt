

PROGRAM FOUR_PROCS
    USE SIZE
    USE MESH
    USE ADAPT
    USE FIELDS
    USE LOAD_PARTITION
    USE MPI
    

    IMPLICIT NONE
    
    !--------------------------------------------------------------------
    CHARACTER(LEN=*), PARAMETER :: MESHFILE = "mesh.rea"
    
    INTEGER :: I
    !-------------------------------------------------------------------
    
    ! START MPI
    CALL START_MPI
    
    ! READ MESH FILE
    CALL READ_MESH(MESHFILE)
    
    ! INITIALIZE ELEMENT POLYNOMIAL ORDER
    CALL HPADAPT_INIT
    
    ! WRITE DATA TO FILE
    CALL WRITE_FIELDS

    DO I=1, 45
        ! H-REFINEMENT
        CALL HADAPT
        
        ! P-REFINEMENT
        CALL PADAPT
        
        CALL DEALLOCATE_GLOBAL
        
        ! WRITE DATA TO FILE
        CALL WRITE_FIELDS
        
        CALL H1_H2_PARTITION_SERIAL
      
        CALL WRITE_FIELDS

    
    ENDDO

    
    CALL MPI_FINALIZE(IERROR)

    IF (RANK == 0) THEN
        PRINT *,  "FINISH WRITING DATA TO FILE"
    ENDIF
    

END PROGRAM FOUR_PROCS


