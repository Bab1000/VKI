!------------------------------------------------------------------------------!
!> This subroutine reads the mesh for 1D and 2D domains
  SUBROUTINE domain_mesh ()

    USE mod_physical_model
    USE mod_numerics_data,               ONLY: flag_stag_line
  
    IMPLICIT NONE
  
    INTEGER :: ndim
  
    ! Number of dimensions 
    ndim = get_nb_dim (physical_model)

    ! Read mesh file 
    SELECT CASE (ndim)

      ! 1D flows 
      CASE(1) 
        CALL read_1D_mesh()

      ! 2D flows
      CASE(2)

        ! Stagnation line flow
        IF (flag_stag_line.EQV..TRUE.) THEN

           CALL read_1D_mesh()

        ELSE

           CALL read_2D_mesh()

           CALL write_2D_mesh()  

        ENDIF

      CASE DEFAULT 
        PRINT*
        WRITE(*,10)'solver_fvmcc_F90:: in mesh.F90, error in number of dimensions ...'
        PRINT*
        STOP

    END SELECT 

10  FORMAT(A)

    CONTAINS
  
      !------------------------------------------------------!
      !> This subroutine reads a 1D mesh 
      SUBROUTINE read_1D_mesh ()

        USE mod_general_data,     ONLY: nb_cells, nb_tot_cells, nb_nodes, nb_ghost_cells, xc, & 
                                      & geom_source, volumes, ios_dir
        USE mod_numerics_data,    ONLY: nb_source_terms, source_terms

        IMPLICIT NONE

        INTEGER, PARAMETER :: in1 = 10, in2 = 11
        INTEGER :: i, j, length
        INTEGER :: ios
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: x

        WRITE(*,20)'solver_fvmcc_f90:: Reading mesh from "mesh.dat" file'
        PRINT*

        ! Opening and reading the mesh file

        IF (ios_dir) THEN
            OPEN(UNIT=in1,FILE='./mesh.dat',STATUS='unknown',IOSTAT=ios)
            CALL system('cp mesh.dat mesh_init.dat')
        ELSE
            OPEN(UNIT=in1,FILE='../input/mesh/mesh.dat',STATUS='unknown',IOSTAT=ios)
            CALL system('cp ../input/mesh/mesh.dat ../input/mesh/mesh_init.dat')
        ENDIF

        IF (ios.NE.0) THEN
           PRINT*
           WRITE(*,20)'solver_fvmcc_f90:: "mesh.dat" file not found'
           PRINT*
           STOP
        ENDIF
        
        ! Number of nodes
        READ(in1,*)nb_nodes        

        ! Number of cells
        nb_cells = nb_nodes - 1

        ! Number of ghost cells (two for each boundary)
        nb_ghost_cells = 4

        ! Total number of cells (physical + ghost)
        nb_tot_cells = nb_cells + nb_ghost_cells

        ALLOCATE(x(nb_nodes), xc(nb_tot_cells), volumes(nb_tot_cells)) 

        ! Initialization
        x  = 0.d0
        xc = 0.d0
        volumes = 0.d0

        ! Node list 
        DO i = 1,nb_nodes
           READ(in1,*)x(i)
        ENDDO

        CLOSE(in1)

        ! Cell centroid location
        DO i = 1,nb_cells
           xc(i + 2)      = 0.5d0*(x(i) + x(i + 1))
           volumes(i + 2) = x(i + 1) - x(i)
        ENDDO

        ! Ghost-cells centroid locations and volumes
        ! Left boundary (wall) modified to account for the stretching next to the wall (A. Turchi)
        volumes(1) = volumes(4) 
        volumes(2) = volumes(3)

        xc(2) = xc(3) - volumes(3) 
        xc(1) = xc(2) - (volumes(2) + volumes(1))*0.5d0

        ! Right boundary
        volumes(nb_cells + 3) = volumes(nb_cells + 2) 
        volumes(nb_cells + 4) = volumes(nb_cells + 2)

        xc(nb_cells + 3) = xc(nb_cells + 2) + volumes(nb_cells + 2) 
        xc(nb_cells + 4) = xc(nb_cells + 3) + volumes(nb_cells + 2)  

        DO i = 1,nb_source_terms

           IF ((source_terms(i).EQ.'quasi1D').OR.(source_terms(i).EQ.'quasi1d')) THEN

              ALLOCATE(geom_source(nb_cells))
 
              length = LEN_TRIM(source_terms(i))
              IF (ios_dir) THEN
                  OPEN(UNIT=in2,FILE='./'//source_terms(i)(1:length)//'.dat',STATUS='unknown',IOSTAT=ios)
              ELSE
              OPEN(UNIT=in2,FILE='../input/mesh/'//source_terms(i)(1:length)//'.dat',STATUS='unknown',IOSTAT=ios)
              ENDIF

              IF (ios.NE.0) THEN
                 PRINT*
                 WRITE(*,20)'solver_fvmcc_f90:: "area.dat" or "geom.dat" file not found'
                 PRINT*
                 STOP
              ENDIF

              DO j = 1,nb_cells
                 READ(in2,*)geom_source(j)
              ENDDO

              CLOSE(in2)

           ENDIF

        ENDDO

        DEALLOCATE(x)

        WRITE(*,20)'solver_fvmcc_f90:: "mesh.dat" file read'
        PRINT*

20    FORMAT(A)

      END SUBROUTINE read_1D_mesh 

      !------------------------------------------------------!
      !> This subroutine reads a 2D mesh
      SUBROUTINE read_2D_mesh ()
  
        PRINT*
        WRITE(*,'(A)')'in mesh.F90, "read_2D_mesh ()" not implemented yet...'
        PRINT*

      END SUBROUTINE read_2D_mesh 
  
      !------------------------------------------------------!
      !> This subroutine writes elements and node location for checking purposes.
      SUBROUTINE write_2D_mesh ()
  
        PRINT*
        WRITE(*,'(A)')'in mesh.F90, "write_2D_mesh ()" not implemented yet...'
        PRINT*

      END SUBROUTINE write_2D_mesh 
  
    END SUBROUTINE domain_mesh
!------------------------------------------------------------------------------!

