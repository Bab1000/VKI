!------------------------------------------------------------------------------!
!> This module contains subrouroutines and functions for dealing with radiation.
  MODULE mod_radiation

  IMPLICIT NONE
  LOGICAL, SAVE :: flag_radiation  = .FALSE.   !< flag indicating if radiative source terms are added to the computation
  LOGICAL, SAVE :: flag_precursor  = .TRUE.   !< flag indicating if radiative source are included in the free-stream to compute precursor
  LOGICAL, SAVE :: flag_photochemistry  = .TRUE.   !< flag indicating if photochemistry source terms are included
  LOGICAL, SAVE :: flag_restart_rad  = .FALSE. !< flag indicating a restart for radiation
  LOGICAL, SAVE :: flag_automatic_rad_mesh = .FALSE. !< flag indicating if an automatic mesh is built for radiation
  INTEGER, SAVE :: coupling_period = 1000      !< coupling period (default value)
  LOGICAL, SAVE :: relax = .FALSE.             !< relaxation option (default value)
  REAL(KIND=8), SAVE :: relax_factor = 1.d0   !< relaxation factor (default value)
  INTEGER, SAVE :: nb_cells_rad, nb_nodes_rad, nb_var  
  INTEGER       :: it_rad, i_ttv=0
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: x_rad, xc_rad, vol_rad !< Radiation mesh data
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: radnm1   !< Radiative source term (flowfield mesh)
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: rad      !< Radiative source term (flowfield mesh)
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: rad_old      !< Radiative source term (flowfield mesh)
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: rad_new     !< Radiative source term (flowfield mesh)
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: rad_     !< Radiative source term (radiation mesh)
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: f, f_    !< Flow field data for radiation
  REAL(KIND=8), SAVE:: wall_Rad, relaxation_parameter

  ! Subroutine and functions
  CONTAINS

    !------------------------------------------------------------------------------!
    !> This subroutine initializes the radiative source terms
    SUBROUTINE initialize_radiation ()

      USE mod_general_data,             ONLY: nb_ns, nb_eq, nb_cells, nb_temp

      IMPLICIT NONE

      INTEGER i, j, pos_var
      REAL temp

      ALLOCATE(rad(nb_eq*nb_cells))
      ALLOCATE(rad_old(nb_eq*nb_cells))
      ALLOCATE(rad_new(nb_eq*nb_cells))
      rad = 0.d0
      rad_old = 0.d0
      rad_new = 0.d0
      relaxation_parameter = 1.0
      
      !ALLOCATE(radnm1(nb_eq*nb_cells))
      !radnm1 = 0.d0
      it_rad = 0

      IF (flag_radiation) THEN
        nb_var = nb_ns + 3 
        !ALLOCATE(f(nb_var*nb_cells))
        IF (nb_temp > 1) i_ttv = 1
!        IF (flag_restart_rad) THEN
!           CALL read_1D_mesh_rad ()
!        ELSE
           !CALL build_1D_regular_mesh_rad(50)
!        ENDIF
        ! Allocate memory for interpolations
        !ALLOCATE(f_(nb_var*nb_cells_rad))
        !ALLOCATE(rad_(nb_eq*nb_cells_rad))
      ENDIF

      !IF (flag_radiation .AND. flag_restart_rad) THEN 
        !Read previous radiation field if restart
       ! OPEN(1,FILE='restart_rad.dat')
       ! DO i=1,nb_cells
       !  pos_var = (i-1)*nb_eq
       !  WRITE(1,'(100ES20.10)') (rad(pos_var+j),j=1,nb_eq)
       !ENDDO
       !CLOSE(1)
     !ENDIF

    END SUBROUTINE initialize_radiation

    !------------------------------------------------------------------------------!
    !> This subroutine gives the radiative source term at the icell index according
    !to the current state of the rad vector
    SUBROUTINE get_rad_source_term (icell, s)

      USE mod_general_data,            ONLY: nb_eq, &
                                           & cell_prop, start_prop_phys, nb_prop, pos_T_cell, &
                                           & shockIndex, shockPosition 

      IMPLICIT NONE

      INTEGER :: pos_rad, j
      !REAL(KIND=8) :: Tcell

      INTEGER, INTENT(IN) :: icell                           !< cell location index
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s           !< source term

      pos_rad = (icell-1)*nb_eq
      !Tcell = cell_prop(start_prop_phys + (icell-1)*nb_prop + pos_T_cell)

      DO j=1,nb_eq
        s(j) = rad(pos_rad+j)
          IF ((flag_precursor .EQV. .FALSE.) .AND.( icell >= shockIndex)) THEN
            s(j) = 0.0;
          ENDIF
      ENDDO
!      IF (s(nb_eq)<0.d0 .AND. Tcell < 200) THEN
!         s(nb_eq) =0.d0
!      ENDIF

    END SUBROUTINE get_rad_source_term

    SUBROUTINE compute_rad_source_term ()

      USE mod_general_data,            ONLY: nb_ns, nb_cells, nb_eq, nb_prop, start_u_phys, start_prop_phys, pos_T_cell, &
                                           & pos_pres_cell, rhoi, xi, xc, u, cell_prop, simulation_id
      USE mod_neq_function_pointer,    ONLY: library_get_molar_fractions

      IMPLICIT NONE

      INTEGER :: pos_phys, pos_u, pos_var, i, j

      !CALL update_1D_regular_mesh_rad()
      !
      it_rad = it_rad + 1
      !radnm1 = rad

      !Get the primitive variables for radiation: T, P, xi
      !f=0.d0
      !DO i = 1, nb_cells
      !  pos_phys = start_prop_phys + (i-1)*nb_prop
      !  pos_u    = start_u_phys + (i-1)*nb_eq
      !  pos_var  = (i-1)*nb_var
      !
      !  f(pos_var+1) = cell_prop(pos_phys+pos_T_cell)     ! Tr
      !  f(pos_var+2) = cell_prop(pos_phys+pos_T_cell+i_ttv)     ! Tv
      !  f(pos_var+3) = cell_prop(pos_phys+pos_pres_cell)  ! pressure

      !  DO j=1,nb_ns
      !    rhoi(j) = u(pos_u+j)
      !  ENDDO
      !  CALL library_get_molar_fractions (rhoi, f(pos_var+4:pos_var+nb_var))
      !ENDDO

      !OPEN(1000+it_rad)
      !DO i = 1, nb_cells
      !  pos_var  = (i-1)*nb_var
      !  WRITE(1000+it_rad,'(100ES20.10)') xc(2+i), (f(pos_var+j),j=1,nb_var)
      !ENDDO
      !CLOSE(1000+it_rad)

      !Interpolation on the radiation mesh and writing of the radiation input file
      !CALL flow_to_rad ()

      !Call the radiation code
      !CALL system('./rte1d radiation_input')
      call system('./update_radiation')

      !Interpolate the radiative source terms on the convective mesh
      CALL rad_to_flow()

      !IF (flag_automatic_rad_mesh) CALL build_1D_adaptive_mesh_rad ()

      !Computation of the radiative source term in case of under relaxation
      !IF (relax) THEN
      !  rad = relax_factor*rad + (1.d0-relax_factor)*radnm1
      !ENDIF

      !OPEN(2000+it_rad)
      !DO i=1,nb_cells
      !  pos_var = (i-1)*nb_eq
      !  WRITE(2000+it_rad,'(100ES20.10)') xc(2+i), (rad(pos_var+j),j=1,nb_eq)
      !ENDDO
      !CLOSE(2000+it_rad)

      !Move the result to the output folder
      CALL system('cp ./rad_source_terms.dat '//TRIM(simulation_id)//'_source_term.dat')

      !Write results for a restart
      OPEN(1,FILE=TRIM(simulation_id)//'_restart_rad.dat')
      DO i=1,nb_cells
        pos_var = (i-1)*nb_eq
        WRITE(1,'(100ES20.10)') (rad(pos_var+j),j=1,nb_eq)
      ENDDO
      CLOSE(1)

    END SUBROUTINE compute_rad_source_term

    SUBROUTINE read_1D_mesh_rad ()

      IMPLICIT NONE

      INTEGER, PARAMETER :: inp = 100
      INTEGER :: ios, i

      WRITE(*,*)'solver_fvmcc_f90:: Reading radiation mesh from "mesh_rad.dat" file'
      PRINT*

      ! Opening and reading the mesh file
      OPEN(UNIT=inp,FILE='mesh_rad.dat',STATUS='unknown',IOSTAT=ios)

      IF (ios.NE.0) THEN
         PRINT*
         WRITE(*,*)'solver_fvmcc_f90:: "mesh_rad.dat" file not found'
         PRINT*
         STOP
      ENDIF

      ! Number of nodes and cells
      READ(inp,*)nb_nodes_rad
      nb_cells_rad = nb_nodes_rad - 1

      ALLOCATE(x_rad(nb_nodes_rad), xc_rad(nb_cells_rad), vol_rad(nb_cells_rad))

      ! Initialization
      x_rad  = 0.d0
      xc_rad = 0.d0
      vol_rad = 0.d0

      ! Node list 
      DO i = 1,nb_nodes_rad
         READ(inp,*)x_rad(i)
      ENDDO

      CLOSE(inp)

      ! Cell centroid location and volume
      DO i = 1,nb_cells_rad
         xc_rad(i)  = 0.5d0*(x_rad(i) + x_rad(i+1))
         vol_rad(i) = x_rad(i+1) - x_rad(i)
      ENDDO

    END SUBROUTINE read_1D_mesh_rad

    SUBROUTINE build_1D_regular_mesh_rad (np)

      USE mod_general_data,     ONLY: nb_cells, xc, volumes

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: np

      INTEGER, PARAMETER :: inp = 100
      INTEGER :: i,j, mesh_ratio, remainder

      nb_cells_rad = np
      nb_nodes_rad = nb_cells_rad + 1
      ALLOCATE(x_rad(nb_nodes_rad), xc_rad(nb_cells_rad), vol_rad(nb_cells_rad))

      mesh_ratio = nb_cells / nb_cells_rad
      remainder = mod(nb_cells,nb_cells_rad)

      x_rad  = 0.d0
      xc_rad = 0.d0
      vol_rad = 0.d0
      
      DO i = 1, nb_cells_rad
        DO j = (i-1)*mesh_ratio+1,i*mesh_ratio 
          vol_rad (i) = vol_rad(i) + volumes(j+2)
        ENDDO
      ENDDO
      IF(remainder .NE. 0) THEN
        DO j = nb_cells_rad*mesh_ratio+1, nb_cells
          vol_rad(nb_cells_rad) = vol_rad(nb_cells_rad) + volumes(j+2)
        ENDDO
      ENDIF

      x_rad(1) = xc(3)-(volumes(3)/2.d0) 
      DO i = 1, nb_cells_rad
        x_rad(i+1) = x_rad(i)+vol_rad(i)
        xc_rad(i) = (x_rad(i+1)+x_rad(i))/2.d0
      ENDDO

      ! Writes the radiation mesh
      CALL system('rm -f mesh_rad.dat')
      OPEN(UNIT=inp,FILE='mesh_rad.dat')
      WRITE(inp,*) nb_nodes_rad
      DO i = 1,nb_nodes_rad
         WRITE(inp,*) x_rad(i)
      ENDDO
      CLOSE(inp)

    END SUBROUTINE

    SUBROUTINE update_1D_regular_mesh_rad ()

      USE mod_general_data,     ONLY: nb_cells, xc, volumes

      IMPLICIT NONE

      INTEGER, PARAMETER :: inp = 100
      INTEGER :: i,j, mesh_ratio, remainder

      mesh_ratio = nb_cells / nb_cells_rad
      remainder = mod(nb_cells,nb_cells_rad)

      x_rad  = 0.d0
      xc_rad = 0.d0
      vol_rad = 0.d0
      
      DO i = 1, nb_cells_rad
        DO j = (i-1)*mesh_ratio+1,i*mesh_ratio 
          vol_rad (i) = vol_rad(i) + volumes(j+2)
        ENDDO
      ENDDO
      IF(remainder .NE. 0) THEN
        DO j = nb_cells_rad*mesh_ratio+1, nb_cells
          vol_rad(nb_cells_rad) = vol_rad(nb_cells_rad) + volumes(j+2)
        ENDDO
      ENDIF

      x_rad(1) = xc(3)-(volumes(3)/2.d0) 
      DO i = 1, nb_cells_rad
        x_rad(i+1) = x_rad(i)+vol_rad(i)
        xc_rad(i) = (x_rad(i+1)+x_rad(i))/2.d0
      ENDDO

      ! Writes the radiation mesh
      CALL system('rm -f mesh_rad.dat')
      OPEN(UNIT=inp,FILE='mesh_rad.dat')
      WRITE(inp,*) nb_nodes_rad
      DO i = 1,nb_nodes_rad
         WRITE(inp,*) x_rad(i)
      ENDDO
      CLOSE(inp)

    END SUBROUTINE

    SUBROUTINE build_1D_adaptive_mesh_rad ()

      USE mod_general_data,             ONLY: nb_ns, nb_eq, nb_cells, nb_temp
      USE mod_adapt,                    ONLY: adapt_mesh

      IMPLICIT NONE
      INTEGER, PARAMETER :: inp = 100
      INTEGER :: i

      ! Call the mesh adaptor for optimizing the mesh according to the radiative power
      CALL adapt_mesh(x_rad, rad_(nb_ns+3:nb_eq*nb_cells_rad:nb_eq))

      ! Update cell centers and volumes
      DO i = 1, nb_cells_rad
        xc_rad(i) = (x_rad(i+1)+x_rad(i))/2.d0
        vol_rad(i) = x_rad(i+1) - x_rad(i)
      ENDDO

      ! Writes the radiation mesh
      CALL system('rm -f mesh_rad.dat')
      OPEN(UNIT=inp,FILE='mesh_rad.dat')
      WRITE(inp,*) nb_nodes_rad
      DO i = 1,nb_nodes_rad
         WRITE(inp,*) x_rad(i)
      ENDDO
      CLOSE(inp)

    END SUBROUTINE

!    SUBROUTINE build_1D_adaptive_mesh_rad ()
!
!      USE mod_general_data,     ONLY: nb_eq, nb_cells, xc, volumes
!
!      IMPLICIT NONE
!      
!      INTEGER :: i, j, posT, cells
!      INTEGER, DIMENSION(nb_cells) :: ind_tab
!      REAL(KIND=8) :: Tmax
!      REAL(KIND=8) :: tol_dT, tol_d2T 
!      INTEGER :: tolC = 14
!      INTEGER, PARAMETER :: inp = 100
!
!      IF (ALLOCATED(x_rad)) DEALLOCATE(x_rad)
!      IF (ALLOCATED(xc_rad)) DEALLOCATE(xc_rad)
!      IF (ALLOCATED(vol_rad)) DEALLOCATE(vol_rad)
!      IF (ALLOCATED(f_)) DEALLOCATE(f_)
!      IF (ALLOCATED(rad_)) DEALLOCATE(rad_)
!
!      Tmax = maxval(f(1:nb_var*nb_cells:nb_var))
!      tol_dT = Tmax / 10.d0
!      tol_d2T = Tmax / 500.d0
!
!      ! Build a table indices for merging cells according to a temperature delta
!      ! and a maximum number of cells 
!      ind_tab(1) = 1
!      j = 1
!      cells = 1
!      DO i = 2, nb_cells-1
!        posT  = (i-1)*nb_var+1
!        !IF (f(posT).EQ.Tmax) tol_d2T = Tmax / 100.d0
!        IF ( (abs(f(posT)-f(posT-nb_var)) .GE. tol_dT) &
!           &  .OR. (abs(f(posT+nb_var) + f(posT-nb_var) - 2*f(posT)) .GE. tol_d2T) &
!           &  .OR. (cells .GE. tolC) ) THEN
!          j = j+1
!          cells =1
!        ELSE
!          cells = cells+1
!        ENDIF
!        ind_tab(i) = j
!        !PRINT *, "TEST AUTOMATIC MESH ", i, ind_tab(i)
!      ENDDO
!      ind_tab(nb_cells) = ind_tab(nb_cells-1)
!
!      ! Build the radiation mesh according to the table indices
!      nb_cells_rad = ind_tab(nb_cells)
!      nb_nodes_rad = nb_cells_rad + 1
!      ALLOCATE(x_rad(nb_nodes_rad), xc_rad(nb_cells_rad), vol_rad(nb_cells_rad))
!      ALLOCATE(f_(nb_var*nb_cells_rad), rad_(nb_eq*nb_cells_rad))
!      x_rad  = 0.d0
!      xc_rad = 0.d0
!      vol_rad = 0.d0
!
!      DO j = 1, nb_cells
!        vol_rad(ind_tab(j)) = vol_rad(ind_tab(j)) + volumes(j+2)
!      ENDDO
!
!      x_rad(1) = xc(3)-(volumes(3)/2.d0) 
!      DO i = 1, nb_cells_rad
!        x_rad(i+1) = x_rad(i)+vol_rad(i)
!        xc_rad(i) = (x_rad(i+1)+x_rad(i))/2.d0
!      ENDDO
!
!      ! Writes the radiation mesh
!      OPEN(UNIT=inp,FILE='../input/mesh/mesh_rad.dat')
!      WRITE(inp,*) nb_nodes_rad
!      DO i = 1,nb_nodes_rad
!         WRITE(inp,*) x_rad(i)
!      ENDDO
!
!    END SUBROUTINE

    SUBROUTINE flow_to_rad()

      USE mod_general_data,     ONLY: xc, volumes, nb_ns
      use mod_neq_function_pointer, ONLY : library_species_name

      IMPLICIT NONE
      INTEGER :: i, j, k, pos_var, pos_var_
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: tmp
      CHARACTER(len=20) :: species_names(nb_ns)
      
      !Interpolation on the coarse radiation mesh by volume averaging
      !All the points of the radiation mesh belongs to the flow mesh
      ALLOCATE(tmp(nb_var))
      j=1
      pos_var = 0
      tmp = 0.d0
      !DO i=1,nb_cells_rad
      i = 1
      DO WHILE (i .LE. nb_cells_rad)
        DO WHILE (x_rad(i+1) > xc(j+2)+0.5*volumes(j+2))
          tmp(:) = tmp(:) + f(pos_var+1:pos_var+nb_var)*volumes(j+2)
          j=j+1
          pos_var = (j-1)*nb_var
        ENDDO
        ! Treatment of flow cells that cross a radiation cell
        !Left contribution
        tmp(:) = tmp(:) + f(pos_var+1:pos_var+nb_var)*(x_rad(i+1)-(xc(j+2)-0.5*volumes(j+2)))
        pos_var_ = (i-1)*nb_var
        f_(pos_var_+1:pos_var_+nb_var) = tmp(:) / vol_rad(i)

        i = i+1
        ! Radiation cells contained within a flow cell
        DO WHILE ((x_rad(i+1) < xc(j+2)+0.5*volumes(j+2)) .AND. (i.LE.nb_cells_rad))
          pos_var_ = (i-1)*nb_var
          f_(pos_var_+1:pos_var_+nb_var) = f(pos_var+1:pos_var+nb_var)
          i = i+1
        ENDDO

        !Right contribution
        tmp = 0.d0
        tmp(:) = tmp(:) + f(pos_var+1:pos_var+nb_var)*(xc(j+2)+0.5*volumes(j+2)-x_rad(i))
        j=j+1
        pos_var = (j-1)*nb_var
        
      ENDDO
      DEALLOCATE(tmp)
      
      f_ = max(f_, 0.0)
      

      !Writing of the radiation input file: heading + flowfield + system_list
      !WARNING length in cm in the radiation solver
      CALL system('rm -f radiation_input')

      OPEN(5,FILE='radiation_input')
      WRITE(5,*) 'Number of points'
      WRITE(5,*) nb_nodes_rad
      CLOSE(5)

      CALL system('cat radiation_heading >> radiation_input')

      OPEN(5,FILE='radiation_input',POSITION='append')
      write(5,'(A20)',advance='no') "#Loc"
      write(5,'(A20)',advance='no') "Tr"
      write(5,'(A20)',advance='no') "Tv"
      write(5,'(A20)',advance='no') "P"
      call library_species_name(species_names)
      DO i = 1, nb_ns
         write(5,'(A20)',advance='no') trim(species_names(i))
      ENDDO
      write(5,*)
      DO i = 1, nb_cells_rad
        pos_var = (i-1)*nb_var
        WRITE(5,'(100ES20.10)') x_rad(i)*1.d2, (f_(pos_var+j),j=1,nb_var)
      ENDDO
      WRITE(5,'(ES20.10)') x_rad(nb_nodes_rad)*1.d2
      WRITE(5,*)
      CLOSE(5)

      CALL system('cat radiation_systems >> radiation_input')

    END SUBROUTINE flow_to_rad

    SUBROUTINE rad_to_flow ()

      USE mod_general_data,     ONLY: nb_ns, nb_eq, nb_cells, xc, cell_prop, start_prop_phys, pos_T_cell


      IMPLICIT NONE
      INTEGER :: posS, posT, posjm1, posj, posjp1, posjp2, posi, i, j, ios
      REAL(KIND=8) :: temp
      REAL(KIND=8) :: Twall_rad, flux_Surface, x_radFlux, emiss, wall_reRad
      REAL(KIND=8), PARAMETER :: sigma = 5.670373D-8
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: a1, a2, a3, a4
      !Open the flux file to obtain the radiation into the surface
      OPEN(3,FILE='flux.dat', IOSTAT = ios)
      IF (ios /= 0) THEN
        WRITE(*,*) "File flux does not exist"
        STOP
      ENDIF

      READ(3,*) x_radFlux, flux_Surface

      CLOSE(3)
      emiss = 0.85D0
      x_radFlux = x_radFlux/100.0
      flux_Surface = flux_Surface*1.0D4
      Twall_rad = cell_prop(start_prop_phys + pos_T_cell)
      wall_reRad = emiss * sigma *Twall_rad* Twall_rad* Twall_rad* Twall_rad
      wall_Rad = abs(flux_Surface - wall_reRad)/emiss

      ! First count the number of lines in the file
      OPEN(2,FILE='rad_source_terms.dat')
      nb_cells_rad = 0
      DO WHILE (.TRUE.)
          read(2,*, iostat=ios) temp
          if (ios /= 0) exit
          nb_cells_rad = nb_cells_rad + 1
      END DO
      REWIND(2)
      
      IF (allocated(xc_rad)) deallocate(xc_rad)
      IF (allocated(rad_))   deallocate(rad_)
      
      allocate(xc_rad(nb_cells_rad))
      allocate(rad_(nb_cells_rad*nb_eq))
 
      rad_ = 0.d0
      rad  = 0.d0
      ! Read the results from the radiation code
      ! Reminder: no radiative source term in the momentum equations
      IF (i_ttv .EQ. 1) THEN
        DO i=1,nb_cells_rad
          posS = (i-1)*nb_eq
          posT = posS + nb_ns + 2
          READ(2,*) xc_rad(i), (rad_(posS+j),j=1,nb_ns), rad_(posT+1), rad_(posT+2)

          ! Temporarily make photochemistry source terms zero
          IF (flag_photochemistry .EQV. .FALSE.) rad_(posS+1:posS+nb_ns) = 0.0;
      
        ENDDO
      ELSE
        DO i=1,nb_cells_rad
          posS = (i-1)*nb_eq
          posT = posS + nb_ns + 2
          READ(2,*) xc_rad(i), (rad_(posS+j),j=1,nb_ns), rad_(posT+1)
        ENDDO
      ENDIF
      CLOSE(2)

      ! Unit conversion (cm-3 +> m-3)
      rad_ = rad_ * 1.0d6
      
      ! Unit conversion (cm -> m)
      xc_rad = xc_rad / 100.0d0

!      ! Lagrange third order interpolation on the flow mesh
!      ALLOCATE(a1(nb_eq),a2(nb_eq),a3(nb_eq),a4(nb_eq))
!      j=2
!      posjm1 = (j-2)*nb_eq
!      posj   = (j-1)*nb_eq
!      posjp1 =  j   *nb_eq
!      posjp2 = (j+1)*nb_eq
!      DO i=1,nb_cells
!        posi = (i-1)*nb_eq
!        !Interpolation coefficients
!        a1(:)=rad_(posjm1+1:posjm1+nb_eq)&
!          &/(xc_rad(j-1)-xc_rad(j  ))/(xc_rad(j-1)-xc_rad(j+1))/(xc_rad(j-1)-xc_rad(j+2))
!        a2(:)=rad_(posj  +1:posj  +nb_eq)&
!          &/(xc_rad(j  )-xc_rad(j-1))/(xc_rad(j  )-xc_rad(j+1))/(xc_rad(j  )-xc_rad(j+2))
!        a3(:)=rad_(posjp1+1:posjp1+nb_eq)&
!          &/(xc_rad(j+1)-xc_rad(j-1))/(xc_rad(j+1)-xc_rad(j  ))/(xc_rad(j+1)-xc_rad(j+2))
!        a4(:)=rad_(posjp2+1:posjp2+nb_eq)&
!          &/(xc_rad(j+2)-xc_rad(j-1))/(xc_rad(j+2)-xc_rad(j  ))/(xc_rad(j+2)-xc_rad(j+1))
!
!        rad(posi+1:posi+nb_eq) = &
!          &  a1(:)*(xc(i+2)-xc_rad(j  ))*(xc(i+2)-xc_rad(j+1))*(xc(i+2)-xc_rad(j+2)) &
!          &+ a2(:)*(xc(i+2)-xc_rad(j-1))*(xc(i+2)-xc_rad(j+1))*(xc(i+2)-xc_rad(j+2)) &
!          &+ a3(:)*(xc(i+2)-xc_rad(j-1))*(xc(i+2)-xc_rad(j  ))*(xc(i+2)-xc_rad(j+2)) &
!          &+ a4(:)*(xc(i+2)-xc_rad(j-1))*(xc(i+2)-xc_rad(j  ))*(xc(i+2)-xc_rad(j+1)) 
!
!        !Change interpolation interval
!        IF ( (xc(i+3) > xc_rad(j+1)).AND.(j<nb_cells_rad-2) ) THEN
!          j=j+1
!          posjm1 = (j-2)*nb_eq
!          posj   = (j-1)*nb_eq
!          posjp1 =  j   *nb_eq
!          posjp2 = (j+1)*nb_eq
!        ENDIF
!      ENDDO
!      DEALLOCATE(a1,a2,a3,a4)

      ! Linear interpolation on the flow mesh
      ALLOCATE(a1(nb_eq),a2(nb_eq))
      j=1
      posj   = (j-1)*nb_eq
      posjp1 =  j   *nb_eq
      DO i=1,nb_cells
        !Change interpolation interval
        DO WHILE ( (xc(i+2) > xc_rad(j+1)).AND.(j<nb_cells_rad-1) )
          j=j+1
          posj   = (j-1)*nb_eq
          posjp1 =  j   *nb_eq
        ENDDO

        posi = (i-1)*nb_eq
        a1(:)=rad_(posj  +1:posj  +nb_eq)/(xc_rad(j+1)-xc_rad(j))
        a2(:)=rad_(posjp1+1:posjp1+nb_eq)/(xc_rad(j+1)-xc_rad(j))

        rad_new(posi+1:posi+nb_eq) = a2(:)*(xc(i+2)-xc_rad(j  )) &
                             & - a1(:)*(xc(i+2)-xc_rad(j+1))
      ENDDO
      rad= (1.0-relaxation_parameter)*rad_old+relaxation_parameter*rad_new
      DEALLOCATE(a1,a2)
      rad_old = rad
    END SUBROUTINE rad_to_flow

  END MODULE mod_radiation
