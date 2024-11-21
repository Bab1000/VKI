!> This subroutine reads the initial field starting from a previously computed solution.
!! For sake of simplicity only primitive variables are read.
  SUBROUTINE read_solution() 

    USE mod_general_data,             ONLY: nb_ns, nb_temp, nb_dim, nb_eq, nb_cells, start_u_phys, R_gas, & 
                                          & p, simulation_id, ios_dir
    USE mod_numerics_data,            ONLY: flag_restart_pg, flag_restart_eq
    USE mod_neq_function_pointer,     ONLY: library_compute_eq_composition

    IMPLICIT NONE

    INTEGER :: ios
    INTEGER :: i, i1, is
    INTEGER, PARAMETER :: in1 = 10
    REAL(KIND=8), PARAMETER :: rhoi_tol = 1.d-15
    REAL(KIND=8), PARAMETER :: T_min = 500.d0
    REAL(KIND=8) :: pres, dens, temp
    REAL(KIND=8), DIMENSION(nb_dim) :: u_vec
    REAL(KIND=8), DIMENSION(nb_ns) :: rhoi
    REAL(KIND=8), DIMENSION(nb_eq) :: prim_var

    ! Opening the restart.dat file
    IF (ios_dir) THEN
        OPEN(UNIT=in1,FILE='./'//TRIM(simulation_id)//'_restart.dat',STATUS='old',IOSTAT=ios)
    ELSE
    	OPEN(UNIT=in1,FILE='../output/restart.dat',STATUS='old',IOSTAT=ios)
    ENDIF

    IF (ios.NE.0) THEN 

       WRITE(*,5)'solver_fvmcc_F90:: restart.dat file not found'
       PRINT*
       STOP

    ENDIF

    ! Index initialization
    i1 = start_u_phys

    ! Solution initialization
    p = 0.d0

    ! Storing primitive variable solution vector
    IF (flag_restart_pg.EQV..TRUE.) THEN 

       ! Starting nonequilibrium flow solution from a calorically perfect gas solution
       DO i = 1,nb_cells 

           READ(in1,*)dens,u_vec,pres

           temp = MAX(pres/(dens*r_gas),T_min)
           
           ! Equilibrium chemical composition 
           CALL library_compute_eq_composition (pres, temp, rhoi)

           ! Species densities
           DO is = 1,nb_ns 
              prim_var(is) = MAX(rhoi(is),rhoi_tol)
           ENDDO

           ! Velocity components
           DO is = 1,nb_dim
              prim_var(nb_ns + is) = u_vec(is)
           ENDDO

           ! Temperatures
           DO is = 1,nb_temp
              prim_var(nb_ns + nb_dim + is) = temp
           ENDDO

           DO is = 1,nb_eq 
              p(i1 + is) = prim_var(is)
           ENDDO

           ! Index update
           i1 = i1 + nb_eq

       ENDDO

    ELSEIF (flag_restart_eq.EQV..TRUE.) THEN 

       ! Starting nonequilibrium flow solution from an equilibrium flow solution
       DO i = 1,nb_cells 

           READ(in1,*)pres,u_vec,temp

           temp = MAX(temp,T_min)
           
           ! Equilibrium chemical composition 
           CALL library_compute_eq_composition (pres, temp, rhoi)

           ! Species densities
           DO is = 1,nb_ns
              prim_var(is) = MAX(rhoi(is),rhoi_tol)
           ENDDO

           ! Velocity components
           DO is = 1,nb_dim
              prim_var(nb_ns + is) = u_vec(is)
           ENDDO

           ! Temperatures
           DO is = 1,nb_temp
              prim_var(nb_ns + nb_dim + is) = temp
           ENDDO

           DO is = 1,nb_eq 
              p(i1 + is) = prim_var(is)
           ENDDO

           ! Index update
           i1 = i1 + nb_eq

       ENDDO

    ELSE

       ! Restarting normally from a previous solution (calorically perfect gas and nonequilibrium flows)
       DO i = 1,nb_cells
       
          READ(in1,*)p(i1 + 1:i1 + nb_eq)

          ! Index update
          i1 = i1 + nb_eq

       ENDDO
  
    ENDIF

    CLOSE(in1)

    ! Initialization of ghost states
    DO is = 1,nb_eq
       prim_var(is) = p(2*nb_eq + is)
    ENDDO

    ! Left boundary ghost states
    p(1:nb_eq) = prim_var
    p(nb_eq + 1:2*nb_eq) = prim_var

    ! Right boundary ghost states
    p((nb_cells + 2)*nb_eq + 1:(nb_cells + 3)*nb_eq) = prim_var
    p((nb_cells + 3)*nb_eq + 1:(nb_cells + 4)*nb_eq) = prim_var

5 FORMAT(A)

  END SUBROUTINE read_solution  
!------------------------------------------------------------------------------!
