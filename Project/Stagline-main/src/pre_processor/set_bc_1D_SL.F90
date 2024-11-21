!------------------------------------------------------------------------------!
!> This subroutine associates pointers for boundary conditions (1D stagnation line calorically perfect gas flows).
!! The free-strem pressure is also selected since it is used for the geometrical inviscid source terms 
  SUBROUTINE set_bc_1D_SL ()

    USE mod_general_data,           ONLY: nb_eq, pos_u, p_inf, V_inf
    USE mod_numerics_data,          ONLY: poly_rec, flag_full_Impl, bc_Jac
    USE mod_domain_boundary,        ONLY: bc, boundary_data, get_boundary_pin, get_boundary_inlet
    USE mod_function_pointer

    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8), DIMENSION(nb_eq) :: inlet_data

    ! Explicit boundary conditions 
    IF (flag_full_Impl.EQV..FALSE.) THEN 

       ! First order accuracy in space
       IF (poly_rec.EQ.'constant') THEN

          ALLOCATE(get_ghost_state_Expl_1D_SL_1st(2))

          ! Loop over domain boundaries
          DO i = 1,2

             SELECT CASE (bc(i))
              
               CASE ('slip_wall')
                 get_ghost_state_Expl_1D_SL_1st(i)%bc_procedure => slip_wall_1D_SL_Expl

               CASE ('sup_in')
                 p_inf = get_boundary_pin (boundary_data(i))
                 inlet_data = get_boundary_inlet(nb_eq, boundary_data(i))
                 V_inf = ABS(inlet_data(pos_u))
                 get_ghost_state_Expl_1D_SL_1st(i)%bc_procedure => sup_in_1D_SL_Expl 

               CASE('no_slip_iso_Twall')
                 get_ghost_state_Expl_1D_SL_1st(i)%bc_procedure => no_slip_iso_Twall_1D_SL_Expl

               CASE DEFAULT 
                 WRITE(*,10)'In "set_bc_1D_SL.F90", bc not implemented yet...'
                 WRITE(*,10)bc(i)
                 PRINT*
                 STOP 

             END SELECT 

          ENDDO

       ! Second order accuracy in space
       ELSE 

          ALLOCATE(get_ghost_state_Expl_1D_SL_2nd(2))

          ! Loop over domain boundaries
          DO i = 1,2

             SELECT CASE (bc(i))

               CASE ('slip_wall')
                 get_ghost_state_Expl_1D_SL_2nd(i)%bc_procedure => slip_wall_1D_SL_Expl_2nd  

               CASE ('sup_in')
                 p_inf = get_boundary_pin (boundary_data(i))
                 inlet_data = get_boundary_inlet(nb_eq, boundary_data(i))
                 V_inf = ABS(inlet_data(pos_u))
                 get_ghost_state_Expl_1D_SL_2nd(i)%bc_procedure => sup_in_1D_SL_Expl_2nd

               CASE('no_slip_iso_Twall')
                 get_ghost_state_Expl_1D_SL_2nd(i)%bc_procedure => no_slip_iso_Twall_1D_SL_Expl_2nd 

               CASE DEFAULT 
                 WRITE(*,10)'In "set_bc_1D_SL.F90", bc not implemented yet...'
                 WRITE(*,10)bc(i)
                 PRINT*
                 STOP 

             END SELECT 

          ENDDO

       ENDIF

    ! Implicit boundary conditions
    ELSE 

        ! First order accuracy in space
       IF (poly_rec.EQ.'constant') THEN

          ALLOCATE(get_ghost_state_Expl_1D_SL_1st(2))
          ALLOCATE(get_ghost_state_Impl_1D_SL_1st(2))

          ! Loop over domain boundaries
          DO i = 1,2
             
             SELECT CASE (bc(i))

               CASE ('slip_wall')
                 IF (bc_Jac(i).EQ.'numerical') THEN
                    get_ghost_state_Expl_1D_SL_1st(i)%bc_procedure => slip_wall_1D_SL_Expl 
                 ELSEIF (bc_Jac(i).EQ.'analytical') THEN
                    get_ghost_state_Impl_1D_SL_1st(i)%bc_procedure => slip_wall_1D_SL_Impl
                 ELSE 
                    WRITE(*,10)'In " set_bc_1D_SL.F90", error in setting the option for the boundary condition Jacobian...'
                    PRINT*
                    STOP
                 ENDIF

               CASE ('sup_in')
                 p_inf = get_boundary_pin (boundary_data(i))
                 inlet_data = get_boundary_inlet(nb_eq, boundary_data(i))
                 V_inf = ABS(inlet_data(pos_u)) 
                 IF (bc_Jac(i).EQ.'numerical') THEN
                    get_ghost_state_Expl_1D_SL_1st(i)%bc_procedure => sup_in_1D_SL_Expl
                 ELSEIF (bc_Jac(i).EQ.'analytical') THEN
                    get_ghost_state_Impl_1D_SL_1st(i)%bc_procedure => sup_in_1D_SL_Impl
                 ELSE 
                    WRITE(*,10)'In " set_bc_1D_SL.F90", error in setting the option for the boundary condition Jacobian...'
                    PRINT*
                    STOP
                 ENDIF

               CASE('no_slip_iso_Twall')
                 IF (bc_Jac(i).EQ.'numerical') THEN
                    get_ghost_state_Expl_1D_SL_1st(i)%bc_procedure => no_slip_iso_Twall_1D_SL_Expl
                 ELSEIF (bc_Jac(i).EQ.'analytical') THEN
                    get_ghost_state_Impl_1D_SL_1st(i)%bc_procedure => no_slip_iso_Twall_1D_SL_Impl
                 ELSE 
                    WRITE(*,10)'In " set_bc_1D_SL.F90", error in setting the option for the boundary condition Jacobian...'
                    PRINT*
                    STOP
                 ENDIF

               CASE DEFAULT 
                 WRITE(*,10)'In "set_bc_1D_SL.F90", bc not implemented yet...'
                 WRITE(*,10)bc(i)
                 PRINT*
                 STOP 

             END SELECT 

          ENDDO

       ! Second order accuracy in space
       ELSE 

          ALLOCATE(get_ghost_state_Expl_1D_SL_2nd(2))
          ALLOCATE(get_ghost_state_Impl_1D_SL_2nd(2))
       
          ! Loop over domain boundaries
          DO i = 1,2
       
             SELECT CASE (bc(i))

               CASE ('slip_wall')
                 IF (bc_Jac(i).EQ.'numerical') THEN
                    get_ghost_state_Expl_1D_SL_2nd(i)%bc_procedure => slip_wall_1D_SL_Expl_2nd
                 ELSEIF (bc_Jac(i).EQ.'analytical') THEN
                    get_ghost_state_Impl_1D_SL_2nd(i)%bc_procedure => slip_wall_1D_SL_Impl_2nd
                 ELSE 
                    WRITE(*,10)'In " set_bc_1D_SL.F90", error in setting the option for the boundary condition Jacobian...'
                    PRINT*
                    STOP
                 ENDIF

               CASE ('sup_in')
                 p_inf = get_boundary_pin (boundary_data(i))
                 inlet_data = get_boundary_inlet(nb_eq, boundary_data(i))
                 V_inf = ABS(inlet_data(pos_u))
                 IF (bc_Jac(i).EQ.'numerical') THEN
                    get_ghost_state_Expl_1D_SL_2nd(i)%bc_procedure => sup_in_1D_SL_Expl_2nd
                 ELSEIF (bc_Jac(i).EQ.'analytical') THEN
                    get_ghost_state_Impl_1D_SL_2nd(i)%bc_procedure => sup_in_1D_SL_Impl_2nd
                 ELSE 
                    WRITE(*,10)'In " set_bc_1D_SL.F90", error in setting the option for the boundary condition Jacobian...'
                    PRINT*
                    STOP
                 ENDIF 

               CASE('no_slip_iso_Twall')
                 IF (bc_Jac(i).EQ.'numerical') THEN
                    get_ghost_state_Expl_1D_SL_2nd(i)%bc_procedure => no_slip_iso_Twall_1D_SL_Expl_2nd 
                 ELSEIF (bc_Jac(i).EQ.'analytical') THEN
                    get_ghost_state_Impl_1D_SL_2nd(i)%bc_procedure => no_slip_iso_Twall_1D_SL_Impl_2nd
                 ELSE 
                    WRITE(*,10)'In "set_bc_1D_SL.F90", error in setting the option for the boundary condition Jacobian...'
                    PRINT*
                    STOP
                 ENDIF

               CASE DEFAULT 
                 WRITE(*,10)'In "set_bc_1D_SL.F90", bc not implemented yet...'
                 WRITE(*,10)bc(i)
                 PRINT*
                 STOP 

             END SELECT 

          ENDDO

       ENDIF

   ENDIF

10 FORMAT(A)

  END SUBROUTINE set_bc_1D_SL
!------------------------------------------------------------------------------!
