!------------------------------------------------------------------------------!
!> This subroutine associates pointers for boundary conditions (1D stagnation line nonequilibrium flows).
!! The free-strem pressure is also computed since it is used for the geometrical inviscid source terms 
  SUBROUTINE set_bc_neq_1D_SL ()

#include "../config.h"

    USE mod_general_data,            ONLY: nb_ns, nb_temp, p_inf, Ri
    USE mod_numerics_data,           ONLY: poly_rec, flag_full_Impl, bc_Jac
    USE mod_domain_boundary,         ONLY: bc, name_bc,boundary_data, get_boundary_rhoi, get_boundary_Tvec, &
                                         & get_boundary_pin
    USE mod_function_pointer
    USE mod_neq_function_pointer,    ONLY: library_get_pressure

#ifdef CARBONABLA
    USE mod_ablation_pointers,       ONLY: procedure_ablation_initialize
#endif
    
    USE mod_general_ablation,        ONLY: Flag_ablation

    IMPLICIT NONE
    
    EXTERNAL set_ablation_pointers

    INTEGER :: i, j
    REAL(KIND=8) :: T, Te
    REAL(KIND=8), DIMENSION(nb_ns) :: rhoi
    REAL(KIND=8), DIMENSION(nb_temp) :: temp 
    
    ! Explicit boundary conditions 
    IF (flag_full_Impl.EQV..FALSE.) THEN 

       ! First order accuracy in space
       IF (poly_rec.EQ.'constant') THEN

          ALLOCATE(get_ghost_state_Expl_1D_SL_1st(2))

          ! Loop over domain boundaries
          DO i = 1,2

             SELECT CASE (bc(i))
              
               CASE ('slip_wall')
                 get_ghost_state_Expl_1D_SL_1st(i)%bc_procedure => slip_wall_neq_1D_SL_Expl

               CASE ('sup_in')
                 get_ghost_state_Expl_1D_SL_1st(i)%bc_procedure => sup_in_neq_1D_SL_Expl 

                 ! Compute the free-stream pressure 
                 rhoi = get_boundary_rhoi (nb_ns, boundary_data(i))
                 temp = get_boundary_Tvec (nb_temp, boundary_data(i))
                
                 CALL library_get_pressure(rhoi, temp, p_inf)

               CASE ('sub_in_vel_equi_comp')
                 get_ghost_state_Expl_1D_SL_1st(i)%bc_procedure => sub_in_vel_neq_1D_SL_Expl
                 ! Load the free-stream pressure from the input
                 p_inf = get_boundary_pin (boundary_data(i))

               CASE('no_slip_iso_Twall')
                 get_ghost_state_Expl_1D_SL_1st(i)%bc_procedure => no_slip_iso_Twall_neq_1D_SL_Expl
               
               CASE('no_slip_adiabatic')
                 get_ghost_state_Expl_1D_SL_1st(i)%bc_procedure => no_slip_adiabatic_neq_1D_SL_Expl

#ifdef CARBONABLA
               CASE('no_slip_SEB_ablation')
                 CALL set_ablation_pointers () 
                 CALL procedure_ablation_initialize (nb_ns)
                 get_ghost_state_Expl_1D_SL_1st(i)%bc_procedure => no_slip_SEB_ablation_neq_1D_SL_Expl   

               CASE('no_slip_isothermal_ablation_Twall')
                 CALL set_ablation_pointers () 
                 CALL procedure_ablation_initialize (nb_ns)
                 get_ghost_state_Expl_1D_SL_1st(i)%bc_procedure => no_slip_isothermal_ablation_neq_1D_SL_Expl   
#endif

              CASE('no_slip_iso_General_Ablation_Twall')
                 get_ghost_state_Expl_1D_SL_1st(i)%bc_procedure => no_slip_isothermal_general_ablation_neq_1D_SL_Expl
              

               CASE('no_slip_seb_General_Ablation')
                 get_ghost_state_Expl_1D_SL_1st(i)%bc_procedure => no_slip_SEB_general_ablation_neq_1D_SL_Expl

 
               CASE DEFAULT 
                 WRITE(*,10)'In "set_neq_bc_1D_SL.F90", bc not implemented yet...'
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
                 get_ghost_state_Expl_1D_SL_2nd(i)%bc_procedure => slip_wall_neq_1D_SL_Expl_2nd  

               CASE ('sup_in')
                 get_ghost_state_Expl_1D_SL_2nd(i)%bc_procedure => sup_in_neq_1D_SL_Expl_2nd

                 ! Compute the free-stream pressure 
                 rhoi = get_boundary_rhoi (nb_ns, boundary_data(i))
                 temp = get_boundary_Tvec (nb_temp, boundary_data(i))
                
                 CALL library_get_pressure(rhoi, temp, p_inf)

               CASE ('sub_in_vel_equi_comp')
                 get_ghost_state_Expl_1D_SL_2nd(i)%bc_procedure => sub_in_vel_neq_1D_SL_Expl_2nd
                 ! Load the free-stream pressure from the input
                 p_inf = get_boundary_pin (boundary_data(i))

               CASE('no_slip_iso_Twall')
                 get_ghost_state_Expl_1D_SL_2nd(i)%bc_procedure => no_slip_iso_Twall_neq_1D_SL_Expl_2nd 

#ifdef CARBONABLA
               CASE('no_slip_SEB_ablation')
                 CALL set_ablation_pointers () 
                 CALL procedure_ablation_initialize (nb_ns)
                 get_ghost_state_Expl_1D_SL_2nd(i)%bc_procedure => no_slip_SEB_ablation_neq_1D_SL_Expl_2nd   

               CASE('no_slip_isothermal_ablation_Twall')
                 CALL set_ablation_pointers () 
                 CALL procedure_ablation_initialize (nb_ns)
                 get_ghost_state_Expl_1D_SL_2nd(i)%bc_procedure => no_slip_isothermal_ablation_neq_1D_SL_Expl_2nd   
#endif
               CASE('no_slip_iso_General_Ablation_Twall')
                 get_ghost_state_Expl_1D_SL_2nd(i)%bc_procedure => no_slip_isothermal_general_ablation_neq_1D_SL_Expl_2nd


               CASE DEFAULT 
                 WRITE(*,10)'In "set_bc_neq_1D_SL.F90", bc not implemented yet...'
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
                    get_ghost_state_Expl_1D_SL_1st(i)%bc_procedure => slip_wall_neq_1D_SL_Expl 
                 ELSEIF (bc_Jac(i).EQ.'analytical') THEN 
                    WRITE(*,10)'In "set_neq_bc_neq_1D_SL.F90", analytical Jacobian for bc'
                    WRITE(*,10)bc(i)
                    WRITE(*,10)'not implemented yet...'
                    PRINT*
                    STOP
                 ELSE 
                    WRITE(*,10)'In "set_bc_neq_1D_SL.F90", error in setting the option for the boundary condition Jacobian...'
                    PRINT*
                    STOP
                 ENDIF

               CASE ('sup_in')
                 IF (bc_Jac(i).EQ.'numerical') THEN
                    get_ghost_state_Expl_1D_SL_1st(i)%bc_procedure => sup_in_neq_1D_SL_Expl
                 ELSEIF (bc_Jac(i).EQ.'analytical') THEN 
                    get_ghost_state_Impl_1D_SL_1st(i)%bc_procedure => sup_in_neq_1D_SL_Impl
                 ELSE 
                    WRITE(*,10)'In "set_bc_neq_1D_SL.F90", error in setting the option for the boundary condition Jacobian...'
                    PRINT*
                    STOP
                 ENDIF

                 ! Compute the free-stream pressure 
                 rhoi = get_boundary_rhoi (nb_ns, boundary_data(i))
                 temp = get_boundary_Tvec (nb_temp, boundary_data(i))
                
                 CALL library_get_pressure(rhoi, temp, p_inf)

               CASE ('sub_in_vel_equi_comp')
                 IF (bc_Jac(i).EQ.'numerical') THEN
                    get_ghost_state_Expl_1D_SL_1st(i)%bc_procedure => sub_in_vel_neq_1D_SL_Expl
                 ELSE 
                    WRITE(*,10)'In "set_bc_neq_1D_SL.F90", error in setting the option for the boundary condition Jacobian...'
                    PRINT*
                    STOP
                 ENDIF
                 ! Load the free-stream pressure from the input
                 p_inf = get_boundary_pin (boundary_data(i))
               
               CASE('no_slip_iso_Twall')
                 IF (bc_Jac(i).EQ.'numerical') THEN 
                    get_ghost_state_Expl_1D_SL_1st(i)%bc_procedure => no_slip_iso_Twall_neq_1D_SL_Expl
                 ELSEIF (bc_Jac(i).EQ.'analytical') THEN
                    WRITE(*,10)'In "set_neq_bc_neq_1D_SL.F90", analytical Jacobian for bc'
                    WRITE(*,10)bc(i)
                    WRITE(*,10)'not implemented yet...'
                    PRINT*
                    STOP 
                 ELSE 
                    WRITE(*,10)'In "set_bc_neq_1D_SL.F90", error in setting the option for the boundary condition Jacobian...'
                    PRINT*
                    STOP
                 ENDIF
               
               CASE('no_slip_adiabatic')
                 IF (bc_Jac(i).EQ.'numerical') THEN 
                    get_ghost_state_Expl_1D_SL_1st(i)%bc_procedure => no_slip_adiabatic_neq_1D_SL_Expl
                 ELSEIF (bc_Jac(i).EQ.'analytical') THEN
                    WRITE(*,10)'In "set_neq_bc_neq_1D_SL.F90", analytical Jacobian for bc'
                    WRITE(*,10)bc(i)
                    WRITE(*,10)'not implemented yet...'
                    PRINT*
                    STOP 
                 ELSE 
                    WRITE(*,10)'In "set_bc_neq_1D_SL.F90", error in setting the option for the boundary condition Jacobian...'
                    PRINT*
                    STOP
                 ENDIF

#ifdef CARBONABLA
               CASE('no_slip_isothermal_ablation_Twall')
                 CALL set_ablation_pointers () 
                 CALL procedure_ablation_initialize (nb_ns)
                 IF (bc_Jac(i).EQ.'numerical') THEN 
                    get_ghost_state_Expl_1D_SL_1st(i)%bc_procedure => no_slip_isothermal_ablation_neq_1D_SL_Expl   
                 ELSEIF (bc_Jac(i).EQ.'analytical') THEN
                    WRITE(*,10)'In "set_neq_bc_neq_1D_SL.F90", analytical Jacobian for bc'
                    WRITE(*,10)bc(i)
                    WRITE(*,10)'not implemented yet...'
                    PRINT*
                    STOP
                 ELSE 
                    WRITE(*,10)'In "set_bc_neq_1D_SL.F90", error in setting the option for the boundary condition Jacobian...'
                    PRINT*
                    STOP
                 ENDIF

               CASE('no_slip_SEB_ablation')
                 CALL set_ablation_pointers () 
                 CALL procedure_ablation_initialize (nb_ns)
                 IF (bc_Jac(i).EQ.'numerical') THEN 
                 get_ghost_state_Expl_1D_SL_1st(i)%bc_procedure => no_slip_SEB_ablation_neq_1D_SL_Expl   
                 ELSEIF (bc_Jac(i).EQ.'analytical') THEN
                    WRITE(*,10)'In "set_neq_bc_neq_1D_SL.F90", analytical Jacobian for bc'
                    WRITE(*,10)bc(i)
                    WRITE(*,10)'not implemented yet...'
                    PRINT*
                    STOP
                 ELSE 
                    WRITE(*,10)'In "set_bc_neq_1D_SL.F90", error in setting the option for the boundary condition Jacobian...'
                    PRINT*
                    STOP
                 ENDIF
#endif

                  CASE('no_slip_iso_General_Ablation_Twall')
                  IF (bc_Jac(i).EQ.'numerical') THEN 
                    get_ghost_state_Expl_1D_SL_1st(i)%bc_procedure => no_slip_isothermal_general_ablation_neq_1D_SL_Expl
                 ELSEIF (bc_Jac(i).EQ.'analytical') THEN
                    WRITE(*,10)'In "set_neq_bc_neq_1D_SL.F90", analytical Jacobian for bc'
                    WRITE(*,10)bc(i)
                    WRITE(*,10)'not implemented yet...'
                    PRINT*
                    STOP 
                 ELSE 
                    WRITE(*,10)'In "set_bc_neq_1D_SL.F90", error in setting the option for the boundary condition Jacobian...'
                    PRINT*
                    STOP
                 ENDIF

                CASE('no_slip_seb_General_Ablation')
                  IF (bc_Jac(i).EQ.'numerical') THEN 
                    get_ghost_state_Expl_1D_SL_1st(i)%bc_procedure => no_slip_SEB_general_ablation_neq_1D_SL_Expl
                 ELSEIF (bc_Jac(i).EQ.'analytical') THEN
                    WRITE(*,10)'In "set_neq_bc_neq_1D_SL.F90", analytical Jacobian for bc'
                    WRITE(*,10)bc(i)
                    WRITE(*,10)'not implemented yet...'
                    PRINT*
                    STOP 
                 ELSE 
                    WRITE(*,10)'In "set_bc_neq_1D_SL.F90", error in setting the option for the boundary condition Jacobian...'
                    PRINT*
                    STOP
                 ENDIF
                 


               CASE DEFAULT 
                 WRITE(*,10)'In "set_neq_bc_neq_1D_SL.F90", bc not implemented yet...'
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
                    get_ghost_state_Expl_1D_SL_2nd(i)%bc_procedure => slip_wall_neq_1D_SL_Expl_2nd 
                 ELSEIF (bc_Jac(i).EQ.'analytical') THEN 
                    WRITE(*,10)'In "set_neq_bc_neq_1D_SL.F90", analytical Jacobian for bc'
                    WRITE(*,10)bc(i)
                    WRITE(*,10)'not implemented yet...'
                    PRINT*
                    STOP
                 ELSE 
                    WRITE(*,10)'In "set_bc_neq_1D_SL.F90", error in setting the option for the boundary condition Jacobian...'
                    PRINT*
                    STOP
                 ENDIF

               CASE ('sup_in')
                 IF (bc_Jac(i).EQ.'numerical') THEN
                    get_ghost_state_Expl_1D_SL_2nd(i)%bc_procedure => sup_in_neq_1D_SL_Expl_2nd
                 ELSEIF (bc_Jac(i).EQ.'analytical') THEN
                    get_ghost_state_Impl_1D_SL_2nd(i)%bc_procedure => sup_in_neq_1D_SL_Impl_2nd  
                 ELSE 
                    WRITE(*,10)'In "set_bc_neq_1D_SL.F90", error in setting the option for the boundary condition Jacobian...'
                    PRINT*
                    STOP
                 ENDIF

                 ! Compute the free-stream pressure 
                 rhoi = get_boundary_rhoi (nb_ns, boundary_data(i))
                 temp = get_boundary_Tvec (nb_temp, boundary_data(i))
                
                 CALL library_get_pressure(rhoi, temp, p_inf)
               
               CASE ('sub_in_vel_equi_comp')
                 IF (bc_Jac(i).EQ.'numerical') THEN
                    get_ghost_state_Expl_1D_SL_2nd(i)%bc_procedure => sub_in_vel_neq_1D_SL_Expl_2nd
                 ELSE 
                    WRITE(*,10)'In "set_bc_neq_1D_SL.F90", error in setting the option for the boundary condition Jacobian...'
                    PRINT*
                    STOP
                 ENDIF
                 ! Load the free-stream pressure from the input
                 p_inf = get_boundary_pin (boundary_data(i))
                
               CASE('no_slip_iso_Twall')
                 IF (bc_Jac(i).EQ.'numerical') THEN 
                    get_ghost_state_Expl_1D_SL_2nd(i)%bc_procedure => no_slip_iso_Twall_neq_1D_SL_Expl_2nd
                 ELSEIF (bc_Jac(i).EQ.'analytical') THEN
                    WRITE(*,10)'In "set_neq_bc_neq_1D_SL.F90", analytical Jacobian for bc'
                    WRITE(*,10)bc(i)
                    WRITE(*,10)'not implemented yet...'
                    PRINT*
                    STOP
                 ELSE 
                    WRITE(*,10)'In "set_bc_neq_1D_SL.F90", error in setting the option for the boundary condition Jacobian...'
                    PRINT*
                    STOP
                 ENDIF

                 CASE('no_slip_adiabatic')
                 IF (bc_Jac(i).EQ.'numerical') THEN 
                    get_ghost_state_Expl_1D_SL_2nd(i)%bc_procedure => no_slip_adiabatic_neq_1D_SL_Expl_2nd
                 ELSEIF (bc_Jac(i).EQ.'analytical') THEN
                    WRITE(*,10)'In "set_neq_bc_neq_1D_SL.F90", analytical Jacobian for bc'
                    WRITE(*,10)bc(i)
                    WRITE(*,10)'not implemented yet...'
                    PRINT*
                    STOP 
                 ELSE 
                    WRITE(*,10)'In "set_bc_neq_1D_SL.F90", error in setting the option for the boundary condition Jacobian...'
                    PRINT*
                    STOP
                 ENDIF


#ifdef CARBONABLA
               CASE('no_slip_isothermal_ablation_Twall')
                 CALL set_ablation_pointers () 
                 CALL procedure_ablation_initialize (nb_ns)
                 IF (bc_Jac(i).EQ.'numerical') THEN 
                    get_ghost_state_Expl_1D_SL_2nd(i)%bc_procedure => no_slip_isothermal_ablation_neq_1D_SL_Expl_2nd   
                 ELSEIF (bc_Jac(i).EQ.'analytical') THEN
                    WRITE(*,10)'In "set_neq_bc_neq_1D_SL.F90", analytical Jacobian for bc'
                    WRITE(*,10)bc(i)
                    WRITE(*,10)'not implemented yet...'
                    PRINT*
                    STOP
                 ELSE 
                    WRITE(*,10)'In "set_bc_neq_1D_SL.F90", error in setting the option for the boundary condition Jacobian...'
                    PRINT*
                    STOP
                 ENDIF

               CASE('no_slip_SEB_ablation')
                 CALL set_ablation_pointers () 
                 CALL procedure_ablation_initialize (nb_ns)
                 IF (bc_Jac(i).EQ.'numerical') THEN 
                    get_ghost_state_Expl_1D_SL_2nd(i)%bc_procedure => no_slip_SEB_ablation_neq_1D_SL_Expl_2nd
                 ELSEIF (bc_Jac(i).EQ.'analytical') THEN
                    WRITE(*,10)'In "set_neq_bc_neq_1D_SL.F90", analytical Jacobian for bc'
                    WRITE(*,10)bc(i)
                    WRITE(*,10)'not implemented yet...'
                    PRINT*
                    STOP
                 ELSE 
                    WRITE(*,10)'In "set_bc_neq_1D_SL.F90", error in setting the option for the boundary condition Jacobian...'
                    PRINT*
                    STOP
                 ENDIF
               
                CASE('no_slip_iso_General_Ablation_Twall')
                  IF (bc_Jac(i).EQ.'numerical') THEN 
                    get_ghost_state_Expl_1D_SL_2nd(i)%bc_procedure => no_slip_isothermal_general_ablation_neq_1D_SL_Expl_2nd
                 ELSEIF (bc_Jac(i).EQ.'analytical') THEN
                    WRITE(*,10)'In "set_neq_bc_neq_1D_SL.F90", analytical Jacobian for bc'
                    WRITE(*,10)bc(i)
                    WRITE(*,10)'not implemented yet...'
                    PRINT*
                    STOP 
                 ELSE 
                    WRITE(*,10)'In "set_bc_neq_1D_SL.F90", error in setting the option for the boundary condition Jacobian...'
                    PRINT*
                    STOP
                 ENDIF


#endif

               CASE DEFAULT 
                 WRITE(*,10)'In "set_neq_bc_neq_1D_SL.F90", bc not implemented yet...'
                 WRITE(*,10)bc(i)
                 PRINT*
                 STOP 

             END SELECT 

          ENDDO

       ENDIF

   ENDIF

10 FORMAT(A)

  END SUBROUTINE set_bc_neq_1D_SL
!------------------------------------------------------------------------------!
