#include"../config.h"

!------------------------------------------------------------------------------!
!> This subroutine applies an general ablative boundary condition coupled with gsi from mutation++ 
!!assuming an isothermal wall with the temperature given in input for 1D stagnation line 
!! nonequilibrium flows. 

        SUBROUTINE no_slip_isothermal_general_ablation_neq_1D_SL_Expl (id, phys_data1, u_phys1, ghost_data1, u_ghost1)

        USE mod_general_data,          ONLY: nb_ns, nb_temp, nb_eq, pos_u, pos_v, pos_T,  & 
                                           & pos_v_cell, pos_u_cell, pos_T_cell, pos_pres_cell, Ri, &
                                           & pos_em, pos_Te, nb_te, volumes
        USE mod_domain_boundary,       ONLY: boundary, boundary_data, get_boundary_Twall
        USE mod_function_pointer,      ONLY: get_cons_phys_from_prim
        USE mod_neq_function_pointer,  ONLY: library_get_molar_fractions, library_comp_tol, &
                                           & library_set_diff_model, library_set_wall_state, &
                                           & library_solve_surface_balance, library_get_wall_state, &
                                           & library_get_mass_blowing_rate, library_get_surface_production_rates, &
                                           & library_set_cond_heat_flux

        USE mod_general_ablation
     
        IMPLICIT NONE

        INTEGER :: i
        REAL(KIND=8) :: p_ov_T, sum_xi_g
        REAL(KIND=8) :: rho_wall, p_wall
        REAL(KIND=8) :: vol_p, dx
        REAL(KIND=8) :: T1
        REAL(KIND=8),  DIMENSION(nb_temp) :: Twall
        REAL(KIND=8), DIMENSION(nb_ns)   :: xi, sum_xi
        REAL(KIND=8), DIMENSION(nb_ns)   :: rhoi_wall, xi_w
        REAL(KIND=8), DIMENSION(nb_ns)   :: xi_g

        REAL(KIND=8), DIMENSION(nb_eq)   :: prim
        TYPE(boundary) :: bound


        INTEGER, INTENT(IN) :: id                               !< boundary id
        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: u_phys1      !< conservative variables of the first physical cell
        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: phys_data1   !< physical properties of the first physical cell
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost1     !< conservative variables of the first ghost cell
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data1  !< physical properties of the first ghost cell

        ! Wall temperature
        bound = boundary_data(id)

         DO i = 1,nb_temp
                Twall(1 + i - 1) = get_boundary_Twall(bound)
         ENDDO
        

        T1 = phys_data1(pos_T_cell)
        ! Getting the temperature of the two physical state next to the wall

        ! Some definitions
        vol_p = volumes(3)
        dx = (vol_p)/2.0   

        ! Species molar fractions of physical state 
        CALL library_get_molar_fractions(u_phys1(1:nb_ns), xi)
     
       ! DO i = 1, 6
            
      !  xi(i) = 1.D-16
       ! ENDDO
 
        CALL library_comp_tol(xi) 
  
                

        ! Wall species densities computed based on the molar fractions and pressureo
        ! of physical state and first guess wall temperature
        p_ov_T = phys_data1(pos_pres_cell)/Twall(1)
        p_wall = phys_data1(pos_pres_cell)
        DO i = 1,nb_ns 
        ! Constant extrapolation of mole fractions
                rhoi_wall(i) = p_ov_T*xi(i)/Ri(i)
        ENDDO
        
        ! Pass mole fractions of first physical cell and distance from the wall (half volume)
        CALL library_set_diff_model( xi, dx )


        ! Set the initial wall state using the first physical cell desities and the imposed wall temperature
        CALL library_set_wall_state( rhoi_wall, Twall)
        
        ! Solve the surface mass balance
        CALL library_solve_surface_balance()


        ! Get the species densities at the wall obtained from the resolution of the surface mass balance
        CALL library_get_wall_state(rhoi_wall, Twall)

        ! Get the surface mass blowing flux (kg/m2-s)
        CALL library_get_mass_blowing_rate(mdot)

        CALL library_get_surface_production_rates(mdot_i)

        mdot_i = -mdot_i

        ! Compute the velocity normal to the wall as the ratio between mdot and the wall gas density
        ! Careful with the sign! Gas should enter the domain. 
        ! It might be necessary to add a minus depending on the CFD code convention 
        mdot = -mdot

        IF (mdot > 35.0D0) THEN
        mdot = 35.0D0
        ENDIF

        uwall=mdot/sum(rhoi_wall)
       
        IF (uwall < 0.0D0) THEN
        uwall = 0.0D0
        ENDIF

       ! Get mole fractions of the wall interface
       CALL library_get_molar_fractions(rhoi_wall(1:nb_ns), xi_w)

       xi_w = abs(xi_w)
       
       CALL library_comp_tol(xi_w) 
 
       sum_xi_g = 0.d0
       DO i = 1,nb_ns 
          xi_g(i) = xi_w(i) - (xi(i)-xi_w(i))
          sum_xi_g = sum_xi_g + xi_g(i)
       ENDDO
       xi_g = abs(xi_g)
       CALL library_comp_tol(xi_g) 

   
       ! Normalization of mole fractions 

      ! xi_g = xi_w
       ! Linear extrapolation of temperature
       prim(pos_T) = MAX(Twall(1) - (T1-Twall(1)), 0.5d0 * Twall(1)) 
       !prim(pos_T) = Twall(1) 
       ! Temperature(s)
       DO i = 1,nb_temp
          prim(pos_T + i - 1) = prim(pos_T)
       ENDDO

       ! Ghost pressure set equal to physical one
       IF (nb_te.GE.1) THEN
          p_ov_T = p_wall/(prim(pos_T)+xi(pos_em)*(prim(pos_Te)-prim(pos_T)))
       ELSE
          p_ov_T = p_wall / prim(pos_T)
       ENDIF

       ! Computation of ghost species densities
       DO i = 1,nb_ns 
          prim(i) = p_ov_T*xi_g(i)/Ri(i)
       ENDDO

        ! Extrapolation of ghost velocities
        prim(pos_u) = uwall - (phys_data1(pos_u_cell)-uwall)
        !prim(pos_u) = uwall
        prim(pos_v) = -phys_data1(pos_v_cell)

        Twall_translation = Twall(1)
        Pwall_ablation= p_wall

        ! Compute the ghost state conservative variables and physical properties from primitive variables 
        CALL get_cons_phys_from_prim(prim, u_ghost1, ghost_data1)

        END SUBROUTINE no_slip_isothermal_general_ablation_neq_1D_SL_Expl
!------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!> This subroutine applies an general ablative boundary condition coupled with gsi from mutation++ 
!! computing the wall temperature using a surface energy balance 1D stagnation line 
!! nonequilibrium flows. 

        SUBROUTINE no_slip_SEB_general_ablation_neq_1D_SL_Expl (id, phys_data1, u_phys1, ghost_data1, u_ghost1)

        USE mod_general_data,          ONLY: nb_ns, nb_temp, nb_eq, pos_u, pos_v, pos_T,  & 
                                           & pos_v_cell, pos_u_cell, pos_T_cell, pos_pres_cell, Ri, &
                                           & pos_em, pos_Te, nb_te, volumes, start_prop_phys, nb_prop, &
                                           & cell_prop
        USE mod_domain_boundary,       ONLY: boundary, boundary_data, get_boundary_Twall
        USE mod_function_pointer,      ONLY: get_cons_phys_from_prim
        USE mod_neq_function_pointer,  ONLY: library_get_molar_fractions, library_comp_tol, &
                                           & library_set_diff_model, library_set_wall_state, &
                                           & library_solve_surface_balance, library_get_wall_state, &
                                           & library_get_mass_blowing_rate, library_get_surface_production_rates, &
                                           & library_set_cond_heat_flux, library_set_wall_radiation

        USE mod_general_ablation
        USE mod_radiation,              ONLY:  wall_Rad
     
        IMPLICIT NONE

        INTEGER :: i
        REAL(KIND=8) :: p_ov_T, sum_xi_g
        REAL(KIND=8) :: T2
        REAL(KIND=8) :: rho_wall, p_wall
        REAL(KIND=8) :: vol_p, dx
        REAL(KIND=8),  DIMENSION(nb_temp) :: Twall, T1
        REAL(KIND=8), DIMENSION(nb_ns)   :: xi, sum_xi
        REAL(KIND=8), DIMENSION(nb_ns)   :: rhoi_wall, xi_w
        REAL(KIND=8), DIMENSION(nb_ns)   :: xi_g

        REAL(KIND=8), DIMENSION(nb_eq)   :: prim
        TYPE(boundary) :: bound


        INTEGER, INTENT(IN) :: id                               !< boundary id
        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: u_phys1      !< conservative variables of the first physical cell
        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: phys_data1   !< physical properties of the first physical cell
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost1     !< conservative variables of the first ghost cell
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data1  !< physical properties of the first ghost cell


        DO i = 1,nb_temp
           T1(i) = phys_data1(pos_T_cell + i - 1)
        ENDDO

        T2 = cell_prop(start_prop_phys + nb_prop + pos_T_cell)

        DO i = 1,nb_temp
         Twall(1 + i - 1) = T1(1)-(T2-T1(1))*0.5d0
       ENDDO
      
        ! Getting the temperature of the two physical state next to the wall

        ! Some definitions
        vol_p = volumes(3)
        dx = (vol_p)/2.0   

        ! Species molar fractions of physical state 
        CALL library_get_molar_fractions(u_phys1(1:nb_ns), xi)
     
       ! DO i = 1, 6
            
      !  xi(i) = 1.D-16
       ! ENDDO
 
        CALL library_comp_tol(xi) 
  
                

        ! Wall species densities computed based on the molar fractions and pressureo
        ! of physical state and first guess wall temperature
        p_ov_T = phys_data1(pos_pres_cell)/Twall(1)
        p_wall = phys_data1(pos_pres_cell)
        DO i = 1,nb_ns 
        ! Constant extrapolation of mole fractions
                rhoi_wall(i) = p_ov_T*xi(i)/Ri(i)
        ENDDO
        
        ! Pass mole fractions of first physical cell and distance from the wall (half volume)
        CALL library_set_diff_model( xi, dx )

        
        CALL library_set_cond_heat_flux( T1, dx )

        ! Set the initial wall state using the first physical cell desities and the imposed wall temperature
        CALL library_set_wall_state( rhoi_wall, Twall)
         

        CALL library_set_wall_radiation(wall_Rad)
        
        ! Solve the surface mass balance
        CALL library_solve_surface_balance()


        ! Get the species densities at the wall obtained from the resolution of the surface mass balance
        CALL library_get_wall_state(rhoi_wall, Twall)
       
        ! Get the surface mass blowing flux (kg/m2-s)
        CALL library_get_mass_blowing_rate(mdot)

        CALL library_get_surface_production_rates(mdot_i)

        mdot_i = -mdot_i

        ! Compute the velocity normal to the wall as the ratio between mdot and the wall gas density
        ! Careful with the sign! Gas should enter the domain. 
        ! It might be necessary to add a minus depending on the CFD code convention 
        mdot = -mdot
        IF (mdot > 50.0D0) THEN
        mdot = 50.0D0
        ENDIF

        uwall=mdot/sum(rhoi_wall)
       
        IF (uwall < 0.0D0) THEN
        uwall = 0.0D0
        ENDIF

       ! Get mole fractions of the wall interface
       CALL library_get_molar_fractions(rhoi_wall(1:nb_ns), xi_w)

       xi_w = abs(xi_w)
       
       CALL library_comp_tol(xi_w) 
 
       sum_xi_g = 0.d0
       DO i = 1,nb_ns 
          xi_g(i) = xi_w(i) - (xi(i)-xi_w(i))
          sum_xi_g = sum_xi_g + xi_g(i)
       ENDDO
       xi_g = abs(xi_g)
       CALL library_comp_tol(xi_g) 


       ! Normalization of mole fractions 

      ! xi_g = xi_w
       ! Linear extrapolation of temperature
       prim(pos_T) = MAX(Twall(1) - (T1(1)-Twall(1)), 0.5d0 * Twall(1)) 
       !prim(pos_T) = Twall(1) 
       ! Temperature(s)
       DO i = 1,nb_temp
          prim(pos_T + i - 1) = prim(pos_T)
       ENDDO

       ! Ghost pressure set equal to physical one
       IF (nb_te.GE.1) THEN
          p_ov_T = p_wall/(prim(pos_T)+xi(pos_em)*(prim(pos_Te)-prim(pos_T)))
       ELSE
          p_ov_T = p_wall / prim(pos_T)
       ENDIF

       ! Computation of ghost species densities
       DO i = 1,nb_ns 
          prim(i) = p_ov_T*xi_g(i)/Ri(i)
       ENDDO

        ! Extrapolation of ghost velocities
        prim(pos_u) = uwall - (phys_data1(pos_u_cell)-uwall)
        !prim(pos_u) = uwall
        prim(pos_v) = -phys_data1(pos_v_cell)

        Twall_translation = Twall(1)
        Pwall_ablation= p_wall
      

        ! Compute the ghost state conservative variables and physical properties from primitive variables 
        CALL get_cons_phys_from_prim(prim, u_ghost1, ghost_data1)

        END SUBROUTINE no_slip_SEB_general_ablation_neq_1D_SL_Expl
!------------------------------------------------------------------------------------------------------------------------------!
