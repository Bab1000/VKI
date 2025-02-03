#include"../config.h"

#ifdef CARBONABLA
!------------------------------------------------------------------------------!
!> This subroutine applies an ablative boundary condition in case of
!! carbon-based (i.e. graphite, carbon preform, carbon-phenolic
!carbon-carbon) wall by solving the Surface Mass Balance (SMB) and
!assuming an isothermal wall with the temperature given in input for 1D stagnation line 
!! nonequilibrium flows. 

        SUBROUTINE no_slip_isothermal_ablation_neq_1D_SL_Expl (id, phys_data1, u_phys1, ghost_data1, u_ghost1)

        USE mod_general_data,          ONLY: nb_ns, nb_temp, nb_eq, pos_u, pos_v, pos_T,  & 
                                           & pos_v_cell, pos_u_cell, pos_T_cell, pos_pres_cell, Ri, &
                                           & pos_em, pos_Te, nb_te
        USE mod_domain_boundary,       ONLY: boundary, boundary_data, get_boundary_Twall
        USE mod_function_pointer,      ONLY: get_cons_phys_from_prim
        USE mod_neq_function_pointer,  ONLY: library_get_molar_fractions, library_comp_tol
        USE mod_ablation_data,         ONLY: uwall, Twall 
        USE mod_ablation_pointers,     ONLY: procedure_ablation_SMB
        IMPLICIT NONE

        INTEGER :: i
        REAL(KIND=8) :: p_ov_T, sum_xi_g
        REAL(KIND=8) :: T1
        REAL(KIND=8) :: rho_r, rho_l, u_r, v_r, T_r
        REAL(KIND=8) :: rho_wall, p_wall
        REAL(KIND=8) :: frac, frac_inv
        REAL(KIND=8), DIMENSION(nb_ns)   :: xi
        REAL(KIND=8), DIMENSION(nb_ns)   :: rhoi_wall, xi_w
        REAL(KIND=8), DIMENSION(nb_ns)   :: rhoi_r
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
        Twall = get_boundary_Twall(bound)

        ! Getting the temperature of the two physical state next to the wall
        T1 = phys_data1(pos_T_cell)
            

        ! Species molar fractions of physical state 
        CALL library_get_molar_fractions(u_phys1(1:nb_ns), xi)
        CALL library_comp_tol(xi) 


        ! Wall species densities computed based on the molar fractions and pressureo
        ! of physical state and first guess wall temperature
        p_ov_T = phys_data1(pos_pres_cell)/Twall
        p_wall = phys_data1(pos_pres_cell)
        DO i = 1,nb_ns 
        ! Constant extrapolation of mole fractions
                rhoi_wall(i) = p_ov_T*xi(i)/Ri(i)
        ENDDO

        ! Call surface mass balance subroutine
        CALL  procedure_ablation_SMB(u_phys1,p_wall,rhoi_wall)

        ! The definition of the Roe's average state is used here to extrapolate the
        ! ghost cell values. Uncomment below to use it.
        ! NOTE: The extrapolation of the mole fractions seems to work better
        ! (better pressure profile next to the wall)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Right (wall/phisical) and average state densities at the wall interface are computed
        !rho_r = 0.d0
        !rho_wall = 0.d0
        !DO i = 1, nb_ns
        !   rhoi_r(i) = u_phys1(i)
        !   rho_r = rho_r + rhoi_r(i)
        !   rho_wall = rho_wall + rhoi_wall(i)
        !ENDDO
        !
        !! Left state (ghost/wall) of the wall interface.
        !! It is calculated usign the definition of the Roe average density
        !rho_l =  (rho_wall * rho_wall)/ rho_r
        !
        !! Calculation of some useful quantities
        !frac     = 1.d0/(1.d0+sqrt(rho_r/rho_l))
        !frac_inv = 1.d0/frac
        !
        !! Species densities
        !DO i = 1,nb_ns 
        !   ! Roe average state extrapolation
        !   prim(i)=(rhoi_wall(i) / rho_wall - rhoi_r(i)  * frac * sqrt(rho_r/rho_l) / rho_r) * rho_l * (frac_inv )
        !   !Check to avoid negative species densities
        !   if  (prim(i).lt.1.d-30) prim(i) = 1.d-30
        !ENDDO
        !
        ! Roe average state extrapolation
        !u_r      =  phys_data1(pos_u_cell)
        !v_r      =  phys_data1(pos_v_cell)
        !prim(pos_u) = (uwall - frac * sqrt(rho_r/rho_l) * u_r) * frac_inv
        !prim(pos_v) = - v_r
        !
        !! Temperature(s)
        !DO i = 1,nb_temp
        !   ! Roe average state extrapolation
        !   T_r = phys_data1(pos_T_cell + i - 1)
        !   prim(pos_T + i - 1) = MAX((Twall - frac * sqrt(rho_r/rho_l) * T_r) * frac_inv, 0.5d0 * Twall)
        !ENDDO
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! Extrapolation of the first ghost state through mole fractions computed
       ! in the SMB subroutine. This ensure the that the same diffusive flux
       ! will result when the ghost state will be used to compute it at the wall
       ! interface (centered difference using ghost and physical state)

       ! Get mole fractions of the wall interface
       CALL library_get_molar_fractions(rhoi_wall(1:nb_ns), xi_w)
       CALL library_comp_tol(xi_w) 

       ! Linear extrapolation of ghost state mole fraction
       sum_xi_g = 0.d0
       DO i = 1,nb_ns 
          xi_g(i) = xi_w(i) - (xi(i)-xi_w(i))
          sum_xi_g = sum_xi_g + xi_g(i)
       ENDDO

       ! Normalization of mole fractions 
       xi_g = xi_g / sum_xi_g

       ! Linear extrapolation of temperature
       prim(pos_T) = MAX(Twall - (T1-Twall), 0.5d0 * Twall) 

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
        prim(pos_v) = -phys_data1(pos_v_cell)

        ! Compute the ghost state conservative variables and physical properties from primitive variables 
        CALL get_cons_phys_from_prim(prim, u_ghost1, ghost_data1)

        END SUBROUTINE no_slip_isothermal_ablation_neq_1D_SL_Expl
!------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------------!
!> This subroutine applies an ablative boundary condition in case of
!! non-pyrolyzing carbon-based (i.e. graphite, carbon preform, carbon-carbon) wall by solving the Surface Mass Balance (SMB) and
!! the Surface Energy Balance (SEB) for 1D stagnation line 
!! nonequilibrium flows. 

        SUBROUTINE no_slip_SEB_ablation_neq_1D_SL_Expl(id, phys_data1, u_phys1, ghost_data1, u_ghost1)

        USE mod_general_data,          ONLY: nb_ns, nb_temp, nb_eq, pos_u, pos_v, pos_T,  & 
                                           & pos_v_cell, pos_u_cell, pos_T_cell, pos_pres_cell, Ri, &
                                           & pos_em, pos_Te, nb_te, cell_prop, nb_prop, start_prop_phys
        USE mod_function_pointer,      ONLY: get_cons_phys_from_prim
        USE mod_neq_function_pointer,  ONLY: library_get_molar_fractions, library_comp_tol
        USE mod_ablation_data,         ONLY: uwall, Twall, phi_pyro, mdot_wall_tot
        USE mod_ablation_num_param,    ONLY: j_max, eps
        USE mod_ablation_pointers,     ONLY: procedure_ablation_SMB, procedure_ablation_compute_wall_SEB, &
                                             & library_ablation_compute_wall_mass_blowing_rate, &
                                             & procedure_ablation_compute_wall_source_terms
        IMPLICIT NONE

        INTEGER :: i, j
        REAL(KIND=8) :: p_ov_T
        REAL(KIND=8) :: T2
        REAL(KIND=8) :: p_wall
        REAL(KIND=8) :: rhowall 
        REAL(KIND=8) :: seb
        REAL(KIND=8) :: sum_xi_g
        REAL(KIND=8), DIMENSION(nb_ns)   :: xi
        REAL(KIND=8), DIMENSION(nb_ns)   :: xi_w
        REAL(KIND=8), DIMENSION(nb_ns)   :: xi_g
        REAL(KIND=8), DIMENSION(nb_ns)   :: rhoi_wall
        REAL(KIND=8), DIMENSION(nb_eq)   :: prim
        REAL(KIND=8), DIMENSION(nb_temp) :: T1
        REAL(KIND=8), DIMENSION(2)       :: seb_vec, Twall_vec 


        INTEGER, INTENT(IN) :: id                               !< boundary id
        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: u_phys1      !< conservative variables of the first physical cell
        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: phys_data1   !< physical properties of the first physical cell
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost1     !< conservative variables of the first ghost cell
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data1  !< physical properties of the first ghost cell

        ! Getting the temperature of the two physical state next to the wall
        DO i = 1,nb_temp
           T1(i) = phys_data1(pos_T_cell + i - 1)
        ENDDO
        T2 = cell_prop(start_prop_phys + nb_prop + pos_T_cell)
        !T2 = phys_data2(pos_T_cell)
        ! First guess of the wall temperature by linear extrapolation of the temperature (tranlational if multiple temperatures) in the first two physical states
        Twall_vec(1) = T1(1)-(T2-T1(1))*0.5d0
        Twall = Twall_vec(1)
            
        ! Species molar fractions of physical state 
        CALL library_get_molar_fractions(u_phys1(1:nb_ns), xi)
        CALL library_comp_tol(xi) 


        ! Wall species densities computed based on the molar fractions and pressureo
        ! of physical state and first guess wall temperature
        p_ov_T = phys_data1(pos_pres_cell)/Twall
        p_wall = phys_data1(pos_pres_cell)
        DO i = 1,nb_ns 
        ! Constant extrapolation of mole fractions
                rhoi_wall(i) = p_ov_T*xi(i)/Ri(i)
        ENDDO

        ! Call surface mass balance subroutine
        CALL  procedure_ablation_SMB(u_phys1,p_wall,rhoi_wall)

        ! Get mole fractions of the wall interface
        !CALL library_get_molar_fractions(rhoi_wall(1:nb_ns), xi_w)


        !************************************************************
        !LOOP FOR NEWTON-RAPHSON FOR SEB

        DO j = 1,j_max


                ! Library that compute the Surface Energy Balance (SEB)
                CALL procedure_ablation_compute_wall_SEB(rhoi_wall,Twall,T1,seb)  
                seb_vec(1) = seb 

                  
                ! Wall temperature perturbation
                Twall = Twall*(1.d0+eps)
                Twall_vec(2) = Twall

                !p_ov_T = phys_data1(pos_pres_cell)/Twall
                !DO i = 1,nb_ns 
                !   rhoi_wall(i) = p_ov_T*xi_w(i)/Ri(i)
                !ENDDO
                ! Library that compute the total mass flow rate to be used in the evaluation of the wall normal velocity
                !CALL library_ablation_compute_wall_mass_blowing_rate(rhoi_wall,Twall,mdot_wall_tot)
                ! Mass blowing rate limiter to avoid too high values
                !IF(mdot_wall_tot.gt.100.d0) mdot_wall_tot = 100.d0
                ! Library that compute the source terms
                !CALL procedure_ablation_compute_wall_source_terms(rhoi_wall,Twall) 

                ! Call surface mass balance subroutine
                CALL  procedure_ablation_SMB(u_phys1,p_wall,rhoi_wall)
                                    
                ! Library that compute the Surface Energy Balance (SEB)
                CALL procedure_ablation_compute_wall_SEB(rhoi_wall,Twall,T1,seb)  
                seb_vec(2) = seb 

                ! New wall temperature evaluation (Newton)
                Twall_vec(1) = Twall_vec(1)  - seb_vec(1)  / (seb_vec(2) - seb_vec(1)) * (Twall_vec(2) - Twall_vec(1))

                ! Wall temperature limiters to avoid too high/low values
                IF(Twall_vec(1).lt.300.d0)  Twall_vec(1) = 300.d0
                IF(Twall_vec(1).gt.5000.d0) Twall_vec(1) = 5000.d0

                ! "New" wall temperature at the end of the Newton's procedure
                Twall = Twall_vec(1)

        ENDDO
            
        !************************************************************


        ! Call surface mass balance subroutine
        CALL  procedure_ablation_SMB(u_phys1,p_wall,rhoi_wall)
                            
        ! Library that compute the Surface Energy Balance (SEB) only needed for checking purposes
        !CALL procedure_ablation_compute_wall_SEB(rhoi_wall,Twall,T1,seb)  

        !************************************************************
        ! Extrapolation of the first ghost state through mole fractions computed
        ! in the SMB subroutine. This ensure the that the same diffusive flux
        ! will result when the ghost state will be used to compute it at the wall
        ! interface (centered difference using ghost and physical state)
        
        ! Get mole fractions of the wall interface
        CALL library_get_molar_fractions(rhoi_wall(1:nb_ns), xi_w)
        CALL library_comp_tol(xi_w) 
        
        ! Linear extrapolation of ghost state mole fraction
        sum_xi_g = 0.d0
        DO i = 1,nb_ns 
           xi_g(i) = xi_w(i) - (xi(i)-xi_w(i))
           sum_xi_g = sum_xi_g + xi_g(i)
        ENDDO
        
        ! Normalization of mole fractions 
        xi_g = xi_g / sum_xi_g
        
        ! Linear extrapolation of temperature
        DO i = 1,nb_temp
           prim(pos_T + i - 1) = MAX(Twall - (T1(i)-Twall), 0.5d0 * Twall) 
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
        prim(pos_v) = -phys_data1(pos_v_cell)

        ! Compute the ghost state conservative variables and physical properties from primitive variables 
        CALL get_cons_phys_from_prim(prim, u_ghost1, ghost_data1)


        END SUBROUTINE no_slip_SEB_ablation_neq_1D_SL_Expl
!------------------------------------------------------------------------------!
#endif

