#include"../config.h"

#ifdef CARBONABLA
!------------------------------------------------------------------------------!
!> This subroutine applies an ablative boundary condition in case of
!! carbon-based (i.e. graphite, carbon preform, carbon-phenolic
!carbon-carbon) wall by solving the Surface Mass Balance (SMB) and
!assuming an isothermal wall with the temperature given in input for 1D stagnation line 
!! nonequilibrium flows. 

        SUBROUTINE no_slip_isothermal_ablation_neq_1D_SL_Expl_2nd (id, phys_data1, u_phys1, phys_data2, u_phys2,   & 
        &  ghost_data1, u_ghost1, ghost_data2, u_ghost2)

        USE mod_general_data,          ONLY: nb_ns, nb_temp, nb_eq, pos_u, pos_v, pos_T,  & 
                                           & pos_v_cell, pos_u_cell, pos_T_cell, pos_pres_cell, Ri 
        USE mod_domain_boundary,       ONLY: boundary, boundary_data, get_boundary_Twall
        USE mod_function_pointer,      ONLY: get_cons_phys_from_prim
        USE mod_neq_function_pointer,  ONLY: library_get_molar_fractions
        USE mod_ablation_data,         ONLY: uwall, Twall
        USE mod_ablation_pointers,     ONLY: procedure_ablation_SMB
        IMPLICIT NONE

        INTEGER :: i
        REAL(KIND=8) :: p_ov_T, sum_xi_g
        REAL(KIND=8) :: T1
        REAL(KIND=8) :: p_wall
        REAL(KIND=8), DIMENSION(nb_ns)   :: xi
        REAL(KIND=8), DIMENSION(nb_ns)   :: xi2
        REAL(KIND=8), DIMENSION(nb_ns)   :: rhoi_wall, xi_w
        REAL(KIND=8), DIMENSION(nb_ns)   :: xi_g
        REAL(KIND=8), DIMENSION(nb_eq)   :: prim
        TYPE(boundary) :: bound


        INTEGER, INTENT(IN) :: id                               !< boundary id
        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: u_phys1      !< conservative variables of the first physical cell
        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: u_phys2      !< conservative variables of the second physical cell
        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: phys_data1   !< physical properties of the first physical cell
        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: phys_data2   !< physical properties of the second physical cell
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost1     !< conservative variables of the first ghost cell
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost2     !< conservative variables of the second ghost cell
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data1  !< physical properties of the first ghost cell
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data2  !< physical properties of the second ghost cell

        ! Wall temperature
        bound = boundary_data(id)
        Twall = get_boundary_Twall(bound)

        ! Getting the temperature of the physical state next to the wall
        T1 = phys_data1(pos_T_cell)
            
        ! Species molar fractions of physical state 
        CALL library_get_molar_fractions(u_phys1(1:nb_ns), xi)
        CALL library_get_molar_fractions(u_phys2(1:nb_ns), xi2)
        
        ! Parabolic extrapolation (zero derivative at wall) of wall pressure. 
        ! Doesn't seem to do any better than the constant (below)
        !p_wall = (9.d0 * phys_data1(pos_pres_cell) - phys_data2(pos_pres_cell)) / 8.d0

        ! Constant extrapolation of wall pressure
        p_wall = phys_data1(pos_pres_cell)

        ! Wall species densities computed based on the linearly extrapolated 
        !molar fractions and wall pressure and temperature
        p_ov_T = p_wall/Twall
        DO i = 1,nb_ns 
           rhoi_wall(i) = p_ov_T*0.5d0*(3.d0*xi(i)-xi2(i))/Ri(i)
        ENDDO

        ! Call surface mass balance subroutine
        CALL  procedure_ablation_SMB(u_phys1,p_wall,rhoi_wall)

        ! Ghost state 1
        ! Extrapolation of the first ghost state through mole fractions computed
        ! in the SMB subroutine. This ensure the that the same diffusive flux
        ! will result when the ghost state will be used to compute it at the wall
        ! interface (centered difference using ghost and physical state)

        ! Get molar fraction of the wall interface
        CALL library_get_molar_fractions(rhoi_wall(1:nb_ns), xi_w)

        ! Linear extrapolation of the first ghost state molar fraction
        sum_xi_g = 0.d0
        DO i = 1,nb_ns 
           xi_g(i) = xi_w(i) - (xi(i)-xi_w(i))
           sum_xi_g = sum_xi_g + xi_g(i)
        ENDDO

        ! Normalization of molar fractions 
        xi_g = xi_g / sum_xi_g


        ! Linear extrapolation of temperature
        prim(pos_T) = MAX(Twall - (T1-Twall), 0.5d0 * Twall) 

        ! Ghost pressure set equal to physical one
        p_ov_T = phys_data1(pos_pres_cell) / prim(pos_T)

        ! Computation of ghost species densities
        DO i = 1,nb_ns 
           prim(i) = p_ov_T*xi_g(i)/Ri(i)
        ENDDO

        ! Extrapolation of ghost velocities
        prim(pos_u) = uwall - (phys_data1(pos_u_cell)-uwall)
        prim(pos_v) = -phys_data1(pos_v_cell)

        ! Compute the ghost state conservative variables and physical properties from primitive variables 
        CALL get_cons_phys_from_prim(prim, u_ghost1, ghost_data1)

        ! Ghost state 2 

        ! Linear extrapolation of the second ghost state
        DO i = 1,nb_ns 
                prim(i) = prim(i) - (rhoi_wall(i)-prim(i)) * 2.d0
                !Check to avoid negative species densities
                if  (prim(i).lt.1.d-30) prim(i) = 1.d-30
        ENDDO
        prim(pos_T) = MAX(prim(pos_T) - (Twall-prim(pos_T)) * 2.d0, 0.5d0 * Twall) 
        prim(pos_u) =  prim(pos_u) - (uwall-prim(pos_u)) * 2.d0
        prim(pos_v) = -phys_data2(pos_v_cell)

        ! Compute the ghost state conservative variables and physical properties from primitive variables 
        CALL get_cons_phys_from_prim(prim, u_ghost2, ghost_data2)

        END SUBROUTINE no_slip_isothermal_ablation_neq_1D_SL_Expl_2nd
!------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------------!
!> This subroutine applies an ablative boundary condition in case of
!! non-pyrolyzing carbon-based (i.e. graphite, carbon preform, carbon-carbon) wall by solving the Surface Mass Balance (SMB) and
!! the Surface Energy Balance (SEB) for 1D stagnation line 
!! nonequilibrium flows. 

        SUBROUTINE no_slip_SEB_ablation_neq_1D_SL_Expl_2nd (id, phys_data1, u_phys1, phys_data2, u_phys2,   & 
        &  ghost_data1, u_ghost1, ghost_data2, u_ghost2)

        USE mod_general_data,          ONLY: nb_ns, nb_temp, nb_eq, pos_u, pos_v, pos_T,  & 
                                           & pos_v_cell, pos_u_cell, pos_T_cell, pos_pres_cell, Ri 
        USE mod_function_pointer,      ONLY: get_cons_phys_from_prim
        USE mod_neq_function_pointer,  ONLY: library_get_molar_fractions, library_comp_tol
        USE mod_ablation_data,         ONLY: uwall, Twall, phi_pyro
        USE mod_ablation_num_param,    ONLY: j_max, eps
        USE mod_ablation_pointers,     ONLY: procedure_ablation_SMB, procedure_ablation_compute_wall_SEB
        IMPLICIT NONE

        INTEGER :: i, j
        REAL(KIND=8) :: p_ov_T, sum_xi_g
        REAL(KIND=8) :: T2
        REAL(KIND=8) :: mdot_wall_tot
        REAL(KIND=8) :: p_wall
        REAL(KIND=8) :: rhowall 
        REAL(KIND=8) :: seb
        REAL(KIND=8), DIMENSION(nb_ns)   :: xi
        REAL(KIND=8), DIMENSION(nb_ns)   :: xi2
        REAL(KIND=8), DIMENSION(nb_ns)   :: xi_g
        REAL(KIND=8), DIMENSION(nb_ns)   :: xi_w
        REAL(KIND=8), DIMENSION(nb_ns)   :: rhoi_wall
        REAL(KIND=8), DIMENSION(nb_eq)   :: prim
        REAL(KIND=8), DIMENSION(nb_temp) :: T1
        REAL(KIND=8), DIMENSION(2)       :: seb_vec, Twall_vec 


        INTEGER, INTENT(IN) :: id                               !< boundary id
        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: u_phys1      !< conservative variables of the first physical cell
        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: u_phys2      !< conservative variables of the second physical cell
        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: phys_data1   !< physical properties of the first physical cell
        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: phys_data2   !< physical properties of the second physical cell
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost1     !< conservative variables of the first ghost cell
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost2     !< conservative variables of the second ghost cell
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data1  !< physical properties of the first ghost cell
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data2  !< physical properties of the second ghost cell


        ! Species molar fractions of physical state 
        CALL library_get_molar_fractions(u_phys1(1:nb_ns), xi)
        CALL library_comp_tol(xi) 
        CALL library_get_molar_fractions(u_phys2(1:nb_ns), xi2)
        CALL library_comp_tol(xi2) 
        
        ! Getting the temperature of the two physical state next to the wall
        ! Getting the temperature of the two physical state next to the wall
        DO i = 1,nb_temp
           T1(i) = phys_data1(pos_T_cell + i - 1)
        ENDDO
        T2 = phys_data2(pos_T_cell)

        ! First guess of the wall temperature by linear extrapolation of the temperature (tranlational if multiple temperatures) in the first two physical states
        Twall_vec(1) = T1(1)-(T2-T1(1))*0.5d0
        Twall = Twall_vec(1)

        ! Constant extrapolation of wall pressure
        p_wall = phys_data1(pos_pres_cell)

        ! Wall species densities computed based on the linearly extrapolated 
        !molar fractions and wall pressure and temperature
        p_ov_T = p_wall/Twall
        DO i = 1,nb_ns 
           rhoi_wall(i) = p_ov_T*0.5d0*(3.d0*xi(i)-xi2(i))/Ri(i)
        ENDDO

        ! Call surface mass balance subroutine
        CALL  procedure_ablation_SMB(u_phys1,p_wall,rhoi_wall)

        !************************************************************
        !LOOP FOR NEWTON-RAPHSON FOR SEB

        DO j = 1,j_max


                ! Library that compute the Surface Energy Balance (SEB)
                CALL procedure_ablation_compute_wall_SEB(rhoi_wall,Twall,T1,seb)  
                seb_vec(1) = seb 

                  
                ! Wall temperature perturbation
                Twall = Twall*(1.d0+eps)
                Twall_vec(2) = Twall

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

        ! Ghost state 1
        ! Extrapolation of the first ghost state through mole fractions computed
        ! in the SMB subroutine. This ensure the that the same diffusive flux
        ! will result when the ghost state will be used to compute it at the wall
        ! interface (centered difference using ghost and physical state)

        ! Get molar fraction of the wall interface
        CALL library_get_molar_fractions(rhoi_wall(1:nb_ns), xi_w)
        CALL library_comp_tol(xi_w) 

        ! Linear extrapolation of the first ghost state molar fraction
        sum_xi_g = 0.d0
        DO i = 1,nb_ns 
           xi_g(i) = xi_w(i) - (xi(i)-xi_w(i))
           sum_xi_g = sum_xi_g + xi_g(i)
        ENDDO

        ! Normalization of molar fractions 
        xi_g = xi_g / sum_xi_g


        ! Linear extrapolation of temperature
        DO i = 1,nb_temp
           prim(pos_T + i - 1) = MAX(Twall - (T1(i)-Twall), 0.5d0 * Twall) 
        ENDDO

        ! Ghost pressure set equal to physical one
        p_ov_T = phys_data1(pos_pres_cell) / prim(pos_T)

        ! Computation of ghost species densities
        DO i = 1,nb_ns 
           prim(i) = p_ov_T*xi_g(i)/Ri(i)
        ENDDO

        ! Extrapolation of ghost velocities
        prim(pos_u) = uwall - (phys_data1(pos_u_cell)-uwall)
        prim(pos_v) = -phys_data1(pos_v_cell)

        ! Compute the ghost state conservative variables and physical properties from primitive variables 
        CALL get_cons_phys_from_prim(prim, u_ghost1, ghost_data1)

        ! Ghost state 2 

        ! Linear extrapolation of the second ghost state
        DO i = 1,nb_ns 
                prim(i) = prim(i) - (rhoi_wall(i)-prim(i)) * 2.d0
                !Check to avoid negative species densities
                if  (prim(i).lt.1.d-30) prim(i) = 1.d-30
        ENDDO
        prim(pos_T) = MAX(prim(pos_T) - (Twall-prim(pos_T)) * 2.d0, 0.5d0 * Twall) 
        prim(pos_u) =  prim(pos_u) - (uwall-prim(pos_u)) * 2.d0
        prim(pos_v) = -phys_data2(pos_v_cell)

        ! Compute the ghost state conservative variables and physical properties from primitive variables 
        CALL get_cons_phys_from_prim(prim, u_ghost2, ghost_data2)

        END SUBROUTINE no_slip_SEB_ablation_neq_1D_SL_Expl_2nd
!------------------------------------------------------------------------------!
#endif
