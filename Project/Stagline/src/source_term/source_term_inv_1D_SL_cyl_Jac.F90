!-----------------------------------------------------------------------------!
!> This subroutine computes the inviscid source term and the related Jacobian for 1D stagnation line 
!! calorically perfect gas flows (cylinder case).
  SUBROUTINE source_term_inv_1D_SL_cyl_Jac (cell_id, source_id, r, prop, cons, s, js)

    USE mod_general_data,        ONLY: pos_u_cell, pos_v_cell, pos_pres_cell, pos_h0_cell,   & 
                                     & pos_rho_cell, pos_rho, pos_rhou, pos_rhov, pos_rhoE,  & 
                                     & p_inf, gamma

    IMPLICIT NONE

    REAL(KIND=8) :: ov_r
    REAL(KIND=8) :: rho, h0, p, u, u2, v
    REAL(KIND=8) :: rho_vel_sum_ov_r, vel_sum, vel_sum_ov_r, vel_sum_ov_r_u, vel_sum_ov_r_v, & 
                  & u_ov_r, v_ov_r, h0_ov_r
    REAL(KIND=8) :: gamma_minus1, gamma_minus1_u2

    INTEGER, INTENT(IN) :: cell_id                   !< cell identifier
    INTEGER, INTENT(IN) :: source_id                 !< source term identifier
    REAL(KIND=8), INTENT(IN) :: r                    !< cell centroid radial location
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prop   !< physical properties
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons   !< conservative variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s     !< source term
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: js  !< source term Jacobian   
 
    ! Density, velocity components, pressure and specific total enthalpy
    rho = prop(pos_rho_cell)
    u   = prop(pos_u_cell)
    v   = prop(pos_v_cell)
    p   = prop(pos_pres_cell)
    h0  = prop(pos_h0_cell)    

    ! Common factors
    u2 = u**2
    ov_r = 1.d0/r
    u_ov_r = u*ov_r
    v_ov_r = v*ov_r
    h0_ov_r = h0*ov_r
    vel_sum = u + v
    vel_sum_ov_r = vel_sum*ov_r
    vel_sum_ov_r_u = vel_sum_ov_r*u
    vel_sum_ov_r_v = vel_sum_ov_r*v
    rho_vel_sum_ov_r = rho*vel_sum_ov_r
    gamma_minus1 = gamma - 1.d0
    gamma_minus1_u2 = gamma_minus1*u2

    ! Source term
    s(pos_rho)  = - rho_vel_sum_ov_r
    s(pos_rhou) = - rho_vel_sum_ov_r*u
    s(pos_rhov) = - 2.d0*(rho_vel_sum_ov_r*v - (p - p_inf)*ov_r)
    s(pos_rhoE) = - rho_vel_sum_ov_r*h0

    ! Source term Jacobian
    ! First column
    js(pos_rho,pos_rho)  = 0.d0
    js(pos_rhou,pos_rho) = vel_sum_ov_r_u
    js(pos_rhov,pos_rho) = (2.d0*vel_sum_ov_r_v + gamma_minus1_u2*ov_r)
    js(pos_rhoE,pos_rho) = - (0.5d0*gamma_minus1_u2 - h0)*vel_sum_ov_r

    ! Second column 
    js(pos_rho,pos_rhou)  = - ov_r
    js(pos_rhou,pos_rhou) = - (vel_sum_ov_r + u_ov_r)
    js(pos_rhov,pos_rhou) = - 2.d0*(gamma_minus1*u_ov_r + v_ov_r)
    js(pos_rhoE,pos_rhou) = - (h0_ov_r - gamma_minus1*vel_sum_ov_r_u)

    ! Third column
    js(pos_rho,pos_rhov)  = - ov_r
    js(pos_rhou,pos_rhov) = - u_ov_r
    js(pos_rhov,pos_rhov) = - 2.d0*(vel_sum_ov_r + v_ov_r)
    js(pos_rhoE,pos_rhov) = - h0_ov_r

    ! Fourth column
    js(pos_rho,pos_rhoE)  = 0.d0
    js(pos_rhou,pos_rhoE) = 0.d0
    js(pos_rhov,pos_rhoE) = 2.d0*gamma_minus1*ov_r
    js(pos_rhoE,pos_rhoE) = - gamma*vel_sum_ov_r   

  END SUBROUTINE source_term_inv_1D_SL_cyl_Jac
!-----------------------------------------------------------------------------!
