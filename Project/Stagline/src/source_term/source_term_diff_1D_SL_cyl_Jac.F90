!------------------------------------------------------------------------------!
!> This subroutine computes the diffusive source term and the related Jacobians for 1D stagnation line flows 
!! for calorically perfect gas (cylinder case).
  SUBROUTINE source_term_diff_1D_SL_cyl_Jac(r_c, vol_l, vol_c, vol_r, prop_left, prop_cell, prop_right, & 
                                          & u_left, u_cell, u_right, s, js_left, js_cell, js_right)

    
    USE mod_general_data,           ONLY: pos_u_cell, pos_v_cell, pos_T_cell, pos_mu_cell, pos_lambda_cell, &
                                        & pos_rho_cell, pos_rho, pos_rhou, pos_rhov, pos_rhoE, gamma, R_gas
    USE mod_numerics_data,          ONLY: Ad, Bd
    USE mod_function_pointer,       ONLY: get_stress_tensor_1D_SL

    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8), PARAMETER :: coeff0 = 1.d0/3.d0
    REAL(KIND=8), PARAMETER :: coeff1 = 2.d0*coeff0
    REAL(KIND=8), PARAMETER :: coeff2 = 5*coeff1
    REAL(KIND=8) :: ov_dr, ov_rc
    REAL(KIND=8) :: mu, lambda
    REAL(KIND=8) :: rho, u, v, uv, us, vs, vel_sum, T, du_dr, dv_dr, dT_dr
    REAL(KIND=8) :: ov_rho, mu_ov_rho, m2_mu_ov_rho, m4_mu_ov_rho, m3_mu_ov_rho, lam_ov_rho_T, lam_ov_rhoR_gm1
    REAL(KIND=8) :: coeff0_mu_ov_rho, coeff1_mu_ov_rho, coeff2_mu_ov_rho
    REAL(KIND=8) :: tau_rr, tau_rt, tau_tt
    REAL(KIND=8) :: ul, ur, vl, vr, Tl, Tr

    REAL(KIND=8), INTENT(IN) :: r_c                        !< radial location of cell centroid
    REAL(KIND=8), INTENT(IN) :: vol_l                      !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_c                      !< volume of central state
    REAL(KIND=8), INTENT(IN) :: vol_r                      !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left       !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_cell       !< conservative variables of central state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right      !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prop_left    !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prop_cell    !< physical properties of central state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prop_right   !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s           !< source term
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: js_left   !< source term Jacobian with respect to the left state
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: js_cell   !< source term Jacobian with respect to the central state
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: js_right  !< source term Jacobian with respect to the right state 

    ! Velocity components, temperature, dynamic viscosity and thermal conductivity of left and right states
    ! Left state
    ul = prop_left(pos_u_cell)
    vl = prop_left(pos_v_cell)
    Tl = prop_left(pos_T_cell)

    ! Center state
    u   = prop_cell(pos_u_cell)
    v   = prop_cell(pos_v_cell)
    rho = prop_cell(pos_rho_cell)
    T   = prop_cell(pos_T_cell)
    mu  = prop_cell(pos_mu_cell)
    lambda = prop_cell(pos_lambda_cell)

    ! Right state
    ur = prop_right(pos_u_cell)
    vr = prop_right(pos_v_cell)
    Tr = prop_right(pos_T_cell) 

    ! Velocity and temperature gradients 
    ov_dr = 2.d0/(vol_l + 2.d0*vol_c + vol_r)
    du_dr = ov_dr*(ur - ul)
    dv_dr = ov_dr*(vr - vl)
    dT_dr = ov_dr*(Tr - Tl)

    ! Stress tensor
    CALL get_stress_tensor_1D_SL (mu, r_c, u, v, du_dr, dv_dr, tau_rr, tau_rt, tau_tt)

    ! Diffusive source term
    ov_rc = 1.d0/r_c
    s(pos_rho)  = 0.d0
    s(pos_rhou) = (tau_rr + tau_rt - tau_tt)*ov_rc
    s(pos_rhov) = (2.d0*tau_rt - tau_tt)*ov_rc
    s(pos_rhoE) = (u*(tau_rr + tau_rt) + v*tau_tt + lambda*dT_dr)*ov_rc

    ! Diffusive source term Jacobians
    ! Common factors
    us = u**2
    vs = v**2
    uv = u*v
    vel_sum = u + v
    ov_rho = ov_dr*ov_rc/rho 
    mu_ov_rho = mu*ov_rho
    m2_mu_ov_rho = 2.d0*mu_ov_rho
    coeff1_mu_ov_rho = coeff1*mu_ov_rho
    lam_ov_rho_T     = lambda*ov_rho*T
    lam_ov_rhoR_gm1  = lambda*ov_rho*(gamma - 1.d0)/R_gas

    ! Ad matrix
    ! First column
    Ad(pos_rho,pos_rho)  = 0.d0
    Ad(pos_rhou,pos_rho) = - mu_ov_rho*(vel_sum + u)
    Ad(pos_rhov,pos_rho) = - coeff1_mu_ov_rho*(2.d0*u + 3.d0*v)
    Ad(pos_rhoE,pos_rho) = - coeff0*mu_ov_rho*u*(vel_sum + 3.d0*u) + lam_ov_rhoR_gm1*us - lam_ov_rho_T 

    ! Second column 
    Ad(pos_rho,pos_rhou)  = 0.d0
    Ad(pos_rhou,pos_rhou) = m2_mu_ov_rho + m2_mu_ov_rho
    Ad(pos_rhov,pos_rhou) = coeff1_mu_ov_rho
    Ad(pos_rhoE,pos_rhou) = coeff1_mu_ov_rho*(2.d0*u - v) - lam_ov_rhoR_gm1*u

    ! Third column 
    Ad(pos_rho,pos_rhov)  = 0.d0
    Ad(pos_rhou,pos_rhov) = mu_ov_rho
    Ad(pos_rhov,pos_rhov) = m2_mu_ov_rho
    Ad(pos_rhoE,pos_rhov) = mu_ov_rho*u

    ! Fourth colum
    Ad(pos_rho,pos_rhoE)  = 0.d0
    Ad(pos_rhou,pos_rhoE) = 0.d0
    Ad(pos_rhov,pos_rhoE) = 0.d0
    Ad(pos_rhoE,pos_rhoE) = lam_ov_rhoR_gm1

    ! Bd matrix
    ! Common factors
    mu_ov_rho = (mu/rho)*ov_rc**2
    m3_mu_ov_rho = 3.d0*mu_ov_rho
    coeff0_mu_ov_rho = coeff0*mu_ov_rho 
    coeff1_mu_ov_rho = coeff1*mu_ov_rho
    coeff2_mu_ov_rho = coeff2*mu_ov_rho
 
    ! First column
    Bd(pos_rho,pos_rho)  = 0.d0
    Bd(pos_rhou,pos_rho) = m3_mu_ov_rho*vel_sum
    Bd(pos_rhov,pos_rho) = coeff2_mu_ov_rho*vel_sum
    Bd(pos_rhoE,pos_rho) = coeff1_mu_ov_rho*(5.d0*us + uv - 4.d0*vs)

    ! Second column 
    Bd(pos_rho,pos_rhou)  = 0.d0
    Bd(pos_rhou,pos_rhou) = - m3_mu_ov_rho
    Bd(pos_rhov,pos_rhou) = - coeff2_mu_ov_rho
    Bd(pos_rhoE,pos_rhou) = coeff0_mu_ov_rho*(10.d0*u - v)

    ! Third column
    Bd(pos_rho,pos_rhov)  = 0.d0
    Bd(pos_rhou,pos_rhov) = - m3_mu_ov_rho
    Bd(pos_rhov,pos_rhov) = - coeff2_mu_ov_rho
    Bd(pos_rhoE,pos_rhov) = coeff0_mu_ov_rho*(8.d0*v - u)

    ! Fourth column 
    Bd(pos_rho,pos_rhoE)  = 0.d0
    Bd(pos_rhou,pos_rhoE) = 0.d0
    Bd(pos_rhov,pos_rhoE) = 0.d0
    Bd(pos_rhoE,pos_rhoE) = 0.d0
    
    ! Left, central and right state diffusive source term Jacobians
    js_left  = - Ad
    js_cell  = Bd
    js_right = Ad
    
  END SUBROUTINE source_term_diff_1D_SL_cyl_Jac
!------------------------------------------------------------------------------!
