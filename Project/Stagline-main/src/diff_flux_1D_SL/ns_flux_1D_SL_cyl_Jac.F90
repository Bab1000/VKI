!------------------------------------------------------------------------------!
!> This subroutine computes the diffusive flux and the related Jacobians when solving the 1D stagnation line equations 
!! for a calorically perfect gas (cylinder case).
  SUBROUTINE ns_flux_1D_SL_cyl_Jac (r_l, r_r, vol_l, vol_r, left_data, right_data, u_left, u_right, fd, jfdl, jfdr)

    USE mod_general_data,               ONLY: pos_u_cell, pos_v_cell, pos_T_cell, pos_mu_cell, pos_lambda_cell, & 
                                            & pos_rho_cell, pos_rho, pos_rhou, pos_rhov, pos_rhoE, gamma, R_gas
    USE mod_numerics_data,              ONLY: Ad, Bd
    USE mod_function_pointer,           ONLY: get_stress_tensor_1D_SL

    IMPLICIT NONE

    REAL(KIND=8), PARAMETER :: coeff1 = 2.d0/3.d0
    REAL(KIND=8), PARAMETER :: coeff2 = 2.d0*coeff1
    REAL(KIND=8) :: mu, lambda
    REAL(KIND=8) :: r, rho, u, v, vel_sum, us, T, du_dr, dv_dr, dT_dr
    REAL(KIND=8) :: ov_dr, ov_r, ov_rho, mu_ov_rho, u_mu_ov_rho, v_mu_ov_rho, us_mu_ov_rho, lam_ov_rho_T,  & 
                  & lam_ov_rhoR_gm1, coeff_mu_ov_rho
    REAL(KIND=8) :: tau_rr, tau_rt, tau_tt
    REAL(KIND=8) :: rhol, rhor, ul, ur, vl, vr, Tl, Tr, mul, mur, lambdal, lambdar

    REAL(KIND=8), INTENT(IN) :: r_l                       !< radial position of left state
    REAL(KIND=8), INTENT(IN) :: r_r                       !< radial position of rigth state
    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state 
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: fd         !< diffusive flux
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfdl     !< diffusive flux Jacobian with respect to the left state
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfdr     !< diffusive flux Jacobian with respect to the right state

    ! Cell interface position
    r = 0.5d0*(r_l + r_r)

    ! Velocity components, temperature, dynamic viscosity and thermal conductivity of left and right states
    ! Left state
    ul   = left_data(pos_u_cell)
    vl   = left_data(pos_v_cell)
    Tl   = left_data(pos_T_cell)
    rhol = left_data(pos_rho_cell)
    mul  = left_data(pos_mu_cell)
    lambdal = left_data(pos_lambda_cell)

    ! Right state
    ur   = right_data(pos_u_cell)
    vr   = right_data(pos_v_cell)
    Tr   = right_data(pos_T_cell)
    rhor = right_data(pos_rho_cell) 
    mur  = right_data(pos_mu_cell)
    lambdar = right_data(pos_lambda_cell)

    ! Velocity and temperature gradients
    ov_dr = 2.d0/(vol_l + vol_r)
    du_dr = ov_dr*(ur - ul)
    dv_dr = ov_dr*(vr - vl)
    dT_dr = ov_dr*(Tr - Tl)

    ! Velocity components, temperature, density, viscosity and thermal conductivity at cell interface
    u   = 0.5d0*(ul + ur)
    v   = 0.5d0*(vl + vr)
    T   = 0.5d0*(Tl + Tr)
    rho = 0.5d0*(rhol + rhor)
    mu  = 0.5d0*(mul + mur) 
    lambda = 0.5d0*(lambdal + lambdar)

    ! Stress tensor
    CALL get_stress_tensor_1D_SL (mu, r, u, v, du_dr, dv_dr, tau_rr, tau_rt, tau_tt)

    ! Diffusive flux  
    fd(pos_rho)  = 0.d0
    fd(pos_rhou) = tau_rr
    fd(pos_rhov) = tau_rt
    fd(pos_rhoE) = tau_rr*u + lambda*dT_dr

    ! Diffusive flux Jacobians 
    ! Ad matrix
    ! Common factors
    ! For the evaluation of the Ad and Bd matrices the radial velocity at the cell interface 
    ! is set equal to the left state velocity for sake of robustness 
    u  = ul
    us = u**2 
    ov_rho = ov_dr/rho 
    mu_ov_rho = mu*ov_rho
    u_mu_ov_rho  = u*mu_ov_rho
    us_mu_ov_rho = u*u_mu_ov_rho
    v_mu_ov_rho  = v*mu_ov_rho
    lam_ov_rho_T = lambda*ov_rho*T
    lam_ov_rhoR_gm1 = lambda*ov_rho*(gamma - 1.d0)/R_gas
   
    ! First column 
    Ad(pos_rho,pos_rho)  = 0.d0
    Ad(pos_rhou,pos_rho) = - coeff2*u_mu_ov_rho 
    Ad(pos_rhov,pos_rho) = - v_mu_ov_rho
    Ad(pos_rhoE,pos_rho) = - coeff2*us_mu_ov_rho + 0.5d0*lam_ov_rhoR_gm1*us - lam_ov_rho_T
 
    ! Second column 
    Ad(pos_rho,pos_rhou)  = 0.d0
    Ad(pos_rhou,pos_rhou) = coeff2*mu_ov_rho
    Ad(pos_rhov,pos_rhou) = 0.d0
    Ad(pos_rhoE,pos_rhou) = coeff2*u_mu_ov_rho - lam_ov_rhoR_gm1*u

    ! Third column 
    Ad(pos_rho,pos_rhov)  = 0.d0
    Ad(pos_rhou,pos_rhov) = 0.d0
    Ad(pos_rhov,pos_rhov) = mu_ov_rho
    Ad(pos_rhoE,pos_rhov) = 0.d0

    ! Fourth column 
    Ad(pos_rho,pos_rhoE)  = 0.d0
    Ad(pos_rhou,pos_rhoE) = 0.d0
    Ad(pos_rhov,pos_rhoE) = 0.d0
    Ad(pos_rhoE,pos_rhoE) = lam_ov_rhoR_gm1

    ! Bd matrix
    ! Common factors
    vel_sum   = (u + v)
    mu_ov_rho = mu/(rho*r) 
    coeff_mu_ov_rho = coeff1*mu_ov_rho

    ! First column 
    Bd(pos_rho,pos_rho)  = 0.d0
    Bd(pos_rhou,pos_rho) = coeff_mu_ov_rho*vel_sum
    Bd(pos_rhov,pos_rho) = mu_ov_rho*vel_sum 
    Bd(pos_rhoE,pos_rho) = 2.d0*coeff_mu_ov_rho*vel_sum*u

    ! Second column 
    Bd(pos_rho,pos_rhou)  = 0.d0
    Bd(pos_rhou,pos_rhou) = - coeff_mu_ov_rho
    Bd(pos_rhov,pos_rhou) = - mu_ov_rho
    Bd(pos_rhoE,pos_rhou) = - coeff_mu_ov_rho*(u + vel_sum)

    ! Third column 
    Bd(pos_rho,pos_rhov)  = 0.d0
    Bd(pos_rhou,pos_rhov) = - coeff_mu_ov_rho
    Bd(pos_rhov,pos_rhov) = - mu_ov_rho
    Bd(pos_rhoE,pos_rhov) = - coeff_mu_ov_rho*u

    ! Fourth column 
    Bd(pos_rho,pos_rhoE)  = 0.d0
    Bd(pos_rhou,pos_rhoE) = 0.d0
    Bd(pos_rhov,pos_rhoE) = 0.d0
    Bd(pos_rhoE,pos_rhoE) = 0.d0 
  
    ! Left and right state diffusive flux Jacobians
    jfdr = Ad + Bd
    jfdl = - Ad + Bd
    
  END SUBROUTINE ns_flux_1D_SL_cyl_Jac
!------------------------------------------------------------------------------!

