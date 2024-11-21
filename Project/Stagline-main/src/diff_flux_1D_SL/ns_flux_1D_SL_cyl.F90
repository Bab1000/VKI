!------------------------------------------------------------------------------!
!> This subroutine computes the diffusive flux when solving the 1D stagnation line equations for a calorically perfect gas.
!! (cylinder case)
  SUBROUTINE ns_flux_1D_SL_cyl (r_l, r_r, vol_l, vol_r, left_data, right_data, u_left, u_right, fd)

    USE mod_general_data,               ONLY: pos_u_cell, pos_v_cell, pos_T_cell, pos_mu_cell, pos_lambda_cell, & 
                                            & pos_rho, pos_rhou, pos_rhov, pos_rhoE
    USE mod_function_pointer,           ONLY: get_stress_tensor_1D_SL

    IMPLICIT NONE

    REAL(KIND=8) :: ov_dr
    REAL(KIND=8) :: mu, lambda
    REAL(KIND=8) :: r, u, v, du_dr, dv_dr, dT_dr
    REAL(KIND=8) :: tau_rr, tau_rt, tau_tt
    REAL(KIND=8) :: ul, ur, vl, vr, Tl, Tr, mul, mur, lambdal, lambdar

    REAL(KIND=8), INTENT(IN) :: r_l                       !< radial position of left state
    REAL(KIND=8), INTENT(IN) :: r_r                       !< radial position of rigth state
    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state 
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: fd         !< diffusive flux

    ! Cell interface position
    r = 0.5d0*(r_l + r_r)

    ! Velocity components, temperature, dynamic viscosity and thermal conductivity of left and right states
    ! Left state
    ul  = left_data(pos_u_cell)
    vl  = left_data(pos_v_cell)
    Tl  = left_data(pos_T_cell)
    mul = left_data(pos_mu_cell)
    lambdal = left_data(pos_lambda_cell)

    ! Right state
    ur  = right_data(pos_u_cell)
    vr  = right_data(pos_v_cell)
    Tr  = right_data(pos_T_cell)
    mur = right_data(pos_mu_cell)
    lambdar = right_data(pos_lambda_cell)

    ! Velocity and temperature gradients
    ov_dr = 2.d0/(vol_l + vol_r)
    du_dr = ov_dr*(ur - ul)
    dv_dr = ov_dr*(vr - vl)
    dT_dr = ov_dr*(Tr - Tl)

    ! Velocity components, viscosity and thermal conductivity at cell interface
    u  = 0.5d0*(ul + ur)
    v  = 0.5d0*(vl + vr)
    mu = 0.5d0*(mul + mur) 
    lambda = 0.5d0*(lambdal + lambdar)

    ! Stress tensor
    CALL get_stress_tensor_1D_SL (mu, r, u, v, du_dr, dv_dr, tau_rr, tau_rt, tau_tt)

    ! Diffusive flux  
    fd(pos_rho)  = 0.d0
    fd(pos_rhou) = tau_rr
    fd(pos_rhov) = tau_rt
    fd(pos_rhoE) = tau_rr*u + lambda*dT_dr 

  END SUBROUTINE ns_flux_1D_SL_cyl
!------------------------------------------------------------------------------!

