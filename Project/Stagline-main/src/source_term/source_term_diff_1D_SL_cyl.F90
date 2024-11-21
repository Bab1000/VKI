!------------------------------------------------------------------------------!
!> This subroutine computes the diffusive source term for 1D stagnation line calorically perfect gas flows (cylinder case).
  SUBROUTINE source_term_diff_1D_SL_cyl (r_c, vol_l, vol_c, vol_r, prop_left, prop_cell, prop_right, & 
                                       & u_left, u_cell, u_right, s)

    USE mod_general_data,           ONLY: pos_u_cell, pos_v_cell, pos_T_cell, pos_mu_cell, pos_lambda_cell, &
                                        & pos_rho, pos_rhou, pos_rhov, pos_rhoE
    USE mod_function_pointer,       ONLY: get_stress_tensor_1D_SL

    IMPLICIT NONE

    REAL(KIND=8) :: ov_dr, ov_rc
    REAL(KIND=8) :: mu, lambda
    REAL(KIND=8) :: u, v, du_dr, dv_dr, dT_dr
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

    ! Velocity components, temperature, dynamic viscosity and thermal conductivity of left and right states
    ! Left state
    ul = prop_left(pos_u_cell)
    vl = prop_left(pos_v_cell)
    Tl = prop_left(pos_T_cell)

    ! Center state
    u  = prop_cell(pos_u_cell)
    v  = prop_cell(pos_v_cell)
    mu = prop_cell(pos_mu_cell)
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

  END SUBROUTINE source_term_diff_1D_SL_cyl
!------------------------------------------------------------------------------!
