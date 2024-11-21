!------------------------------------------------------------------------------!
!> This subrotine provides the numerical flux Jacobian for 1D stagnation line calorically perfect gas flows
!! according to the positive and negative split of the inviscid flux Jacobian.
  SUBROUTINE Pos_Neg_split_Jac_1D_SL (nx, vol_l, vol_r, left_data, right_data, u_left, u_right, jfl, jfr)

    USE mod_general_data,               ONLY: nb_eq, pos_u_cell, pos_v_cell, pos_c_cell, pos_ek_cell, pos_h0_cell, gamma

    IMPLICIT NONE

    INTEGER :: i, j
    REAL(KIND=8) :: c, ek, h0, nx_ov_c, ov_c, ov_c2, u, v, v2, vn, vn_ov_c
    REAL(KIND=8) :: l1, l2, l3, l1p, l2p, l3p, l1m, l2m, l3m, l1p_ek, l1m_ek, nx2
    REAL(KIND=8) :: eig_diff, eig_sum,  eig_sum_m_l1p,  eig_sum_m_l1m
    REAL(KIND=8) :: gamma_minus1, gamma_minus1_ov_c, gamma_minus1_ov_c2, gamma_minus1_ov_c2_ek, & 
                  & gamma_minus1_ov_c2_u, gamma_minus1_ov_c2_v2, gamma_minus1_ov_c2_ek_u,       &  
                  & t_m_gamma

    REAL(KIND=8), INTENT(IN) :: nx                        !< normal to the cell interface
    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfl      !< numerical flux Jacobian with respect to the left state 
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfr      !< numerical flux Jacobian with respect to the right state 

    ! Useful quantities
    nx2 = nx**2
    gamma_minus1 = gamma - 1.d0 
    t_m_gamma    = 2.d0 - gamma

    ! Left state
    ! A^+ matrix (positive split Jacobian) 
    c   = left_data(pos_c_cell)
    ek  = left_data(pos_ek_cell)
    h0  = left_data(pos_h0_cell)
    u   = left_data(pos_u_cell)
    v   = left_data(pos_v_cell)
    v2  = u**2
    vn  = u*nx

    ! Eigenvalues 
    l1 = vn 
    l2 = vn - c
    l3 = vn + c

    ! Positive and negative eigenvalues 
    l1p = 0.5d0*(l1 + ABS(l1))
    l2p = 0.5d0*(l2 + ABS(l2))
    l3p = 0.5d0*(l3 + ABS(l3))

    ! Common factors
    ov_c  = 1.d0/c
    ov_c2 = ov_c/c
    nx_ov_c = nx*ov_c
    vn_ov_c = vn*ov_c
    l1p_ek  = l1p*ek
    gamma_minus1_ov_c       = gamma_minus1*ov_c
    gamma_minus1_ov_c2      = gamma_minus1_ov_c*ov_c
    gamma_minus1_ov_c2_ek   = gamma_minus1_ov_c2*ek
    gamma_minus1_ov_c2_ek_u = gamma_minus1_ov_c2_ek*u
    gamma_minus1_ov_c2_u    = gamma_minus1_ov_c2*u
    gamma_minus1_ov_c2_v2   = gamma_minus1_ov_c2*v2

    eig_sum  = 0.5d0*(l3p + l2p)
    eig_diff = 0.5d0*(l3p - l2p) 
    eig_sum_m_l1p = eig_sum - l1p  

    ! Fill A^+ matrix
    ! First column 
    jfl(1,1) = l1p + eig_sum_m_l1p*gamma_minus1_ov_c2_ek - eig_diff*vn_ov_c
    jfl(2,1) = l1p*u + eig_sum_m_l1p*gamma_minus1_ov_c2_ek_u - eig_sum*nx*vn + & 
             & eig_diff*(gamma_minus1_ov_c*ek*nx - u*vn_ov_c)
    jfl(3,1) = gamma_minus1_ov_c2_ek*v*eig_sum_m_l1p - v*vn_ov_c*eig_diff
    jfl(4,1) = l1p_ek + (eig_sum*h0 - l1p_ek)*gamma_minus1_ov_c2_ek - &
             & eig_sum*vn**2 + eig_diff*(gamma_minus1_ov_c*ek*vn - h0*vn_ov_c) 

    ! Second column 
    jfl(1,2) = - eig_sum_m_l1p*gamma_minus1_ov_c2_u + eig_diff*nx_ov_c
    jfl(2,2) = l1p*gamma_minus1_ov_c2_v2 + eig_sum*(nx2 - gamma_minus1_ov_c2_v2) + & 
             & eig_diff*t_m_gamma*vn_ov_c
    jfl(3,2) = - gamma_minus1_ov_c2_u*v*eig_sum_m_l1p + v*nx_ov_c*eig_diff
    jfl(4,2) = l1p*gamma_minus1_ov_c2_ek_u + eig_sum*(vn*nx - gamma_minus1_ov_c2_u*h0) + &
             & eig_diff*(h0*nx_ov_c - gamma_minus1_ov_c*u*vn)

    ! Third column
    jfl(1,3) = 0.d0
    jfl(2,3) = 0.d0
    jfl(3,3) = l1p
    jfl(3,4) = 0.d0

    ! Fourth column 
    jfl(1,4) = eig_sum_m_l1p*gamma_minus1_ov_c2
    jfl(2,4) = eig_sum_m_l1p*gamma_minus1_ov_c2_u + eig_diff*gamma_minus1*nx_ov_c
    jfl(3,4) = gamma_minus1_ov_c2*v*eig_sum_m_l1p
    jfl(4,4) = - l1p*gamma_minus1_ov_c2_ek + eig_sum*gamma_minus1_ov_c2*h0 + & 
             & eig_diff*gamma_minus1_ov_c*vn


    ! Right state
    ! A^- matrix (negative split Jacobian) 
    c   = right_data(pos_c_cell)
    ek  = right_data(pos_ek_cell)
    h0  = right_data(pos_h0_cell)
    u   = right_data(pos_u_cell)
    v   = right_data(pos_v_cell)
    v2  = u**2
    vn  = u*nx

    ! Eigenvalues 
    l1 = vn 
    l2 = vn - c
    l3 = vn + c

    ! Positive and negative eigenvalues 
    l1m = 0.5d0*(l1 - ABS(l1))
    l2m = 0.5d0*(l2 - ABS(l2))
    l3m = 0.5d0*(l3 - ABS(l3))

    ! Common factors
    ov_c  = 1.d0/c
    ov_c2 = ov_c/c
    nx_ov_c = nx*ov_c
    vn_ov_c = vn*ov_c
    l1m_ek  = l1m*ek
    gamma_minus1_ov_c       = gamma_minus1*ov_c
    gamma_minus1_ov_c2      = gamma_minus1_ov_c*ov_c
    gamma_minus1_ov_c2_ek   = gamma_minus1_ov_c2*ek
    gamma_minus1_ov_c2_ek_u = gamma_minus1_ov_c2_ek*u
    gamma_minus1_ov_c2_u    = gamma_minus1_ov_c2*u
    gamma_minus1_ov_c2_v2   = gamma_minus1_ov_c2*v2

    eig_sum  = 0.5d0*(l3m + l2m)
    eig_diff = 0.5d0*(l3m - l2m) 
    eig_sum_m_l1m = eig_sum - l1m 

    ! Fill A^- matrix
    ! First column
    jfr(1,1) = l1m + eig_sum_m_l1m*gamma_minus1_ov_c2_ek - eig_diff*vn_ov_c
    jfr(2,1) = l1m*u + eig_sum_m_l1m*gamma_minus1_ov_c2_ek_u - eig_sum*nx*vn + & 
             & eig_diff*(gamma_minus1_ov_c*ek*nx - u*vn_ov_c)
    jfr(3,1) = gamma_minus1_ov_c2_ek*v*eig_sum_m_l1m - v*vn_ov_c*eig_diff
    jfr(4,1) = l1m_ek + (eig_sum*h0 - l1m_ek)*gamma_minus1_ov_c2_ek - &
             & eig_sum*vn**2 + eig_diff*(gamma_minus1_ov_c*ek*vn - h0*vn_ov_c) 

    ! Second column 
    jfr(1,2) = - eig_sum_m_l1m*gamma_minus1_ov_c2_u + eig_diff*nx_ov_c
    jfr(2,2) = l1m*gamma_minus1_ov_c2_v2 + eig_sum*(nx2 - gamma_minus1_ov_c2_v2) + & 
             & eig_diff*t_m_gamma*vn_ov_c
    jfr(3,2) = - gamma_minus1_ov_c2_u*v*eig_sum_m_l1m + v*nx_ov_c*eig_diff
    jfr(4,2) = l1m*gamma_minus1_ov_c2_ek_u + eig_sum*(vn*nx - gamma_minus1_ov_c2_u*h0) + &
             & eig_diff*(h0*nx_ov_c - gamma_minus1_ov_c*u*vn)

     ! Third column
    jfr(1,3) = 0.d0
    jfr(2,3) = 0.d0
    jfr(3,3) = l1m
    jfr(3,4) = 0.d0

    ! Fourth column 
    jfr(1,4) = eig_sum_m_l1m*gamma_minus1_ov_c2
    jfr(2,4) = eig_sum_m_l1m*gamma_minus1_ov_c2_u + eig_diff*gamma_minus1*nx_ov_c
    jfr(3,4) = gamma_minus1_ov_c2*v*eig_sum_m_l1m
    jfr(4,4) = - l1m*gamma_minus1_ov_c2_ek + eig_sum*gamma_minus1_ov_c2*h0 + & 
             & eig_diff*gamma_minus1_ov_c*vn

  END SUBROUTINE Pos_Neg_split_Jac_1D_SL
!------------------------------------------------------------------------------!
