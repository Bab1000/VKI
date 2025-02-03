!------------------------------------------------------------------------------!
!> This subroutine computes the convective flux by means of fvs proposed by Steger and Warming for a 
!! 1D stagnation line calorically perfect gas flow.
  SUBROUTINE Steger_Warming_1D_SL (nx, vol_l, vol_r, left_data, right_data, u_left, u_right, f)
  
    USE mod_general_data,                ONLY: nb_eq, pos_u_cell, pos_v_cell, pos_c_cell, pos_rho_cell, pos_ek_cell,  & 
                                             & gamma
    USE mod_numerics_data,               ONLY: fl, fr

    IMPLICIT NONE

    REAL(KIND=8) :: l1, l2, l3, l1p, l2p, l3p, l1m, l2m, l3m
    REAL(KIND=8) :: gamma_minus1, ov_gamma, ov_gamma_minus1, alpha, diff
    REAL(KIND=8) :: cl, cr, ekl, ekr, rhol, rhor, ul, ur, vl, vr, unl, unr
 
    REAL(KIND=8), INTENT(IN) :: nx                        !< normal to the cell interface
    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: f          !< numerical flux
 
    ! Common factors
    ov_gamma = 1.d0/gamma
    gamma_minus1    = gamma - 1.d0 
    ov_gamma_minus1 = 1.d0/gamma_minus1
  
    ! Data of left and right states
    ul   = left_data(pos_u_cell)
    vl   = left_data(pos_v_cell)
    cl   = left_data(pos_c_cell)
    ekl  = left_data(pos_ek_cell)
    rhol = left_data(pos_rho_cell)
    unl  = ul*nx     

    ur   = right_data(pos_u_cell)  
    vr   = right_data(pos_v_cell)  
    cr   = right_data(pos_c_cell)
    ekr  = right_data(pos_ek_cell)
    rhor = right_data(pos_rho_cell)
    unr  = ur*nx

    ! Left state positive split-flux
    ! Eigenvalues
    l1   = unl 
    l2   = unl - cl
    l3   = unl + cl 
    l1p  = 0.5d0*(l1 + ABS(l1))
    l2p  = 0.5d0*(l2 + ABS(l2))
    l3p  = 0.5d0*(l3 + ABS(l3))

    ! Common factors
    alpha = 2*gamma_minus1*l1p + l2p + l3p
    diff  = (l3p - l2p)*cl*nx

    fl(1) = alpha 
    fl(2) = alpha*ul  + diff
    fl(3) = alpha*vl
    fl(4) = alpha*ekl + (l3p + l2p)*(cl**2)*ov_gamma_minus1 + diff*ul

    ! Right state negative split-flux
    ! Eigenvalues
    l1   = unr 
    l2   = unr - cr
    l3   = unr + cr  
    l1m  = 0.5d0*(l1 - ABS(l1))
    l2m  = 0.5d0*(l2 - ABS(l2))
    l3m  = 0.5d0*(l3 - ABS(l3))

    ! Common factors
    alpha = 2*gamma_minus1*l1m + l2m + l3m
    diff  = (l3m - l2m)*cr*nx

    fr(1) = alpha 
    fr(2) = alpha*ur  + diff
    fr(3) = alpha*vr
    fr(4) = alpha*ekr + (l3m + l2m)*(cr**2)*ov_gamma_minus1 + diff*ur 

    ! Sum of left and right split-fluxes
    f = 0.5d0*ov_gamma*(rhol*fl + rhor*fr)

  END SUBROUTINE Steger_Warming_1D_SL
!------------------------------------------------------------------------------!
