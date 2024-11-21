!------------------------------------------------------------------------------!
!> This subroutine computes the diffusive flux for the 1D Navier-Stokes equations for calorically perfect gas flows.
  SUBROUTINE ns_flux_1D (vol_l, vol_r, left_data, right_data, u_left, u_right, fd)

    USE mod_general_data,       ONLY: pos_u_cell, pos_mu_cell, pos_kappa_cell, pos_lambda_cell, pos_T_cell,  & 
                                    & pos_rho, pos_rhou, pos_rhoE, R_gas

    IMPLICIT NONE

    REAL(KIND=8), PARAMETER :: coeff = 4.d0/3.d0
    REAL(KIND=8) :: tmp
    REAL(KIND=8) :: ov_dx
    REAL(KIND=8) :: kappa, mu, lambda, u
    REAL(KIND=8) :: dT_dx, du_dx
    REAL(KIND=8) :: lambdal, lambdar, kappal, kappar, mul, mur, Tl, Tr, ul, ur

    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: fd         !< diffusive flux

    ! Temperature, velocity, dynamic viscosity, bulk viscosity and thermal conductivity of left and right states
    ! Left state
    ul      = left_data(pos_u_cell)
    Tl      = left_data(pos_T_cell) 
    mul     = left_data(pos_mu_cell)
    kappal  = left_data(pos_kappa_cell)
    lambdal = left_data(pos_lambda_cell) 

    ! Right state
    ur      = right_data(pos_u_cell)
    Tr      = right_data(pos_T_cell)
    mur     = right_data(pos_mu_cell)
    kappar  = right_data(pos_kappa_cell)
    lambdar = right_data(pos_lambda_cell)
    
    ! Velocity and temperature gradients
    ov_dx = 2.d0/(vol_l + vol_r)
    du_dx  = (ur - ul)*ov_dx  
    dT_dx  = (Tr - Tl)*ov_dx

    ! Velocity, dynamic viscosity, bulk viscosity and thermal conductivity at cell interface
    u      = 0.5d0*(ul + ur)
    mu     = 0.5d0*(mul + mur)
    kappa  = 0.5d0*(kappal + kappar)
    lambda = 0.5d0*(lambdal + lambdar) 

    ! Diffusive flux
    tmp   = (coeff*mu + kappa)*du_dx
    fd(pos_rho)  = 0.d0
    fd(pos_rhou) = tmp
    fd(pos_rhoE) = tmp*u + lambda*dT_dx

  END SUBROUTINE ns_flux_1D
!------------------------------------------------------------------------------!
