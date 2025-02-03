!------------------------------------------------------------------------------!
!> This subroutine computes the diffusive flux and diffusive flux Jacobians for the 1D Navier-Stokes equations 
!! for calorically perfect gas flows.
  SUBROUTINE ns_flux_1D_Jac (vol_l, vol_r, left_data, right_data, u_left, u_right, fd, jfdl, jfdr)

    USE mod_general_data,       ONLY: pos_u_cell, pos_rho_cell, pos_mu_cell, pos_kappa_cell, pos_lambda_cell, & 
                                    & pos_T_cell, pos_rho, pos_rhou, pos_rhoE, gamma, R_gas

    IMPLICIT NONE

    REAL(KIND=8), PARAMETER :: coeff = 4.d0/3.d0
    REAL(KIND=8) :: tmp1, tmp2, tmp3, tmp4, tmp5
    REAL(KIND=8) :: ov_dx, gamma_minus1
    REAL(KIND=8) :: kappa, mu, lambda, rho, T, u, v2 
    REAL(KIND=8) :: dT_dx, du_dx
    REAL(KIND=8) :: lambdal, lambdar, kappal, kappar, mul, mur, rhol, rhor, Tl, Tr, ul, ur

    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: fd         !< diffusive flux
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfdl     !< diffusive flux Jacobian with respect to the left state
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfdr     !< diffusive flux Jacobian with respect to the right state

    ! Density, temperature, velocity, dynamic viscosity, bulk viscosity and thermal conductivity of left and right states
    ! Left state
    ul      = left_data(pos_u_cell)
    rhol    = left_data(pos_rho_cell)
    Tl      = left_data(pos_T_cell)
    mul     = left_data(pos_mu_cell)
    kappal  = left_data(pos_kappa_cell)
    lambdal = left_data(pos_lambda_cell) 

    ! Right state
    ur      = right_data(pos_u_cell)
    rhor    = right_data(pos_rho_cell)
    Tr      = right_data(pos_T_cell)
    mur     = right_data(pos_mu_cell)
    kappar  = right_data(pos_kappa_cell)
    lambdar = right_data(pos_lambda_cell)
    
    ! Velocity and temperature gradients
    ov_dx = 2.d0/(vol_l + vol_r)
    du_dx = (ur - ul)*ov_dx  
    dT_dx = (Tr - Tl)*ov_dx

    ! Velocity, temperature, density, dynamic viscosity, bulk viscosity and thermal conductivity at cell interface
    u   = 0.5d0*(ul + ur)
    v2  = u**2
    T   = 0.5d0*(Tl + Tr) 
    rho = 0.5d0*(rhol + rhor)
    mu     = 0.5d0*(mul + mur)
    kappa  = 0.5d0*(kappal + kappar)
    lambda = 0.5d0*(lambdal + lambdar)
 
    ! Diffusive flux and diffusive flux Jacobians 
    ! Common factors
    gamma_minus1 = gamma - 1.d0
    tmp1 = coeff*mu + kappa
    tmp2 = tmp1*ov_dx/rho 
    tmp3 = tmp2*u
    tmp4 = lambda/rho*ov_dx
    tmp5 = tmp4/R_gas*gamma_minus1    

    ! Diffusive flux
    tmp1  = tmp1*du_dx
    fd(pos_rho)  = 0.d0
    fd(pos_rhou) = tmp1
    fd(pos_rhoE) = tmp1*u + lambda*dT_dx

    ! Right state diffusive flux Jacobian 
    ! First column 
    jfdr(pos_rho,pos_rho)  = 0.d0
    jfdr(pos_rhou,pos_rho) = - tmp3    
    jfdr(pos_rhoe,pos_rho) = tmp4*(0.5d0*v2*gamma_minus1/R_gas - T) - tmp2*v2

    ! Second column
    jfdr(pos_rho,pos_rhou)  = 0.d0
    jfdr(pos_rhou,pos_rhou) = tmp2  
    jfdr(pos_rhoE,pos_rhou) = - tmp5*u + tmp3

    ! Third column
    jfdr(pos_rho,pos_rhoE)  = 0.d0
    jfdr(pos_rhou,pos_rhoE) = 0.d0    
    jfdr(pos_rhoE,pos_rhoE) = tmp5

    ! Left state diffusive flux Jacobian
    jfdl = - jfdr

  END SUBROUTINE ns_flux_1D_Jac
!------------------------------------------------------------------------------!
