!------------------------------------------------------------------------------!
!> Thus subroutine provides the inviscid flux Jacobian spectral radius for the 1D Euler equations.
  SUBROUTINE inv_spectral_radius_1D (phys_prop, sp)

    USE mod_general_data,      ONLY: pos_u_cell, pos_c_cell

    IMPLICIT NONE

    REAL(KIND=8) :: u, c

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_prop  !< physical properties
    REAL(KIND=8), INTENT(OUT) :: sp                      !< spectral radius

    ! Velocity and speed of sound 
    u = phys_prop(pos_u_cell)
    c = phys_prop(pos_c_cell)

    ! Inviscid flux Jacobian spectral radius
    sp = ABS(u) + c

  END SUBROUTINE inv_spectral_radius_1D
!------------------------------------------------------------------------------!
!> This subroutine sets the diffusive flux Jacobian spectral radius to zero when
!! solving the 1D Euler equations.
  SUBROUTINE null_visc_spectral_radius_1D (phys_prop, sp)

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_prop  !< physical properties
    REAL(KIND=8), INTENT(OUT) :: sp                      !< spectral radius

    sp = 0.d0

  END SUBROUTINE null_visc_spectral_radius_1D
!------------------------------------------------------------------------------!
!> This subroutine provides the diffusive flux Jacobian spectral radius for the 1D Navier-Stokes equations
!! for a calorically perfect gas.
  SUBROUTINE visc_spectral_radius_1D (phys_prop, sp)

    USE mod_general_data,      ONLY: gamma, R_gas, pos_mu_cell, pos_rho_cell, pos_lambda_cell 

    IMPLICIT NONE

    REAL(KIND=8), PARAMETER :: coeff = 4.d0/3.d0
    REAL(KIND=8) :: rho, mu, lambda
    REAL(KIND=8) :: l1, l2

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_prop  !< physical properties
    REAL(KIND=8), INTENT(OUT) :: sp                      !< spectral radius

    ! Density, dynamic viscosity and thermal conductivity 
    rho    = phys_prop(pos_rho_cell)
    mu     = phys_prop(pos_mu_cell) 
    lambda = phys_prop(pos_lambda_cell) 

    ! Viscous eigenvalues associated to the global energy and momentum equations
    l1 = coeff*mu/rho
    l2 = lambda*(gamma - 1.d0)/(rho*R_gas)

    ! Spectral radius (the maximum between the two eigenvalues is taken) 
    sp = MAX(l1,l2)

  END SUBROUTINE visc_spectral_radius_1D
!------------------------------------------------------------------------------!
!> This subroutine provides the diffusive flux Jacobian spectral radius for the 1D Navier-Stokes equations
!! for a nonequilibrium flow.
  SUBROUTINE visc_spectral_radius_neq_1D (phys_prop, sp)

    USE mod_general_data,      ONLY: pos_mu_cell, pos_lambda_cell, pos_rho_cell, pos_beta_cell

    IMPLICIT NONE

    REAL(KIND=8), PARAMETER :: coeff = 4.d0/3.d0
    REAL(KIND=8) :: rho, mu, lambda, beta
    REAL(KIND=8) :: l1, l2
   
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_prop  !< physical properties
    REAL(KIND=8), INTENT(OUT) :: sp                      !< spectral radius

    ! Density, dynamic viscosity, thermal conducticity and beta factor  
    rho    = phys_prop(pos_rho_cell)
    mu     = phys_prop(pos_mu_cell) 
    lambda = phys_prop(pos_lambda_cell) 
    beta   = phys_prop(pos_beta_cell)

    ! Viscous eigenvalues associated to the global energy and momentum equations
    l1 = coeff*mu/rho
    l2 = lambda/beta

    ! Spectral radius (the maximum between the two eigenvalues is taken) 
    sp = MAX(l1,l2)

  END SUBROUTINE visc_spectral_radius_neq_1D
!------------------------------------------------------------------------------!

