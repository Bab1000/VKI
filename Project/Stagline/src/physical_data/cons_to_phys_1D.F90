!------------------------------------------------------------------------------!
!> This subroutine computes the physical properties from conservative variables for 1D calorically perfect gas flows
!! when solving the Euler equations. 
  SUBROUTINE cons_to_phys_1D_Eu (cons, phys_prop)

    USE mod_general_data,        ONLY: pos_u_cell, pos_h0_cell, pos_c_cell, pos_ek_cell, pos_pres_cell,  & 
                                     & pos_rho_cell, pos_T_cell, gamma, R_gas

    IMPLICIT NONE

    REAL(KIND=8) :: rho, rhoE
    REAL(KIND=8) :: ov_rho, rhoE_ov_rho, p_ov_rho
    REAL(KIND=8) :: c, ek, p, u

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons        !< vector of conservative variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: phys_prop  !< vector of physical properties 

    ! Density
    rho = cons(1)
    ov_rho = 1.d0/rho

    ! Velocity and kinetic energy per unit mass
    u  = cons(2)*ov_rho 
    ek = 0.5d0*u**2

    ! Pressure and speed of sound
    rhoE = cons(3) 
    rhoE_ov_rho = rhoE*ov_rho

    p  = (gamma - 1.d0)*(rhoE - rho*ek)
    p_ov_rho = p*ov_rho
    c        = DSQRT(gamma*p_ov_rho)

    ! Fill vector of physical properties 
    phys_prop(pos_u_cell)    = u 
    phys_prop(pos_h0_cell)   = rhoE_ov_rho + p_ov_rho
    phys_prop(pos_c_cell)    = c
    phys_prop(pos_ek_cell)   = ek
    phys_prop(pos_pres_cell) = p
    phys_prop(pos_rho_cell)  = rho
    phys_prop(pos_T_cell)    = p_ov_rho/R_gas
    
  END SUBROUTINE cons_to_phys_1D_Eu
!------------------------------------------------------------------------------!
!> This subroutine computes the physical properties from conservative variables for 1D calorically perfect gas flows
!! when solving the Navier-Stokes equations. 
  SUBROUTINE cons_to_phys_1D_Ns (cons, phys_prop)

    USE mod_general_data,        ONLY: pos_u_cell, pos_h0_cell, pos_c_cell, pos_ek_cell, pos_pres_cell,  & 
                                     & pos_rho_cell, pos_T_cell, pos_mu_cell, pos_kappa_cell,            & 
                                     & pos_lambda_cell, R_gas, gamma
    USE mod_function_pointer,    ONLY: get_transpCoeff

    IMPLICIT NONE

    REAL(KIND=8) :: rho, rhoE
    REAL(KIND=8) :: ov_rho, rhoE_ov_rho, p_ov_rho
    REAL(KIND=8) :: c, lambda, ek, mu, p, T, u

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons         !< vector of conservative variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: phys_prop   !< vector of physical properties

    ! Density
    rho = cons(1)
    ov_rho = 1.d0/rho

    ! Velocity and kinetic energy per unit mass
    u  = cons(2)*ov_rho 
    ek = 0.5d0*u**2

    ! Pressure, temperature and speed of sound
    rhoE = cons(3) 
    rhoE_ov_rho = rhoE*ov_rho

    p  = (gamma - 1.d0)*(rhoE - rho*ek)
    p_ov_rho = p*ov_rho
    T        = p_ov_rho/R_gas
    c        = DSQRT(gamma*p_ov_rho)

    ! Compute transport properties
    CALL get_transpCoeff (T, mu, lambda)

    ! Fill vector of physical properties 
    phys_prop(pos_u_cell)      = u 
    phys_prop(pos_h0_cell)     = rhoE_ov_rho + p_ov_rho
    phys_prop(pos_c_cell)      = c
    phys_prop(pos_ek_cell)     = ek
    phys_prop(pos_pres_cell)   = p
    phys_prop(pos_rho_cell)    = rho
    phys_prop(pos_T_cell)      = T
    phys_prop(pos_mu_cell)     = mu
    phys_prop(pos_kappa_cell)  = 0.d0
    phys_prop(pos_lambda_cell) = lambda
 
  END SUBROUTINE cons_to_phys_1D_Ns
!------------------------------------------------------------------------------!
