!------------------------------------------------------------------------------!
!> This subroutine computes the conservative and physical variables from the primitive variables for 
!! 1D stagnation line calorically perfect gas flows when solving the Euler equations
  SUBROUTINE prim_to_cons_phys_1D_SL_Eu (prim, cons, phys_prop) 

    USE mod_general_data,             ONLY: pos_u_cell, pos_v_cell, pos_h0_cell, pos_pres_cell, pos_c_cell, & 
                                          & pos_ek_cell, pos_rho_cell, pos_T_cell, gamma, R_gas

    IMPLICIT NONE

    REAL(KIND=8) :: ov_gamma_minus1 
    REAL(KIND=8) :: rho, p, u, v
    REAL(KIND=8) :: ek, p_ov_rho, ov_rho, rhoE

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prim        !< vector of primitive variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cons       !< vector of conservative variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: phys_prop  !< vector of physical variables

    ! Density, velocity components and pressure
    rho = prim(1)
    u   = prim(2)
    v   = prim(3)
    p   = prim(4)

    ! Useful quantity
    ov_gamma_minus1 = 1.d0/(gamma - 1.d0)

    ek       = 0.5d0*u*u
    ov_rho   = 1.d0/rho
    p_ov_rho = p*ov_rho
    rhoE     = p*ov_gamma_minus1 + rho*ek

    ! Fill the vector of conservative variables
    cons(1) = rho
    cons(2) = rho*u
    cons(3) = rho*v
    cons(4) = rhoE

    ! Fill the vector of physical propeties
    phys_prop(pos_u_cell)    = u
    phys_prop(pos_v_cell)    = v
    phys_prop(pos_h0_cell)   = rhoE*ov_rho + p_ov_rho
    phys_prop(pos_c_cell)    = DSQRT(gamma*p_ov_rho)
    phys_prop(pos_ek_cell)   = ek 
    phys_prop(pos_pres_cell) = p
    phys_prop(pos_rho_cell)  = rho 
    phys_prop(pos_T_cell)    = p_ov_rho/R_gas

  END SUBROUTINE prim_to_cons_phys_1D_SL_Eu
!------------------------------------------------------------------------------!
!> This subroutine computes the conservative and physical variables from the primitive variables for 
!! 1D stagnation line calorically perfect gas flows when solving the Navier-Stokes equations
  SUBROUTINE prim_to_cons_phys_1D_SL_Ns (prim, cons, phys_prop) 

    USE mod_general_data,             ONLY: pos_u_cell, pos_v_cell, pos_h0_cell, pos_pres_cell, pos_c_cell,      &   
                                          & pos_ek_cell, pos_rho_cell, pos_T_cell, pos_mu_cell, pos_lambda_cell, & 
                                          & pos_kappa_cell, R_gas, gamma
    USE mod_function_pointer,         ONLY: get_transpCoeff

    IMPLICIT NONE

    REAL(KIND=8) :: ov_gamma_minus1 
    REAL(KIND=8) :: mu, lambda
    REAL(KIND=8) :: rho, p, u, v
    REAL(KIND=8) :: ek, p_ov_rho, ov_rho, rhoE, T

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prim        !< vector of primitive variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cons       !< vector of conservative variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: phys_prop  !< vector of physical properties

    ! Density velocity components and pressure
    rho = prim(1)
    u   = prim(2)
    v   = prim(3)
    p   = prim(4)

    ! Useful quantity
    ov_gamma_minus1 = 1.d0/(gamma - 1.d0)

    ek       = 0.5d0*u*u
    ov_rho   = 1.d0/rho
    p_ov_rho = p*ov_rho
    rhoE     = p*ov_gamma_minus1 + rho*ek
    T = p_ov_rho/R_gas

    ! Dynamic viscosity and thermal conductivity 
    CALL get_transpCoeff (T, mu, lambda)

    ! Fill the vector of conservative variables
    cons(1) = rho
    cons(2) = rho*u
    cons(3) = rho*v
    cons(4) = rhoE

    ! Fill the vector of physical propeties
    phys_prop(pos_u_cell)    = u
    phys_prop(pos_v_cell)    = v
    phys_prop(pos_h0_cell)   = rhoE*ov_rho + p_ov_rho
    phys_prop(pos_c_cell)    = DSQRT(gamma*p_ov_rho)
    phys_prop(pos_ek_cell)   = ek 
    phys_prop(pos_pres_cell) = p
    phys_prop(pos_rho_cell)  = rho 
    phys_prop(pos_T_cell)    = T
    phys_prop(pos_mu_cell)     = mu
    phys_prop(pos_kappa_cell)  = 0.d0 
    phys_prop(pos_lambda_cell) = lambda 

  END SUBROUTINE prim_to_cons_phys_1D_SL_Ns
!------------------------------------------------------------------------------!
