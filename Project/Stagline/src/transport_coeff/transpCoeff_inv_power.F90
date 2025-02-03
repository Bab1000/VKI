!------------------------------------------------------------------------------!
! This subroutine computes the dynamic viscosity for an inverse power interaction potential. 
! The thermal conductivity is then computed according to Kinetic Theory of gases
! (this subroutine shuold be used only for mono-atomic gases).
  SUBROUTINE transpCoeff_inv_power (T, mu, lambda)

    USE mod_general_data,            ONLY: T_ref, mu_ref, exp_visc, R_gas

    IMPLICIT NONE

    REAL(KIND=8), PARAMETER :: fac = 15.d0/4.d0

    REAL(KIND=8), INTENT(IN) :: T
    REAL(KIND=8), INTENT(OUT) :: mu, lambda

    mu     = mu_ref*(T/T_ref)**exp_visc
    lambda = fac*mu*R_gas

  END SUBROUTINE transpCoeff_inv_power
!------------------------------------------------------------------------------!
