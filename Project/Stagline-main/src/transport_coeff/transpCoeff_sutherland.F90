!------------------------------------------------------------------------------!
! This subroutine computes the dynamic viscosity according to Sutherlans's law. 
! The thermal conductivity is computed by assuming a constant Prandtl number.
  SUBROUTINE transpCoeff_sutherland (T, mu, lambda)

     USE mod_general_data,            ONLY: T_ref, mu_ref, cp_gas, R_gas, Pr

     IMPLICIT NONE

     REAL(KIND=8), INTENT(IN) :: T
     REAL(KIND=8), INTENT(OUT) :: mu, lambda

     mu     = mu_ref*(T**1.5/(T + T_ref))
     lambda = mu*cp_gas/Pr

  END SUBROUTINE transpCoeff_sutherland
!------------------------------------------------------------------------------!

