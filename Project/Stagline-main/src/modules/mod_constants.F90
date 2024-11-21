!------------------------------------------------------------------------------!
!> This module provides some useful physical constants.
  MODULE mod_constants
    REAL(KIND=8), PARAMETER :: ukb = 1.380658d-23   !< Boltzmann's constant [J/K]
    REAL(KIND=8), PARAMETER :: una = 6.0221367d23   !< Avogadros' number [1/mol]
    REAL(KIND=8), PARAMETER :: urg = ukb*una        !< Universal gas constant [J/mol/K]
    REAL(KIND=8), PARAMETER :: ue  = 1.602191d-19   !< Electron charge [C]
    REAL(KIND=8), PARAMETER :: uh  = 6.626075d-34   !< Planck's constant [J*s]
    REAL(KIND=8), PARAMETER :: upi = 3.14159265d0   !< Pi-greek 
  END MODULE mod_constants
!------------------------------------------------------------------------------!
