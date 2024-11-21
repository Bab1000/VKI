!------------------------------------------------------------------------------!
! This module provides pointer definition for the Argon CR library.
  MODULE mod_argon_CR_function_pointer

    IMPLICIT NONE

    ABSTRACT INTERFACE

      ! Interface for subroutine computing species energies and constant volume specific heats
      SUBROUTINE get_energy_cv (temp, e, cv)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e, cv
      END SUBROUTINE get_energy_cv

      ! Interface for subroutine computing species energies and constant volume specific heats
      SUBROUTINE get_energy (temp, e)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e
      END SUBROUTINE get_energy

      ! Interface for subroutine computing the temperatures from the energy
      ! densities 
      SUBROUTINE get_temp (rhoi, rho_e, temp)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, rho_e
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp
      END SUBROUTINE get_temp

      ! Interface for subroutine computing mass production terms  
      SUBROUTINE get_source (rhoi, temp, omega)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: omega
      END SUBROUTINE get_source

    END INTERFACE

    PROCEDURE(get_energy), POINTER, SAVE :: get_species_energy
    PROCEDURE(get_energy) :: energy_CR, energy_MT

    PROCEDURE(get_energy_cv), POINTER, SAVE :: get_species_energy_cv
    PROCEDURE(get_energy_cv) :: energy_cv_CR, energy_cv_MT 

    PROCEDURE(get_temp), POINTER, SAVE :: get_temperatures
    PROCEDURE(get_temp) :: get_temperatures_CR, get_temperature_MT

    PROCEDURE(get_source), POINTER, SAVE :: get_mass_prod_terms
    PROCEDURE(get_source) :: source_CR, source_MT 

  END MODULE mod_argon_CR_function_pointer
!------------------------------------------------------------------------------!
