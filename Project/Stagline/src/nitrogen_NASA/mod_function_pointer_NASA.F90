!------------------------------------------------------------------------------!
! This module provides pointer definition for the Nitrogen NASA library.
! It also provides some useful fitting functions for vTm and vv processes (N2 + N2 system)
  MODULE mod_function_pointer_NASA

    IMPLICIT NONE

    ! Interface for the subroutine enabling temperature computation  
    ABSTRACT INTERFACE 
      SUBROUTINE get_temperatures (rhoi, rho_eint, temp)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, rho_eint
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp
      END SUBROUTINE get_temperatures 

      ! Interface for the subroutine computing the frozen specific heat
      SUBROUTINE get_frozen_cv (T, cv)
        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv
      END SUBROUTINE get_frozen_cv

      ! Interface for subroutine computing species energy 
      SUBROUTINE get_energy (temp, e)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e
      END SUBROUTINE get_energy

      ! Interface for subroutine computing species energy and constant volume specific heat
      SUBROUTINE get_energy_cv (T, e, cv)
        REAL(KIND=8), INTENT(IN) :: T 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e, cv
      END SUBROUTINE get_energy_cv

      ! Interface for subroutine computing mass production terms  
      SUBROUTINE get_source (rhoi, temp, omega)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: omega
      END SUBROUTINE get_source

      ! Interface for subroutine computing mass production terms  
      SUBROUTINE get_source_Jac (rhoi, temp, omega, js_omega)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: omega, js_omega
      END SUBROUTINE get_source_Jac

      ! Interface for subrutine computing the internal temperature
      SUBROUTINE get_internal_temperature (ni, temp_in, temp_out) 
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp_in, ni 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp_out
      END SUBROUTINE get_internal_temperature

      ! Interface for subroutine computing the internal component of thermal conductivity
      SUBROUTINE get_internal_therm_cond(nb, T, xi, xig, Dij, lambda_int)
        REAL(KIND=8), INTENT(IN) :: nb, T
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, xig, Dij
        REAL(KIND=8), INTENT(OUT) :: lambda_int
      END SUBROUTINE get_internal_therm_cond

    END INTERFACE

    PROCEDURE(get_temperatures), POINTER, SAVE :: get_temp_NASA
    PROCEDURE(get_temperatures) :: compute_T_VC, compute_T_VC_rr, compute_T_BRVC, compute_T_RVC, & 
                                 & compute_temp_MT_TTint

    PROCEDURE(get_frozen_cv), POINTER, SAVE :: get_species_cv_NASA
    PROCEDURE(get_frozen_cv) :: compute_frozen_cv_VC, compute_frozen_cv_VC_rr, compute_frozen_cv_BRVC, compute_frozen_cv_RVC, & 
                              & compute_frozen_cv_MT_TTint  

    PROCEDURE(get_energy), POINTER, SAVE :: get_species_energy_NASA 
    PROCEDURE(get_energy) :: energy_VC, energy_VC_rr, energy_BRVC, energy_RVC, energy_MT_TTint

    PROCEDURE(get_energy_cv), POINTER, SAVE :: get_species_energy_cv_NASA 
    PROCEDURE(get_energy_cv) :: energy_cv_VC, energy_cv_VC_rr, energy_cv_BRVC, energy_cv_RVC, energy_cv_MT_TTint

    PROCEDURE(get_source), POINTER, SAVE :: get_mass_prod_terms_NASA
    PROCEDURE(get_source), POINTER :: source_VC, source_VC_rr, source_BRVC, source_RVC, source_MT_TTint

    PROCEDURE(get_source_Jac), POINTER, SAVE :: get_mass_prod_terms_Jac_NASA
    PROCEDURE(get_source_Jac), POINTER, SAVE :: source_VC_Jac, source_VC_rr_Jac, source_BRVC_Jac, source_RVC_Jac, & 
                                              & source_MT_TTint_Jac

    PROCEDURE(get_internal_temperature), POINTER, SAVE :: get_Tint_NASA 
    PROCEDURE(get_internal_temperature) :: compute_Tint_VC_rr, compute_Tint_NonUnif

    PROCEDURE(get_internal_therm_cond), POINTER, SAVE :: get_lambda_int
    PROCEDURE(get_internal_therm_cond) :: lambda_int_VC_rr, lambda_int_BRVC

  END MODULE mod_function_pointer_NASA
!------------------------------------------------------------------------------!
