!------------------------------------------------------------------------------!
! This module provides pointer definition for the Nitrogen Park library.
  MODULE mod_function_pointer_Park

    ! Interface for the subroutine enabling temperature computation  
    ABSTRACT INTERFACE 
      SUBROUTINE get_temperatures (rhoi, rho_eint, temp)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, rho_eint
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp
      END SUBROUTINE get_temperatures 
    END INTERFACE 

    ! Interface for the subroutine computing the frozen specific heats
    ABSTRACT INTERFACE 
      SUBROUTINE get_frozen_cv (T, cv)
        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv
      END SUBROUTINE get_frozen_cv
    END INTERFACE

    ! Interface for subroutine computing frozen specific heats and energy species energy components in 
    ! thermal equilibrium with translation
    ABSTRACT INTERFACE 
      SUBROUTINE get_frozen_energy_cv (T, cv, e)
        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv, e
      END SUBROUTINE get_frozen_energy_cv
    END INTERFACE

    ! Interface for subroutine computing the source term vector
    ABSTRACT INTERFACE 
      SUBROUTINE get_source (rhoi, temp, omega)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: omega
      END SUBROUTINE get_source
    END INTERFACE

    ! Interface for subroutine computing the source term vector and the related
    ! Jacobian
    ABSTRACT INTERFACE 
      SUBROUTINE get_source_Jac (rhoi, temp, omega, js_omega)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: omega, js_omega
      END SUBROUTINE get_source_Jac
    END INTERFACE

    PROCEDURE(get_temperatures), POINTER, SAVE :: get_temp_Park
    PROCEDURE(get_temperatures) :: compute_T, compute_T_Tv, compute_T_Tr_Tv

    PROCEDURE(get_frozen_cv), POINTER, SAVE :: get_species_cv_Park
    PROCEDURE(get_frozen_cv) :: compute_frozen_cv_T, compute_frozen_cv_T_Tv, compute_frozen_cv_T_Tr_Tv   

    PROCEDURE(get_frozen_energy_cv), POINTER, SAVE :: get_species_energy_cv_Park
    PROCEDURE(get_frozen_energy_cv) :: compute_frozen_energy_cv_T, compute_frozen_energy_cv_T_Tv,  & 
                                     & compute_frozen_energy_cv_T_Tr_Tv 

    PROCEDURE(get_source), POINTER, SAVE :: get_source_term
    PROCEDURE(get_source) :: source, VT_RT_transfer     

    PROCEDURE(get_source_Jac), POINTER, SAVE :: get_source_term_Jac
    PROCEDURE(get_source_Jac) :: source_Jac, VT_RT_transf_Jac

  END MODULE mod_function_pointer_Park
!-------------------------------------------------------------------------------!
