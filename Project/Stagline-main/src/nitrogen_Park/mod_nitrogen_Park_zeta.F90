!------------------------------------------------------------------------------!
! This module contains the subroutine initializing pointers.
  MODULE mod_nitrogen_Park_zeta

    USE mod_nitrogen_Park_initialize_CFD,        ONLY: nb_ns, nb_tvib, nb_trot, solver 
    USE mod_nitrogen_Park_CFD_prop
    USE mod_nitrogen_Park_CFD_source
    USE mod_nitrogen_Park_CFD_source_VT_RT 
    USE mod_function_pointer_Park,               ONLY: get_temp_Park, get_species_cv_Park, get_species_energy_cv_Park, & 
                                                    &  get_source_term, get_source_term_Jac

    IMPLICIT NONE

    CONTAINS

      !----------------------------------------------------!
      ! Subroutine for pointer initialization
      SUBROUTINE set_pointer ()

        INTEGER :: nb_temp

        nb_temp = 1 + nb_trot + nb_tvib

        ! Subroutine selection for source term and related Jacobian
        IF (nb_ns.EQ.1) THEN

          get_source_term => VT_RT_transf
          get_source_term_Jac => VT_RT_transf_Jac

        ELSE

          get_source_term => source
          get_source_term_Jac => source_Jac

        ENDIF

        ! Subroutine selection for thermodynamics
        SELECT CASE (nb_temp)

          ! Thermal equilibrium
          CASE(1) 
            get_temp_Park => compute_T
            get_species_cv_Park => compute_frozen_cv_T
            get_species_energy_cv_Park => compute_frozen_energy_cv_T

          ! Thermal nonequilibrium
          CASE DEFAULT 

            IF (nb_tvib.EQ.1.AND.nb_trot.EQ.0) THEN

               get_temp_Park => compute_T_Tv
               get_species_cv_Park => compute_frozen_cv_T_Tv
               get_species_energy_cv_Park => compute_frozen_energy_cv_T_Tv

            ELSEIF (nb_tvib.EQ.1.AND.nb_trot.EQ.1) THEN

               get_temp_Park => compute_T_Tr_Tv
               get_species_cv_Park => compute_frozen_cv_T_Tr_Tv
               get_species_energy_cv_Park => compute_frozen_energy_cv_T_Tr_Tv 

            ELSE 

              PRINT*
              WRITE(*,10)solver(1:LEN_TRIM(solver)),':: in mod_nitrogen_Park_zeta, error in temperature selection...'
              PRINT*
              STOP

            ENDIF

        END SELECT 

10    FORMAT(A,A)

      END SUBROUTINE set_pointer

  END MODULE mod_nitrogen_Park_zeta
!------------------------------------------------------------------------------!
