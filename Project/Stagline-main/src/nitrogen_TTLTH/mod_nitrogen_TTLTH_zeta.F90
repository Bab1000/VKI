!------------------------------------------------------------------------------!
! This module contains the subroutine initializing pointers.
  MODULE mod_nitrogen_TTLTH_zeta

    USE mod_nitrogen_TTLTH_initialize_CFD,        ONLY: nb_tvib, nb_trot, solver 
    USE mod_nitrogen_TTLTH_CFD_prop 
    USE mod_function_pointer_TTLTH,               ONLY: get_temp_TTLTH, get_species_cv_TTLTH, get_species_energy_cv_TTLTH

    IMPLICIT NONE

    CONTAINS

      !----------------------------------------------------!
      ! Subroutine for pointer initialization
      SUBROUTINE set_pointer ()

        INTEGER :: nb_temp

        nb_temp = 1 + nb_trot + nb_tvib

        ! Subroutine selection
        SELECT CASE (nb_temp)

          ! Thermal equilibrium
          CASE(1) 
            get_temp_TTLTH => compute_T
            get_species_cv_TTLTH => compute_frozen_cv_T
            get_species_energy_cv_TTLTH => compute_frozen_energy_cv_T

          ! Thermal nonequilibrium
          CASE DEFAULT 

            IF (nb_tvib.GE.1.AND.nb_trot.EQ.0) THEN

               get_temp_TTLTH => compute_T_Tv
               get_species_cv_TTLTH => compute_frozen_cv_T_Tv
               get_species_energy_cv_TTLTH => compute_frozen_energy_cv_T_Tv

            ELSEIF (nb_tvib.EQ.1.AND.nb_trot.EQ.1) THEN

               get_temp_TTLTH => compute_T_Tr_Tv
               get_species_cv_TTLTH => compute_frozen_cv_T_Tr_Tv
               get_species_energy_cv_TTLTH => compute_frozen_energy_cv_T_Tr_Tv 

            ELSE 

              PRINT*
              WRITE(*,10)solver(1:LEN_TRIM(solver)),':: in mod_nitrogen_TTLTH_zeta, error in temperature selection...'
              PRINT*
              STOP

            ENDIF

        END SELECT 

10    FORMAT(A,A)

      END SUBROUTINE set_pointer

  END MODULE mod_nitrogen_TTLTH_zeta
!------------------------------------------------------------------------------!
