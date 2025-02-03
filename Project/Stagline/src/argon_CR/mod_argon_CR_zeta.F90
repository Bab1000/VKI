!------------------------------------------------------------------------------!
! This module contains the subroutine initializing pointers.
  MODULE mod_argon_CR_zeta

    USE mod_argon_CR_initialize_CFD,       ONLY: model, solver
    USE mod_argon_CR_CFD_prop
    USE mod_argon_CR_CFD_source_CR
    USE mod_argon_CR_CFD_source_MT
    USE mod_argon_CR_function_pointer,     ONLY: get_species_energy, get_temperatures,      &  
                                               & get_species_energy_cv, get_mass_prod_terms

    IMPLICIT NONE

    CONTAINS 

      !----------------------------------------------------!
      ! Subroutine for pointer initialization
      SUBROUTINE set_pointer 

        INTEGER :: length

        length = LEN_TRIM(solver)

        ! Subroutine selection
        SELECT CASE(model)

          CASE('CR') 
            WRITE(*,20)solver(1:length),':: Argon CR library -> pointer initialization for CR model'
            PRINT*

            get_mass_prod_terms => source_CR
            get_temperatures => get_temperatures_CR
            get_species_energy => energy_CR
            get_species_energy_cv => energy_cv_CR

          CASE('MT')
            WRITE(*,20)solver(1:length),':: Argon CR library -> pointer initialization for MT model'
            PRINT*

            get_mass_prod_terms => source_MT
            get_temperatures => get_temperatures_MT
            get_species_energy => energy_MT 
            get_species_energy_cv => energy_cv_MT

          CASE DEFAULT 
            PRINT*
            WRITE(*,20)solver(1:length),':: error in model selection, in mod_argon_CR_zeta ..'
            PRINT*
            STOP

         END SELECT            

20     FORMAT(A,A)

      END SUBROUTINE set_pointer 

  END MODULE mod_argon_CR_zeta
!------------------------------------------------------------------------------!
