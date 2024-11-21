!------------------------------------------------------------------------------!
! This module contains the subroutine initializing pointers.
  MODULE mod_nitrogen_NASA_zeta

    USE mod_nitrogen_NASA_initialize_CFD,       ONLY: model, solver 
    USE mod_nitrogen_NASA_CFD_prop
    USE mod_nitrogen_NASA_CFD_source_VC
    USE mod_nitrogen_NASA_CFD_source_VC_rr
    USE mod_nitrogen_NASA_CFD_source_BRVC
    USE mod_nitrogen_NASA_CFD_source_RVC
    USE mod_nitrogen_NASA_CFD_source_MT_TTint 
    USE mod_nitrogen_NASA_CFD_transport
    USE mod_function_pointer_NASA,              ONLY: get_temp_NASA, get_species_cv_NASA, get_species_energy_NASA, & 
                                                   &  get_species_energy_cv_NASA, get_Tint_NASA,                   & 
                                                   &  get_mass_prod_terms_NASA, get_mass_prod_terms_Jac_NASA,      & 
                                                   &  get_lambda_int

    IMPLICIT NONE

    CONTAINS 

      !----------------------------------------------------!
      ! Subroutine for pointer initialization
      SUBROUTINE set_pointer ()

        INTEGER :: length

        length = LEN_TRIM(solver)

        ! Subroutine selection
        SELECT CASE(model)

          CASE('RVC') 
            WRITE(*,20)solver(1:length),':: Nitrogen NASA library -> pointer initialization for RVC model'
            PRINT*

            get_temp_NASA => compute_T_RVC
            get_species_cv_NASA => compute_frozen_cv_RVC
            get_species_energy_NASA => energy_RVC
            get_species_energy_cv_NASA => energy_cv_RVC
            get_mass_prod_terms_NASA => source_RVC 
            get_mass_prod_terms_Jac_NASA => source_RVC_Jac
            get_Tint_NASA => compute_Tint_NonUnif 
            get_lambda_int => lambda_int_RVC
 
           CASE('BRVC')
            WRITE(*,20)solver(1:length),':: Nitrogen NASA library -> pointer initialization for BRVC model'
            PRINT*

            get_temp_NASA => compute_T_BRVC
            get_species_cv_NASA => compute_frozen_cv_BRVC
            get_species_energy_NASA => energy_BRVC
            get_species_energy_cv_NASA => energy_cv_BRVC
            get_mass_prod_terms_NASA => source_BRVC
            get_mass_prod_terms_Jac_NASA => source_BRVC_Jac
            get_Tint_NASA => compute_Tint_NonUnif
            get_lambda_int => lambda_int_BRVC

          CASE('VC')
            WRITE(*,20)solver(1:length),':: Nitrogen NASA library -> pointer initialization for VC model'
            PRINT*

            get_temp_NASA => compute_T_VC
            get_species_cv_NASA => compute_frozen_cv_VC
            get_species_energy_NASA => energy_VC
            get_species_energy_cv_NASA => energy_cv_VC
            get_mass_prod_terms_NASA => source_VC
            get_mass_prod_terms_Jac_NASA => source_VC_Jac
            get_Tint_NASA => compute_Tint_NonUnif
            get_lambda_int => lambda_int_VC

          CASE('VC_rr')
            WRITE(*,20)solver(1:length),':: Nitrogen NASA library -> pointer initialization for VC_rr model' 
            PRINT*

            get_temp_NASA => compute_T_VC_rr
            get_species_cv_NASA => compute_frozen_cv_VC_rr
            get_species_energy_NASA => energy_VC_rr
            get_species_energy_cv_NASA => energy_cv_VC_rr
            get_mass_prod_terms_NASA => source_VC_rr
            get_mass_prod_terms_Jac_NASA => source_VC_rr_Jac
            get_Tint_NASA => compute_Tint_VC_rr
            get_lambda_int => lambda_int_VC_rr

          CASE('MT_TTint')
            WRITE(*,20)solver(1:length),':: Nitrogen NASA library -> pointer initialization for MT_TTint model' 
            PRINT*

            get_temp_NASA => compute_temp_MT_TTint
            get_species_cv_NASA => compute_frozen_cv_MT_TTint
            get_species_energy_NASA => energy_RVC
            get_species_energy_cv_NASA => energy_cv_MT_TTint
            get_mass_prod_terms_NASA => source_MT_TTint 
            get_mass_prod_terms_Jac_NASA => source_MT_TTint_Jac
            get_lambda_int => lambda_int_MT_TTint

          CASE('EQ')
            WRITE(*,20)solver(1:length),':: Nitrogen NASA library -> equilibrium- no pointer initialization'
            PRINT*

          CASE DEFAULT 
            PRINT*
            WRITE(*,20)solver(1:length),':: error in model selection, in mod_nitrogen_NASA_zeta ..'
            PRINT*
            STOP

        END SELECT

20     FORMAT(A,A)

      END SUBROUTINE set_pointer

  END MODULE mod_nitrogen_NASA_zeta
!------------------------------------------------------------------------------!
