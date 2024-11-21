!------------------------------------------------------------------------------!
!> This subroutine associates pointers to subroutines of the library interface 
  SUBROUTINE set_library_pointer ()

#include "../config.h"

    USE mod_general_data,              ONLY: library_name, flag_post_CR
    USE mod_neq_function_pointer,      ONLY: library_initialize, library_finalize, library_get_Ri,  &
                                           & library_get_el_data,                                   &
                                           & library_get_density, library_get_pressure,             & 
                                           & library_get_mass_fractions,                            &
                                           & library_get_molar_fractions,                           &
                                           & library_get_mass_fractions_from_molar_fractions,       &   
                                           & library_comp_tol,                                      &                       
                                           & library_compute_eq_composition, library_get_source,    & 
                                           & library_compute_eq_composition_pyro,                   & 
                                           & library_get_source_Jac, library_get_temperatures,      & 
                                           & library_get_prandtl_number,library_get_energy_densities,&
                                           & library_get_data,                                      &   
                                           & library_get_thermodynamic_data,                        & 
                                           & library_get_transpCoeff, library_get_species_DiffFlux, &  
                                           & library_get_thermal_cond,                              &
                                           & library_get_species_enthalpies,                        &  
                                           & library_get_enthalpy,                                  &  
                                           & library_get_species_index,                             &  
                                           & library_get_species_name,                              &  
                                           & library_write_fvmcc_solution,                          &
                                           & library_species_name, library_number_of_elements,      &
                                           & library_compute_eq_composition_meteor,                 &
                                           & library_compute_elemental_mole_fraction,               &
                                           & library_get_molar_fractions_from_mass_fractions,       &  
                                           & library_surface_mass_balance_constraint,               &
                                           & library_get_element_mass_fractions,                    &
                                           & library_set_diff_model,                                &
                                           & library_set_cond_heat_flux,                            &
                                           & library_set_wall_radiation,                            &
                                           & library_set_wall_state,                                &
                                           & library_solve_surface_balance,                         &
                                           & library_get_wall_state,                                &
                                           & library_get_surface_production_rates,                  &
                                           & library_get_mass_blowing_rate,                         &
                                           & library_get_frozen_sound_speed

    IMPLICIT NONE

    ! Library selection 
    SELECT CASE(library_name)

      ! Nitrogen_Park library
#ifdef NITROGEN_PARK
      CASE ('Nitrogen_Park')
        CALL set_pointer_Nitrogen_Park ()
#endif

      ! Nitrogen_DSMC library
#ifdef NITROGEN_DSMC
      CASE ('Nitrogen_DSMC')
        CALL set_pointer_Nitrogen_DSMC ()
#endif

      ! Nitrogen_NASA library 
#ifdef NITROGEN_NASA
      CASE ('Nitrogen_NASA')
        CALL set_pointer_Nitrogen_NASA ()
#endif

      ! Nitrogen_Bari library 
#ifdef NITROGEN_BARI
      CASE ('Nitrogen_Bari')
        CALL set_pointer_Nitrogen_Bari ()
#endif

      ! Nitrogen_FHO library 
#ifdef NITROGEN_FHO
      CASE ('Nitrogen_FHO')
        CALL set_pointer_Nitrogen_FHO ()
#endif

      ! Nitrogen_TTLTH library 
#ifdef NITROGEN_TTLTH
      CASE ('Nitrogen_TTLTH')
        CALL set_pointer_Nitrogen_TTLTH ()
#endif

      ! Argon_CR library 
#ifdef ARGON_CR
      CASE ('Argon_CR')
        CALL set_pointer_Argon_CR ()
#endif

#ifdef MUTATIONPP
      ! Mutation++ library
      CASE('Mutation++','Mutationpp')
        CALL set_pointer_Mutationpp ()
#endif

#ifdef POLYAT_GAS_WCU
      ! Polyat_gas_WCU library
      CASE('Polyat_gas_WCU','polyat_gas_WCU')
        CALL set_pointer_Polyat_gas_WCU ()
#endif
      ! Library not chosen
      CASE('none')
        WRITE(*,10)'In "set_library_pointer.F90", thermodynamic library not correctly set...'
        PRINT*
        STOP

      CASE DEFAULT 
        WRITE(*,10)'In "set_library_pointer.F90", thermodynamic library not implemented yet...'
        WRITE(*,10)library_name
        PRINT*
        STOP

     END SELECT 

10 FORMAT(A)

    CONTAINS

      !----------------------------------------------------!
      !> This is a dummy subroutine 
      SUBROUTINE set_none_library()

      END SUBROUTINE set_none_library

      !----------------------------------------------------!
      !> This subroutine associates pointers to the interface of the Nitrogen_Park thermodynamic library.
#ifdef NITROGEN_PARK
      SUBROUTINE set_pointer_Nitrogen_Park ()

        USE mod_Nitrogen_Park_library

        IMPLICIT NONE

        ! Library inizialization
        library_initialize => Nitrogen_Park_initialize   

        ! Library closing 
        library_finalize => Nitrogen_Park_finalize 

        ! Specific gas constants
        library_get_Ri => Nitrogen_Park_get_Ri

        ! Subroutine for mixture density 
        library_get_density => Nitrogen_Park_get_density

        ! Subroutine for mixture pressure 
        library_get_pressure => Nitrogen_Park_get_pressure

        ! Subroutine for species mass fractions
        library_get_mass_fractions => Nitrogen_Park_get_mass_fractions

        ! Subroutine for species molar fractions
        library_get_molar_fractions => Nitrogen_Park_get_molar_fractions  

        ! Subroutine for computing mass fractions and mixture molar mass from 
        ! molar fractions 
        library_get_mass_fractions_from_molar_fractions => Nitrogen_Park_mass_frac_from_molar_frac  

        ! Subroutine for fix application to species molar fractions
        library_comp_tol => Nitrogen_Park_comp_tol

        ! Subroutine for equilibrium composition computation
        library_compute_eq_composition => Nitrogen_Park_compute_eq_composition

        ! Subroutine for source term 
        library_get_source => Nitrogen_Park_get_source 

        ! Subroutine for source term and source term Jacobian
        library_get_source_Jac => Nitrogen_Park_get_source_Jac 

        ! Subroutine for temperatures 
        library_get_temperatures => Nitrogen_Park_get_temperatures 

        ! Subroutine for computing energy densities
        library_get_energy_densities => Nitrogen_Park_get_energy_densities

        ! Subroutines for computing thermodynamic data
        library_get_thermodynamic_data => Nitrogen_Park_get_thermodynamic_data

        ! Subroutine for computing temperatures and other data
        library_get_data => Nitrogen_Park_get_data

        ! Subroutine for computing transport properties and species mass diffusion fluxes
        library_get_transpCoeff => Nitrogen_Park_compute_transpCoeff
        library_get_species_DiffFlux => Nitrogen_Park_compute_species_DiffFlux

        ! Flag for post-processing
        flag_post_CR = .TRUE.

        ! Subroutine for post-processing
        library_write_fvmcc_solution => Nitrogen_Park_write_fvmcc_solution 

      END SUBROUTINE set_pointer_Nitrogen_Park 
#endif

      !----------------------------------------------------!
      !> This subroutine associates pointers to the interface of the Nitrogen_DSMC thermodynamic library. 
#ifdef NITROGEN_DSMC
      SUBROUTINE set_pointer_Nitrogen_DSMC ()

        USE mod_Nitrogen_DSMC_library

        IMPLICIT NONE

        ! Library inizialization
        library_initialize => Nitrogen_DSMC_initialize   

        ! Library closing 
        library_finalize => Nitrogen_DSMC_finalize 

        ! Specific gas constants
        library_get_Ri => Nitrogen_DSMC_get_Ri

        ! Subroutine for mixture density 
        library_get_density => Nitrogen_DSMC_get_density

        ! Subroutine for mixture pressure 
        library_get_pressure => Nitrogen_DSMC_get_pressure

        ! Subroutine for species mass fractions
        library_get_mass_fractions => Nitrogen_DSMC_get_mass_fractions

        ! Subroutine for species maolar fractions
        library_get_molar_fractions => Nitrogen_DSMC_get_molar_fractions 

        ! Subroutine for species maolar fractions
        library_get_molar_fractions => Nitrogen_DSMC_get_molar_fractions

        ! Subroutine for computing mass fractions and mixture molar mass from 
        ! molar fractions 
        library_get_mass_fractions_from_molar_fractions => Nitrogen_DSMC_mass_frac_from_molar_frac   

        ! Subroutine for fix application to species molar fractions
        library_comp_tol => Nitrogen_DSMC_comp_tol

        ! Subroutine for equilibrium composition computation
        library_compute_eq_composition => Nitrogen_DSMC_compute_eq_composition

        ! Subroutine for source term 
        library_get_source => Nitrogen_DSMC_get_source 

        ! Subroutine for source term and source term Jacobian
        library_get_source_Jac => Nitrogen_DSMC_get_source_Jac 

        ! Subroutine for temperatures 
        library_get_temperatures => Nitrogen_DSMC_get_temperatures 

        ! Subroutine for computing energy densities
        library_get_energy_densities => Nitrogen_DSMC_get_energy_densities

        ! Subroutines for computing thermodynamic data
        library_get_thermodynamic_data => Nitrogen_DSMC_get_thermodynamic_data

        ! Subroutine for computing temperatures and other data
        library_get_data => Nitrogen_DSMC_get_data

        ! Subroutine for computing transport properties and species mass diffusion fluxes
        library_get_transpCoeff => Nitrogen_DSMC_compute_transpCoeff
        library_get_species_DiffFlux => Nitrogen_DSMC_compute_species_DiffFlux

        ! Flag for post-processing
        flag_post_CR = .TRUE.

        ! Subroutine for post-processing
        library_write_fvmcc_solution => Nitrogen_DSMC_write_fvmcc_solution 

      END SUBROUTINE set_pointer_Nitrogen_DSMC
#endif

      !----------------------------------------------------!
      !> This subroutine associates pointers to the interface of the Nitrogen_NASA thermodynamic library. 
#ifdef NITROGEN_NASA
      SUBROUTINE set_pointer_Nitrogen_NASA ()

        USE mod_Nitrogen_NASA_library

        IMPLICIT NONE

        ! Library inizialization
        library_initialize => Nitrogen_NASA_initialize   

        ! Library closing 
        library_finalize => Nitrogen_NASA_finalize 

        ! Specific gas constants
        library_get_Ri => Nitrogen_NASA_get_Ri

        ! Subroutine for mixture density 
        library_get_density => Nitrogen_NASA_get_density

        ! Subroutine for mixture pressure 
        library_get_pressure => Nitrogen_NASA_get_pressure

        ! Subroutine for species mass fractions
        library_get_mass_fractions => Nitrogen_NASA_get_mass_fractions

        ! Subroutine for species maolar fractions
        library_get_molar_fractions => Nitrogen_NASA_get_molar_fractions 

        ! Subroutine for computing mass fractions and mixture molar mass from 
        ! molar fractions 
        library_get_mass_fractions_from_molar_fractions => Nitrogen_NASA_mass_frac_from_molar_frac  

        ! Subroutine for fix application to species molar fractions
        library_comp_tol => Nitrogen_NASA_comp_tol

        ! Subroutine for equilibrium composition computation
        library_compute_eq_composition => Nitrogen_NASA_compute_eq_composition

        ! Subroutine for source term 
        library_get_source => Nitrogen_NASA_get_source 

        ! Subroutine for source term and source term Jacobian
        library_get_source_Jac => Nitrogen_NASA_get_source_Jac 

        ! Subroutine for temperatures 
        library_get_temperatures => Nitrogen_NASA_get_temperatures 

        ! Subroutine for computing energy densities
        library_get_energy_densities => Nitrogen_NASA_get_energy_densities

        ! Subroutines for computing thermodynamic data
        library_get_thermodynamic_data => Nitrogen_NASA_get_thermodynamic_data

        ! Subroutine for computing temperatures and other data
        library_get_data => Nitrogen_NASA_get_data

        ! Subroutine for computing transport properties and species mass diffusion fluxes
        library_get_transpCoeff => Nitrogen_NASA_compute_transpCoeff
        library_get_species_DiffFlux => Nitrogen_NASA_compute_species_DiffFlux

        ! Flag for post-processing
        flag_post_CR = .TRUE.

        ! Subroutine for post-processing
        library_write_fvmcc_solution => Nitrogen_NASA_write_fvmcc_solution

      END SUBROUTINE set_pointer_Nitrogen_NASA
#endif

      !----------------------------------------------------!
      !> This subroutine associates pointers to the interface of the Nitrogen_Bari thermodynamic library. 
#ifdef NITROGEN_BARI
      SUBROUTINE set_pointer_Nitrogen_Bari ()

        USE mod_Nitrogen_Bari_library

        IMPLICIT NONE

        ! Library inizialization
        library_initialize => Nitrogen_Bari_initialize   

        ! Library closing 
        library_finalize => Nitrogen_Bari_finalize 

        ! Specific gas constants
        library_get_Ri => Nitrogen_Bari_get_Ri

        ! Subroutine for mixture density 
        library_get_density => Nitrogen_Bari_get_density

        ! Subroutine for mixture pressure 
        library_get_pressure => Nitrogen_Bari_get_pressure 

        ! Subroutine for species mass fractions
        library_get_mass_fractions => Nitrogen_Bari_get_mass_fractions

        ! Subroutine for species maolar fractions
        library_get_molar_fractions => Nitrogen_Bari_get_molar_fractions 

        ! Subroutine for computing mass fractions and mixture molar mass from 
        ! molar fractions 
        library_get_mass_fractions_from_molar_fractions => Nitrogen_Bari_mass_frac_from_molar_frac   

        ! Subroutine for fix application to species molar fractions
        library_comp_tol => Nitrogen_Bari_comp_tol

        ! Subroutine for equilibrium composition computation
        library_compute_eq_composition => Nitrogen_Bari_compute_eq_composition

        ! Subroutine for source term 
        library_get_source => Nitrogen_Bari_get_source 

        !  Subroutine for source term and source term Jacobian 
        library_get_source_Jac => Nitrogen_Bari_get_source_Jac 

        ! Subroutine for temperatures 
        library_get_temperatures => Nitrogen_Bari_get_temperatures 

        ! Subroutine for computing thermodynamic data
        library_get_thermodynamic_data => Nitrogen_Bari_get_thermodynamic_data

        ! Subroutine for computing energy densities
        library_get_energy_densities => Nitrogen_Bari_get_energy_densities

        ! Subroutine for computing temperatures and other data
        library_get_data => Nitrogen_Bari_get_data

        ! Subroutine for computing transport properties and species mass diffusion fluxes
        library_get_transpCoeff => Nitrogen_Bari_compute_transpCoeff
        library_get_species_DiffFlux => Nitrogen_Bari_compute_species_DiffFlux

        ! Flag for post-processing
        flag_post_CR = .TRUE.

        ! Subroutine for post-processing
        library_write_fvmcc_solution => Nitrogen_Bari_write_fvmcc_solution 

      END SUBROUTINE set_pointer_Nitrogen_Bari
#endif

      !----------------------------------------------------!
      !> This subroutine associates pointers to the interface of the Nitrogen_FHO thermodynamic library.
#ifdef NITROGEN_FHO
      SUBROUTINE set_pointer_Nitrogen_FHO ()

        USE mod_Nitrogen_FHO_library

        IMPLICIT NONE

        ! Library inizialization
        library_initialize => Nitrogen_FHO_initialize   

        ! Library closing 
        library_finalize => Nitrogen_FHO_finalize 

        ! Specific gas constants
        library_get_Ri => Nitrogen_FHO_get_Ri

        ! Subroutine for mixture density 
        library_get_density => Nitrogen_FHO_get_density

        ! Subroutine for mixture density 
        library_get_pressure => Nitrogen_FHO_get_pressure

        ! Subroutine for species mass fractions
        library_get_mass_fractions => Nitrogen_FHO_get_mass_fractions

        ! Subroutine for species maolar fractions
        library_get_molar_fractions => Nitrogen_FHO_get_molar_fractions 

        ! Subroutine for computing mass fractions and mixture molar mass from 
        ! molar fractions 
        library_get_mass_fractions_from_molar_fractions => Nitrogen_FHO_mass_frac_from_molar_frac  

        ! Subroutine for fix application to species molar fractions
        library_comp_tol => Nitrogen_FHO_comp_tol

        ! Subroutine for equilibrium composition computation
        library_compute_eq_composition => Nitrogen_FHO_compute_eq_composition

        ! Subroutine for source term 
        library_get_source => Nitrogen_FHO_get_source 

        ! Subroutine for source term and source term Jacobian 
        library_get_source_Jac => Nitrogen_FHO_get_source_Jac 

        ! Subroutine for temperatures 
        library_get_temperatures => Nitrogen_FHO_get_temperatures 

        ! Subroutine for computing energy densities
        library_get_energy_densities => Nitrogen_FHO_get_energy_densities

        ! Subroutines for computing thermodynamic data
        library_get_thermodynamic_data => Nitrogen_FHO_get_thermodynamic_data 

        ! Subroutine for computing temperatures and other data
        library_get_data => Nitrogen_FHO_get_data

        ! Subroutine for computing transport properties and species diffusion fluxes
        library_get_transpCoeff => Nitrogen_FHO_compute_transpCoeff
        library_get_species_DiffFlux => Nitrogen_FHO_compute_species_DiffFlux

        ! Flag for post-processing
        flag_post_CR = .TRUE.

        ! Subroutine for post-processing
        library_write_fvmcc_solution => Nitrogen_FHO_write_fvmcc_solution 

      END SUBROUTINE set_pointer_Nitrogen_FHO      
#endif

      !----------------------------------------------------!
      !> This subroutine associates pointers to the interface of the Nitrogen_TTLTH thermodynamic library.
#ifdef NITROGEN_TTLTH
      SUBROUTINE set_pointer_Nitrogen_TTLTH ()

        USE mod_Nitrogen_TTLTH_library

        IMPLICIT NONE

        ! Library inizialization
        library_initialize => Nitrogen_TTLTH_initialize   

        ! Library closing 
        library_finalize => Nitrogen_TTLTH_finalize 

        ! Specific gas constants
        library_get_Ri => Nitrogen_TTLTH_get_Ri

        ! Subroutine for mixture density 
        library_get_density => Nitrogen_TTLTH_get_density

        ! Subroutine for mixture density 
        library_get_pressure => Nitrogen_TTLTH_get_pressure 

        ! Subroutine for species mass fractions
        library_get_mass_fractions => Nitrogen_TTLTH_get_mass_fractions

        ! Subroutine for computing mass fractions and mixture molar mass from 
        ! molar fractions 
        library_get_mass_fractions_from_molar_fractions => Nitrogen_TTLTH_mass_frac_from_molar_frac  

        ! Subroutine for fix application to species molar fractions
        library_comp_tol => Nitrogen_TTLTH_comp_tol

        ! Subroutine for equilibrium composition computation
        library_compute_eq_composition => Nitrogen_TTLTH_compute_eq_composition

        ! Subroutine for source term 
        library_get_source => Nitrogen_TTLTH_get_source 

        ! Subroutine for source term and source term Jacobian 
        library_get_source_Jac => Nitrogen_TTLTH_get_source_Jac 

        ! Subroutine for temperatures 
        library_get_temperatures => Nitrogen_TTLTH_get_temperatures 

        ! Subroutine for computing energy densities
        library_get_energy_densities => Nitrogen_TTLTH_get_energy_densities

        ! Subroutines for computing thermodynamic data
        library_get_thermodynamic_data => Nitrogen_TTLTH_get_thermodynamic_data 

        ! Subroutine for computing temperatures and other data
        library_get_data => Nitrogen_TTLTH_get_data

        ! Subroutine for computing transport properties and species mass diffusion fluxes
        library_get_transpCoeff => Nitrogen_TTLTH_compute_transpCoeff
        library_get_species_DiffFlux => Nitrogen_TTLTH_compute_species_DiffFlux

        ! Flag for post-processing
        flag_post_CR = .TRUE.

        ! Subroutine for post-processing
        library_write_fvmcc_solution => Nitrogen_TTLTH_write_fvmcc_solution 

      END SUBROUTINE set_pointer_Nitrogen_TTLTH      
#endif

      !----------------------------------------------------!
      !> This subroutine associates pointers to the interface of the Argon_CR thermodynamic library.
#ifdef ARGON_CR
      SUBROUTINE set_pointer_Argon_CR ()

        USE mod_Argon_CR_library

        IMPLICIT NONE

        ! Library inizialization
        library_initialize => Argon_CR_initialize   

        ! Library closing 
        library_finalize => Argon_CR_finalize 

        ! Specific gas constants
        library_get_Ri => Argon_CR_get_Ri

        ! Free electron position in the species list and free electron specific gas constant
        library_get_el_data => Argon_CR_get_el_data

        ! Subroutine for mixture density 
        library_get_density => Argon_CR_get_density

        ! Subroutine for mixture pressure 
        library_get_pressure => Argon_CR_get_pressure

        ! Subroutine for species mass fractions
        library_get_mass_fractions => Argon_CR_get_mass_fractions

        ! Subroutine for species maolar fractions
        library_get_molar_fractions => Argon_CR_get_molar_fractions 

        ! Subroutine for computing mass fractions and mixture molar mass from 
        ! molar fractions 
        library_get_mass_fractions_from_molar_fractions => Argon_CR_mass_frac_from_molar_frac  

        ! Subroutine for fix application to species molar fractions
        library_comp_tol => Argon_CR_comp_tol

        ! Subroutine for equilibrium composition computation
        library_compute_eq_composition => Argon_CR_compute_eq_composition

        ! Subroutine for source term 
        library_get_source => Argon_CR_get_source 

        ! Subroutine for source term and source term Jacobian 
        library_get_source_Jac => Argon_CR_get_source_Jac 

        ! Subroutine for temperatures 
        library_get_temperatures => Argon_CR_get_temperatures 

        ! Subroutine for computing energy densities
        library_get_energy_densities => Argon_CR_get_energy_densities

        ! Subroutines for computing thermodynamic data
        library_get_thermodynamic_data => Argon_CR_get_thermodynamic_data 

        ! Subroutine for computing temperatures and other data
        library_get_data => Argon_CR_get_data

        ! Subroutine for computing transport properties and species mass diffusion fluxes
        library_get_transpCoeff => Argon_CR_compute_transpCoeff
        library_get_species_DiffFlux => Argon_CR_compute_species_DiffFlux

        ! Flag for post-processing
        flag_post_CR = .TRUE.

        ! Subroutine for post-processing
        library_write_fvmcc_solution => Argon_CR_write_fvmcc_solution 

      END SUBROUTINE set_pointer_Argon_CR  
#endif

      !----------------------------------------------------!
      !> This subroutine associates pointers to the interface of the Mutationpp thermodynamic library.
#ifdef MUTATIONPP
      SUBROUTINE set_pointer_Mutationpp ()

        USE mod_Mutationpp_library

        IMPLICIT NONE

        ! Library inizialization
        library_initialize => Mutationpp_initialize   

        ! Library closing 
        library_finalize => Mutationpp_finalize 

        ! Specific gas constants
        library_get_Ri => Mutationpp_get_Ri

        ! Subroutine for mixture density 
        library_get_density => Mutationpp_get_density

        ! Subroutine for mixture pressure 
        library_get_pressure => Mutationpp_get_pressure

        ! Subroutine for species mass fractions
        library_get_mass_fractions => Mutationpp_get_mass_fractions

        ! Subroutine for species maolar fractions
        library_get_molar_fractions => Mutationpp_get_molar_fractions  

        ! Subroutine for computing mass fractions and mixture molar mass from 
        ! molar fractions 
        library_get_mass_fractions_from_molar_fractions => Mutationpp_mass_frac_from_molar_frac

        ! Subroutine for fix application to species molar fractions
        library_comp_tol => Mutationpp_comp_tol

        ! Subroutine for equilibrium composition computation
        library_compute_eq_composition => Mutationpp_compute_eq_composition

        ! Subroutine for pyrolyis gas equilibrium composition computation
        library_compute_eq_composition_pyro => Mutationpp_compute_eq_composition_pyro

        ! Subroutine for source term 
        library_get_source => Mutationpp_get_source 

        ! Subroutine for source term and source term Jacobian 
        library_get_source_Jac => Mutationpp_get_source_Jac 

        ! Subroutine for temperatures 
        library_get_temperatures => Mutationpp_get_temperatures 

        ! Subroutine for temperatures 
        library_get_prandtl_number => Mutationpp_get_prandtl_number 

        ! Subroutine for computing energy densities
        library_get_energy_densities => Mutationpp_get_energy_densities

        ! Subroutines for computing thermodynamic data
        library_get_thermodynamic_data => Mutationpp_get_thermodynamic_data 

        ! Subroutine for computing temperatures and other data
        library_get_data => Mutationpp_get_data       
 
        ! Subroutine for computing transport properties and species mass diffusion fluxes
        library_get_transpCoeff => Mutationpp_compute_transpCoeff
        library_get_species_DiffFlux => Mutationpp_compute_species_DiffFlux
        library_get_thermal_cond => Mutationpp_compute_thermal_cond 

        ! Subroutine for computing species enthalpies
        library_get_species_enthalpies => Mutationpp_get_species_enthalpy

        ! Subroutine for computing the mixture enthalpy
        library_get_enthalpy => Mutationpp_get_enthalpy

        ! Subroutine for getting the species index
        library_get_species_index => Mutationpp_get_species_index

        ! Subroutine for getting the species index
        library_get_species_name => Mutationpp_get_species_name

        ! Flag for post-processing
        flag_post_CR = .FALSE.

        ! Subroutine for post-processing
        library_write_fvmcc_solution => Mutationpp_write_fvmcc_solution 

        ! Subroutine to obtain the species name
        library_species_name => Mutationpp_species_name 

        ! Subroutine to obtain the elements of a mixture
        library_number_of_elements => Mutationpp_number_of_elements

        ! Subroutine for computing molar fractions and mixture molar mass from
        ! mass fractions 
        library_get_molar_fractions_from_mass_fractions => Mutationpp_molar_frac_from_mass_frac

        ! Subroutine for equilibrium composition computation considering
        ! elements for meteor ablation
        library_compute_eq_composition_meteor => Mutationpp_compute_eq_composition_meteor
       
        ! Subroutine for calculating the elemental mole fraction of a mixture 
        library_compute_elemental_mole_fraction => Mutationpp_compute_elemental_mole_fraction

        ! Subroutine to obtain the surface mass balance with contraints
        library_surface_mass_balance_constraint => Mutationpp_surface_mass_balance_constraint

        ! Subroutine to obtain the elemental mass fractions from species mass fractions
        library_get_element_mass_fractions => Mutationpp_get_element_mass_fractions
        
        ! Subroutine to set the diffusion model for GSI 
        library_set_diff_model => Mutationpp_set_diff_model

        ! Subroutine to set the conductive heat flux for GSI 
        library_set_cond_heat_flux => Mutationpp_set_cond_heat_flux

        library_set_wall_radiation => Mutationpp_set_wall_radiation

        ! Subroutine to set the wall state for GSI 
        library_set_wall_state => Mutationpp_set_wall_state

        ! Subroutine to solve SMB from GSI 
        library_solve_surface_balance => Mutationpp_solve_surface_balance

        ! Subroutine to get the wall state from GSI 
        library_get_wall_state => Mutationpp_get_wall_state

        ! Subroutine to get the production rates from GSI 
        library_get_surface_production_rates => Mutationpp_get_surface_production_rates

        ! Subroutine to get the production rates from GSI 
        library_get_mass_blowing_rate => Mutationpp_get_mass_blowing_rate

        library_get_frozen_sound_speed => Mutationpp_get_frozen_sound_speed
                    
      END SUBROUTINE set_pointer_Mutationpp      
#endif

      !----------------------------------------------------!
      !> This subroutine associates pointers to the interface of the Polyat_gas_WCU thermodynamic library.
#ifdef POLYAT_GAS_WCU
      SUBROUTINE set_pointer_Polyat_gas_WCU ()

        USE mod_Polyat_gas_WCU_library

        IMPLICIT NONE

        ! Library inizialization
        library_initialize => Polyat_gas_WCU_initialize   

        ! Library closing 
        library_finalize => Polyat_gas_WCU_finalize 

        ! Specific gas constants
        library_get_Ri => Polyat_gas_WCU_get_Ri

        ! Subroutine for mixture density 
        library_get_density => Polyat_gas_WCU_get_density

        ! Subroutine for mixture pressure 
        library_get_pressure => Polyat_gas_WCU_get_pressure

        ! Subroutine for species mass fractions
        library_get_mass_fractions => Polyat_gas_WCU_get_mass_fractions

        ! Subroutine for species maolar fractions
        library_get_molar_fractions => Polyat_gas_WCU_get_molar_fractions  

        ! Subroutine for computing mass fractions and mixture molar mass from 
        ! molar fractions 
        library_get_mass_fractions_from_molar_fractions => Polyat_gas_WCU_mass_frac_from_molar_frac

        ! Subroutine for fix application to species molar fractions
        library_comp_tol => Polyat_gas_WCU_comp_tol

        ! Subroutine for equilibrium composition computation
        library_compute_eq_composition => Polyat_gas_WCU_compute_eq_composition

        ! Subroutine for source term 
        library_get_source => Polyat_gas_WCU_get_source 

        ! Subroutine for source term and source term Jacobian 
        library_get_source_Jac => Polyat_gas_WCU_get_source_Jac 

        ! Subroutine for temperatures 
        library_get_temperatures => Polyat_gas_WCU_get_temperatures 

        ! Subroutine for computing energy densities
        library_get_energy_densities => Polyat_gas_WCU_get_energy_densities

        ! Subroutines for computing thermodynamic data
        library_get_thermodynamic_data => Polyat_gas_WCU_get_thermodynamic_data 

        ! Subroutine for computing temperatures and other data
        library_get_data => Polyat_gas_WCU_get_data       
 
        ! Subroutine for computing transport properties and species mass diffusion fluxes
        library_get_transpCoeff => Polyat_gas_WCU_compute_transpCoeff
        library_get_species_DiffFlux => Polyat_gas_WCU_compute_species_DiffFlux

        ! Flag for post-processing
        flag_post_CR = .FALSE.

        ! Subroutine for post-processing
        library_write_fvmcc_solution => Polyat_gas_WCU_write_fvmcc_solution 

      END SUBROUTINE set_pointer_Polyat_gas_WCU      
#endif

   END SUBROUTINE set_library_pointer
!------------------------------------------------------------------------------!
