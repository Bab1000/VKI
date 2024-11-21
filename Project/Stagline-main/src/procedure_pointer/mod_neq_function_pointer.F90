!------------------------------------------------------------------------------!
!> This module provides the support for function and subroutine pointers for thermodynamic
!! libraries to be used for 1D or 2D nonequilibrium flow computations.
  MODULE mod_neq_function_pointer

#include"../config.h"

    IMPLICIT NONE

    !------------------------------------------------------!
    ! Interfaces for function/subroutines and related pointers

    ! Interface for subroutine initializing the thermodynamic library.
    ABSTRACT INTERFACE

      SUBROUTINE lib_in (ns, ntrot, ntvib, nte, ntemp, neq, ndim, solver, mixture, reaction, state, transf, path, xi_tol)
        INTEGER, INTENT(IN) :: ns, ntrot, ntvib, nte, ntemp, neq, ndim
        REAL(KIND=8), INTENT(IN) :: xi_tol
        CHARACTER*(*), INTENT(IN) ::  solver, mixture, reaction, state, transf, path
      END SUBROUTINE lib_in

      ! Interface for subroutine shutting down the thermodynamic library 
      SUBROUTINE lib_fin ()  
      END SUBROUTINE lib_fin 

      ! Inteface for subroutine providing species gas constants
      SUBROUTINE lib_get_Ri (Ri) 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: Ri
      END SUBROUTINE lib_get_Ri

      ! Interface for subroutine providing the position of the electron species
      ! in the species list and the electron gas constant
      SUBROUTINE lib_get_el_data (pos_el, Re)
        INTEGER, INTENT(OUT) :: pos_el
        REAL(KIND=8), INTENT(OUT) :: Re
      END SUBROUTINE lib_get_el_data

      ! Interface for subroutine computing the mixture density
      SUBROUTINE lib_get_density (rhoi, rho)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), INTENT(OUT) :: rho
      END SUBROUTINE lib_get_density

      ! Interface for subroutine computing the mixture density
      SUBROUTINE lib_get_pressure (rhoi, temp, p )
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp
        REAL(KIND=8), INTENT(OUT) :: p
      END SUBROUTINE lib_get_pressure 

      ! Interface for subroutine computing the species mass fractions
      SUBROUTINE lib_get_mass_fractions (rho, rhoi, yi)
        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi
      END SUBROUTINE lib_get_mass_fractions

      ! Interface for subroutine computing the species molar fractions
      SUBROUTINE lib_get_molar_fractions (rhoi, xi)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: xi
      END SUBROUTINE lib_get_molar_fractions

      ! Interface for subroutine computing the species enthalpies 
      SUBROUTINE lib_get_species_enthalpies (T, hi)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: hi
      END SUBROUTINE lib_get_species_enthalpies

      ! Interface for subroutine computing the mixture enthalpy
      SUBROUTINE lib_get_enthalpy (rho, rhoi, T, h )
        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: h
      END SUBROUTINE lib_get_enthalpy

      ! Interface for subroutine computing the species mass fractions and
      ! mixture molar mass from the species mass fractions
      SUBROUTINE lib_get_mass_fractions_from_molar_fractions(xi, mm, yi)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi
        REAL(KIND=8), INTENT(OUT) :: mm
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi
      END SUBROUTINE lib_get_mass_fractions_from_molar_fractions

      ! Interface for subroutine computing the species mole fractions and
      ! mixture molar mass from the species mass fractions
      SUBROUTINE lib_get_molar_fractions_from_mass_fractions(yi, mm, xi)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: yi
        REAL(KIND=8), INTENT(OUT) :: mm
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: xi
      END SUBROUTINE lib_get_molar_fractions_from_mass_fractions


      ! Interface for subroutine apply a small number to the composition respecting the mass constraint
      ! (sum of species mole fractions equal to one). The tolerance must be specified when initializing the library
      SUBROUTINE lib_comp_tol(xi)
        REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: xi
      END SUBROUTINE lib_comp_tol

      ! Interface for subroutine computing source term vector (no Jacobian in output).
      SUBROUTINE lib_compute_source_term_kinetics (rhoi, temp, source) 
        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: source
      END SUBROUTINE lib_compute_source_term_kinetics

      ! Interface for subroutine computing source term vector (Jacobian in output).
      SUBROUTINE lib_compute_source_term_kinetics_Jac (rhoi, temp, source, jacobian)
        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: source, jacobian
      END SUBROUTINE lib_compute_source_term_kinetics_Jac

      ! Interface for subroutine computing temperatures
      SUBROUTINE lib_compute_temperatures (rhoi, rho_eint, temp)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, rho_eint
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp
      END SUBROUTINE lib_compute_temperatures 

      ! Interface for subroutine computing temperatures
      SUBROUTINE lib_get_prandtl_number (rhoi, temp, pr)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), INTENT(OUT) :: pr
      END SUBROUTINE lib_get_prandtl_number

      ! Interface for subroutine computing the internal energy density vector
      SUBROUTINE lib_compute_energy_densities (rhoi, temp, rho_eint)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp     
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: rho_eint
      END SUBROUTINE lib_compute_energy_densities

      ! Interface for subroutine computing thermodynamic data
      SUBROUTINE lib_compute_thermodynamic_data (rho, rhoi, temp, c, gamma, p, alpha, beta_vec, energy_data)
        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), INTENT(OUT) :: c, gamma, p, alpha
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: beta_vec, energy_data
      END SUBROUTINE lib_compute_thermodynamic_data

      ! Interface for subroutine computing temperatures thermodynamic data
      SUBROUTINE lib_compute_data (rho, rhoi, rho_eint, temp, c, gamma, p, alpha, beta_vec, energy_data)
        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, rho_eint
        REAL(KIND=8), INTENT(OUT) :: c, gamma, p, alpha
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp, beta_vec, energy_data
      END SUBROUTINE lib_compute_data

      ! Interface for subroutine computing transport properties
      SUBROUTINE lib_compute_transpCoeff (p, xi, temp, mu, kappa, lambda, Di, chi)
        REAL(KIND=8), INTENT(IN) :: p
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, temp
        REAL(KIND=8), INTENT(OUT) :: mu, kappa
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: Di, chi, lambda
      END SUBROUTINE lib_compute_transpCoeff

      ! Interface for subroutine computing species mass diffusion flux
      SUBROUTINE lib_compute_species_DiffFlux (p, T, Te, xi, diff_driv, Ji)
        REAL(KIND=8), INTENT(IN) :: p, T, Te
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, diff_driv
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: Ji
      END SUBROUTINE lib_compute_species_DiffFlux

      ! Interface for subroutine computing the equilibrium composition
      SUBROUTINE lib_compute_eq_composition (p, T, rhoi)
        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: rhoi
      END SUBROUTINE lib_compute_eq_composition

      ! Interface for subroutine computing the equilibrium composition of the
      ! pyrolysis gas
      SUBROUTINE lib_compute_eq_composition_pyro (p, T, yi)
        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi
      END SUBROUTINE lib_compute_eq_composition_pyro

      ! Interface for subroutine computing thermal conductivity
      SUBROUTINE lib_compute_thermal_cond (rhoi, T, lambda)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: T
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: lambda
      END SUBROUTINE lib_compute_thermal_cond

      ! Interface for subroutine performing solution post-processing 
      SUBROUTINE lib_write_fvmcc_solution (flow_file, pop_file, cell, xold, x, uold, u, rhoi, temp)
        INTEGER, INTENT(IN) :: flow_file, pop_file, cell
        REAL(KIND=8), INTENT(IN) :: xold, x, uold, u
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
      END SUBROUTINE lib_write_fvmcc_solution

      ! Interface for subroutine to obtain species name from Mutation++
      SUBROUTINE lib_compute_species_name (species_name)
        CHARACTER(Len=*), INTENT(OUT), DIMENSION(:) :: species_name
      END SUBROUTINE lib_compute_species_name

      ! Interface for subroutine getting the species index
      SUBROUTINE lib_get_species_index (sp_name, sp_index)
        CHARACTER*(*), INTENT(IN) :: sp_name
        INTEGER, INTENT(OUT)      :: sp_index
      END SUBROUTINE lib_get_species_index

      ! Interface for subroutine getting the species name
      SUBROUTINE lib_get_species_name (sp_index, sp_name)
        INTEGER, INTENT(IN)        :: sp_index
        CHARACTER*(*), INTENT(OUT) :: sp_name
      END SUBROUTINE lib_get_species_name 

     SUBROUTINE lib_number_of_elements(elements) 
        INTEGER, INTENT(INOUT) :: elements
      END SUBROUTINE lib_number_of_elements

           ! Interface for subroutine computing the equilibrium composition for
      ! meteor
      SUBROUTINE lib_compute_eq_composition_meteor (p, T, elements, rhoi)
        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: elements
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: rhoi
      END SUBROUTINE lib_compute_eq_composition_meteor

      ! Interface for subroutine computing the elemental mole fraction for a
      ! mixture
      SUBROUTINE lib_compute_elemental_mole_fraction (x_species,x_elements)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: x_species
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: x_elements    
      END SUBROUTINE lib_compute_elemental_mole_fraction

     ! Interface for subroutine to obtain surface mass balance with contraint from Mutation++
     SUBROUTINE lib_compute_surface_mass_balance_constraint (input_mixture, T, P, xi) 
       
        CHARACTER*(*), INTENT(IN) :: input_mixture
        REAL(KIND=8), INTENT(IN) :: T, P
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: xi

     END SUBROUTINE lib_compute_surface_mass_balance_constraint

     SUBROUTINE lib_get_element_mass_fractions (yi, ye)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: yi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ye
     END SUBROUTINE lib_get_element_mass_fractions

     SUBROUTINE lib_set_diff_model (xi, dx)
       
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi  
        REAL(KIND=8), INTENT(IN) :: dx 

     END SUBROUTINE lib_set_diff_model

     SUBROUTINE lib_set_wall_radiation (rad_flux)
       
        REAL(KIND=8), INTENT(IN) :: rad_flux

     END SUBROUTINE lib_set_wall_radiation
          
     SUBROUTINE lib_set_cond_heat_flux (T, dx)
       
       REAL(KIND=8), DIMENSION(:), INTENT(IN) :: T  
       REAL(KIND=8), INTENT(IN) :: dx 

     END SUBROUTINE lib_set_cond_heat_flux

     SUBROUTINE lib_set_wall_state (rhoi, T)
       
       REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi 
       REAL(KIND=8), DIMENSION(:), INTENT(IN) :: T 

     END SUBROUTINE lib_set_wall_state


     SUBROUTINE lib_solve_surface_balance()
       
      
     END SUBROUTINE lib_solve_surface_balance

   
     SUBROUTINE lib_get_wall_state (rhoi, T)
       
       REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: rhoi 
       REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: T
      
     END SUBROUTINE lib_get_wall_state

     SUBROUTINE lib_get_surface_production_rates (wdoti)

       REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: wdoti
                   
     END SUBROUTINE lib_get_surface_production_rates

     SUBROUTINE lib_get_mass_blowing_rate(mdot)

       REAL(KIND=8), INTENT(OUT) :: mdot 
                   
     END SUBROUTINE lib_get_mass_blowing_rate

     SUBROUTINE lib_get_frozen_sound_speed(rhoi, temp, c)

      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp, rhoi
      REAL(KIND=8), INTENT(OUT) :: c
                   
     END SUBROUTINE lib_get_frozen_sound_speed



    END INTERFACE

    ! Function/subroutine pointer definition
    PROCEDURE(lib_in), POINTER, SAVE :: library_initialize
    PROCEDURE(lib_fin), POINTER, SAVE :: library_finalize
    PROCEDURE(lib_get_Ri), POINTER, SAVE :: library_get_Ri
    PROCEDURE(lib_get_el_data), POINTER, SAVE :: library_get_el_data
    PROCEDURE(lib_get_density), POINTER, SAVE :: library_get_density
    PROCEDURE(lib_get_pressure), POINTER, SAVE :: library_get_pressure 
    PROCEDURE(lib_get_species_enthalpies), POINTER, SAVE :: library_get_species_enthalpies
    PROCEDURE(lib_get_enthalpy), POINTER, SAVE :: library_get_enthalpy
    PROCEDURE(lib_get_mass_fractions), POINTER, SAVE :: library_get_mass_fractions
    PROCEDURE(lib_get_molar_fractions), POINTER, SAVE :: library_get_molar_fractions
    PROCEDURE(lib_get_mass_fractions_from_molar_fractions), POINTER, SAVE :: library_get_mass_fractions_from_molar_fractions
    PROCEDURE(lib_comp_tol), POINTER, SAVE :: library_comp_tol
    PROCEDURE(lib_compute_source_term_kinetics), POINTER, SAVE ::library_get_source
    PROCEDURE(lib_compute_source_term_kinetics_Jac), POINTER, SAVE ::library_get_source_Jac
    PROCEDURE(lib_compute_temperatures), POINTER, SAVE :: library_get_temperatures
    PROCEDURE(lib_get_prandtl_number), POINTER, SAVE :: library_get_prandtl_number
    PROCEDURE(lib_compute_energy_densities), POINTER, SAVE :: library_get_energy_densities
    PROCEDURE(lib_compute_thermodynamic_data), POINTER, SAVE :: library_get_thermodynamic_data
    PROCEDURE(lib_compute_data), POINTER, SAVE :: library_get_data
    PROCEDURE(lib_compute_transpCoeff), POINTER, SAVE :: library_get_transpCoeff
    PROCEDURE(lib_compute_species_DiffFlux), POINTER, SAVE :: library_get_species_DiffFlux
    PROCEDURE(lib_compute_eq_composition), POINTER, SAVE :: library_compute_eq_composition
    PROCEDURE(lib_compute_thermal_cond), POINTER, SAVE :: library_get_thermal_cond
    PROCEDURE(lib_write_fvmcc_solution), POINTER, SAVE :: library_write_fvmcc_solution
    PROCEDURE(lib_compute_eq_composition_pyro), POINTER, SAVE :: library_compute_eq_composition_pyro
    PROCEDURE(lib_get_species_index), POINTER, SAVE :: library_get_species_index
    PROCEDURE(lib_get_species_name), POINTER, SAVE :: library_get_species_name
    PROCEDURE(lib_compute_species_name), POINTER, SAVE :: library_species_name
    PROCEDURE(lib_number_of_elements), POINTER, SAVE :: library_number_of_elements
    PROCEDURE(lib_get_molar_fractions_from_mass_fractions), POINTER, SAVE :: library_get_molar_fractions_from_mass_fractions
    PROCEDURE(lib_compute_eq_composition_meteor), POINTER, SAVE :: library_compute_eq_composition_meteor
    PROCEDURE(lib_compute_elemental_mole_fraction), POINTER, SAVE :: library_compute_elemental_mole_fraction
    PROCEDURE(lib_compute_surface_mass_balance_constraint), POINTER, SAVE :: library_surface_mass_balance_constraint
    PROCEDURE(lib_get_element_mass_fractions), POINTER, SAVE :: library_get_element_mass_fractions
    PROCEDURE(lib_set_diff_model), POINTER, SAVE :: library_set_diff_model
    PROCEDURE(lib_set_cond_heat_flux), POINTER, SAVE :: library_set_cond_heat_flux 
    PROCEDURE(lib_set_wall_radiation), POINTER, SAVE :: library_set_wall_radiation
    PROCEDURE(lib_set_wall_state), POINTER, SAVE :: library_set_wall_state
    PROCEDURE(lib_solve_surface_balance), POINTER, SAVE :: library_solve_surface_balance
    PROCEDURE(lib_get_wall_state), POINTER, SAVE :: library_get_wall_state
    PROCEDURE(lib_get_surface_production_rates), POINTER, SAVE :: library_get_surface_production_rates
    PROCEDURE(lib_get_mass_blowing_rate), POINTER, SAVE :: library_get_mass_blowing_rate
    PROCEDURE(lib_get_frozen_sound_speed), POINTER, SAVE :: library_get_frozen_sound_speed



    ! Subroutines of Nitrogen_Park library interface
#ifdef NITROGEN_PARK
    PROCEDURE(lib_in) :: Nitrogen_Park_initialize
    PROCEDURE(lib_fin) :: Nitrogen_Park_finalize
    PROCEDURE(lib_get_Ri) :: Nitrogen_Park_get_Ri
    PROCEDURE(lib_get_density) :: Nitrogen_Park_get_density 
    PROCEDURE(lib_get_pressure) :: Nitrogen_Park_get_pressure
    PROCEDURE(lib_get_mass_fractions) :: Nitrogen_Park_get_mass_fractions 
    PROCEDURE(lib_get_molar_fractions) :: Nitrogen_Park_get_molar_fractions
    PROCEDURE(lib_get_mass_fractions_from_molar_fractions) :: Nitrogen_Park_mass_frac_from_molar_frac
    PROCEDURE(lib_comp_tol) :: Nitrogen_Park_comp_tol
    PROCEDURE(lib_compute_source_term_kinetics) :: Nitrogen_Park_get_source
    PROCEDURE(lib_compute_source_term_kinetics_Jac) :: Nitrogen_Park_get_source_Jac
    PROCEDURE(lib_compute_temperatures) :: Nitrogen_Park_get_temperatures
    PROCEDURE(lib_compute_energy_densities) :: Nitrogen_Park_get_energy_densities
    PROCEDURE(lib_compute_thermodynamic_data) :: Nitrogen_Park_get_thermodynamic_data
    PROCEDURE(lib_compute_data) :: Nitrogen_Park_get_data
    PROCEDURE(lib_compute_transpCoeff) :: Nitrogen_Park_compute_transpCoeff
    PROCEDURE(lib_compute_species_DiffFlux) :: Nitrogen_Park_compute_species_DiffFlux
    PROCEDURE(lib_compute_eq_composition) :: Nitrogen_Park_compute_eq_composition
    PROCEDURE(lib_write_fvmcc_solution) :: Nitrogen_Park_write_fvmcc_solution
#endif

    ! Subroutines of Nitrogen_DSMC library interface
#ifdef NITROGEN_DSMC
    PROCEDURE(lib_in) :: Nitrogen_DSMC_initialize
    PROCEDURE(lib_fin) :: Nitrogen_DSMC_finalize
    PROCEDURE(lib_get_Ri) :: Nitrogen_DSMC_get_Ri
    PROCEDURE(lib_get_density) :: Nitrogen_DSMC_get_density 
    PROCEDURE(lib_get_pressure) :: Nitrogen_DSMC_get_pressure
    PROCEDURE(lib_get_mass_fractions) :: Nitrogen_DSMC_get_mass_fractions
    PROCEDURE(lib_get_molar_fractions) :: Nitrogen_DSMC_get_molar_fractions
    PROCEDURE(lib_get_mass_fractions_from_molar_fractions) :: Nitrogen_DSMC_mass_frac_from_molar_frac
    PROCEDURE(lib_comp_tol) :: Nitrogen_DSMC_comp_tol
    PROCEDURE(lib_compute_source_term_kinetics) :: Nitrogen_DSMC_get_source
    PROCEDURE(lib_compute_source_term_kinetics_Jac) :: Nitrogen_DSMC_get_source_Jac
    PROCEDURE(lib_compute_temperatures) :: Nitrogen_DSMC_get_temperatures
    PROCEDURE(lib_compute_energy_densities) :: Nitrogen_DSMC_get_energy_densities
    PROCEDURE(lib_compute_thermodynamic_data) :: Nitrogen_DSMC_get_thermodynamic_data
    PROCEDURE(lib_compute_data) :: Nitrogen_DSMC_get_data
    PROCEDURE(lib_compute_transpCoeff) :: Nitrogen_DSMC_compute_transpCoeff
    PROCEDURE(lib_compute_species_DiffFlux) :: Nitrogen_DSMC_compute_species_DiffFlux
    PROCEDURE(lib_compute_eq_composition) :: Nitrogen_DSMC_compute_eq_composition
    PROCEDURE(lib_write_fvmcc_solution) :: Nitrogen_DSMC_write_fvmcc_solution
#endif

    ! Subroutines of Nitrogen_NASA library interface
#ifdef NITROGEN_NASA
    PROCEDURE(lib_in) :: Nitrogen_NASA_initialize
    PROCEDURE(lib_fin) :: Nitrogen_NASA_finalize
    PROCEDURE(lib_get_Ri) :: Nitrogen_NASA_get_Ri
    PROCEDURE(lib_get_density) :: Nitrogen_NASA_get_density 
    PROCEDURE(lib_get_pressure) :: Nitrogen_NASA_get_pressure 
    PROCEDURE(lib_get_mass_fractions) :: Nitrogen_NASA_get_mass_fractions
    PROCEDURE(lib_get_molar_fractions) :: Nitrogen_NASA_get_molar_fractions
    PROCEDURE(lib_get_mass_fractions_from_molar_fractions) :: Nitrogen_NASA_mass_frac_from_molar_frac
    PROCEDURE(lib_comp_tol) :: Nitrogen_NASA_comp_tol 
    PROCEDURE(lib_compute_source_term_kinetics) :: Nitrogen_NASA_get_source
    PROCEDURE(lib_compute_source_term_kinetics_Jac) :: Nitrogen_NASA_get_source_Jac
    PROCEDURE(lib_compute_temperatures) :: Nitrogen_NASA_get_temperatures
    PROCEDURE(lib_compute_energy_densities) :: Nitrogen_NASA_get_energy_densities
    PROCEDURE(lib_compute_thermodynamic_data) :: Nitrogen_NASA_get_thermodynamic_data
    PROCEDURE(lib_compute_data) :: Nitrogen_NASA_get_data
    PROCEDURE(lib_compute_transpCoeff) :: Nitrogen_NASA_compute_transpCoeff
    PROCEDURE(lib_compute_species_DiffFlux) :: Nitrogen_NASA_compute_species_DiffFlux
    PROCEDURE(lib_compute_eq_composition) :: Nitrogen_NASA_compute_eq_composition
    PROCEDURE(lib_write_fvmcc_solution) :: Nitrogen_NASA_write_fvmcc_solution
#endif

    ! Subroutines of Nitrogen_Bari library interface
#ifdef NITROGEN_BARI
    PROCEDURE(lib_in) :: Nitrogen_Bari_initialize
    PROCEDURE(lib_fin) :: Nitrogen_Bari_finalize
    PROCEDURE(lib_get_Ri) :: Nitrogen_Bari_get_Ri
    PROCEDURE(lib_get_density) :: Nitrogen_Bari_get_density 
    PROCEDURE(lib_get_pressure) :: Nitrogen_Bari_get_pressure 
    PROCEDURE(lib_get_mass_fractions) :: Nitrogen_Bari_get_mass_fractions
    PROCEDURE(lib_get_molar_fractions) :: Nitrogen_Bari_get_molar_fractions 
    PROCEDURE(lib_get_mass_fractions_from_molar_fractions) :: Nitrogen_Bari_mass_frac_from_molar_frac 
    PROCEDURE(lib_comp_tol) :: Nitrogen_Bari_comp_tol
    PROCEDURE(lib_compute_source_term_kinetics) :: Nitrogen_Bari_get_source
    PROCEDURE(lib_compute_source_term_kinetics_Jac) :: Nitrogen_Bari_get_source_Jac
    PROCEDURE(lib_compute_temperatures) :: Nitrogen_Bari_get_temperatures
    PROCEDURE(lib_compute_energy_densities) :: Nitrogen_Bari_get_energy_densities
    PROCEDURE(lib_compute_thermodynamic_data) :: Nitrogen_Bari_get_thermodynamic_data
    PROCEDURE(lib_compute_data) :: Nitrogen_Bari_get_data
    PROCEDURE(lib_compute_transpCoeff) :: Nitrogen_Bari_compute_transpCoeff
    PROCEDURE(lib_compute_species_DiffFlux) :: Nitrogen_Bari_compute_species_DiffFlux 
    PROCEDURE(lib_compute_eq_composition) :: Nitrogen_Bari_compute_eq_composition  
    PROCEDURE(lib_write_fvmcc_solution) :: Nitrogen_Bari_write_fvmcc_solution
#endif

    ! Subroutines of Nitrogen_FHO library interface
#ifdef NITROGEN_FHO
    PROCEDURE(lib_in) :: Nitrogen_FHO_initialize
    PROCEDURE(lib_fin) :: Nitrogen_FHO_finalize
    PROCEDURE(lib_get_Ri) :: Nitrogen_FHO_get_Ri
    PROCEDURE(lib_get_density) :: Nitrogen_FHO_get_density 
    PROCEDURE(lib_get_pressure) :: Nitrogen_FHO_get_pressure
    PROCEDURE(lib_get_mass_fractions) :: Nitrogen_FHO_get_mass_fractions
    PROCEDURE(lib_get_molar_fractions) :: Nitrogen_FHO_get_molar_fractions 
    PROCEDURE(lib_get_mass_fractions_from_molar_fractions) :: Nitrogen_FHO_mass_frac_from_molar_frac 
    PROCEDURE(lib_comp_tol) :: Nitrogen_FHO_comp_tol
    PROCEDURE(lib_compute_source_term_kinetics) :: Nitrogen_FHO_get_source
    PROCEDURE(lib_compute_source_term_kinetics_Jac) :: Nitrogen_FHO_get_source_Jac
    PROCEDURE(lib_compute_temperatures) :: Nitrogen_FHO_get_temperatures
    PROCEDURE(lib_compute_energy_densities) :: Nitrogen_FHO_get_energy_densities
    PROCEDURE(lib_compute_thermodynamic_data) :: Nitrogen_FHO_get_thermodynamic_data
    PROCEDURE(lib_compute_data) :: Nitrogen_FHO_get_data
    PROCEDURE(lib_compute_transpCoeff) :: Nitrogen_FHO_compute_transpCoeff
    PROCEDURE(lib_compute_species_DiffFlux) :: Nitrogen_FHO_compute_species_DiffFlux
    PROCEDURE(lib_compute_eq_composition) :: Nitrogen_FHO_compute_eq_composition  
    PROCEDURE(lib_write_fvmcc_solution) :: Nitrogen_FHO_write_fvmcc_solution
#endif

! Subroutines of Nitrogen_TTLTH library interface
#ifdef NITROGEN_TTLTH
    PROCEDURE(lib_in) :: Nitrogen_TTLTH_initialize
    PROCEDURE(lib_fin) :: Nitrogen_TTLTH_finalize
    PROCEDURE(lib_get_Ri) :: Nitrogen_TTLTH_get_Ri
    PROCEDURE(lib_get_density) :: Nitrogen_TTLTH_get_density 
    PROCEDURE(lib_get_pressure) :: Nitrogen_TTLTH_get_pressure
    PROCEDURE(lib_get_mass_fractions) :: Nitrogen_TTLTH_get_mass_fractions
    PROCEDURE(lib_get_molar_fractions) :: Nitrogen_TTLTH_get_molar_fractions 
    PROCEDURE(lib_get_mass_fractions_from_molar_fractions) :: Nitrogen_TTLTH_mass_frac_from_molar_frac
    PROCEDURE(lib_comp_tol) :: Nitrogen_TTLTH_comp_tol
    PROCEDURE(lib_compute_source_term_kinetics) :: Nitrogen_TTLTH_get_source
    PROCEDURE(lib_compute_source_term_kinetics_Jac) :: Nitrogen_TTLTH_get_source_Jac
    PROCEDURE(lib_compute_temperatures) :: Nitrogen_TTLTH_get_temperatures
    PROCEDURE(lib_compute_energy_densities) :: Nitrogen_TTLTH_get_energy_densities
    PROCEDURE(lib_compute_thermodynamic_data) :: Nitrogen_TTLTH_get_thermodynamic_data
    PROCEDURE(lib_compute_data) :: Nitrogen_TTLTH_get_data
    PROCEDURE(lib_compute_transpCoeff) :: Nitrogen_TTLTH_compute_transpCoeff
    PROCEDURE(lib_compute_species_DiffFlux) :: Nitrogen_TTLTH_compute_species_DiffFlux
    PROCEDURE(lib_compute_eq_composition) :: Nitrogen_TTLTH_compute_eq_composition  
    PROCEDURE(lib_write_fvmcc_solution) :: Nitrogen_TTLTH_write_fvmcc_solution
#endif

    ! Subroutines of Argon_CR library interface
#ifdef ARGON_CR
    PROCEDURE(lib_in) :: Argon_CR_initialize
    PROCEDURE(lib_fin) :: Argon_CR_finalize
    PROCEDURE(lib_get_Ri) :: Argon_CR_get_Ri
    PROCEDURE(lib_get_el_data) :: Argon_CR_get_el_data
    PROCEDURE(lib_get_density) :: Argon_CR_get_density 
    PROCEDURE(lib_get_pressure) :: Argon_CR_get_pressure
    PROCEDURE(lib_get_mass_fractions) :: Argon_CR_get_mass_fractions
    PROCEDURE(lib_get_molar_fractions) :: Argon_CR_get_molar_fractions
    PROCEDURE(lib_get_mass_fractions_from_molar_fractions) :: Argon_CR_mass_frac_from_molar_frac
    PROCEDURE(lib_comp_tol) :: Argon_CR_comp_tol 
    PROCEDURE(lib_compute_source_term_kinetics) :: Argon_CR_get_source
    PROCEDURE(lib_compute_source_term_kinetics_Jac) :: Argon_CR_get_source_Jac
    PROCEDURE(lib_compute_temperatures) :: Argon_CR_get_temperatures
    PROCEDURE(lib_compute_energy_densities) :: Argon_CR_get_energy_densities
    PROCEDURE(lib_compute_thermodynamic_data) :: Argon_CR_get_thermodynamic_data
    PROCEDURE(lib_compute_data) :: Argon_CR_get_data
    PROCEDURE(lib_compute_transpCoeff) :: Argon_CR_compute_transpCoeff
    PROCEDURE(lib_compute_species_DiffFlux) :: Argon_CR_compute_species_DiffFlux  
    PROCEDURE(lib_compute_eq_composition) :: Argon_CR_compute_eq_composition  
    PROCEDURE(lib_write_fvmcc_solution) :: Argon_CR_write_fvmcc_solution
#endif 

    ! Subroutines of Mutationpp library interface
#ifdef MUTATIONPP
    PROCEDURE(lib_in) :: Mutationpp_initialize
    PROCEDURE(lib_fin) :: Mutationpp_finalize
    PROCEDURE(lib_get_Ri) :: Mutationpp_get_Ri
    PROCEDURE(lib_get_density) :: Mutationpp_get_density 
    PROCEDURE(lib_get_mass_fractions) :: Mutationpp_get_mass_fractions
    PROCEDURE(lib_get_molar_fractions) :: Mutationpp_get_molar_fractions
    PROCEDURE(lib_get_species_enthalpies) :: Mutationpp_get_species_enthalpy
    PROCEDURE(lib_get_enthalpy) :: Mutationpp_get_enthalpy
    PROCEDURE(lib_get_mass_fractions_from_molar_fractions) :: Mutationpp_mass_frac_from_molar_frac
    PROCEDURE(lib_compute_source_term_kinetics) :: Mutationpp_get_source
    PROCEDURE(lib_compute_source_term_kinetics_Jac) :: Mutationpp_get_source_Jac
    PROCEDURE(lib_compute_temperatures) :: Mutationpp_get_temperatures
    PROCEDURE(lib_get_prandtl_number) :: Mutationpp_get_prandtl_number
    PROCEDURE(lib_compute_energy_densities) :: Mutationpp_get_energy_densities
    PROCEDURE(lib_compute_thermodynamic_data) :: Mutationpp_get_thermodynamic_data
    PROCEDURE(lib_compute_data) :: Mutationpp_get_data
    PROCEDURE(lib_compute_transpCoeff) :: Mutationpp_compute_transpCoeff
    PROCEDURE(lib_compute_species_DiffFlux) :: Mutationpp_compute_species_DiffFlux 
    PROCEDURE(lib_compute_eq_composition) :: Mutationpp_compute_eq_composition  
    PROCEDURE(lib_compute_thermal_cond) :: Mutationpp_compute_thermal_cond
    PROCEDURE(lib_write_fvmcc_solution) :: Mutationpp_write_fvmcc_solution
    PROCEDURE(lib_compute_eq_composition_pyro) :: Mutationpp_compute_eq_composition_pyro
    PROCEDURE(lib_get_species_index) :: Mutationpp_get_species_index
    PROCEDURE(lib_get_species_name) :: Mutationpp_get_species_name
    PROCEDURE(lib_compute_species_name) :: Mutationpp_species_name
    PROCEDURE(lib_number_of_elements) :: Mutationpp_number_of_elements
    PROCEDURE(lib_get_molar_fractions_from_mass_fractions) :: Mutationpp_molar_frac_from_mass_frac 
    PROCEDURE(lib_compute_eq_composition_meteor) :: Mutationpp_compute_eq_composition_meteor 
    PROCEDURE(lib_compute_elemental_mole_fraction) :: Mutationpp_compute_elemental_mole_fraction
    PROCEDURE(lib_compute_surface_mass_balance_constraint) :: Mutationpp_surface_mass_balance_constraint
    PROCEDURE(lib_get_element_mass_fractions) :: Mutationpp_get_element_mass_fractions
    PROCEDURE(lib_set_diff_model) :: Mutationpp_set_diff_model
    PROCEDURE(lib_set_cond_heat_flux) :: Mutationpp_set_cond_heat_flux
    PROCEDURE(lib_set_wall_radiation) :: Mutationpp_set_wall_radiation
    PROCEDURE(lib_set_wall_state) :: Mutationpp_set_wall_state
    PROCEDURE(lib_solve_surface_balance) :: Mutationpp_solve_surface_balance
    PROCEDURE(lib_get_wall_state) :: Mutationpp_get_wall_state
    PROCEDURE(lib_get_surface_production_rates) :: Mutationpp_get_surface_production_rates
    PROCEDURE(lib_get_mass_blowing_rate) :: Mutationpp_get_mass_blowing_rate
    PROCEDURE(lib_get_frozen_sound_speed) :: Mutationpp_get_frozen_sound_speed


#endif

#ifdef POLYAT_GAS_WCU
    PROCEDURE(lib_in) :: Polyat_gas_WCU_initialize
    PROCEDURE(lib_fin) :: Polyat_gas_WCU_finalize
    PROCEDURE(lib_get_Ri) :: Polyat_gas_WCU_get_Ri
    PROCEDURE(lib_get_density) :: Polyat_gas_WCU_get_density 
    PROCEDURE(lib_get_pressure) :: Polyat_gas_WCU_get_pressure
    PROCEDURE(lib_get_mass_fractions) :: Polyat_gas_WCU_get_mass_fractions
    PROCEDURE(lib_get_molar_fractions) :: Polyat_gas_WCU_get_molar_fractions
    PROCEDURE(lib_get_mass_fractions_from_molar_fractions) :: Polyat_gas_WCU_mass_frac_from_molar_frac
    PROCEDURE(lib_comp_tol) :: Polyat_gas_WCU_comp_tol
    PROCEDURE(lib_compute_source_term_kinetics) :: Polyat_gas_WCU_get_source
    PROCEDURE(lib_compute_source_term_kinetics_Jac) :: Polyat_gas_WCU_get_source_Jac
    PROCEDURE(lib_compute_temperatures) :: Polyat_gas_WCU_get_temperatures
    PROCEDURE(lib_compute_energy_densities) :: Polyat_gas_WCU_get_energy_densities
    PROCEDURE(lib_compute_thermodynamic_data) :: Polyat_gas_WCU_get_thermodynamic_data
    PROCEDURE(lib_compute_data) :: Polyat_gas_WCU_get_data
    PROCEDURE(lib_compute_transpCoeff) :: Polyat_gas_WCU_compute_transpCoeff
    PROCEDURE(lib_compute_species_DiffFlux) :: Polyat_gas_WCU_compute_species_DiffFlux
    PROCEDURE(lib_compute_eq_composition) :: Polyat_gas_WCU_compute_eq_composition
    PROCEDURE(lib_write_fvmcc_solution) :: Polyat_gas_WCU_write_fvmcc_solution
#endif

  END MODULE mod_neq_function_pointer
!------------------------------------------------------------------------------!
