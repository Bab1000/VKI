!------------------------------------------------------------------------------!
! This module implements an interface for the Polyat_gas_WCU thermodynamic library. It provides an
! implementation of functions and subroutines used to compute thermodynamic and transport properties and 
! source terms related to collisional and radiative processes.
! The interface is written such to have no data and variable dependency from the flow solver being used. 
! Therefore it can plugged to any code without big efforts. 
  MODULE mod_Polyat_gas_WCU_library 

#include"../config.h"

#ifdef POLYAT_GAS_WCU
    IMPLICIT NONE
  
    INTEGER, SAVE :: nb_ns, nb_tvib, nb_te, nb_trot, nb_temp, nb_int_temp, nb_eq, nb_dim, posTr, posTv 
    REAL(KIND=8), PARAMETER ::  kb = 1.380658d-23
    REAL(KIND=8), PARAMETER ::  na = 6.0221367d23   
    REAL(KIND=8), PARAMETER ::  ru = kb*na
    REAL(KIND=8), SAVE :: xi_tol
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: mass, Ri

    ! Subroutines of Polyat_gas_WCU thermodynamic library interface. In some cases they call 
    ! other subroutines for the model under use (even these are not dependent on flow solver data). 
    CONTAINS
  
      !------------------------------------------------------!
      ! This subroutine initializes the thermodynamic library in use.  
      SUBROUTINE Polyat_gas_WCU_initialize (input_ns, input_ntrot, input_ntvib, input_nte, input_ntemp, input_neq,  & 
                                         & input_ndim, input_solver, input_mixture, input_reaction, state_model,    &
                                         & input_transf, input_path, input_xi_tol) 
  
        USE mod_polyat_gas_wcu_initialize_CFD,            ONLY: initialize
 
        INTEGER :: is, nt, neq, length
        
        INTEGER, INTENT(IN) :: input_ns, input_ntrot, input_ntvib, input_nte, input_ntemp, input_neq, input_ndim
        REAL(KIND=8), INTENT(IN) :: input_xi_tol
        CHARACTER*(*), INTENT(IN) :: input_solver, input_mixture, input_reaction, state_model, input_transf, input_path
         
        ! Parameters passed from the solver
        nb_ns   = input_ns
        nb_trot = input_ntrot  
        nb_tvib = input_ntvib
        nb_te   = input_nte
        nb_temp = input_ntemp
        nb_eq   = input_neq
        nb_dim  = input_ndim

        nb_int_temp = nb_temp - (nb_te + 1)

        length  = LEN_TRIM(input_solver)
 
        ! Consistency checking on number of temperatures and conservation equations
        nt = 1 + input_ntrot + input_ntvib + input_nte
        IF (nt.NE.nb_temp) THEN
           PRINT*
           WRITE(*,5)input_solver(1:length),'::Polyat_gas_WCU_library interface'   
           PRINT*
           WRITE(*,5)input_solver(1:length),'::Error in number of temperatures'     
           PRINT*
           STOP
        ENDIF

        neq = input_ns + 1 + input_ntrot + input_ntvib + input_nte + input_ndim
        IF (neq.NE.nb_eq) THEN
           PRINT*
           WRITE(*,5)input_solver(1:length),'::Polyat_gas_WCU_library interface'   
           PRINT*
           WRITE(*,5)input_solver(1:length),'::Error in number of temperatures'     
           PRINT*
           STOP
        ENDIF
 
        ! Tolerance on species molar fractions
        xi_tol = input_xi_tol

        ! Allocation of molecular mass and specific gas constant vectors
        ALLOCATE(mass(nb_ns), Ri(nb_ns))

        ! Polyat_gas_WCU library is initialized
        CALL initialize (input_solver, input_mixture, input_transf, input_reaction, input_path, & 
                       & nb_ns, nb_trot, nb_tvib, nb_te, nb_eq, nb_dim, mass)
  
        ! Speficic gas constant of single species
        DO is = 1,nb_ns 
           Ri(is) = 1.d0/mass(is)
        ENDDO 
        Ri = ru*Ri
       
5     FORMAT(A,A)

      END SUBROUTINE Polyat_gas_WCU_initialize

      !------------------------------------------------------!
      ! This subroutine finalizes the Nitrogen DSMC thermodynamic library
      SUBROUTINE Polyat_gas_WCU_finalize

         USE mod_polyat_gas_wcu_initialize_CFD,            ONLY: finalize

         IF (ALLOCATED(mass))       DEALLOCATE(mass)
         IF (ALLOCATED(Ri))         DEALLOCATE(Ri)

         CALL finalize ()

      END SUBROUTINE Polyat_gas_WCU_finalize

      !------------------------------------------------------!
      ! This subroutine computes post-shock equilibrium conditions
      SUBROUTINE Polyat_gas_WCU_post_shock_eq (p1, u1, T1, p2, u2, T2, yi, rhoi, xN)

        USE mod_polyat_gas_wcu_initialize_CFD,            ONLY: Rgas
        USE mod_polyat_gas_wcu_CFD_prop,                  ONLY: get_post_shock

        REAL(KIND=8), INTENT(IN) :: p1, u1, T1
        REAL(KIND=8), INTENT(IN), OPTIONAL :: xN
        REAL(KIND=8), INTENT(OUT) :: p2, u2, T2
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi, rhoi

        ! Compute post-shock conditions
        CALL get_post_shock (p1, u1, T1, p2, u2, T2)

        ! Species mass fractions and densities 
        yi = 1.d0
        rhoi = p2/(Rgas*T2)

      END SUBROUTINE Polyat_gas_WCU_post_shock_eq

      !------------------------------------------------------!
      ! This subroutine write down the solution when solving the ODE system for
      ! studying the nonequilibrium flow behind a shock wave
      SUBROUTINE Polyat_gas_WCU_write_ODE_solution (flow_file, pop_file, k, xold, x, uold, m_dot, z)

        INTEGER :: is
        REAL(KIND=8) :: mm
        REAL(KIND=8) :: u, p, rho
        REAL(KIND=8), SAVE :: time
        REAL(KIND=8), DIMENSION(nb_ns) :: xi
        REAL(KIND=8), DIMENSION(nb_temp) :: temp

        INTEGER, INTENT(IN) :: flow_file, pop_file, k
        REAL(KIND=8), INTENT(IN) :: xold, x, uold, m_dot
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: z

        ! Species molar fractions
        mm = 0.d0
        DO is = 1,nb_ns 
           mm = mm + z(is)/mass(is)
        ENDDO 
        mm = 1.d0/mm

        DO is = 1,nb_ns 
           xi(is) = z(is)/mass(is) 
        ENDDO
        xi = xi*mm

        ! Density, velocity and temperature
        u   = z(nb_ns + 1)
        rho = m_dot/u

        DO is = 1,nb_temp
           temp(is) = z(nb_ns + 1 + is)
        ENDDO

        ! Pressure 
        p = 0.d0
        DO is = 1,nb_ns 
           p = p + z(is)/mass(is)
        ENDDO
        p = p*ru*rho*temp(1)

        ! Lagrangian time 
        IF (k.EQ.1) THEN 
           time = 0.d0
        ELSE 
           time = time + 0.5d0*(x - xold)*(1.d0/uold + 1.d0/u)
        ENDIF

        ! Flowfield output file
        WRITE(flow_file,10)x,time,temp,xi,u,rho,p
        CALL flush(flow_file)

10    FORMAT(100E20.10)

      END SUBROUTINE Polyat_gas_WCU_write_ODE_solution 

      !------------------------------------------------------!
      ! This subroutine writes down the solution when solving the flow equations 
      ! by means of Finite Volume method
      SUBROUTINE Polyat_gas_WCU_write_fvmcc_solution (flow_file, pop_file, k, xold, x, uold, u, rhoi, temp)

        INTEGER :: is
        REAL(KIND=8) :: mm, fac
        REAL(KIND=8) :: p, rho
        REAL(KIND=8), SAVE :: time
        REAL(KIND=8), DIMENSION(nb_ns) :: xi, yi

        INTEGER, INTENT(IN) :: flow_file, pop_file, k
        REAL(KIND=8), INTENT(IN) :: xold, x, uold, u
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp

        ! Density 
        CALL Polyat_gas_WCU_get_density (rhoi, rho)

        ! Species mass fractions
        CALL Polyat_gas_WCU_get_mass_fractions (rho, rhoi, yi)
 
        ! Species molar fractions 
        mm = 0.d0
        DO is = 1,nb_ns 
           mm = mm + yi(is)/mass(is)
        ENDDO 
        mm = 1.d0/mm

        DO is = 1,nb_ns 
           xi(is) = yi(is)/mass(is) 
        ENDDO
        xi = xi*mm

        ! Pressure 
        p = 0.d0
        DO is = 1,nb_ns 
           p = p + yi(is)/mass(is)
        ENDDO
        p = p*ru*rho*temp(1)

        ! Lagrangian time 
        IF (k.EQ.1) THEN 
           time = 0.d0
        ELSE 
           time = time + 0.5d0*(x - xold)*(1.d0/uold + 1.d0/u)
        ENDIF

        ! Flowfield output file
        WRITE(flow_file,10)x,time,temp,xi,u,rho,p
        CALL flush(flow_file)

10    FORMAT(100E20.10)

      END SUBROUTINE Polyat_gas_WCU_write_fvmcc_solution   

      !------------------------------------------------------!
      ! This subroutine computes the equilibrium composition given pressure and temperature
      SUBROUTINE Polyat_gas_WCU_compute_eq_composition (p, T, rhoi)

        USE mod_polyat_gas_wcu_initialize_CFD,    ONLY: Rgas

        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: rhoi

        rhoi = p/(Rgas*T)

      END SUBROUTINE Polyat_gas_WCU_compute_eq_composition

      !------------------------------------------------------!
      ! This subroutine computes the N-N2 equilibrium internal energy per unit mass 
      SUBROUTINE Polyat_gas_WCU_compute_eq_species_energy (T, e)

        USE mod_polyat_gas_wcu_CFD_prop,          ONLY: get_energy

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: e

        CALL get_energy(T, e(1))

      END SUBROUTINE Polyat_gas_WCU_compute_eq_species_energy 

      !------------------------------------------------------!
      ! This subroutine computes the N-N2 equilibrium enthalpy per unit mass 
      SUBROUTINE Polyat_gas_WCU_compute_eq_species_enthalpy (T, h)

        USE mod_polyat_gas_wcu_initialize_CFD,    ONLY: Rgas
        USE mod_polyat_gas_wcu_CFD_prop,          ONLY: get_energy

        REAL(KIND=8) :: e

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: h

        CALL get_energy(T, e)
        h = e + Rgas*T

      END SUBROUTINE Polyat_gas_WCU_compute_eq_species_enthalpy 

      !------------------------------------------------------!
      ! This subroutine computes the temperatures and thermodynamic data,
      ! additional data for Jacobians are given in output
      SUBROUTINE Polyat_gas_WCU_get_data (rho, rhoi, rho_eint, temp, c, gamma, p, alpha, beta, ei)

        USE mod_polyat_gas_wcu_initialize_CFD,    ONLY: Rgas
        USE mod_polyat_gas_wcu_CFD_prop,          ONLY: get_cv, get_temperature
 
        IMPLICIT NONE
  
        REAL(KIND=8) :: cv, e, T

        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, rho_eint
        REAL(KIND=8), INTENT(OUT) :: c, gamma, p, alpha
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp, beta, ei
  
        ! Temperature computation
        e = rho_eint(1)/rho
        CALL get_temperature(e, T)
        temp = T
       
        ! Constant volume specific heat
        CALL get_cv (T, cv)  
        
        gamma = (cv + Rgas)/cv
        c     = DSQRT(gamma*Rgas*T)
        p     = rho*Rgas*T
        alpha = gamma - 1.d0
        beta  = rho*cv
        ei    = e
       
      END SUBROUTINE Polyat_gas_WCU_get_data
  
      !------------------------------------------------------!
      ! This subroutine computes thermodynamic data needed to be stored
      ! for nonequilibrium gas flow computations. 
      SUBROUTINE Polyat_gas_WCU_get_thermodynamic_data (rho, rhoi, temp, c, gamma, p, alpha, beta, ei)
  
        USE mod_polyat_gas_wcu_initialize_CFD,    ONLY: Rgas
        USE mod_polyat_gas_wcu_CFD_prop,          ONLY: get_cv, get_energy 

        REAL(KIND=8) :: cv, T
  
        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), INTENT(OUT) :: c, gamma, p, alpha
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: beta, ei

        ! Temperature   
        T = temp(1)
        
        ! Specific energy and constant volume specific heat
        CALL get_energy(T, ei(1))
        CALL get_cv (T, cv)  
        
        gamma = (cv + Rgas)/cv
        c     = DSQRT(gamma*Rgas*T)
        p     = rho*Rgas*T
        alpha = gamma - 1.d0
        beta  = rho*cv
        
      END SUBROUTINE Polyat_gas_WCU_get_thermodynamic_data
 
      !------------------------------------------------------!
      ! This subroutine computes the energy density for total internal,
      ! rotational, vibrational and electronic components.
      SUBROUTINE Polyat_gas_WCU_get_energy_densities (rhoi, temp, rho_e)
 
        USE mod_polyat_gas_wcu_CFD_prop,          ONLY: get_energy

        REAL(KIND=8) :: e

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: rho_e

        CALL get_energy(temp(1), e)
        rho_e = rhoi(1)*e 
        
      END SUBROUTINE Polyat_gas_WCU_get_energy_densities
  
      !------------------------------------------------------!
      ! This subroutine computes the total internal, rotational, vibrational and electronic energy 
      ! per unit mass of the mixture.
      SUBROUTINE Polyat_gas_WCU_get_energy (rho, rhoi, temp, e)
 
        USE mod_polyat_gas_wcu_CFD_prop,          ONLY: get_energy

        REAL(KIND=8) :: T

        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e
       
        ! Temperature
        T = temp(1)

        CALL get_energy(T, e(1)) 

      END SUBROUTINE Polyat_gas_WCU_get_energy
  
      !------------------------------------------------------!
      ! This subroutine computes the total internal, rotational, vibrational and electronic enthalpy 
      ! per unit mass of the mixture.
      SUBROUTINE Polyat_gas_WCU_get_enthalpy (rho, rhoi, temp, h)
  
        USE mod_polyat_gas_wcu_initialize_CFD,    ONLY: Rgas
        USE mod_polyat_gas_wcu_CFD_prop,          ONLY: get_energy 
  
        REAL(KIND=8) :: e, T
  
        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: h
 
        ! Temperature
        T = temp(1)

        CALL get_energy(T, e) 
        h = e + Rgas*T
       
      END SUBROUTINE Polyat_gas_WCU_get_enthalpy 
 
      !----------------------------------------------------!
      ! This subroutine computes the species internal energies and constant 
      ! volume specific heat per unit mass
      SUBROUTINE Polyat_gas_WCU_get_species_energy (temp, e) 
 
        USE mod_polyat_gas_wcu_CFD_prop,          ONLY: get_energy
        
        REAL(KIND=8) :: T

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e

        ! Temperature
        T = temp(1)

        CALL get_energy(T, e(1)) 

      END SUBROUTINE Polyat_gas_WCU_get_species_energy

      !----------------------------------------------------!
      ! This subroutine computes the species internal energies heat per unit mass
      SUBROUTINE Polyat_gas_WCU_get_species_energy_cv (temp, e, cv) 
 
        USE mod_polyat_gas_wcu_CFD_prop,          ONLY: get_cv, get_energy
        
        REAL(KIND=8) :: T

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e, cv

        ! Temperature
        T = temp(1)

        CALL get_energy(T, e(1))  
        CALL get_cv(T, cv(1))

      END SUBROUTINE Polyat_gas_WCU_get_species_energy_cv

      !----------------------------------------------------!
      ! This subroutine computes the species internal energies and constant 
      ! volume specific heat per unit mass
      SUBROUTINE Polyat_gas_WCU_get_species_enthalpy (temp, h) 
 
        USE mod_polyat_gas_wcu_initialize_CFD,    ONLY: Rgas
        USE mod_polyat_gas_wcu_CFD_prop,          ONLY: get_energy 

        REAL(KIND=8) :: e, T

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: h

        ! Temperature
        T = temp(1)

        ! Specific energy and constant volume specific heat
        CALL get_energy(T, e)  

        h = e + Rgas*T

      END SUBROUTINE Polyat_gas_WCU_get_species_enthalpy

      !----------------------------------------------------!
      ! This subroutine computes the species enthalpies and constant 
      ! pressure specific heat per unit mass
      SUBROUTINE Polyat_gas_WCU_get_species_enthalpy_cp (temp, h, cp) 

        USE mod_polyat_gas_wcu_initialize_CFD,    ONLY: Rgas
        USE mod_polyat_gas_wcu_CFD_prop,          ONLY: get_cv, get_energy 

        REAL(KIND=8) :: cv, e, T

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: h, cp

        ! Temperature
        T = temp(1)

        ! Specific energy and constant volume specific heat
        CALL get_energy(T, e)  
        CALL get_cv(T, cv)

        h  = e  + Rgas*T
        cp = cv + Rgas

      END SUBROUTINE Polyat_gas_WCU_get_species_enthalpy_cp

      !------------------------------------------------------!
      ! This subroutine computes the frozen specific heats at constant volume of each species
      SUBROUTINE Polyat_gas_WCU_get_species_frozen_cv (T, cv)  
 
        USE mod_polyat_gas_wcu_CFD_prop,          ONLY: get_cv 

        REAL(KIND=8), INTENT(IN) :: T 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv
  
        CALL get_cv(T, cv(1))

      END SUBROUTINE Polyat_gas_WCU_get_species_frozen_cv 
  
      !------------------------------------------------------!
      ! This subroutine computes the frozen specific heats at constant pressure of each species
      SUBROUTINE Polyat_gas_WCU_get_species_frozen_cp (T, cp)  
  
        USE mod_polyat_gas_wcu_initialize_CFD,    ONLY: Rgas
        USE mod_polyat_gas_wcu_CFD_prop,          ONLY: get_cv 

        REAL(KIND=8), INTENT(IN) :: T 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cp
  
        CALL get_cv(T, cp(1))
        cp = cp + Rgas

      END SUBROUTINE Polyat_gas_WCU_get_species_frozen_cp

      !------------------------------------------------------!
      ! This subroutine computes the frozen specific heats at constant volume of each species and 
      ! the species energy components in thermal equilibrium with translation.
      SUBROUTINE Polyat_gas_WCU_get_species_frozen_energy_cv (temp, e, cv)  
 
        USE mod_polyat_gas_wcu_CFD_prop,          ONLY: get_cv, get_energy   
 
        REAL(KIND=8) :: T

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv, e
  
        ! Temperature
        T = temp(1)

        ! Specific energy and constant volume specific heat
        CALL get_energy(T, e(1))  
        CALL get_cv(T, cv(1))       

      END SUBROUTINE Polyat_gas_WCU_get_species_frozen_energy_cv 
  
      !------------------------------------------------------!
      ! This subroutine computes the frozen specific heats at constant pressure of each species and 
      ! the species enthalpy components in thermal equilibrium with translation.
      SUBROUTINE Polyat_gas_WCU_get_species_frozen_enthalpy_cp (temp, h, cp)  
  
        USE mod_polyat_gas_wcu_initialize_CFD,    ONLY: Rgas
        USE mod_polyat_gas_wcu_CFD_prop,          ONLY: get_cv, get_energy
        
        REAL(KIND=8) :: cv, e, T 

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cp, h
 
        ! Temperature 
        T = temp(1)

        ! Specific energy and constant volume specific heat
        CALL get_energy(T, e)  
        CALL get_cv(T, cv)

        h  = e  + Rgas*T
        cp = cv + Rgas 

      END SUBROUTINE Polyat_gas_WCU_get_species_frozen_enthalpy_cp

      !------------------------------------------------------!
      ! This subroutine computes the temperatures given the energy densities.
      SUBROUTINE Polyat_gas_WCU_get_temperatures (rhoi, rho_eint, temp)
  
        USE mod_function_pointer_DSMC,    ONLY: get_temp_DSMC
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rho_eint
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp
  
        CALL get_temp_DSMC (rhoi, rho_eint, temp) 
        
      END SUBROUTINE Polyat_gas_WCU_get_temperatures
 
      !------------------------------------------------------!
      ! This subroutine computes the equilibrium speed of sound
      SUBROUTINE Polyat_gas_WCU_get_eq_sound_speed (rhoi, p, T, c)
 
        USE mod_polyat_gas_wcu_initialize_CFD,    ONLY: Rgas
        USE mod_polyat_gas_wcu_CFD_prop,          ONLY: get_cv

        REAL(KIND=8) :: cv, gamma
  
        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), INTENT(IN), DIMENSION(:) :: rhoi
        REAL(KIND=8), INTENT(OUT) :: c
 
        ! Frozen specific heat (constant volume)
        CALL get_cv (T, cv)
        
        gamma = (cv + Rgas)/cv
        c     = DSQRT(gamma*Rgas*T) 

      END SUBROUTINE Polyat_gas_WCU_get_eq_sound_speed 

      !------------------------------------------------------!
      ! This subroutine computes the equilibrium specific heat ratio
      SUBROUTINE Polyat_gas_WCU_get_eq_gamma (rhoi, p, T, gamma)

        USE mod_polyat_gas_wcu_initialize_CFD,    ONLY: Rgas
        USE mod_polyat_gas_wcu_CFD_prop,          ONLY: get_cv

        REAL(KIND=8) :: cv, c
  
        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), INTENT(IN), DIMENSION(:) :: rhoi
        REAL(KIND=8), INTENT(OUT) :: gamma
 
        ! Frozen specific heat (constant volume)
        CALL get_cv (T, cv)
        
        gamma = (cv + Rgas)/cv
        c     = DSQRT(gamma*Rgas*T) 

      END SUBROUTINE Polyat_gas_WCU_get_eq_gamma 

      !------------------------------------------------------!
      ! This subroutine computes the equilibrium specific heat ratio
      SUBROUTINE Polyat_gas_WCU_get_eq_gamma_sound_speed (rhoi, p, T, gamma, c)
 
        USE mod_polyat_gas_wcu_initialize_CFD,    ONLY: Rgas
        USE mod_polyat_gas_wcu_CFD_prop,          ONLY: get_cv

        REAL(KIND=8) :: cv
  
        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), INTENT(IN), DIMENSION(:) :: rhoi
        REAL(KIND=8), INTENT(OUT) :: gamma, c
 
        ! Frozen specific heat (constant volume)
        CALL get_cv (T, cv)
        
        gamma = (cv + Rgas)/cv
        c     = DSQRT(gamma*Rgas*T)

      END SUBROUTINE Polyat_gas_WCU_get_eq_gamma_sound_speed 
 
      !------------------------------------------------------!
      ! This subroutine computes the frozen speed of sound
      SUBROUTINE Polyat_gas_WCU_get_frozen_sound_speed (rhoi, temp, c)
 
        USE mod_polyat_gas_wcu_initialize_CFD,    ONLY: Rgas
        USE mod_polyat_gas_wcu_CFD_prop,          ONLY: get_cv

        REAL(KIND=8) :: cv, gamma, T
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), INTENT(OUT) :: c 
  
        ! Temperature 
        T = temp(1)

        ! Frozen specific heat (constant volume)
        CALL get_cv (T, cv)
        
        gamma = (cv + Rgas)/cv
        c     = DSQRT(gamma*Rgas*T) 
 
      END SUBROUTINE Polyat_gas_WCU_get_frozen_sound_speed
 
      !------------------------------------------------------!
      ! This subroutine gets the frozen speed specific heat ratio.
      SUBROUTINE Polyat_gas_WCU_get_frozen_gamma (rhoi, temp, gamma)
 
        USE mod_polyat_gas_wcu_initialize_CFD,    ONLY: Rgas
        USE mod_polyat_gas_wcu_CFD_prop,          ONLY: get_cv

        REAL(KIND=8) :: cv, T
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), INTENT(OUT) :: gamma
  
        ! Temperature 
        T = temp(1)

        ! Frozen specific heat (constant volume)
        CALL get_cv (T, cv)
        
        gamma = (cv + Rgas)/cv

      END SUBROUTINE Polyat_gas_WCU_get_frozen_gamma
 
      !------------------------------------------------------!
      ! This subroutine gets the frozen speed of sound and specific heat ratio.
      SUBROUTINE Polyat_gas_WCU_get_frozen_gamma_sound_speed (rhoi, temp, gamma, c)
  
        USE mod_polyat_gas_wcu_initialize_CFD,    ONLY: Rgas
        USE mod_polyat_gas_wcu_CFD_prop,          ONLY: get_cv

        REAL(KIND=8) :: cv, T
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), INTENT(OUT) :: gamma, c 
  
        ! Temperature 
        T = temp(1)

        ! Frozen specific heat (constant volume)
        CALL get_cv (T, cv)
        
        gamma = (cv + Rgas)/cv
        c     = DSQRT(gamma*Rgas*T)
  
      END SUBROUTINE Polyat_gas_WCU_get_frozen_gamma_sound_speed
  
      !------------------------------------------------------!
      ! This subroutine computes the mixture static pressure (presence of free electrons is accounted for)
      SUBROUTINE Polyat_gas_WCU_get_pressure (rhoi, temp, p)
   
        INTEGER :: is
        REAL(KIND=8) :: p_h, p_el, T, Te
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), INTENT(OUT) :: p
       
        ! Temperatures
        T  = temp(1)
        Te = temp(SIZE(temp))

        ! Free electron contribution
        p_el = 0.d0
        IF (nb_te.EQ.1) p_el = rhoi(1)*Ri(1)*Te
  
        ! Heavy particle contribution
        p_h = 0.d0   
        DO is = nb_te + 1,nb_ns 
           p_h = p_h + rhoi(is)*Ri(is)
        ENDDO
        p_h = p_h*T
  
        ! Mixture pressure
        p = p_el + p_h
  
      END SUBROUTINE Polyat_gas_WCU_get_pressure
  
      !------------------------------------------------------!
      ! This subroutine computes the mixture density.
      SUBROUTINE Polyat_gas_WCU_get_density (rhoi, rho)
  
        INTEGER :: is
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), INTENT(OUT) :: rho
  
        ! Mixture density
        rho = 0.d0
        DO is = 1,nb_ns 
           rho = rho + rhoi(is)
        ENDDO
  
      END SUBROUTINE Polyat_gas_WCU_get_density
  
      !------------------------------------------------------!
      ! This subroutine computes the specific gas constant R.
      SUBROUTINE Polyat_gas_WCU_get_Rgas (rhoi, R)
  
        INTEGER :: is 
        REAL(KIND=8) :: rho
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), INTENT(OUT) :: R
  
        ! Inizialization
        rho = 0.d0   
        R = 0.d0
        DO is = 1,nb_ns 
           R = R + rhoi(is)*Ri(is)
           rho = rho + rhoi(is) 
        ENDDO
        R = R/rho
  
      END SUBROUTINE Polyat_gas_WCU_get_Rgas
  
      !------------------------------------------------------!
      ! This subroutine gives the species gas constants
      SUBROUTINE Polyat_gas_WCU_get_Ri (out_Ri)
  
        INTEGER :: is
  
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: out_Ri
  
        DO is = 1,nb_ns 
           out_Ri(is) = Ri(is) 
        ENDDO
  
      END SUBROUTINE Polyat_gas_WCU_get_Ri 
 
      !------------------------------------------------------!
      ! This subroutine gives the species molar masses
      SUBROUTINE Polyat_gas_WCU_get_mi (out_mi)
  
        INTEGER :: is
  
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: out_mi
  
        DO is = 1,nb_ns 
           out_mi(is) = mass(is) 
        ENDDO
  
      END SUBROUTINE Polyat_gas_WCU_get_mi

      !------------------------------------------------------!
      ! This subroutine computes the mixture number density 
      SUBROUTINE Polyat_gas_WCU_get_mix_numb_dens (p, T, Te, xi, nb) 

        REAL(KIND=8), INTENT(IN) :: p, T, Te
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi
        REAL(KIND=8) :: nb

        nb = p/(kb*T)

      END SUBROUTINE Polyat_gas_WCU_get_mix_numb_dens 

      !------------------------------------------------------!
      ! This subroutine computes the species mass fractions
      SUBROUTINE Polyat_gas_WCU_get_mass_fractions (rho, rhoi, yi)
  
        INTEGER :: is
  
        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi
  
        DO is = 1,nb_ns 
           yi(is) = rhoi(is)/rho
        ENDDO
  
      END SUBROUTINE Polyat_gas_WCU_get_mass_fractions
  
      !------------------------------------------------------!
      ! This subroutine computes the species molar fractions
      SUBROUTINE Polyat_gas_WCU_get_molar_fractions (rhoi, xi)

        INTEGER :: is
        REAL(KIND=8) :: sum  

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: xi

        sum = 0.d0
        DO is = 1,nb_ns 
           sum = sum + rhoi(is)/mass(is)
        ENDDO
        sum = 1.d0/sum

        DO is = 1,nb_ns 
           xi(is) = (rhoi(is)/mass(is))*sum
        ENDDO 

      END SUBROUTINE Polyat_gas_WCU_get_molar_fractions 

      !------------------------------------------------------!
      ! This subroutine computes the molar fractions from mass fractions. 
      ! The mixture molar mass is also computed and given in output. 
      SUBROUTINE Nitrogen_DMSC_molar_frac_from_mass_frac (yi, mm, xi)

        INTEGER :: is

        REAL(KIND=8), INTENT(OUT) :: mm
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: yi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: xi

        ! Mixture molar mass
        mm = 0.d0
        DO is = 1,nb_ns 
           mm = mm + yi(is)/mass(is) 
        ENDDO
        mm = 1.d0/mm
  
        DO is = 1,nb_ns 
           xi(is) = yi(is)*mm/mass(is)
        ENDDO

      END SUBROUTINE Nitrogen_DMSC_molar_frac_from_mass_frac 
     
      !------------------------------------------------------!
      ! This subroutine computes the mass fractions from molar fractions. 
      ! The mixture molar mass is also computed and given in output.
      SUBROUTINE Polyat_gas_WCU_mass_frac_from_molar_frac (xi, mm, yi)

        INTEGER :: is
        REAL(KIND=8) :: ov_mm

        REAL(KIND=8), INTENT(OUT) :: mm
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi

        ! Mixture molar mass
        mm = 0.d0
        DO is = 1,nb_ns 
           mm = mm + xi(is)*mass(is)
        ENDDO

        ov_mm = 1.d0/mm
        DO is = 1,nb_ns 
           yi(is) = xi(is)*mass(is)*ov_mm
        ENDDO

      END SUBROUTINE Polyat_gas_WCU_mass_frac_from_molar_frac

      !------------------------------------------------------!
      ! This subroutine adds a small number to the composition respecting the mass constraint
      ! (sum of species molar fractions equal to one)
      SUBROUTINE Polyat_gas_WCU_comp_tol(xi)

        INTEGER :: i
        REAL(KIND=8) :: tmp
        REAL(KIND=8) :: sum_xi

        REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: xi

        ! Initialization
        sum_xi = 0.d0
        DO i = 1,nb_ns
           tmp = xi(i) + xi_tol 
           xi(i)  = tmp
           sum_xi = sum_xi + tmp
        ENDDO

        sum_xi = 1.d0/sum_xi
        DO i = 1,nb_ns 
           xi(i) = xi(i)*sum_xi
        ENDDO

      END SUBROUTINE Polyat_gas_WCU_comp_tol

      !------------------------------------------------------!
      ! This subroutine computes the source term due to collisional and radiative processes.
      SUBROUTINE Polyat_gas_WCU_get_source (rhoi, temp, s) 
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s
  
        s = 0.d0
  
      END SUBROUTINE Polyat_gas_WCU_get_source 
  
      !------------------------------------------------------!
      ! This subroutine computes the source term due to collisional and radiative processes 
      ! and its Jacobian (with respect to primitive variables).
      SUBROUTINE Polyat_gas_WCU_get_source_Jac (rhoi, temp, s, js) 
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s, js
       
        s  = 0.d0
        js = 0.d0

      END SUBROUTINE Polyat_gas_WCU_get_source_Jac
 
      !------------------------------------------------------!
      ! This subroutine computes the species mass diffusion flux 
      SUBROUTINE Polyat_gas_WCU_compute_species_DiffFlux (p, T, Te, xi, diff_driv, Ji)

        REAL(KIND=8), INTENT(IN) :: p, T, Te
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, diff_driv        
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: Ji
      
        ! No mass diffusion for a single species gas
        Ji = 0.d0
   
      END SUBROUTINE Polyat_gas_WCU_compute_species_DiffFlux
 
      !------------------------------------------------------!
      ! This subrouotine computes the transport coefficients 
      SUBROUTINE Polyat_gas_WCU_compute_transpCoeff (p, xi, temp, mu, kappa, lambda, Di, chi)

        USE mod_polyat_gas_wcu_CFD_prop,       ONLY: get_transport_coeff 
  
        REAL(KIND=8), INTENT(IN) :: p
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, temp
        REAL(KIND=8), INTENT(OUT) :: kappa, mu
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: chi, Di, lambda

        ! Compute trasport coefficients
        Di  = 0.d0
        chi = 0.d0
        CALL get_transport_coeff (temp(1), mu, kappa, lambda(1)) 

      END SUBROUTINE Polyat_gas_WCU_compute_transpCoeff
#endif
  
  END MODULE mod_Polyat_gas_WCU_library
!------------------------------------------------------------------------------!  
