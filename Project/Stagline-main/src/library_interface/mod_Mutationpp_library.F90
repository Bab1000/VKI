!------------------------------------------------------------------------------!
! This module implements an interface for the Mutation++ thermodynamic library. It provides an
! implementation of functions and subroutines used to compute thermodynamic and transport properties and 
! source terms related to collisional and radiative processes.
! The interface is written such to have no data and variable dependency from the flow solver being used. 
! Therefore it can plugged to any code without big efforts. 
  MODULE mod_Mutationpp_library 

#include "../config.h"

#ifdef MUTATIONPP
    USE mutationpp

    IMPLICIT NONE
  
    INTEGER, SAVE :: nb_ns, nb_tvib, nb_te, nb_trot, nb_temp, nb_int_temp, nb_eq, nb_dim 
    REAL(KIND=8), PARAMETER :: kb = 1.380658d-23
    REAL(KIND=8), PARAMETER :: na = 6.0221367d23   
    REAL(KIND=8), PARAMETER :: ru = kb*na
    REAL(KIND=8), SAVE :: xi_tol
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: mass, Ri
    CHARACTER*80, SAVE :: solver

    ! Subroutines of Mutation++ thermodynamic library interface. In some cases they call 
    ! other subroutines for the model under use (even these are not dependent on flow solver data). 
    CONTAINS
  
      !------------------------------------------------------!
      ! This subroutine initializes the thermodynamic library in use.  
      SUBROUTINE Mutationpp_initialize (input_ns, input_ntrot, input_ntvib, input_nte, input_ntemp, input_neq,        &
                                      & input_ndim, input_solver, input_mixture, input_reaction, input_state_model,   &
                                      &  input_transf, input_path, input_xi_tol) 
  

        INTEGER :: is, nt, neq, length

        INTEGER, INTENT(IN) :: input_ns, input_ntrot, input_ntvib, input_nte, input_ntemp, input_neq, input_ndim
        REAL(KIND=8), INTENT(IN) :: input_xi_tol
        CHARACTER*(*), INTENT(IN) :: input_solver, input_mixture, input_reaction, input_state_model, input_transf, input_path
        
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

        CALL mpp_initialize(input_mixture, input_state_model)

        ! Consistency checking on number of temperatures and conservation equations
        nt = 1 + input_ntrot + input_ntvib + input_nte
        IF (nt.NE.nb_temp) THEN
           PRINT*
           WRITE(*,5)input_solver(1:length),'::Mutationpp_library interface'   
           PRINT*
           WRITE(*,5)input_solver(1:length),'::Error in number of temperatures'     
           PRINT*
           STOP
        ENDIF

        neq = input_ns + 1 + input_ntrot + input_ntvib + input_nte + input_ndim
        IF (neq.NE.nb_eq) THEN
           PRINT*
           WRITE(*,5)input_solver(1:length),'::Mutationpp_library interface'  
           PRINT* 
           WRITE(*,5)input_solver(1:length),'::Error in number of equations'     
           PRINT*
           STOP
        ENDIF

        ! Tolerance on species molar fractions
        xi_tol = input_xi_tol

        ! Allocation of molecular mass and specific gas constant vectors
        ALLOCATE(mass(nb_ns), Ri(nb_ns))
        
        CALL mpp_species_mw(mass)

        ! Speficic gas constant of single species
        DO is = 1,nb_ns 
           Ri(is) = 1.d0/mass(is)
        ENDDO 
        Ri = ru*Ri
       
        print*, "Call set_pointer"
        print*, "Useful indices"
        ! Solver name
        solver = input_solver

5     FORMAT(A,A)

      END SUBROUTINE Mutationpp_initialize

      !------------------------------------------------------!
      ! This subroutine finalizes the Nitrogen Park thermodynamic library
      SUBROUTINE Mutationpp_finalize ()

         IF (ALLOCATED(mass))       DEALLOCATE(mass)
         IF (ALLOCATED(Ri))         DEALLOCATE(Ri)

         CALL mpp_destroy()

      END SUBROUTINE Mutationpp_finalize

      !------------------------------------------------------!
      ! This subroutine computes post-shock nonequilibrium conditions
      !SUBROUTINE Mutationpp_post_shock_neq (p1, u1, T1, p2, u2, T2, yi, extra)

      !INTEGER :: is
      !INTEGER :: state_var = 1

      !REAL(KIND=8), PARAMETER :: tol = 1.d-8
      !REAL(KIND=8) :: tmp
      !REAL(KIND=8) :: f, fp, resR, resT, ratio, rhs, ratio_old, T_old
      !REAL(KIND=8) :: cv, cp, R, rho
      !REAL(KIND=8) :: gamma, gp1, gm1, h1, h2, rho1, rho2, c1, g, M1, M1s, m_dot
      !REAL(KIND=8), DIMENSION(nb_ns) :: cpi, hi, rhoi
      !REAL(KIND=8), DIMENSION(3) :: left, right, res

      !REAL(KIND=8), INTENT(IN) :: p1, u1, T1
      !REAL(KIND=8), INTENT(IN), OPTIONAL :: extra
      !REAL(KIND=8), INTENT(OUT) :: p2, u2, T2
      !REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi 

      ! Not used. Can be removed in the future

      ! Compute the equilibrium composition
      !mpp_equilibrate(T1, p1)
    
      ! Compute the mixture density 
      !rho = mpp_density()
    
      ! Compute mixture constant (R)
      !R = RU/mpp_mixture_mw()

      ! Compute the species mass fractions
      !CALL mpp_y(yi)

      ! Thermal nonequilibrium case (Rankine-Hugoniot jump relations are applied)
      !IF (nb_temp.GT.1) THEN

      ! Compute the mixture frozen gamma (low temperature before the shock) 
      !gamma = mpp_mixture_frozen_gamma()      

      ! Useful factors
      !gp1 = gamma + 1.d0
      !gm1 = gamma - 1.d0
  
      !c1  = SQRT(gamma*R*T1)
      !M1  = u1/c1
      !M1s = M1**2

      ! Frozen post-shock conditions (Rankine-Hugoniot relations are applied)
      ! Pressure, velocity and temperature after the shock
      !p2 = p1*(2.d0*gamma*M1s - gamma + 1.d0)/gp1
      !u2 = u1 - c1*2.d0/gp1*(M1 - 1.d0/M1)
      !T2 = T1*(2.d0*gamma*M1s - gamma + 1.d0)*(gamma - 1.d0 + 2.d0/M1s)/(gp1**2) 

      ! Thermal equilibrium case
      !ELSE

      ! Mass, momentum and energy fluxes (pre-shock conditions)
      !rho1  = p1/(R*T1)
      !m_dot = rho1*u1

      !rhoi = rho1*yi
      !CALL mpp_set_state(rhoi, T1, state_var)

      ! Species specific enthalpies and constant pressure specific heats
      !CALL mpp_species_h_mass(hi)
      !CALL mpp_species_cp_mass(cpi)

      !h1   = 0.d0
      !DO is = 1,nb_ns 
      !    h1 = h1 + yi(is)*hi(is)
      !ENDDO

      !left(1) = m_dot 
      !left(2) = m_dot*u1 + p1 
      !left(3) = h1 + 0.5d0*u1**2
   
      ! Guess value for the density ratio 
      !ratio = 0.1d0 
      !resR  = 1.d0
      !T2    = T1 
      !DO WHILE (resR.GT.tol) 

      !    p2   = p1 + m_dot*u1*(1.d0 - ratio)
      !    rho2 = rho1/ratio
        
      !    rhs  = h1 + 0.5d0*u1**2*(1.d0 - ratio**2)
        
      !    ! Newton loop for post-shock temperature
      !    resT = 1.d0
      !    DO WHILE (resT.GT.tol) 

      !        f  = 0.d0
      !        fp = 0.d0

      !        rhoi = rho2*yi
      !        CALL mpp_set_state(rhoi, T2, state_var)

      !        ! Species specific enthalpies and constant pressure specific heats
      !        CALL mpp_species_h_mass(hi)
      !        CALL mpp_species_cp_mass(cpi)

      !        DO is = 1,nb_ns
      !            tmp = yi(is)
      !            f   = f  + tmp*hi(is)
      !            fp  = fp + tmp*cpi(is)
      !        ENDDO

      !        T_old = T2 

      !        ! Post-shock temperature residual
      !        f  = f - rhs
      !        T2 = T2 - f/fp
      !        resT = ABS(T2 - T_old)/T_old

      !    ENDDO

      !    ratio_old = ratio
        
      !    rho2 = p2/(R*T2)
        
      !    ! Density ratio update and residual
      !    ratio = rho1/rho2
      !    resR = ABS(ratio - ratio_old)/ratio

      !ENDDO

      !u2 = u1*ratio

      ! Mass momentum and energy flux (post-shock)
      !m_dot = rho2*u2
    
      !rhoi = rho2*yi
      !CALL mpp_set_state(rhoi, T2, state_var)

      ! Species specific enthalpies and constant pressure specific heats
      !CALL mpp_species_h_mass(hi)
      !CALL mpp_species_cp_mass(cpi)

      !h2   = 0.d0
      !DO is = 1,nb_ns 
      !    h2 = h2 + yi(is)*hi(is)
      !ENDDO   
 
      !right(1) = m_dot
      !right(2) = m_dot*u2 + p2 
      !right(3) = h2 + 0.5d0*u2**2

      !DO is = 1,3 
      !    res(is) = ABS(right(is) - left(is))/left(is)*100.d0
      !ENDDO
    
      !WRITE(*,5)solver(1:LEN_TRIM(solver)),':: Mutation++ -> thermal equilibrium post-shock conditions'
      !    PRINT*
      !    WRITE(*,10)'Residual on mass, momemtum and energy fluxes:'
      !    PRINT*
      !    WRITE(*,15)'Mass    ',res(1),' [%]'
      !    PRINT*
      !    WRITE(*,15)'Momentum',res(2),' [%]'
      !    PRINT*
      !    WRITE(*,15)'Energy  ',res(3),' [%]'
      !    PRINT* 
          
      !  ENDIF

!5   FORMAT(A,A)
!10  FORMAT(A)
!15  FORMAT(A,E14.6,A)

!      END SUBROUTINE Mutationpp_post_shock_neq

      !------------------------------------------------------!
      ! This subroutine computes post-shock equilibrium conditions
      SUBROUTINE Mutationpp_post_shock_eq (p1, u1, T1, p2, u2, T2, yi, rhoi, extra)

        REAL(KIND=8), INTENT(IN) :: p1, u1, T1
        REAL(KIND=8), INTENT(IN), OPTIONAL :: extra
        REAL(KIND=8), INTENT(OUT) :: p2, u2, T2
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi, rhoi

        p2   = 0.d0
        u2   = 0.d0
        T2   = 0.d0
        yi   = 0.d0
        rhoi = 0.d0

        PRINT*,'in "Mutationpp_post_shock_eq", not implemented yet'
        STOP

      END SUBROUTINE Mutationpp_post_shock_eq

      !------------------------------------------------------!
      ! This subroutine write down the solution when solving the ODE system for
      ! studying the nonequilibrium flow behind a shock wave
      SUBROUTINE Mutationpp_write_ODE_solution (flow_file, pop_file, k, xold, x, uold, m_dot, z)

        INTEGER, INTENT(IN) :: flow_file, pop_file, k
        REAL(KIND=8), INTENT(IN) :: xold, x, uold, m_dot
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: z

        WRITE(flow_file,10)x,z

10    FORMAT(100E14.6)

      END SUBROUTINE Mutationpp_write_ODE_solution 

      !----------------------------------------------------!
      ! This subroutine post-process the ODE solution for the flow 
      ! behind a normal shock wave
      SUBROUTINE Mutationpp_post_process_ODE_solution (post_file, x, ni, temp)  

        INTEGER, INTENT(IN) :: post_file
        REAL(KIND=8), INTENT(IN) :: x
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: ni, temp

        PRINT*,'in "Mutationpp_post_process_ODE_solution", not implementedyet..'
        STOP

      END SUBROUTINE Mutationpp_post_process_ODE_solution

      !------------------------------------------------------!
      ! This subroutine writes down the solution when solving the flow equations 
      ! by means of Finite Volume method
      SUBROUTINE Mutationpp_write_fvmcc_solution (flow_file, pop_file, k, xold, x, uold, u, rhoi, temp)

        INTEGER, INTENT(IN) :: flow_file, pop_file, k
        REAL(KIND=8), INTENT(IN) :: xold, x, uold, u
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp

        PRINT*,'in "Mutationpp_write_fvmcc_solution", not implementedyet..'
        STOP

      END SUBROUTINE Mutationpp_write_fvmcc_solution   

      !------------------------------------------------------!
      ! This subroutine computes the equilibrium composition given pressure and temperature
      SUBROUTINE Mutationpp_compute_eq_composition (p, T, rhoi)

        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: rhoi

        ! Not used. Can be removed in the future

        ! Compute the equilibrium composition 
        !mpp_equilibrate(T, p) 
        
        CALL mpp_species_densities(rhoi)

      END SUBROUTINE Mutationpp_compute_eq_composition
            !------------------------------------------------------!
     ! This subroutine computes the equilibrium composition given pressure and temperature for a given
     ! elemental composition. The output is mole fraction
      SUBROUTINE Mutationpp_compute_eq_composition_meteor (p, T,elements, rhoi)

        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: elements
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: rhoi 

        CALL mpp_pyro_equilibrium_composition(T,p,elements,rhoi)

      END SUBROUTINE Mutationpp_compute_eq_composition_meteor

      !------------------------------------------------------!
     ! This subroutine computes  moles of each element that exist in a mixture with
     ! the given set of species moles.

      SUBROUTINE Mutationpp_compute_elemental_mole_fraction (x_species,x_elements)

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: x_species
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: x_elements

       CALL mpp_convert_x_to_y (x_species, x_species)
       CALL mpp_convert_ys_to_ye(x_species, x_elements)
       CALL mpp_convert_ye_to_xe (x_elements, x_elements)

      END SUBROUTINE Mutationpp_compute_elemental_mole_fraction
      !------------------------------------------------------!
      ! This subroutine computes the species specific equilibrium energy per unit mass 
      SUBROUTINE Mutationpp_compute_eq_species_energy (T, e)

        INTEGER :: is

        REAL(KIND=8), INTENT(IN),  DIMENSION(:) :: T
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: e

        CALL Mutationpp_compute_eq_species_enthalpy(T, e)

        DO is = 1,nb_ns 
           e(is) = e(is) - Ri(is)*T(1)
        ENDDO

      END SUBROUTINE Mutationpp_compute_eq_species_energy 

      !------------------------------------------------------!
      ! This subroutine computes the species specific equilibrium enthalpy per unit mass 
      SUBROUTINE Mutationpp_compute_eq_species_enthalpy (T, h)

        INTEGER :: state_var = 1
        REAL(KIND=8), DIMENSION(nb_ns) :: rhoi

        REAL(KIND=8), INTENT(IN),  DIMENSION(:) :: T
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: h

        ! Set state
        rhoi = 1.d0
        CALL mpp_set_state(rhoi, T, state_var)

        CALL mpp_species_h_mass(h)

      END SUBROUTINE Mutationpp_compute_eq_species_enthalpy 

      !------------------------------------------------------!
      ! This subroutine computes the species specific equilibrium entropy per unit mass. The
      ! entropy of mixing contributions is computed separately. 
      SUBROUTINE Mutationpp_compute_eq_species_entropy (p, T, yi, xi, smix, s)

        REAL(KIND=8), INTENT(IN) :: p
        REAL(KIND=8), INTENT(IN), DIMENSION(:) :: T, xi, yi
        REAL(KIND=8), INTENT(OUT) :: smix
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: s

        s    = 0.d0
        smix = 0.d0

        PRINT*,'in "Mutationpp_compute_eq_species_entropy", not implemented yet..'
        STOP

      END SUBROUTINE Mutationpp_compute_eq_species_entropy    

      !------------------------------------------------------!
      ! This subroutine computes the temperatures and thermodynamic data,
      ! additional data for Jacobians are given in output
      SUBROUTINE Mutationpp_get_data (rho, rhoi, rho_eint, temp, c, gamma, p, alpha, beta, ei)
  
        INTEGER :: state_var, i, is
        REAL(KIND=8), DIMENSION(nb_ns*nb_temp) :: cvi

        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, rho_eint
        REAL(KIND=8), INTENT(OUT) :: c, gamma, p, alpha
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp, beta, ei

        ! Set the state based on species densities and energy density
        state_var = 0
        CALL mpp_set_state(rhoi, rho_eint, state_var)
        
        IF (nb_temp == 1) THEN
          temp    = mpp_mixture_t()
        ELSE
          CALL mpp_get_temperatures(temp)
        ENDIF
        p       = mpp_pressure()
        gamma   = mpp_mixture_frozen_gamma()
        c       = sqrt(gamma*p/rho)
        
        CALL mpp_species_cv_mass(cvi)
        DO i = 1, nb_temp
          beta(i) = sum(rhoi*cvi((i-1)*nb_ns+1:i*nb_ns))
        ENDDO
        alpha   = sum(rhoi*Ri)/beta(1)
        
        CALL mpp_species_e_mass(ei)
        
      END SUBROUTINE Mutationpp_get_data 
  
      !------------------------------------------------------!
      ! This subroutine computes thermodynamic data needed to be stored
      ! for nonequilibrium gas flow computations. 
      SUBROUTINE Mutationpp_get_thermodynamic_data (rho, rhoi, temp, c, gamma, p, alpha, beta, ei)
 
        INTEGER :: state_var, i, is
        REAL(KIND=8), DIMENSION(nb_ns*nb_temp) :: cvi

        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), INTENT(OUT) :: c, gamma, p, alpha
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: beta, ei

        ! Set state 
        state_var = 1
        CALL mpp_set_state(rhoi, temp, state_var)

        p       = mpp_pressure()
        gamma   = mpp_mixture_frozen_gamma()
        c       = sqrt(gamma*p/rho)
        
        CALL mpp_species_cv_mass(cvi)
        DO i = 1, nb_temp
          beta(i) = sum(rhoi*cvi((i-1)*nb_ns+1:i*nb_ns))
        ENDDO
        alpha   = sum(rhoi*Ri)/beta(1)
        
        CALL mpp_species_e_mass(ei)

      END SUBROUTINE Mutationpp_get_thermodynamic_data
 
      !------------------------------------------------------!
      ! This subroutine computes the energy density for total internal,
      ! rotational, vibrational and electronic components.
      SUBROUTINE Mutationpp_get_energy_densities (rhoi, temp, rho_e)

        INTEGER :: state_var, itemp
        REAL(KIND=8), DIMENSION(nb_ns*nb_temp) :: ei

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: rho_e

        ! Set state of the mixture and then get rho*e
        state_var = 1
        CALL mpp_set_state(rhoi, temp, state_var)
        CALL mpp_species_e_mass(ei)
        DO itemp = 1, nb_temp
          rho_e(itemp) = sum(rhoi*ei((itemp-1)*nb_ns+1:itemp*nb_ns))
        ENDDO

      END SUBROUTINE Mutationpp_get_energy_densities
  
      !------------------------------------------------------!
      ! This subroutine computes the total internal, rotational, vibrational and electronic energy 
      ! per unit mass of the mixture.
      SUBROUTINE Mutationpp_get_energy (rho, rhoi, temp, e)
  
        INTEGER :: state_var

        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e

        ! Set state 
        state_var = 1
        CALL mpp_set_state(rhoi, temp, state_var)
        e = mpp_mixture_e_mass()  

      END SUBROUTINE Mutationpp_get_energy
  
      !------------------------------------------------------!
      ! This subroutine computes the total internal, rotational, vibrational and electronic enthalpy 
      ! per unit mass of the mixture.
      SUBROUTINE Mutationpp_get_enthalpy (rho, rhoi, temp, h)
 
        INTEGER :: state_var

        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: h

        ! Set state 
        state_var = 1
        CALL mpp_set_state(rhoi, temp, state_var)
        h = mpp_mixture_h_mass()       
 
      END SUBROUTINE Mutationpp_get_enthalpy 
 
      !----------------------------------------------------!
      ! This subroutine computes the species internal energies per unit mass
      SUBROUTINE Mutationpp_get_species_energy (temp, e) 

        INTEGER :: state_var, is
        REAL(KIND=8), DIMENSION(nb_ns) :: rhoi

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e

        ! Set state 
        rhoi = 1.d0
        state_var = 1
        CALL mpp_set_state(rhoi, temp, state_var)

        ! Species specific enthalpies and constant pressure specific heats
        CALL mpp_species_e_mass(e)

      END SUBROUTINE Mutationpp_get_species_energy

      !----------------------------------------------------!
      ! This subroutine computes the species internal energies and constant volume specific heats 
      ! per unit mass
      SUBROUTINE Mutationpp_get_species_energy_cv (temp, e, cv) 

        INTEGER :: state_var, is
        REAL(KIND=8), DIMENSION(nb_ns) :: rhoi

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e, cv

        ! Set state 
        rhoi = 1.d0
        state_var = 1
        CALL mpp_set_state(rhoi, temp, state_var)

        ! Species specific energies and constant volume specific heats
        CALL mpp_species_e_mass(e)
        CALL mpp_species_cv_mass(cv)

      END SUBROUTINE Mutationpp_get_species_energy_cv

      !----------------------------------------------------!
      ! This subroutine computes the species enthalpies per unit mass
      SUBROUTINE Mutationpp_get_species_enthalpy (temp, h) 

        INTEGER :: state_var 
        REAL(KIND=8), DIMENSION(nb_ns) :: rhoi

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: h

        ! Set state
        rhoi = 1.d0
        state_var = 1
        CALL mpp_set_state(rhoi, temp, state_var)
        
        CALL mpp_species_h_mass(h)

      END SUBROUTINE Mutationpp_get_species_enthalpy

      !----------------------------------------------------!
      ! This subroutine computes the species enthalpies and constant 
      ! pressure specific heat per unit mass
      SUBROUTINE Mutationpp_get_species_enthalpy_cp (temp, h, cp) 

        INTEGER :: state_var

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: h, cp
        
        REAL(KIND=8), DIMENSION(nb_ns) :: rhoi

        ! Set state
        rhoi = 1.0 
        CALL mpp_set_state(rhoi, temp, state_var)
        
        CALL mpp_species_h_mass(h)
        CALL mpp_species_cp_mass(cp)

      END SUBROUTINE Mutationpp_get_species_enthalpy_cp

      !------------------------------------------------------!
      ! This subroutine computes the frozen specific heats at constant volume of each species
      SUBROUTINE Mutationpp_get_species_frozen_cv (T, cv)  
 
        REAL(KIND=8), INTENT(IN) :: T 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv
  
        cv = 0.d0

        PRINT*,'in " Mutationpp_get_species_frozen_cv", not implemented yet...'
        STOP

      END SUBROUTINE Mutationpp_get_species_frozen_cv 
  
      !------------------------------------------------------!
      ! This subroutine computes the frozen specific heats at constant pressure of each species
      SUBROUTINE Mutationpp_get_species_frozen_cp (T, cp)  
  
        REAL(KIND=8), INTENT(IN) :: T 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cp
  
        cp = 0.d0

        PRINT*,'in " Mutationpp_get_species_frozen_cv", not implemented yet...'
        STOP

      END SUBROUTINE Mutationpp_get_species_frozen_cp

      !------------------------------------------------------!
      ! This subroutine computes the frozen specific heats at constant volume of each species and 
      ! the species energy components in thermal equilibrium with translation.
      SUBROUTINE Mutationpp_get_species_frozen_energy_cv (T, cv, e)  
 
        REAL(KIND=8), INTENT(IN) :: T 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv, e
 
        cv = 0.d0
        e  = 0.d0 

        PRINT*,'in " Mutationpp_get_species_frozen_energy_cv", not implemented yet...'
        STOP

      END SUBROUTINE Mutationpp_get_species_frozen_energy_cv 
  
      !------------------------------------------------------!
      ! This subroutine computes the frozen specific heats at constant pressure of each species and 
      ! the species enthalpy components in thermal equilibrium with translation.
      SUBROUTINE Mutationpp_get_species_frozen_enthalpy_cp (T, cp, h)  
  
        REAL(KIND=8), INTENT(IN) :: T 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cp, h

        cp = 0.d0
        h  = 0.d0  

        PRINT*,'in " Mutationpp_get_species_frozen_enthalpy_cp", not implemented yet...'
        STOP 

      END SUBROUTINE Mutationpp_get_species_frozen_enthalpy_cp

      !------------------------------------------------------!
      ! This subroutine computes the temperatures given the energy densities.
      SUBROUTINE Mutationpp_get_temperatures (rhoi, rho_eint, temp)
  
        INTEGER :: state_var

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rho_eint
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp

        state_var = 0
        call mpp_set_state(rhoi, rho_eint, state_var)
        call mpp_get_temperatures(temp) 

      END SUBROUTINE Mutationpp_get_temperatures
 
      !------------------------------------------------------!
      ! This subroutine computes the equilibrium speed of sound
      SUBROUTINE Mutationpp_get_eq_sound_speed (rhoi, p, T, c)
 
        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), INTENT(IN), DIMENSION(:) :: rhoi
        REAL(KIND=8), INTENT(OUT) :: c

        c = 0.d0

        PRINT*,'in "Mutationpp_get_eq_sound_speed ", not implemented yet..'
        STOP

      END SUBROUTINE Mutationpp_get_eq_sound_speed 

      !------------------------------------------------------!
      ! This subroutine computes the equilibrium specific heat ratio
      SUBROUTINE Mutationpp_get_eq_gamma (rhoi, p, T, gamma)
 
        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), INTENT(IN), DIMENSION(:) :: rhoi
        REAL(KIND=8), INTENT(OUT) :: gamma

        gamma = 0.d0

        PRINT*,'in "Mutationpp_get_eq_gamma", not implemented yet..'
        STOP

      END SUBROUTINE Mutationpp_get_eq_gamma 

      !------------------------------------------------------!
      ! This subroutine computes the equilibrium specific heat ratio
      SUBROUTINE Mutationpp_get_eq_gamma_sound_speed (rhoi, p, T, gamma, c)
 
        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), INTENT(IN), DIMENSION(:) :: rhoi
        REAL(KIND=8), INTENT(OUT) :: gamma, c
 
        gamma = 0.d0
        c     = 0.d0

        PRINT*,'in "Mutationpp_get_eq_gamma_sound_speed",  not implemented yet..'
        STOP

      END SUBROUTINE Mutationpp_get_eq_gamma_sound_speed 

      !------------------------------------------------------!
      ! This subroutine computes the frozen speed of sound
      SUBROUTINE Mutationpp_get_frozen_sound_speed (rhoi, temp, c)
 
        INTEGER :: state_var 

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp, rhoi
        REAL(KIND=8), INTENT(OUT) :: c 

        ! Set state
        state_var = 1
        CALL mpp_set_state(rhoi, temp, state_var)
        c = mpp_mixture_frozen_sound_speed()
 
      END SUBROUTINE Mutationpp_get_frozen_sound_speed
 
      !------------------------------------------------------!
      ! This subroutine gets the frozen speed specific heat ratio.
      SUBROUTINE Mutationpp_get_frozen_gamma (rhoi, temp, gamma)
 
        INTEGER :: state_var 

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp, rhoi
        REAL(KIND=8), INTENT(OUT) :: gamma 

        ! Set state
        state_var = 1
        CALL mpp_set_state(rhoi, temp, state_var)
        gamma = mpp_mixture_frozen_gamma()

      END SUBROUTINE Mutationpp_get_frozen_gamma
            !------------------------------------------------------!
      ! This subroutine gets the frozen speed of sound and specific heat ratio.
      SUBROUTINE Mutationpp_get_prandtl_number (rhoi, temp, pr)
  
        INTEGER :: state_var, is

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp, rhoi
        REAL(KIND=8), INTENT(OUT) :: pr

        REAL(KIND=8) :: mu, cp
        REAL(KIND=8), DIMENSION(nb_temp) :: lambda
        REAL(KIND=8), DIMENSION(nb_ns) :: cpi
        ! Set state
        state_var = 1

        CALL mpp_set_state(rhoi, temp, state_var)
       ! Mixture dynamic viscosity 
        mu = mpp_viscosity()

        ! Mixture total thermal conductivity 
        CALL mpp_frozen_thermal_conductivity(lambda)

        CALL mpp_species_cp_mass(cpi)
        cp = 0.0D0
        DO is = 1,nb_ns 
           cp = cp + cpi(is)*(rhoi(is)/sum(rhoi))
        ENDDO
      
        pr = cp*mu/lambda(1)

      END SUBROUTINE Mutationpp_get_prandtl_number
  

 
      !------------------------------------------------------!
      ! This subroutine gets the frozen speed of sound and specific heat ratio.
      SUBROUTINE Mutationpp_get_frozen_gamma_sound_speed (rhoi, temp, gamma, c)
  
        INTEGER :: state_var 

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp, rhoi
        REAL(KIND=8), INTENT(OUT) :: gamma, c 

        ! Set state
        state_var = 1
        CALL mpp_set_state(rhoi, temp, state_var)
        gamma = mpp_mixture_frozen_gamma()
        c = sqrt(gamma * mpp_pressure() / mpp_density())

      END SUBROUTINE Mutationpp_get_frozen_gamma_sound_speed
  
      !------------------------------------------------------!
      ! This subroutine computes the mixture static pressure (presence of free electrons is accounted for)
      SUBROUTINE Mutationpp_get_pressure (rhoi, temp, p)
   
        INTEGER :: state_var 

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp, rhoi
        REAL(KIND=8), INTENT(OUT) :: p

        ! Set state
        state_var = 1
        CALL mpp_set_state(rhoi, temp, state_var)
        p = mpp_pressure()
  
      END SUBROUTINE Mutationpp_get_pressure
  
      !------------------------------------------------------!
      ! This subroutine computes the mixture density.
      SUBROUTINE Mutationpp_get_density (rhoi, rho)
  
        INTEGER :: is
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), INTENT(OUT) :: rho
  
        ! Mixture density
        rho = 0.d0
        DO is = 1,nb_ns 
           rho = rho + rhoi(is)
        ENDDO
  
      END SUBROUTINE Mutationpp_get_density
  
      !------------------------------------------------------!
      ! This subroutine computes the specific gas constant R.
      SUBROUTINE Mutationpp_get_Rgas (rhoi, R)
  
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
  
      END SUBROUTINE Mutationpp_get_Rgas
  
      !------------------------------------------------------!
      ! This subroutine gives the species gas constants
      SUBROUTINE Mutationpp_get_Ri (out_Ri)

        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: out_Ri
  
        out_Ri = Ri
 
      END SUBROUTINE Mutationpp_get_Ri 
 
      !------------------------------------------------------!
      ! This subroutine gives the species molar masses
      SUBROUTINE Mutationpp_get_mi (out_mi)

        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: out_mi
  
        out_mi = mass
  
      END SUBROUTINE Mutationpp_get_mi
 
      !------------------------------------------------------!
      ! This subroutine computes the species mass fraction
      SUBROUTINE Mutationpp_get_mass_fractions (rho, rhoi, yi)

        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi
  
        yi = rhoi/rho
  
      END SUBROUTINE Mutationpp_get_mass_fractions

      !------------------------------------------------------!
      ! This subroutine computes the species molar fractions
      SUBROUTINE Mutationpp_get_molar_fractions (rhoi, xi)

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: xi

        call mpp_convert_rho_to_x(rhoi, xi)

      END SUBROUTINE Mutationpp_get_molar_fractions
 
      !------------------------------------------------------!
      ! This subroutine computes the species mass fractions to elements mass fractions
      SUBROUTINE Mutationpp_get_element_mass_fractions (yi, ye)

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: yi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ye
       

       call mpp_convert_ys_to_ye(yi, ye)
  
      END SUBROUTINE Mutationpp_get_element_mass_fractions 
 
      !------------------------------------------------------!

      ! This subroutine computes the molar fractions from mass fractions. 
      ! The mixture molar mass is also computed and given in output.
      SUBROUTINE Mutationpp_molar_frac_from_mass_frac (yi, mm, xi)

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

      END SUBROUTINE Mutationpp_molar_frac_from_mass_frac 
     
      !------------------------------------------------------!
      ! This subroutine computes the mass fractions from molar fractions.
      ! The mixture molar mass is also computed and given in output.
      SUBROUTINE Mutationpp_mass_frac_from_molar_frac (xi, mm, yi)

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

      END SUBROUTINE Mutationpp_mass_frac_from_molar_frac 

      !------------------------------------------------------!
      ! This subroutine adds a small number to the composition respecting the mass constraint
      ! (sum of species molar fractions equal to one)
      SUBROUTINE Mutationpp_comp_tol(xi)

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

      END SUBROUTINE Mutationpp_comp_tol
 
      !------------------------------------------------------!
      ! This subroutine computes the source term due to collisional and radiative processes.
      SUBROUTINE Mutationpp_get_source (rhoi, temp, s) 

        INTEGER :: state_var, i
 
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s

        ! Set the current state of the mixture
        state_var = 1
        CALL mpp_set_state(rhoi, temp, state_var)
        
        ! Compute the species production terms
        s = 0.d0
        CALL mpp_net_production_rates(s(1:nb_ns))

        ! Transfer source terms for internal energy equation
        CALL mpp_source_energy_transfer(s(nb_ns+nb_dim+2:nb_eq)) 
 
      END SUBROUTINE Mutationpp_get_source 
  
      !------------------------------------------------------!
      ! This subroutine computes the source term due to collisional and radiative processes 
      ! and its Jacobian (with respect to primitiva variables).
      SUBROUTINE Mutationpp_get_source_Jac (rhoi, temp, s, js)   
   
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s, js
       
        s  = 0.d0
        js = 0.d0

        PRINT*,'in "Mutationpp_get_source_Jac", not implemented yet...'
        STOP

      END SUBROUTINE Mutationpp_get_source_Jac
 
      !------------------------------------------------------!
      ! This subroutine computes the species mass diffusion flux 
      SUBROUTINE Mutationpp_compute_species_DiffFlux (p, T, Te, xi, diff_driv, Ji)

        INTEGER :: d, l, i, r
        INTEGER :: state_var = 1

        REAL(KIND=8) :: ambE, p_ov_T
        REAL(KIND=8), DIMENSION(nb_ns) :: rhoi
        REAL(KIND=8), DIMENSION(nb_temp) :: temp

        REAL(KIND=8), INTENT(IN) :: p, T, Te
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, diff_driv        
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: Ji

        ! Compute the species densities 
        p_ov_T = p/T
        rhoi(1) = p*xi(1)/(Ri(1)*Te)
        DO i = 2,nb_ns 
           rhoi(i) = p_ov_T*xi(i)/Ri(i)
        ENDDO

        temp = T
        temp(nb_temp) = Te
        ! Set the current state of the mixture
        CALL mpp_set_state(rhoi, temp, state_var)

        ! Compute the species mass diffusion flux 
        DO d = 1,nb_dim      

           l = (d - 1)*nb_ns + 1 
           r = l + nb_ns - 1 
           CALL mpp_stefan_maxwell(diff_driv(l:r), Ji(l:r), ambE)

           DO i = 1,nb_ns 
              Ji(l + i - 1) = rhoi(i)*Ji(l + i - 1) 
           ENDDO
        ENDDO       
       
      END SUBROUTINE Mutationpp_compute_species_DiffFlux
 
      !------------------------------------------------------!
      ! This subrouotine computes transport coefficients 
      SUBROUTINE Mutationpp_compute_transpCoeff (p, xi, temp, mu, kappa, lambda, Di, chi) 

        INTEGER :: i
        INTEGER :: state_var

        REAL(KIND=8) :: T, p_ov_T
        REAL(KIND=8), DIMENSION(nb_ns) :: xi_tol, rhoi

        REAL(KIND=8), INTENT(IN) :: p
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, temp
        REAL(KIND=8), INTENT(OUT) :: mu, kappa
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: chi, Di, lambda
        
        ! Apply tolerance on species molar fractions
        xi_tol = xi
        CALL Mutationpp_comp_tol(xi_tol)

        ! Temperature 
        T = temp(1)
       
        ! Compute species densities
        p_ov_T = p/T
        DO i = 1,nb_ns 
           rhoi(i) = p_ov_T*xi(i)/Ri(i)
        ENDDO

        ! Set the current state of the mixture
        state_var = 1
        CALL mpp_set_state(rhoi, temp, state_var)
  
        ! Mixture dynamic viscosity 
        mu = mpp_viscosity()

        ! Mixture bulk viscosity (set to zero)
        kappa = 0.d0

        ! Mixture total thermal conductivity 
        CALL mpp_frozen_thermal_conductivity(lambda)
        CALL mpp_average_diffusion_coeffs(Di) 

        ! Species thermal diffusion ratios (set to zero)
        chi = 0.d0

      END SUBROUTINE Mutationpp_compute_transpCoeff  

      !------------------------------------------------------!
      ! This subroutine computes the thermal conductivity
      SUBROUTINE Mutationpp_compute_thermal_cond (rhoi, T, lambda) 

        INTEGER :: state_var = 1

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: T
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: lambda
        
         
        ! Set the current state of the mixture
        CALL mpp_set_state(rhoi, T, state_var)


        ! Mixture total thermal conductivity 
        CALL mpp_frozen_thermal_conductivity(lambda)        

      END SUBROUTINE Mutationpp_compute_thermal_cond

      !------------------------------------------------------!
      ! This subroutine get the species names
      SUBROUTINE Mutationpp_species_name(species_name)
     
        INTEGER :: i
        CHARACTER(Len=*), INTENT(OUT), DIMENSION(:) :: species_name
 
        DO i= 1,nb_ns
          CALL mpp_species_name(i,species_name(i)) 
        ENDDO

      END SUBROUTINE Mutationpp_species_name


     !------------------------------------------------------!
     ! This subroutine computes the species mass diffusion flux 
      SUBROUTINE Mutationpp_compute_ambipolar_Electric_field (p, T, Te, xi, diff_driv, Ji, ambE)

        INTEGER :: d, l, i, r
        INTEGER :: state_var = 1

        REAL(KIND=8) :: p_ov_T
        REAL(KIND=8), INTENT(OUT) :: ambE 
        REAL(KIND=8), DIMENSION(nb_ns) :: rhoi
        REAL(KIND=8), DIMENSION(nb_temp) :: temp

        REAL(KIND=8), INTENT(IN) :: p, T, Te
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, diff_driv        
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: Ji


        ! Compute the species densities 
        p_ov_T = p/T
        DO i = 1,nb_ns 
           rhoi(i) = p_ov_T*xi(i)/Ri(i)
        ENDDO

        ! WARNING quick debug for the 2T case without electron
        ! To be fixed more correctly
        temp = T
        ! Set the current state of the mixture
        CALL mpp_set_state(rhoi, temp, state_var)

                
        nb_dim = 1
        ! Compute the species mass diffusion flux 
        DO d = 1,nb_dim      

           l = (d - 1)*nb_ns + 1 
           r = l + nb_ns - 1 
           
           CALL mpp_stefan_maxwell(diff_driv(l:r), Ji(l:r), ambE)
           DO i = 1,nb_ns 
              Ji(l + i - 1) = rhoi(i)*Ji(l + i - 1) 
           ENDDO
        ENDDO  
 
      END SUBROUTINE Mutationpp_compute_ambipolar_Electric_field
      

      !------------------------------------------------------!
! This subroutine computes the surface mass balance with equilibrium and species constraint in the condensed phase

      SUBROUTINE Mutationpp_surface_mass_balance_constraint (input_mixture, T, P, xi) 
       
        integer :: i
        REAL(KIND=8), INTENT(IN) :: T, P
        CHARACTER*(*), INTENT(IN) :: input_mixture

        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: xi
                
        !CALL mpp_surface_mass_balance_stag_line (input_mixture, T, P, xi)
        
      WRITE(*,*) "surface mass balance not implemented"
      STOP

      END SUBROUTINE Mutationpp_surface_mass_balance_constraint
  
      !------------------------------------------------------!
! This subroutine computes the number of elements in mixture

      SUBROUTINE Mutationpp_number_of_elements(elements) 

        INTEGER, INTENT(INOUT) :: elements


        elements = mpp_nelements()

      END SUBROUTINE Mutationpp_number_of_elements
  
      !------------------------------------------------------!


      !------------------------------------------------------!
      ! This subroutine get the species index from the species name
      SUBROUTINE Mutationpp_get_species_index (sp_name, sp_index) 

        CHARACTER*(*), INTENT(IN) :: sp_name
        INTEGER, INTENT(OUT)      :: sp_index
        
        ! Get the index of the  of the mixture
        sp_index= mpp_species_index(sp_name)

      END SUBROUTINE Mutationpp_get_species_index  

      !------------------------------------------------------!
      ! This subroutine get the species name from the species index
      SUBROUTINE Mutationpp_get_species_name (sp_index, sp_name) 

        INTEGER, INTENT(IN)        :: sp_index
        CHARACTER*(*), INTENT(OUT) :: sp_name

        ! Get the species name
        CALL mpp_species_name(sp_index, sp_name)

      END SUBROUTINE Mutationpp_get_species_name

      !------------------------------------------------------!
      ! This subroutine computes the equilibrium composition given pressure and temperature
      ! and a particular elemental composition (To be generalized)
      SUBROUTINE Mutationpp_compute_eq_composition_pyro (p, T, yi)

        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi
        REAL(KIND=8), DIMENSION(5) :: el
        !REAL(KIND=8), DIMENSION(4) :: el

        ! Pyrolysis gas elemental mole fractions in case of e,C,H,N,O (order) mixture
        el(1) = 0.d0
        el(2) = 0.22930d0
        el(3) = 0.66059d0
        el(4) = 0.d0 
        el(5) = 0.11011d0
        ! Pyrolysis gas elemental compostion in case of C,H,N,O (order) mixture
        !el(1) = 0.22930d0 
        !el(2) = 0.66059d0
        !el(3) = 0.d0
        !el(4) = 0.11011d0

        ! Compute the equilibrium composition giving mole fractions as output
        CALL mpp_pyro_equilibrium_composition(T, p, el, yi)
        CALL mpp_convert_x_to_y(yi, yi)

      END SUBROUTINE Mutationpp_compute_eq_composition_pyro
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Below are the subroutines to be added for the GSI in M++ !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !------------------------------------------------------!
      ! This subroutine set the diffusion model
      SUBROUTINE Mutationpp_set_diff_model (xi, dx) 

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi !<mole fraction of first physical cell 
        REAL(KIND=8), INTENT(IN) :: dx !< distance wall-cell_center

        ! Set diffusion model
        CALL mpp_set_diffusion_model(xi, dx)

     END SUBROUTINE Mutationpp_set_diff_model

     !------------------------------------------------------!
     ! This subroutine set the conductive heat flux calculation for the SEB
     SUBROUTINE Mutationpp_set_cond_heat_flux (T, dx) 

       REAL(KIND=8),DIMENSION(:), INTENT(IN) :: T  !< Temperature of first physical cell
       REAL(KIND=8), INTENT(IN) :: dx !< distance wall-cell_center

       ! Set diffusion model
        CALL mpp_set_cond_heat_flux(T, dx)

     END SUBROUTINE Mutationpp_set_cond_heat_flux

     !------------------------------------------------------!
     ! This subroutine set the wall state (species densities and temperature)
     SUBROUTINE Mutationpp_set_wall_state (rhoi, T) 

       REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi !< wall species densities
       REAL(KIND=8), DIMENSION(:), INTENT(IN) :: T !< wall temperature
       INTEGER, PARAMETER :: set_state_rhoi_T = 1 !< set state (rhoi-T) flag

       ! Set the state of the wall
        CALL mpp_set_surface_state(rhoi, T, set_state_rhoi_T)

     END SUBROUTINE Mutationpp_set_wall_state

     !------------------------------------------------------!
      ! This subroutine set the radiation flux into the wall
     SUBROUTINE Mutationpp_set_wall_radiation (rad_flux) 

       REAL(KIND=8), INTENT(IN) :: rad_flux !< wall radiation flux

        !CALL mpp_set_wall_radiation(rad_flux)

        PRINT*,'in "Mutationpp", not implemented yet'
        STOP

     END SUBROUTINE Mutationpp_set_wall_radiation

     !------------------------------------------------------!
     ! This subroutine call M++ to solve the surface mass balance
     SUBROUTINE Mutationpp_solve_surface_balance () 

       ! Solve the surface mass balance
        CALL mpp_solve_surface_balance()

     END SUBROUTINE Mutationpp_solve_surface_balance

     !------------------------------------------------------!
     ! This subroutine returns the wall state as stored in M++
     SUBROUTINE Mutationpp_get_wall_state (rhoi, T) 

       REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: rhoi !< wall species densities
       REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: T !< wall temperature
       INTEGER, PARAMETER :: set_state_rhoi_T = 1 !< set state (rhoi-T) flag

       ! Solve the surface mass balance
        CALL mpp_get_surface_state(rhoi, T, set_state_rhoi_T)

     END SUBROUTINE Mutationpp_get_wall_state

     !------------------------------------------------------!
     ! This subroutine returns the surface source terms (kg/m2-s)
     SUBROUTINE Mutationpp_get_surface_production_rates (wdoti) 

       REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: wdoti !< wall source terms

       ! Get the computed wall source terms
        CALL mpp_surface_production_rates(wdoti)

     END SUBROUTINE Mutationpp_get_surface_production_rates

     !------------------------------------------------------!
     ! This subroutine returns the total mass blowing flux (kg/m2-s)
     SUBROUTINE Mutationpp_get_mass_blowing_rate(mdot) 

       REAL(KIND=8), INTENT(OUT) :: mdot !< mass blowing gas char + pyrolysis (if present)

       ! Get the computed wall source terms
        CALL mpp_mass_blowing_rate(mdot)

      END SUBROUTINE Mutationpp_get_mass_blowing_rate

#endif  
  END MODULE mod_Mutationpp_library
!------------------------------------------------------------------------------!  
