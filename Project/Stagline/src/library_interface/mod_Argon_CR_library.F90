!------------------------------------------------------------------------------!
! This module implements an interface for the Argon_CR thermodynamic library. It provides an
! implementation of functions and subroutines used to compute thermodynamic and transport properties and 
! source terms related to collisional and radiative processes.
! The interface is written such to have no data and variable dependency from the flow solver being used. 
! Therefore it can plugged to any code without big efforts. 
  MODULE mod_Argon_CR_library 

#include"../config.h"

#ifdef ARGON_CR
    IMPLICIT NONE
  
    INTEGER, SAVE :: nb_ns, nb_tvib, nb_te, nb_trot, nb_temp, nb_int_temp, nb_eq, nb_dim 
    REAL(KIND=8), PARAMETER :: kb = 1.380658d-23
    REAL(KIND=8), PARAMETER :: na = 6.0221367d23  
    REAL(KIND=8), PARAMETER :: ue = 1.602191d-19 
    REAL(KIND=8), PARAMETER :: ru = kb*na
    REAL(KIND=8), PARAMETER :: gamma_e    = 5.d0/3.d0
    REAL(KIND=8), PARAMETER :: gamma_e_m1 = gamma_e - 1.d0
    REAL(KIND=8), SAVE :: Re
    REAL(KIND=8), SAVE :: xi_tol
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: mass, Ri

    ! Subroutines of Argon_CR thermodynamic library interface. In some cases they call 
    ! other subroutines for the model under use (even these are not dependent on flow solver data). 
    CONTAINS
  
      !------------------------------------------------------!
      ! This subroutine initializes the thermodynamic library in use.  
      SUBROUTINE Argon_CR_initialize (input_ns, input_ntrot, input_ntvib, input_nte, input_ntemp, input_neq,  & 
                                    & input_ndim, input_solver, input_mixture, input_reaction, state_model,   &
                                    & input_transf, input_path, input_xi_tol) 
  
        USE mod_argon_CR_initialize_CFD,            ONLY: initialize, pos_em
        USE mod_argon_CR_zeta,                      ONLY: set_pointer

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
           WRITE(*,5)input_solver(1:length),'::Argon_CR_library interface'   
           PRINT*
           WRITE(*,5)input_solver(1:length),'::Error in number of temperatures'     
           PRINT*
           STOP
        ENDIF

        neq = input_ns + 1 + input_ntrot + input_ntvib + input_nte + input_ndim
        IF (neq.NE.nb_eq) THEN
           PRINT*
           WRITE(*,5)input_solver(1:length),'::Argon_CR_library interface'  
           PRINT* 
           WRITE(*,5)input_solver(1:length),'::Error in number of equations'     
           PRINT*
           STOP
        ENDIF
 
        ! Tolerance on species molar fractions
        xi_tol = input_xi_tol

        ! Allocation of molecular mass and specific gas constant vectors
        ALLOCATE(mass(nb_ns), Ri(nb_ns))
        
        ! Argon_CR library is initialized
        CALL initialize (input_solver, input_mixture, input_transf, input_reaction, input_path, & 
                       & nb_ns, nb_trot, nb_tvib, nb_te, nb_eq, nb_dim, mass)

        ! Speficic gas constant of single species
        DO is = 1,nb_ns 
           Ri(is) = 1.d0/mass(is)
        ENDDO 
        Ri = ru*Ri

        ! Free electron gas constant
        Re = Ri(1)

        ! Pointer initialization
        CALL set_pointer()

5     FORMAT(A,A)

      END SUBROUTINE Argon_CR_initialize

      !------------------------------------------------------!
      ! This subroutine finalizes the Argon_CR thermodynamic library
      SUBROUTINE Argon_CR_finalize ()

         USE mod_argon_CR_initialize_CFD,           ONLY: finalize

         IF (ALLOCATED(mass))       DEALLOCATE(mass)
         IF (ALLOCATED(Ri))         DEALLOCATE(Ri)

         CALL finalize ()

      END SUBROUTINE Argon_CR_finalize

      !------------------------------------------------------!
      ! This subroutine provides the position of the electron in the species list 
      ! and the electron gas constant
      SUBROUTINE Argon_CR_get_el_data (pos, R)

        USE mod_argon_CR_initialize_CFD,            ONLY: initialize, pos_em

        INTEGER, INTENT(OUT) :: pos
        REAL(KIND=8), INTENT(OUT) :: R

        pos = pos_em
        R   = Re

      END SUBROUTINE Argon_CR_get_el_data

      !------------------------------------------------------!
      ! This subroutine computes the free electron pressure 
      SUBROUTINE Argon_CR_get_el_pressure (rho, Te, pe)

        REAL(KIND=8), INTENT(IN) :: rho, Te
        REAL(KIND=8), INTENT(OUT) :: pe

        pe = rho*Re*Te

      END SUBROUTINE Argon_CR_get_el_pressure

      !------------------------------------------------------!
      ! This subroutine computes post-shock nonequilibrium conditions
      SUBROUTINE Argon_CR_post_shock_neq (p1, u1, T1, p2, u2, T2, yi, extra)

        USE mod_argon_CR_CFD_prop,                  ONLY: post_shock_neq
        USE mod_argon_CR_CFD_eq,                    ONLY: eq_composition_CR
 
        REAL(KIND=8), INTENT(IN) :: p1, u1, T1
        REAL(KIND=8), INTENT(IN), OPTIONAL :: extra
        REAL(KIND=8), INTENT(OUT) :: p2, u2, T2
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi 

        CALL post_shock_neq (p1, u1, T1, p2, u2, T2, yi)

      END SUBROUTINE Argon_CR_post_shock_neq

      !------------------------------------------------------!
      ! This subroutine computes post-shock equilibrium conditions
      SUBROUTINE Argon_CR_post_shock_eq (p1, u1, T1, p2, u2, T2, yi, rhoi, extra)

        REAL(KIND=8), INTENT(IN) :: p1, u1, T1
        REAL(KIND=8), INTENT(IN), OPTIONAL :: extra
        REAL(KIND=8), INTENT(OUT) :: p2, u2, T2
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi, rhoi

        p2   = 0.d0
        u2   = 0.d0
        T2   = 0.d0
        yi   = 0.d0
        rhoi = 0.d0

        PRINT*,'in "Argon_CR_post_shock_eq", not implemented yet..'
        STOP 

      END SUBROUTINE Argon_CR_post_shock_eq

      !------------------------------------------------------!
      ! This subroutine write down the solution when solving the ODE system for
      ! studying the nonequilibrium flow behind a shock wave
      SUBROUTINE Argon_CR_write_ODE_solution (flow_file, pop_file, k, xold, x, uold, m_dot, z)

        USE mod_argon_CR_initialize_CFD,            ONLY: model, lev_Ar, lev_Arp, pos_em, pos_Ar, pos_Arp, & 
                                                        & EiJ_Ar, gi_Ar
        USE mod_argon_CR_CFD_prop,                  ONLY: Q_int

        INTEGER :: is
        REAL(KIND=8) :: fac, Qint
        REAL(KIND=8) :: alpha, u, p_e, p_h, p, rho, y_em, y_Ar, y_Arp, T, Te
        REAL(KIND=8), SAVE :: time

        INTEGER, INTENT(IN) :: flow_file, pop_file, k
        REAL(KIND=8), INTENT(IN) :: xold, x, uold, m_dot
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: z

        ! Density, velocity and temperatures
        u   = z(nb_ns + 1)
        rho = m_dot/u
        T   = z(nb_ns + 2)
        Te  = z(nb_ns + 3)

        ! Mass fractions of em, Ar and Arp 
        y_em  = z(1)

        y_Ar = 0.d0
        DO is = pos_Ar,pos_Ar + lev_Ar - 1 
           y_Ar = y_Ar + z(is)
        ENDDO

        y_Arp = 0.d0
        DO is = pos_Arp,pos_Arp + lev_Arp - 1
           y_Arp = y_Arp + z(is)
        ENDDO
        
        ! Ionization fraction
        alpha = y_Arp/(y_Ar + y_Arp)

        ! Pressure 
        p_e = rho*z(pos_em)*Ri(pos_em)*Te
        p_h = 0.d0
        DO is = pos_em + 1,nb_ns 
           p_h = p_h + z(is)*Ri(is)
        ENDDO
        p_h = p_h*rho*T

        p = p_e + p_h

        ! Lagrangian time 
        IF (k.EQ.1) THEN 
           time = 0.d0
           WRITE(pop_file,5)lev_Ar
        ELSE 
           time = time + 0.5d0*(x - xold)*(1.d0/uold + 1.d0/u)
        ENDIF

        ! Flowfield output file
        WRITE(flow_file,10)x,time,y_em,y_Ar,y_Arp,alpha,rho,u,p,T,Te
        CALL flush(flow_file)

        ! Internal parition function
        CALL Q_int (Te, 'Ar', Qint)

        ! Population file 
        fac = rho*na/(mass(pos_Ar)*Qint)
        DO is = 1,lev_Ar 
           WRITE(pop_file,10)EiJ_Ar(is)/ue,fac*z(pos_Ar +is - 1)/gi_Ar(is)
           CALL flush(pop_file)  
        ENDDO 
        WRITE(pop_file,20)'& ',k,time,x

5     FORMAT(I5)
10    FORMAT(100E20.10)
20    FORMAT(A,I6,100E20.10)

      END SUBROUTINE Argon_CR_write_ODE_solution 

      !----------------------------------------------------!
      ! This subroutine post-process the ODE solution for the flow 
      ! behind a normal shock wave
      SUBROUTINE Argon_CR_post_process_ODE_solution (post_file, x, ni, temp)  

        INTEGER, INTENT(IN) :: post_file
        REAL(KIND=8), INTENT(IN) :: x
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: ni, temp

        PRINT*,'in "Argon_CR_post_process_ODE_solution", not implemented yet..'
        STOP 

      END SUBROUTINE Argon_CR_post_process_ODE_solution

      !------------------------------------------------------!
      ! This subroutine writes down the solution when solving the flow equations 
      ! by means of Finite Volume method
      SUBROUTINE Argon_CR_write_fvmcc_solution (flow_file, pop_file, k, xold, x, uold, u, rhoi, temp)

        USE mod_argon_CR_initialize_CFD,            ONLY: model, lev_Ar, lev_Arp, pos_em, pos_Te, pos_Ar, pos_Arp, & 
                                                        & EiJ_Ar, gi_Ar
        USE mod_argon_CR_CFD_prop,                  ONLY: Q_int

        INTEGER :: is
        REAL(KIND=8) :: alpha, x_em, x_Ar, x_Arp, y_Ar, y_Arp 
        REAL(KIND=8) :: p, ph, pe, rho, T, Te 
        REAL(KIND=8) :: fac, Qint, mm
        REAL(KIND=8), SAVE :: time
        REAL(KIND=8), DIMENSION(nb_ns) :: yi, xi

        INTEGER, INTENT(IN) :: flow_file, pop_file, k
        REAL(KIND=8), INTENT(IN) :: xold, x, uold, u
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp

        ! Temperatures
        T  = temp(1)
        Te = temp(pos_Te)

        ! Compute mixture density 
        CALL Argon_CR_get_density (rhoi, rho)

        ! Compute species mass fractions
        CALL Argon_CR_get_mass_fractions (rho, rhoi, yi)
        
        ! Compute species molar fractions
        CALL Argon_CR_molar_frac_from_mass_frac (yi, mm, xi)

        ! Molar fractions of em, Ar and Arp and mass fractions of Ar and Arp
        x_em = xi(pos_em)
        x_Ar = 0.d0
        y_Ar = 0.d0
        DO is = pos_Ar,pos_Ar + lev_Ar - 1 
           x_Ar = x_Ar + xi(is)
           y_Ar = y_Ar + yi(is)
        ENDDO

        x_Arp = 0.d0
        y_Arp = 0.d0
        DO is = pos_Arp,pos_Arp + lev_Arp - 1
           x_Arp = x_Arp + xi(is)
           y_Arp = y_Arp + yi(is)
        ENDDO

        ! Ionization fraction
        alpha = y_Arp/(y_Ar + y_Arp)
       
        ! Heavy-particle pressure pressure 
        ph = 0.d0 
        DO is = pos_em + 1,nb_ns 
           ph = ph + rhoi(is)*Ri(is)
        ENDDO
        ph = ph*T 
 
        ! Free-electron pressure 
        CALL Argon_CR_get_el_pressure (rhoi(pos_em), Te, pe)

        ! Mixture pressure 
        p = ph + pe

        WRITE(flow_file,10)x,x_em,x_Ar,x_Arp,alpha,rho,p,u,temp

        ! Write population file in case of use of the CR model
        IF (model.EQ.'CR') THEN

           ! Lagrangian time 
           IF (k.EQ.1) THEN 
              time = 0.d0
              WRITE(pop_file,5)lev_Ar
           ELSE 
              time = time + 0.5d0*(x - xold)*(1.d0/uold + 1.d0/u)
           ENDIF

           ! Internal parition function
           CALL Q_int (Te, 'Ar', Qint)

           fac = na/(mass(pos_Ar)*Qint)
           DO is = 1,lev_Ar 
              WRITE(pop_file,10)EiJ_Ar(is)/ue,fac*rhoi(pos_Ar + is - 1)/gi_Ar(is)
              CALL flush(pop_file)  
           ENDDO 
           WRITE(pop_file,20)'& ',k,time,x

        ENDIF

5     FORMAT(I5)
10    FORMAT(100E20.10)
20    FORMAT(A,I6,100E20.10)

      END SUBROUTINE Argon_CR_write_fvmcc_solution   

      !------------------------------------------------------!
      ! This subroutine computes the equilibrium composition given pressure and temperature
      SUBROUTINE Argon_CR_compute_eq_composition (p, T, rhoi)

        USE mod_argon_CR_CFD_eq,         ONLY: eq_composition_CR 

        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: rhoi

        CALL eq_composition_CR (p, T, rhoi)

      END SUBROUTINE Argon_CR_compute_eq_composition

      !------------------------------------------------------!
      ! This subroutine computes the species equilibrium internal energy per unit mass 
      SUBROUTINE Argon_CR_compute_eq_species_energy (T, e)

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: e

        e = 0.d0

        PRINT*,'in "Argon_CR_compute_eq_species_energy", not implemented yet..'
        STOP 

      END SUBROUTINE Argon_CR_compute_eq_species_energy 

      !------------------------------------------------------!
      ! This subroutine computes the species equilibrium enthalpy per unit mass 
      SUBROUTINE Argon_CR_compute_eq_species_enthalpy (T, h)

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: h

        h = 0.d0

        PRINT*,'in "Argon_CR_compute_eq_species_enthalpy", not implemented yet..'
        STOP 

      END SUBROUTINE Argon_CR_compute_eq_species_enthalpy 

      !------------------------------------------------------!
      ! This subroutine computes the species equilibrium entropy per unit mass. The
      ! entropy of mixing contributions is computed separately. 
      SUBROUTINE Argon_CR_compute_eq_species_entropy (p, T, yi, xi, smix, s)

        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), INTENT(IN), DIMENSION(:) :: xi, yi
        REAL(KIND=8), INTENT(OUT) :: smix
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: s

        s    = 0.d0
        smix = 0.d0

        PRINT*,'in "Argon_CR_compute_eq_species_entropy", not implemented yet..'
        STOP 

      END SUBROUTINE Argon_CR_compute_eq_species_entropy    

      !------------------------------------------------------!
      ! This subroutine computes the temperatures and thermodynamic data,
      ! additional data for Jacobians are given in output
      SUBROUTINE Argon_CR_get_data (rho, rhoi, rho_eint, temp, c, gamma, p, alpha, beta, ei)
  
        USE mod_argon_CR_initialize_CFD,            ONLY: pos_em, pos_Te
        USE mod_argon_CR_function_pointer,          ONLY: get_temperatures, get_species_energy_cv        

        INTEGER :: is
        REAL(KIND=8) :: tmp, sum1, sum2
        REAL(KIND=8) :: ki, T, Te, pe, ph
        REAL(KIND=8), DIMENSION(nb_ns) :: cvi

        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, rho_eint
        REAL(KIND=8), INTENT(OUT) :: c, gamma, p, alpha
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp, beta, ei
  
        ! Compute the temperatures
        CALL get_temperatures (rhoi, rho_eint, temp)

        T  = temp(1)
        Te = temp(pos_Te) 

        ! Compute the free electron pressure 
        CALL Argon_CR_get_el_pressure (rhoi(pos_em), Te, pe)

        ! Compute species energies and constant volume specific heats
        CALL get_species_energy_cv (temp, ei, cvi)

        sum1 = 0.d0   
        sum2 = 0.d0      
        DO is = pos_em + 1,nb_ns 
           tmp  = rhoi(is)
           sum1 = sum1 + tmp*Ri(is)
           sum2 = sum2 + tmp*cvi(is)
        ENDDO
  
        ph = sum1*T
        p  = ph + pe
        ki = sum1/sum2
        gamma   = 1.d0 + ki 
        alpha   = ki 
        beta(1) = sum2
        c       = DSQRT((gamma*ph + gamma_e*pe)/rho)

      END SUBROUTINE Argon_CR_get_data 
  
      !------------------------------------------------------!
      ! This subroutine computes thermodynamic data needed to be stored
      ! for nonequilibrium gas flow computations. 
      SUBROUTINE Argon_CR_get_thermodynamic_data (rho, rhoi, temp, c, gamma, p, alpha, beta, ei)
  
        USE mod_argon_CR_initialize_CFD,            ONLY: pos_em, pos_Te
        USE mod_argon_CR_function_pointer,          ONLY: get_temperatures, get_species_energy_cv

        INTEGER :: is
        REAL(KIND=8) :: tmp, sum1, sum2
        REAL(KIND=8) :: ki, T, Te, pe, ph
        REAL(KIND=8), DIMENSION(nb_ns) :: cvi

        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), INTENT(OUT) :: c, gamma, p, alpha
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: beta, ei

        ! Temperatures
        T  = temp(1)
        Te = temp(pos_Te) 

        ! Compute the free electron pressure 
        CALL Argon_CR_get_el_pressure (rhoi(pos_em), Te, pe)

        ! Compute species energies and constant volume specific heats
        CALL get_species_energy_cv (temp, ei, cvi)

        sum1 = 0.d0   
        sum2 = 0.d0      
        DO is = pos_em + 1,nb_ns 
           tmp  = rhoi(is)
           sum1 = sum1 + tmp*Ri(is)
           sum2 = sum2 + tmp*cvi(is)
        ENDDO
  
        ph = sum1*T
        p  = ph + pe
        ki = sum1/sum2
        gamma   = 1.d0 + ki 
        alpha   = ki 
        beta(1) = sum2
        c       = DSQRT((gamma*ph + gamma_e*pe)/rho)

      END SUBROUTINE Argon_CR_get_thermodynamic_data
 
      !------------------------------------------------------!
      ! This subroutine computes the energy density for total internal,
      ! rotational, vibrational and electronic components. In case of free
      ! electrons what is computed is the free electron pseudo-entropy density
      SUBROUTINE Argon_CR_get_energy_densities (rhoi, temp, rho_e)
 
        USE mod_argon_CR_initialize_CFD,            ONLY: pos_em, pos_Te
        USE mod_argon_CR_function_pointer,          ONLY: get_species_energy
 
        INTEGER :: is
        REAL(KIND=8) :: tmp1, tmp2, rho, pe
        REAL(KIND=8), DIMENSION(nb_ns) :: ei

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: rho_e

        ! Compute the species specific energies
        CALL get_species_energy (temp, ei) 

        tmp1 = 0.d0
        rho  = 0.d0
        DO is = 1,nb_ns 
           tmp2 = rhoi(is)
           tmp1 = tmp1 + tmp2*ei(is)
           rho  = rho + tmp2
        ENDDO
        
        ! Free electron pressure 
        CALL Argon_CR_get_el_pressure (rhoi(pos_em), temp(pos_Te), pe)

        rho_e(1) = tmp1
        rho_e(2) = pe*(rho**(-gamma_e_m1))
        
      END SUBROUTINE Argon_CR_get_energy_densities
  
      !------------------------------------------------------!
      ! This subroutine computes the total internal, rotational, vibrational and electronic energy 
      ! per unit mass of the mixture.
      SUBROUTINE Argon_CR_get_energy (rho, rhoi, temp, e)
  
        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e
      
        e = 0.d0
 
        PRINT*,'in "Argon_CR_get_species_energy", not implemented yet..'
        STOP 

      END SUBROUTINE Argon_CR_get_energy
  
      !------------------------------------------------------!
      ! This subroutine computes the total internal, rotational, vibrational and electronic enthalpy 
      ! per unit mass of the mixture.
      SUBROUTINE Argon_CR_get_enthalpy (rho, rhoi, temp, h)
  
        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: h
      
        h = 0.d0
 
        PRINT*,'in "Argon_CR_get_species_energy", not implemented yet..'
        STOP 

      END SUBROUTINE Argon_CR_get_enthalpy 
 
      !----------------------------------------------------!
      ! This subroutine computes the species internal energies and constant 
      ! volume specific heat per unit mass
      SUBROUTINE Argon_CR_get_species_energy (temp, e) 

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e

        e = 0.d0

        PRINT*,'in "Argon_CR_get_species_energy", not implemented yet..'
        STOP 

      END SUBROUTINE Argon_CR_get_species_energy

      !----------------------------------------------------!
      ! This subroutine computes the species internal energies heat per unit mass
      SUBROUTINE Argon_CR_get_species_energy_cv (temp, e, cv) 

        USE mod_argon_CR_function_pointer,              ONLY: get_species_energy_cv
 
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e, cv

        CALL get_species_energy_cv (temp, e, cv)

      END SUBROUTINE Argon_CR_get_species_energy_cv

      !----------------------------------------------------!
      ! This subroutine computes the species internal energies and constant 
      ! volume specific heat per unit mass
      SUBROUTINE Argon_CR_get_species_enthalpy (temp, h) 
 
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: h

        h = 0.d0 

        PRINT*,'in "Argon_CR_get_species_enthalpy", not implemented yet..'
        STOP 

      END SUBROUTINE Argon_CR_get_species_enthalpy

      !----------------------------------------------------!
      ! This subroutine computes the species enthalpies and constant 
      ! pressure specific heat per unit mass
      SUBROUTINE Argon_CR_get_species_enthalpy_cp (temp, h, cp) 

        USE mod_argon_CR_initialize_CFD,            ONLY: pos_em, pos_Te
        USE mod_argon_CR_function_pointer,          ONLY: get_species_energy_cv

        INTEGER :: is
        REAL(KIND=8) :: tmp
        REAL(KIND=8) :: T, Te
 
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: h, cp

        CALL get_species_energy_cv (temp, h, cp)
 
        ! Temperatures 
        T  = temp(1)
        Te = temp(nb_temp)

        ! em
        tmp = Re
        h(pos_em)  = h(pos_em)  + tmp*Te
        cp(pos_em) = cp(pos_em) + tmp 

        ! Electronic levels of Ar and Arp
        DO is = 2,nb_ns 
           tmp = Ri(is)
           h(is)  = h(is)  + tmp*T
           cp(is) = cp(is) + tmp
        ENDDO         

      END SUBROUTINE Argon_CR_get_species_enthalpy_cp

      !------------------------------------------------------!
      ! This subroutine computes the frozen specific heats at constant volume of each species
      SUBROUTINE Argon_CR_get_species_frozen_cv (T, cv)  
 
        REAL(KIND=8), INTENT(IN) :: T 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv
  
        cv = 0.d0

        PRINT*,'in "Argon_CR_get_species_frozen_cv", not implemented yet..'
        STOP 

      END SUBROUTINE Argon_CR_get_species_frozen_cv 
  
      !------------------------------------------------------!
      ! This subroutine computes the frozen specific heats at constant pressure of each species
      SUBROUTINE Argon_CR_get_species_frozen_cp (T, cp)  
  
        REAL(KIND=8), INTENT(IN) :: T 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cp
  
        cp = 0.d0

        PRINT*,'in "Argon_CR_get_species_frozen_cp", not implemented yet..'
        STOP 

      END SUBROUTINE Argon_CR_get_species_frozen_cp

      !------------------------------------------------------!
      ! This subroutine computes the frozen specific heats at constant volume of each species and 
      ! the species energy components in thermal equilibrium with translation.
      SUBROUTINE Argon_CR_get_species_frozen_energy_cv (T, cv, e)  
 
        REAL(KIND=8), INTENT(IN) :: T 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv, e
 
        cv = 0.d0
        e  = 0.d0 

        PRINT*,'in "Argon_CR_get_species_frozen_energy_cv", not implemented yet..'
        STOP 

      END SUBROUTINE Argon_CR_get_species_frozen_energy_cv 
  
      !------------------------------------------------------!
      ! This subroutine computes the frozen specific heats at constant pressure of each species and 
      ! the species enthalpy components in thermal equilibrium with translation.
      SUBROUTINE Argon_CR_get_species_frozen_enthalpy_cp (T, cp, h)  
  
        REAL(KIND=8), INTENT(IN) :: T 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cp, h

        cp = 0.d0
        h  = 0.d0  

        PRINT*,'in "Argon_CR_get_species_frozen_enthalpy_cp", not implemented yet..'
        STOP 

      END SUBROUTINE Argon_CR_get_species_frozen_enthalpy_cp

      !------------------------------------------------------!
      ! This subroutine computes the temperatures given the energy densities.
      SUBROUTINE Argon_CR_get_temperatures (rhoi, rho_eint, temp)
  
        USE mod_argon_CR_function_pointer,          ONLY: get_temperatures

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rho_eint
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp
  
        CALL get_temperatures (rhoi, rho_eint, temp)

      END SUBROUTINE Argon_CR_get_temperatures
 
      !------------------------------------------------------!
      ! This subroutine computes the equilibrium speed of sound
      SUBROUTINE Argon_CR_get_eq_sound_speed (rhoi, p, T, c)
 
        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), INTENT(IN), DIMENSION(:) :: rhoi
        REAL(KIND=8), INTENT(OUT) :: c

        c = 0.d0

        PRINT*,'in "Argon_CR_get_eq_sound_speed", not implemented yet..'
        STOP 

      END SUBROUTINE Argon_CR_get_eq_sound_speed 

      !------------------------------------------------------!
      ! This subroutine computes the equilibrium specific heat ratio
      SUBROUTINE Argon_CR_get_eq_gamma (rhoi, p, T, gamma)
 
        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), INTENT(IN), DIMENSION(:) :: rhoi
        REAL(KIND=8), INTENT(OUT) :: gamma

        gamma = 0.d0

        PRINT*,'in "Argon_CR_get_eq_gamma", not implemented yet..'
        STOP 

      END SUBROUTINE Argon_CR_get_eq_gamma 

      !------------------------------------------------------!
      ! This subroutine computes the equilibrium specific heat ratio
      SUBROUTINE Argon_CR_get_eq_gamma_sound_speed (rhoi, p, T, gamma, c)
 
        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), INTENT(IN), DIMENSION(:) :: rhoi
        REAL(KIND=8), INTENT(OUT) :: gamma, c
 
        gamma = 0.d0
        c     = 0.d0

        PRINT*,'in "Argon_CR_get_eq_gamma_sound_speed", not implemented yet..'
        STOP  

      END SUBROUTINE Argon_CR_get_eq_gamma_sound_speed 

      !------------------------------------------------------!
      ! This subroutine computes the frozen speed of sound
      SUBROUTINE Argon_CR_get_frozen_sound_speed (rhoi, temp, c)
 
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp, rhoi
        REAL(KIND=8), INTENT(OUT) :: c 

        c = 0.d0 

        PRINT*,'in "Argon_CR_get_frozen_sound_speed", not implemented yet..'
        STOP  

      END SUBROUTINE Argon_CR_get_frozen_sound_speed
 
      !------------------------------------------------------!
      ! This subroutine gets the frozen speed specific heat ratio.
      SUBROUTINE Argon_CR_get_frozen_gamma (rhoi, temp, gamma)
 
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), INTENT(OUT) :: gamma 
  
        gamma = 0.d0 
 
        PRINT*,'in "Argon_CR_get_frozen_gamma", not implemented yet..'
        STOP 

      END SUBROUTINE Argon_CR_get_frozen_gamma
 
      !------------------------------------------------------!
      ! This subroutine gets the frozen speed of sound and specific heat ratio.
      SUBROUTINE Argon_CR_get_frozen_gamma_sound_speed (rhoi, temp, gamma, c)
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), INTENT(OUT) :: gamma, c 
  
        gamma = 0.d0
        c     = 0.d0

        PRINT*,'in "Argon_CR_get_frozen_gamma_sound_speed", not implemented yet..'
        STOP 

      END SUBROUTINE Argon_CR_get_frozen_gamma_sound_speed
  
      !------------------------------------------------------!
      ! This subroutine computes the mixture static pressure (presence of free electrons is accounted for)
      SUBROUTINE Argon_CR_get_pressure (rhoi, temp, p)
  
        USE mod_argon_CR_initialize_CFD,            ONLY: pos_em, pos_Te 
 
        INTEGER :: is
        REAL(KIND=8) :: p_h, p_el, T, Te
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), INTENT(OUT) :: p
       
        ! Temperatures
        T  = temp(1)
        Te = temp(pos_Te)

        ! Free electron contribution
        p_el = 0.d0
        IF (nb_te.EQ.1) p_el = rhoi(pos_em)*Re*Te
  
        ! Heavy particle contribution
        p_h = 0.d0   
        DO is = 2,nb_ns 
           p_h = p_h + rhoi(is)*Ri(is)
        ENDDO
        p_h = p_h*T
  
        ! Mixture pressure
        p = p_el + p_h
  
      END SUBROUTINE Argon_CR_get_pressure
  
      !------------------------------------------------------!
      ! This subroutine computes the mixture density.
      SUBROUTINE Argon_CR_get_density (rhoi, rho)
  
        INTEGER :: is
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), INTENT(OUT) :: rho
  
        ! Mixture density
        rho = 0.d0
        DO is = 1,nb_ns 
           rho = rho + rhoi(is)
        ENDDO
  
      END SUBROUTINE Argon_CR_get_density
  
      !------------------------------------------------------!
      ! This subroutine computes the specific gas constant R.
      SUBROUTINE Argon_CR_get_Rgas (rhoi, R)
  
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
  
      END SUBROUTINE Argon_CR_get_Rgas
  
      !------------------------------------------------------!
      ! This subroutine gives the species gas constants
      SUBROUTINE Argon_CR_get_Ri (out_Ri)
  
        INTEGER :: is
  
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: out_Ri
  
        DO is = 1,nb_ns 
           out_Ri(is) = Ri(is) 
        ENDDO
  
      END SUBROUTINE Argon_CR_get_Ri 
 
      !------------------------------------------------------!
      ! This subroutine gives the species molar masses
      SUBROUTINE Argon_CR_get_mi (out_mi)
  
        INTEGER :: is
  
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: out_mi
  
        DO is = 1,nb_ns 
           out_mi(is) = mass(is) 
        ENDDO
  
      END SUBROUTINE Argon_CR_get_mi

      !------------------------------------------------------!
      ! This subroutine computes the mixture number density 
      SUBROUTINE Argon_CR_get_mix_numb_dens (p, T, Te, xi, nb) 

        REAL(KIND=8), INTENT(IN) :: p, T, Te
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi
        REAL(KIND=8) :: nb

        nb = p/(kb*T*(1.d0 + xi(1)*(Te/T - 1.d0))) 

      END SUBROUTINE Argon_CR_get_mix_numb_dens
 
      !------------------------------------------------------!
      ! This subroutine computes the species mass fraction
      SUBROUTINE Argon_CR_get_mass_fractions (rho, rhoi, yi)
  
        INTEGER :: is
  
        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi
  
        DO is = 1,nb_ns 
           yi(is) = rhoi(is)/rho
        ENDDO
  
      END SUBROUTINE Argon_CR_get_mass_fractions

      !------------------------------------------------------!
      ! This subroutine computes the species molar fractions
      SUBROUTINE Argon_CR_get_molar_fractions (rhoi, xi)

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

      END SUBROUTINE Argon_CR_get_molar_fractions
 
      !------------------------------------------------------!
      ! This subroutine computes the molar fractions from mass fractions. 
      ! The mixture molar mass is also computed and given in output. 
      SUBROUTINE Argon_CR_molar_frac_from_mass_frac (yi, mm, xi)

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

      END SUBROUTINE Argon_CR_molar_frac_from_mass_frac 
     
      !------------------------------------------------------!
      ! This subroutine computes the mass fractions from molar fractions. 
      ! The mixture molar mass is also computed and given in output.
      SUBROUTINE Argon_CR_mass_frac_from_molar_frac (xi, mm, yi)

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

      END SUBROUTINE Argon_CR_mass_frac_from_molar_frac 

      !------------------------------------------------------!
      ! This subroutine adds a small number to the composition respecting the mass constraint
      ! (sum of species molar fractions equal to one)
      SUBROUTINE Argon_CR_comp_tol(xi)

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

      END SUBROUTINE Argon_CR_comp_tol
 
      !------------------------------------------------------!
      ! This subroutine computes the source term due to collisional and radiative processes.
      SUBROUTINE Argon_CR_get_source (rhoi, temp, s) 

        USE mod_argon_CR_function_pointer,              ONLY: get_mass_prod_terms

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s
  
        CALL get_mass_prod_terms (rhoi, temp, s)

      END SUBROUTINE Argon_CR_get_source 
  
      !------------------------------------------------------!
      ! This subroutine computes the source term due to collisional and radiative processes 
      ! and its Jacobian (with respect to primitive variables).
      SUBROUTINE Argon_CR_get_source_Jac (rhoi, temp, s, js)   
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s, js
       
        s  = 0.d0
        js = 0.d0

        PRINT*,'in "Argon_CR_get_source_Jac", not implemented yet..'
        STOP 

      END SUBROUTINE Argon_CR_get_source_Jac

      !------------------------------------------------------!
      ! This subroutine computes the species mass diffusion flux 
      SUBROUTINE Argon_CR_compute_species_DiffFlux (p, T, Te, xi, diff_driv, Ji) 
 
        REAL(KIND=8), INTENT(IN) :: p, T, Te
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, diff_driv        
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: Ji

        Ji = 0.d0

        PRINT*,'in "Argon_CR_compute_species_DiffFlux", not implemented yet..'
        STOP 

      END SUBROUTINE Argon_CR_compute_species_DiffFlux
  
      !------------------------------------------------------!
      ! This subrouotine computes the transport coefficients 
      SUBROUTINE Argon_CR_compute_transpCoeff (p, xi, temp, mu, kappa, lambda, Di, chi) 
                         
        REAL(KIND=8), INTENT(IN) :: p
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, temp
        REAL(KIND=8), INTENT(OUT) :: mu, kappa
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: chi, Di, lambda

        mu  = 0.d0
        kappa = 0.d0
        chi = 0.d0
        Di  = 0.d0
        lambda = 0.d0 

        PRINT*,'in "Argon_CR_compute_transpCoeff", not implemented yet..'
        STOP 

      END SUBROUTINE  Argon_CR_compute_transpCoeff   
#endif
  
  END MODULE mod_Argon_CR_library
!------------------------------------------------------------------------------!  
