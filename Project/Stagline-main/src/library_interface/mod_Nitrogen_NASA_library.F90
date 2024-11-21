!------------------------------------------------------------------------------!
! This module implements an interface for the Nitrogen_NASA thermodynamic library. It provides an
! implementation of functions and subroutines used to compute thermodynamic and transport properties and 
! source terms related to collisional and radiative processes.
! The interface is written such to have no data and variable dependency from the flow solver being used. 
! Therefore it can plugged to any code without big efforts. 
  MODULE mod_Nitrogen_NASA_library 

#include"../config.h"

#ifdef NITROGEN_NASA
    IMPLICIT NONE
  
    INTEGER, SAVE :: nb_ns, nb_bins, nb_tvib, nb_te, nb_trot, nb_temp, nb_post_temp, nb_eq, nb_dim 
    REAL(KIND=8), PARAMETER :: kb = 1.380658d-23
    REAL(KIND=8), PARAMETER :: na = 6.0221367d23
    REAL(KIND=8), PARAMETER :: ue = 1.602191d-19  
    REAL(KIND=8), PARAMETER :: ru = kb*na
    REAL(KIND=8), SAVE :: xi_tol
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: mass, Ri, Ek
  
    ! Subroutines of Nitrogen_NASA thermodynamic library interface. In some cases they call 
    ! other subroutines for the model under use (even these are not dependent on flow solver data). 
    CONTAINS
  
      !------------------------------------------------------!
      ! This subroutine initializes the thermodynamic library in use.  
      SUBROUTINE Nitrogen_NASA_initialize (input_ns, input_ntrot, input_ntvib, input_nte, input_ntemp, input_neq,  & 
                                         & input_ndim, input_solver, input_mixture, input_reaction, state_model,   &
                                         & input_transf, input_path, input_xi_tol) 

        USE mod_nitrogen_NASA_initialize_CFD,      ONLY: initialize, EkJ
        USE mod_nitrogen_NASA_zeta,                ONLY: set_pointer 
        
        INTEGER :: is, nt, neq, length

        INTEGER, INTENT(IN) :: input_ns, input_ntrot, input_ntvib, input_nte, input_ntemp, input_neq, input_ndim
        REAL(KIND=8), INTENT(IN) :: input_xi_tol
        CHARACTER*(*), INTENT(IN) :: input_solver, input_mixture, input_reaction, state_model, input_transf, input_path
        
        ! Parameters passed from the solver
        nb_ns   = input_ns
        nb_bins = input_ns - 1
        nb_trot = input_ntrot  
        nb_tvib = input_ntvib
        nb_te   = input_nte
        nb_temp = input_ntemp
        nb_eq   = input_neq
        nb_dim  = input_ndim
 
        length  = LEN_TRIM(input_solver)

        ! Consistency checking on number of temperatures and conservation equations
        nt = 1 + input_ntrot + input_ntvib + input_nte
        IF (nt.NE.nb_temp) THEN
           PRINT*
           WRITE(*,5)input_solver(1:length),'::Nitrogen_NASA_library interface'   
           PRINT*
           WRITE(*,5)input_solver(1:length),'::Error in number of temperatures'     
           PRINT*
           STOP
        ENDIF

        neq = input_ns + 1 + input_ntrot + input_ntvib + input_nte + input_ndim
        IF (neq.NE.nb_eq) THEN
           PRINT*
           WRITE(*,5)input_solver(1:length),'::Nitrogen_NASA_library interface'   
           PRINT*
           WRITE(*,5)input_solver(1:length),'::Error in number of equations'     
           PRINT*
           STOP
        ENDIF
 
        ! Tolerance on species molar fractions
        xi_tol = input_xi_tol

        ! Number of post-processing temperatures (the minimum is always one)
        nb_post_temp = MAX(1,nb_trot + nb_tvib + nb_te)

        ! Allocation of molecular mass and specific gas constatn vectors
        ALLOCATE(mass(nb_ns), Ri(nb_ns))
        
        ! Nitrogen_NASA library is initialized
        CALL initialize (input_solver, input_mixture, input_reaction, input_path, nb_ns, nb_trot, nb_tvib, & 
       &                 nb_dim, nb_eq, mass)
        
        ! Speficic gas constant of single species
        DO is = 1,nb_ns 
           Ri(is) = 1.d0/mass(is)
        ENDDO 
       
        Ri = ru*Ri

        ! Bin energy levels [ev]
        IF(ALLOCATED(EkJ)) THEN 

           ALLOCATE(Ek(nb_bins))

           DO is = 1,nb_bins
              Ek(is) = EkJ(is) 
           ENDDO
           Ek = EkJ/ue

        ENDIF

        ! Pointer initialization
        CALL set_pointer ()

5     FORMAT(A,A)
 
      END SUBROUTINE Nitrogen_NASA_initialize
 
      !------------------------------------------------------!
      ! This subroutine finalizes the Nitrogen NASA thermodynamic library
      SUBROUTINE Nitrogen_NASA_finalize

         USE mod_nitrogen_NASA_initialize_CFD,            ONLY: finalize

         IF (ALLOCATED(mass))         DEALLOCATE(mass)
         IF (ALLOCATED(Ri))           DEALLOCATE(Ri)
         IF (ALLOCATED(Ek))           DEALLOCATE(Ek)

         CALL finalize ()

      END SUBROUTINE Nitrogen_NASA_finalize

      !------------------------------------------------------!
      ! This subroutine computes post-shock nonequilibrium conditions
      SUBROUTINE Nitrogen_NASA_post_shock_neq (p1, u1, T1, p2, u2, T2, yi, xN)

        USE mod_nitrogen_NASA_CFD_prop,                  ONLY: post_shock_neq

        REAL(KIND=8), INTENT(IN) :: p1, u1, T1
        REAL(KIND=8), INTENT(IN), OPTIONAL :: xN
        REAL(KIND=8), INTENT(OUT) :: p2, u2, T2
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi 

        CALL post_shock_neq (p1, u1, T1, p2, u2, T2, yi, xN)

      END SUBROUTINE Nitrogen_NASA_post_shock_neq

      !------------------------------------------------------!
      ! This subroutine computes post-shock equilibrium conditions
      SUBROUTINE Nitrogen_NASA_post_shock_eq (p1, u1, T1, p2, u2, T2, yi, rhoi, xN)

        USE mod_nitrogen_NASA_CFD_eq,                  ONLY: post_shock_eq

        REAL(KIND=8), INTENT(IN) :: p1, u1, T1
        REAL(KIND=8), INTENT(IN), OPTIONAL :: xN
        REAL(KIND=8), INTENT(OUT) :: p2, u2, T2
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi, rhoi

        CALL post_shock_eq (p1, u1, T1, p2, u2, T2, yi, rhoi, xN)

      END SUBROUTINE Nitrogen_NASA_post_shock_eq
 
      !----------------------------------------------------!
      ! This subroutine computes common factors for the ODE system in 
      ! case of thermal nonequilibrium    
      SUBROUTINE Nitrogen_NASA_compute_ODE_factors_tcneq (temp, omega, yi, fac_num, fac_den)

        USE mod_nitrogen_NASA_CFD_prop,                ONLY: get_neq_energy_cv_MT_TTint

        REAL(KIND=8) :: Tint
        REAL(KIND=8) :: cv_int, e_int

        REAL(KIND=8), INTENT(IN), DIMENSION(:) :: temp, omega, yi
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: fac_num, fac_den

        ! Internal temperature
        Tint = temp(2)
 
        ! Nonequilibrium internal energy and constant volume specific heat
        CALL get_neq_energy_cv_MT_TTint (Tint, e_int, cv_int)

        fac_num = omega(2)*e_int
        fac_den = yi(2)*cv_int

      END SUBROUTINE Nitrogen_NASA_compute_ODE_factors_tcneq

      !------------------------------------------------------!
      ! This subroutine writes down the solution when solving the ODE system for
      ! studying the nonequilibrium flow behind a shock wave
      SUBROUTINE Nitrogen_NASA_write_ODE_solution (flow_file, pop_file, k, xold, x, uold, m_dot, z)
 
        USE mod_nitrogen_NASA_CFD_prop,            ONLY: Q_int_bins 

        INTEGER :: is
        REAL(KIND=8) :: mm, fac
        REAL(KIND=8) :: u, p, rho
        REAL(KIND=8), SAVE :: time
        REAL(KIND=8), DIMENSION(2) :: xi
        REAL(KIND=8), DIMENSION(nb_temp) :: temp
        REAL(KIND=8), DIMENSION(nb_bins) :: Qk

        INTEGER, INTENT(IN) :: flow_file, pop_file, k
        REAL(KIND=8), INTENT(IN) :: xold, x, m_dot, uold
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: z

        ! Species molar fractions (N-N2)
        mm = 0.d0
        DO is = 1,nb_ns 
           mm = mm + z(is)/mass(is)
        ENDDO 
        mm = 1.d0/mm

        xi(1) = z(1)/mass(1)
        xi(2) = 0.d0
        DO is = 2,nb_ns 
           xi(2) = xi(2) + z(is)/mass(is) 
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
           WRITE(pop_file,5)nb_ns - 1
        ELSE 
           time = time + 0.5d0*(x - xold)*(1.d0/uold + 1.d0/u)
        ENDIF

        ! Flowfield output file
        WRITE(flow_file,10)x,time,temp,xi,u,rho,p
        CALL flush(flow_file)

        ! Bin internal partition functions
        CALL Q_int_bins (temp(1), Qk)

        ! Population file 
        fac = rho*na/mass(2)
        DO is = 1,nb_bins 
           WRITE(pop_file,10)Ek(is),fac*z(is + 1)/Qk(is)
           CALL flush(pop_file)  
        ENDDO 
        WRITE(pop_file,20)'& ',k,time,x

5     FORMAT(I5)
10    FORMAT(100E20.10)
20    FORMAT(A,I6,100E20.10)

      END SUBROUTINE Nitrogen_NASA_write_ODE_solution  

      !------------------------------------------------------!
      ! This subroutine post-process the ODE solution for the flow 
      ! behind a normal shock wave
      SUBROUTINE Nitrogen_NASA_post_process_ODE_solution (post_file, x, ni, temp) 

        USE mod_nitrogen_NASA_initialize_CFD,      ONLY: nb_bins
        USE mod_function_pointer_NASA,             ONLY: get_Tint_NASA
        USE mod_nitrogen_NASA_CFD_prop,            ONLY: Q_int_bins

        INTEGER :: is
        REAL(KIND=8) :: T
        REAL(KIND=8), DIMENSION(nb_bins) :: Qk, dens
        REAL(KIND=8), DIMENSION(nb_post_temp + 1) :: temp_out

        INTEGER, INTENT(IN) :: post_file
        REAL(KIND=8), INTENT(IN) :: x
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: ni, temp
        
        ! Translational temperature 
        T = temp(1)

        CALL Q_int_bins (T, Qk)

        ! Population input data is normalized (ni/Qi)
        DO is = 1,nb_bins 
           dens(is) = ni(is)*Qk(is)
        ENDDO

        ! Internal temperature(s)
        CALL get_Tint_NASA (dens, temp, temp_out)

        ! Post-processing file file
        WRITE(post_file,10)x,temp_out
        CALL flush(post_file)

10    FORMAT(100E20.10)

      END SUBROUTINE Nitrogen_NASA_post_process_ODE_solution 

      !------------------------------------------------------!
      ! This subroutine writes down the solution when solving the flow equations 
      ! by means of Finite Volume method
      SUBROUTINE Nitrogen_NASA_write_fvmcc_solution (flow_file, pop_file, k, xold, x, uold, u, rhoi, temp)

        USE mod_nitrogen_NASA_initialize_CFD,      ONLY: model
        USE mod_function_pointer_NASA,             ONLY: get_Tint_NASA
        USE mod_nitrogen_NASA_CFD_prop,            ONLY: Q_int_bins 

        INTEGER :: is
        REAL(KIND=8) :: mm, fac
        REAL(KIND=8) :: p, rho, T
        REAL(KIND=8), SAVE :: time
        REAL(KIND=8), DIMENSION(2) :: xi
        REAL(KIND=8), DIMENSION(nb_post_temp + 1) :: temp_out
        REAL(KIND=8), DIMENSION(nb_ns) :: yi
        REAL(KIND=8), DIMENSION(nb_bins) :: ni, Qi

        INTEGER, INTENT(IN) :: flow_file, pop_file, k
        REAL(KIND=8), INTENT(IN) :: xold, x, uold, u
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp

        ! Density 
        CALL Nitrogen_NASA_get_density (rhoi, rho)

        ! Species mass fractions
        CALL Nitrogen_NASA_get_mass_fractions (rho, rhoi, yi)
 
        ! Species molar fractions (N-N2)
        mm = 0.d0
        DO is = 1,nb_ns 
           mm = mm + yi(is)/mass(is)
        ENDDO 
        mm = 1.d0/mm

        xi(1) = yi(1)/mass(1)
        xi(2) = 0.d0
        DO is = 2,nb_ns 
           xi(2) = xi(2) + yi(is)/mass(is) 
        ENDDO
        xi = xi*mm

        ! Translational temperature
        T = temp(1)

        ! Bin internal partition functions
        CALL Q_int_bins (T, Qi)

        ! Bin number densities
        fac = rho*na/mass(2)
        DO is = 1,nb_bins 
           ni(is) = fac*yi(is + 1)
        ENDDO

        ! Internal temperature(s)
        CALL get_Tint_NASA (ni, temp, temp_out)

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
        WRITE(flow_file,10)x,time,temp_out,xi,u,rho,p
        CALL flush(flow_file)

        ! Population file
        DO is = 1,nb_bins 
           WRITE(pop_file,10)Ek(is),ni(is)/Qi(is)
           CALL flush(pop_file)  
        ENDDO 
        WRITE(pop_file,20)'& ',k,time,x

10    FORMAT(100E20.10)
20    FORMAT(A,I6,100E20.10)

      END SUBROUTINE Nitrogen_NASA_write_fvmcc_solution  

      !------------------------------------------------------!
      ! This subroutine computes the equilibrium composition given pressure and temperature
      SUBROUTINE Nitrogen_NASA_compute_eq_composition (p, T, rhoi)

        USE mod_nitrogen_NASA_CFD_eq,           ONLY: eq_composition_bins

        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: rhoi

        CALL eq_composition_bins (p, T, rhoi)

      END SUBROUTINE Nitrogen_NASA_compute_eq_composition

      !------------------------------------------------------!
      ! This subroutine computes the macroscopic equilibrium composition
      SUBROUTINE Nitrogen_NASA_compute_eq_composition_macro (p, T, rhoi)

        USE mod_nitrogen_NASA_CFD_eq,           ONLY: eq_composition

        INTEGER :: is
        REAL(KIND=8) :: xn, xn2
        REAL(KIND=8) :: mm, R, rho
        REAL(KIND=8), DIMENSION(nb_ns) :: yi

        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: rhoi

        CALL eq_composition (p, T, xn, xn2) 

        ! Mixture molecular mass
        mm = xn*mass(1) + xn2*mass(2)

        R = ru/mm

        yi(1) = xn*mass(1)/mm
        yi(2) = xn2*mass(2)/mm

        rho = p/(R*T)

        ! Species densities
        DO is = 1,nb_ns 
           rhoi(is) = rho*yi(is)
        ENDDO

      END SUBROUTINE Nitrogen_NASA_compute_eq_composition_macro 

      !------------------------------------------------------!
      ! This subroutine computes the N-N2 equilibrium internal energy per unit mass 
      SUBROUTINE Nitrogen_NASA_compute_eq_species_energy (T, e)

        USE mod_nitrogen_NASA_CFD_eq,           ONLY: energy

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: e

         CALL energy (T, e)

      END SUBROUTINE Nitrogen_NASA_compute_eq_species_energy 

      !------------------------------------------------------!
      ! This subroutine computes the N-N2 equilibrium enthalpy per unit mass 
      SUBROUTINE Nitrogen_NASA_compute_eq_species_enthalpy (T, h)

        USE mod_nitrogen_NASA_CFD_eq,           ONLY: enthalpy

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: h

        CALL enthalpy (T, h)

      END SUBROUTINE Nitrogen_NASA_compute_eq_species_enthalpy 

      !------------------------------------------------------!
      ! This subroutine computes the N-N2 equilibrium entropy per unit mass. The
      ! entropy of mixing contributions is computed separately. 
      SUBROUTINE Nitrogen_NASA_compute_eq_species_entropy (p, T, yi, xi, smix, s)

        USE mod_nitrogen_NASA_CFD_eq,           ONLY: entropy

        REAL(KIND=8) :: tmp1
        REAL(KIND=8) :: QtN, QtN2

        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), INTENT(IN), DIMENSION(:) :: xi, yi
        REAL(KIND=8), INTENT(OUT) :: smix
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: s

        CALL entropy (p, T, yi, xi, smix, s)

      END SUBROUTINE Nitrogen_NASA_compute_eq_species_entropy   

      !------------------------------------------------------!
      ! This subroutine computes the temperatures and thermodynamic data,
      ! additional data for Jacobians are given in output
      SUBROUTINE Nitrogen_NASA_get_data (rho, rhoi, rho_eint, temp, c, gamma, p, alpha, beta, ei)
  
        USE mod_function_pointer_NASA,        ONLY: get_temp_NASA, get_species_cv_NASA, & 
                                                  & get_species_energy_NASA  
 
        IMPLICIT NONE
  
        INTEGER :: is 
        REAL(KIND=8) :: sum1, sum2,  ki, T, tmp       
        REAL(KIND=8), DIMENSION(nb_ns) :: cv  
  
        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, rho_eint
        REAL(KIND=8), INTENT(OUT) :: c, gamma, p, alpha
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp, ei, beta
  
        ! Temperature computation
        CALL get_temp_NASA (rhoi, rho_eint, temp)
        T = temp(1)
        CALL get_species_cv_NASA (T, cv)   
 
        ! Initialization
        sum1 = 0.d0   
        sum2 = 0.d0      
        DO is = 1,nb_ns 
           tmp  = rhoi(is)
           sum1 = sum1 + tmp*Ri(is)
           sum2 = sum2 + tmp*cv(is)
        ENDDO
  
        p  = sum1*T 
        ki = sum1/sum2
  
        gamma = 1.d0 + ki
        c     = DSQRT(gamma*p/rho)
  
        alpha   = ki 
        beta(1) = sum2
      
        CALL get_species_energy_NASA (temp, ei)
  
      END SUBROUTINE Nitrogen_NASA_get_data
  
      !------------------------------------------------------!
      ! This subroutine computes thermodynamic data needed to be stored
      ! for nonequilibrium gas flow computations. 
      SUBROUTINE Nitrogen_NASA_get_thermodynamic_data (rho, rhoi, temp, c, gamma, p, alpha, beta, ei)

        USE mod_function_pointer_NASA,        ONLY: get_species_cv_NASA, get_species_energy_NASA
  
        INTEGER :: is 
        REAL(KIND=8) :: sum1, sum2,  ki, T, tmp 
        REAL(KIND=8), DIMENSION(nb_ns) :: cv  
  
        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), INTENT(OUT) :: c, gamma, p, alpha
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ei, beta
  
        ! Temperature 
        T = temp(1)
  
        ! Frozen specific heats
        CALL get_species_cv_NASA (T, cv)

        ! Initialization
        sum1 = 0.d0   
        sum2 = 0.d0      
        DO is = 1,nb_ns 
           tmp  = rhoi(is)
           sum1 = sum1 + tmp*Ri(is)
           sum2 = sum2 + tmp*cv(is)
        ENDDO
  
        p  = sum1*T 
        ki = sum1/sum2
  
        gamma = 1.d0 + ki
        c     = DSQRT(gamma*p/rho)
  
        alpha   = ki 
        beta(1) = sum2
      
        CALL get_species_energy_NASA (temp, ei)

      END SUBROUTINE Nitrogen_NASA_get_thermodynamic_data
 
      !------------------------------------------------------!
      ! This subroutine computes the energy density for total internal,
      ! rotational, vibrational and electronic components.
      SUBROUTINE Nitrogen_NASA_get_energy_densities (rhoi, temp, rho_e)
  
        USE mod_function_pointer_NASA,     ONLY: get_species_energy_NASA
  
        INTEGER :: is
        REAL(KIND=8) :: tmp
        REAL(KIND=8), DIMENSION(nb_ns) :: eint
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: rho_e

        ! Temperature         
        CALL get_species_energy_NASA (temp, eint)
        
        ! Internal energy per unit volume of the mixture
        tmp = 0.d0
        DO is = 1,nb_ns 
           tmp = tmp + rhoi(is)*eint(is)
        ENDDO
        rho_e(1) = tmp 
 
      END SUBROUTINE Nitrogen_NASA_get_energy_densities
  
      !------------------------------------------------------!
      ! This subroutine computes the total internal, rotational, vibrational and electronic energy 
      ! per unit mass of the mixture.
      SUBROUTINE Nitrogen_NASA_get_energy (rho, rhoi, temp, e)
  
        USE mod_function_pointer_NASA,     ONLY: get_species_energy_NASA  
  
        INTEGER :: is
        REAL(KIND=8) :: tmp
        REAL(KIND=8), DIMENSION(nb_ns) :: eint, yi
  
        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e
       
        ! Species mass fractions
        CALL Nitrogen_NASA_get_mass_fractions (rho, rhoi, yi)
  
        CALL get_species_energy_NASA (temp, eint)
  
        ! Energy per unit mass of the mixture
        tmp = 0.d0
        DO is = 1,nb_ns 
           tmp = tmp + yi(is)*eint(is)
        ENDDO
        e(1) = tmp  

      END SUBROUTINE Nitrogen_NASA_get_energy
  
      !------------------------------------------------------!
      ! This subroutine computes the total internal, rotational, vibrational and electronic enthalpy 
      ! per unit mass of the mixture.
      SUBROUTINE Nitrogen_NASA_get_enthalpy (rho, rhoi, temp, h)
  
        USE mod_function_pointer_NASA,     ONLY: get_species_energy_NASA   
  
        INTEGER :: is
        REAL(KIND=8) :: tmp, T
        REAL(KIND=8), DIMENSION(nb_ns) :: yi, eint
  
        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: h
       
        ! Species mass fractions
        CALL Nitrogen_NASA_get_mass_fractions (rho, rhoi, yi)
  
        T = temp(1)
        CALL get_species_energy_NASA (temp, eint)
  
        ! Static enthalpy per unit mass of the mixture
        tmp = 0.d0
        DO is = 1,nb_ns 
           tmp = tmp + yi(is)*(eint(is) + Ri(is)*T)
        ENDDO
        h(1) = tmp
  
      END SUBROUTINE Nitrogen_NASA_get_enthalpy 
  
      !------------------------------------------------------!
      ! This subroutine computes the species internal energy per unit mass  
      SUBROUTINE Nitrogen_NASA_get_species_energy (temp, e)
  
        USE mod_function_pointer_NASA,     ONLY: get_species_energy_NASA
           
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e
 
        CALL get_species_energy_NASA (temp, e)
  
      END SUBROUTINE Nitrogen_NASA_get_species_energy

      !------------------------------------------------------!
      ! This subroutine computes the species enthalpy per unit mass 
      SUBROUTINE Nitrogen_NASA_get_species_enthalpy (temp, h)
  
        USE mod_function_pointer_NASA,     ONLY: get_species_energy_NASA
           
        INTEGER :: is
        REAL(KIND=8) :: T
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: h
 
        T = temp(1)
        CALL get_species_energy_NASA (temp, h)

        DO is = 1,nb_ns 
           h(is) = h(is) + Ri(is)*T
        ENDDO
  
      END SUBROUTINE Nitrogen_NASA_get_species_enthalpy
  
      !------------------------------------------------------!
      ! This subroutine computes the species internal energy and constant 
      ! volume specific heat per unit mass  
      SUBROUTINE Nitrogen_NASA_get_species_energy_cv (temp, e, cv)
  
        USE mod_function_pointer_NASA,     ONLY: get_species_energy_cv_NASA
           
        INTEGER :: is
        REAL(KIND=8) :: T
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e, cv
 
        T = temp(1)
        CALL get_species_energy_cv_NASA (T, e, cv)
  
      END SUBROUTINE Nitrogen_NASA_get_species_energy_cv

      !------------------------------------------------------!
      ! This subroutine computes the species enthalpy and constant pressure
      ! specific heat per unit mass
      SUBROUTINE Nitrogen_NASA_get_species_enthalpy_cp (temp, h, cp)
  
        USE mod_function_pointer_NASA,     ONLY: get_species_energy_cv_NASA
           
        INTEGER :: is
        REAL(KIND=8) :: T, tmp1
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: h, cp
 
        T = temp(1)
        CALL get_species_energy_cv_NASA (T, h, cp)
  
        DO is = 1,nb_ns 
           tmp1   = Ri(is)
           h(is)  = h(is) + tmp1*T
           cp(is) = cp(is) + tmp1
        ENDDO

      END SUBROUTINE Nitrogen_NASA_get_species_enthalpy_cp

      !------------------------------------------------------!
      ! This subroutine computes the frozen specific heats at constant volume of each species
      SUBROUTINE Nitrogen_NASA_get_frozen_cv (T, cv)  
  
        USE mod_function_pointer_NASA,     ONLY: get_species_cv_NASA

        INTEGER :: is

        REAL(KIND=8), INTENT(IN) :: T 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv
  
        CALL get_species_cv_NASA (T, cv)   
  
      END SUBROUTINE Nitrogen_NASA_get_frozen_cv 
  
      !------------------------------------------------------!
      ! This subroutine computes the frozen specific heats at constant pressure of each species
      SUBROUTINE Nitrogen_NASA_get_frozen_cp (T, cp)  
  
        USE mod_function_pointer_NASA,     ONLY: get_species_cv_NASA 

        INTEGER :: is

        REAL(KIND=8), INTENT(IN) :: T 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cp
  
        CALL get_species_cv_NASA (T, cp)
  
        DO is = 1,nb_ns 
           cp(is) = cp(is) + Ri(is)
        ENDDO   
  
      END SUBROUTINE Nitrogen_NASA_get_frozen_cp

      !------------------------------------------------------!
      ! This subroutine computes the frozen specific heats at constant volume of each species and 
      ! the species energy components in thermal equilibrium with translation.
      SUBROUTINE Nitrogen_NASA_get_species_frozen_energy_cv (temp, e, cv)  
 
        USE mod_function_pointer_NASA,      ONLY: get_species_energy_cv_NASA   

        REAL(KIND=8) :: T
 
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv, e
  
        T = temp(1)
        CALL get_species_energy_cv_NASA (T, e, cv)

      END SUBROUTINE Nitrogen_NASA_get_species_frozen_energy_cv 
  
      !------------------------------------------------------!
      ! This subroutine computes the frozen specific heats at constant pressure of each species and 
      ! the species enthalpy components in thermal equilibrium with translation.
      SUBROUTINE Nitrogen_NASA_get_species_frozen_enthalpy_cp (temp, h, cp)  
  
        USE mod_function_pointer_NASA,      ONLY: get_species_energy_cv_NASA   

        INTEGER :: i

        REAL(KIND=8) :: T
        REAL(KIND=8) :: tmp

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cp, h

        T = temp(1)
        CALL get_species_energy_cv_NASA (T, h, cp)

        DO i = 1,nb_ns 
           tmp   = Ri(i)
           h(i)  = h(i) + tmp*T
           cp(i) = cp(i) + tmp
        ENDDO 

      END SUBROUTINE Nitrogen_NASA_get_species_frozen_enthalpy_cp

      !------------------------------------------------------!
      ! This subroutine computes the equilibrium speed of sound
      SUBROUTINE Nitrogen_NASA_get_eq_sound_speed (rhoi, p, T, c)
 
        USE mod_Nitrogen_NASA_CFD_eq,    ONLY: energy

        INTEGER :: is
        REAL(KIND=8), PARAMETER :: eps = 1.d-5
        REAL(KIND=8) :: tmp1, tmp2
        REAL(KIND=8) :: e, e_pert, h, h_pert
        REAL(KIND=8) :: rho, rho_pert, T_pert, p_pert
        REAL(KIND=8) :: drho_dp
        REAL(KIND=8) :: gamma
        REAL(KIND=8), DIMENSION(nb_ns) :: yi, ei, dens

        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), INTENT(IN), DIMENSION(:) :: rhoi
        REAL(KIND=8), INTENT(OUT) :: c

        ! Density and mass fractions
        CALL Nitrogen_NASA_get_density (rhoi, rho)
        CALL Nitrogen_NASA_get_mass_fractions (rho, rhoi, yi)

        CALL energy (T, ei)

        e = 0.d0
        h = 0.d0
        DO is = 1,nb_ns
           tmp1 = yi(is) 
           tmp2 = ei(is)
           e    = e + tmp1*tmp2
           h    = h + tmp1*(tmp2 + Ri(is)*T)
        ENDDO

        ! Derivative of density with respect to pressure 
        p_pert = p*(1.d0 + eps)
        
        CALL Nitrogen_NASA_compute_eq_composition_macro (p_pert, T, dens)
        CALL Nitrogen_NASA_get_density (dens, rho_pert)

        drho_dp = (rho_pert - rho)/(eps*p)

        ! Equilibrium specific heat ratio
        T_pert = T*(1.d0 + eps)
        
        CALL Nitrogen_NASA_compute_eq_composition_macro (p, T_pert, dens)
        CALL Nitrogen_NASA_get_density (dens, rho_pert)
        CALL Nitrogen_NASA_get_mass_fractions (rho_pert, dens, yi)

        CALL energy (T_pert, ei)

        e_pert = 0.d0
        h_pert = 0.d0
        DO is = 1,nb_ns
           tmp1   = yi(is) 
           tmp2   = ei(is)
           e_pert = e_pert + tmp1*tmp2
           h_pert = h_pert + tmp1*(tmp2 + Ri(is)*T_pert)
        ENDDO

        gamma = (h_pert - h)/(e_pert - e)
        c     = DSQRT(gamma/drho_dp)

      END SUBROUTINE Nitrogen_NASA_get_eq_sound_speed 

      !------------------------------------------------------!
      ! This subroutine computes the equilibrium specific heat ratio
      SUBROUTINE Nitrogen_NASA_get_eq_gamma (rhoi, p, T, gamma)
 
        USE mod_Nitrogen_NASA_CFD_eq,    ONLY: energy

        INTEGER :: is
        REAL(KIND=8), PARAMETER :: eps = 1.d-5
        REAL(KIND=8) :: tmp1, tmp2
        REAL(KIND=8) :: e, e_pert, h, h_pert
        REAL(KIND=8) :: rho, rho_pert, T_pert
        REAL(KIND=8), DIMENSION(nb_ns) :: yi, ei, dens

        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), INTENT(IN), DIMENSION(:) :: rhoi
        REAL(KIND=8), INTENT(OUT) :: gamma

        ! Density and mass fractions
        CALL Nitrogen_NASA_get_density (rhoi, rho)
        CALL Nitrogen_NASA_get_mass_fractions (rho, rhoi, yi)

        CALL energy (T, ei)

        e = 0.d0
        h = 0.d0
        DO is = 1,nb_ns
           tmp1 = yi(is) 
           tmp2 = ei(is)
           e    = e + tmp1*tmp2
           h    = h + tmp1*(tmp2 + Ri(is)*T)
        ENDDO

        ! Equilibrium specific heat ratio
        T_pert = T*(1.d0 + eps)
        
        CALL Nitrogen_NASA_compute_eq_composition_macro (p, T_pert, dens)
        CALL Nitrogen_NASA_get_density (dens, rho_pert)
        CALL Nitrogen_NASA_get_mass_fractions (rho_pert, dens, yi)

        CALL energy (T_pert, ei)

        e_pert = 0.d0
        h_pert = 0.d0
        DO is = 1,nb_ns
           tmp1   = yi(is) 
           tmp2   = ei(is)
           e_pert = e_pert + tmp1*tmp2
           h_pert = h_pert + tmp1*(tmp2 + Ri(is)*T_pert)
        ENDDO

        gamma = (h_pert - h)/(e_pert - e)

      END SUBROUTINE Nitrogen_NASA_get_eq_gamma 

      !------------------------------------------------------!
      ! This subroutine computes the equilibrium specific heat ratio
      SUBROUTINE Nitrogen_NASA_get_eq_gamma_sound_speed (rhoi, p, T, gamma, c)
 
         USE mod_Nitrogen_NASA_CFD_eq,    ONLY: energy

        INTEGER :: is
        REAL(KIND=8), PARAMETER :: eps = 1.d-5
        REAL(KIND=8) :: tmp1, tmp2
        REAL(KIND=8) :: e, e_pert, h, h_pert
        REAL(KIND=8) :: rho, rho_pert, T_pert, p_pert
        REAL(KIND=8) :: drho_dp
        REAL(KIND=8), DIMENSION(nb_ns) :: yi, ei, dens

        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), INTENT(IN), DIMENSION(:) :: rhoi
        REAL(KIND=8), INTENT(OUT) :: gamma, c

        ! Density and mass fractions
        CALL Nitrogen_NASA_get_density (rhoi, rho)
        CALL Nitrogen_NASA_get_mass_fractions (rho, rhoi, yi)

        CALL energy (T, ei)

        e = 0.d0
        h = 0.d0
        DO is = 1,nb_ns
           tmp1 = yi(is) 
           tmp2 = ei(is)
           e    = e + tmp1*tmp2
           h    = h + tmp1*(tmp2 + Ri(is)*T)
        ENDDO

        ! Derivative of density with respect to pressure 
        p_pert = p*(1.d0 + eps)
        
        CALL Nitrogen_NASA_compute_eq_composition_macro (p_pert, T, dens)
        CALL Nitrogen_NASA_get_density (dens, rho_pert)

        drho_dp = (rho_pert - rho)/(eps*p)

        ! Equilibrium specific heat ratio
        T_pert = T*(1.d0 + eps)
        
        CALL Nitrogen_NASA_compute_eq_composition_macro (p, T_pert, dens)
        CALL Nitrogen_NASA_get_density (dens, rho_pert)
        CALL Nitrogen_NASA_get_mass_fractions (rho_pert, dens, yi)

        CALL energy (T_pert, ei)

        e_pert = 0.d0
        h_pert = 0.d0
        DO is = 1,nb_ns
           tmp1   = yi(is) 
           tmp2   = ei(is)
           e_pert = e_pert + tmp1*tmp2
           h_pert = h_pert + tmp1*(tmp2 + Ri(is)*T_pert)
        ENDDO

        gamma = (h_pert - h)/(e_pert - e)
        c     = DSQRT(gamma/drho_dp)

      END SUBROUTINE Nitrogen_NASA_get_eq_gamma_sound_speed 

      !------------------------------------------------------!
      ! This subroutine computes the temperatures given the energy densities
      SUBROUTINE Nitrogen_NASA_get_temperatures (rhoi, eint, temp)
  
        USE mod_function_pointer_NASA,     ONLY: get_temp_NASA
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: eint
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp
  
        CALL get_temp_NASA (rhoi, eint, temp) 

      END SUBROUTINE Nitrogen_NASA_get_temperatures
  
      !------------------------------------------------------!
      ! This subroutine computes the frozen speed of sound
      SUBROUTINE Nitrogen_NASA_get_frozen_sound_speed (rhoi, temp, c)
 
        USE mod_function_pointer_NASA,     ONLY: get_species_cv_NASA
 
        INTEGER :: is
        REAL(KIND=8) :: ki, sum1, sum2, R, gamma, rho, tmp, T
        REAL(KIND=8), DIMENSION(nb_ns) :: cv  
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp, rhoi
        REAL(KIND=8), INTENT(OUT) :: c 
  
        ! Temperature 
        T = temp(1)

        ! Frozen specific heats
        CALL get_species_cv_NASA (T, cv)
        CALL Nitrogen_NASA_get_density (rhoi, rho)
  
        ! Initialization
        sum1 = 0.d0   
        sum2 = 0.d0   
        DO is = 1,nb_ns 
           tmp  = rhoi(is)
           sum1 = sum1 + tmp*Ri(is)
           sum2 = sum2 + tmp*cv(is)
        ENDDO
        R = sum1/rho
  
        ki    = sum1/sum2
        gamma = 1.d0 + ki
        c     = DSQRT(gamma*R*T)
  
      END SUBROUTINE Nitrogen_NASA_get_frozen_sound_speed
 
      !------------------------------------------------------!
      ! This subroutine gets the frozen speed specific heat ratio.
      SUBROUTINE Nitrogen_NASA_get_frozen_gamma (rhoi, temp, gamma)
  
        USE mod_function_pointer_NASA,     ONLY:  get_species_cv_NASA

        INTEGER :: is
        REAL(KIND=8) :: ki, sum1, sum2, rho, tmp, T
        REAL(KIND=8), DIMENSION(nb_ns) :: cv  
 
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), INTENT(OUT) :: gamma 
  
        ! Temperature 
        T = temp(1)

        ! Frozen specific heats 
        CALL get_species_cv_NASA (T, cv)
   
        ! Inizialization
        sum1 = 0.d0   
        sum2 = 0.d0 
        rho  = 0.d0 
        DO is = 1,nb_ns 
           tmp  = rhoi(is)
           sum1 = sum1 + tmp*Ri(is)
           sum2 = sum2 + tmp*cv(is)
        ENDDO
  
        ki    = sum1/sum2
        gamma = 1.d0 + ki
  
      END SUBROUTINE Nitrogen_NASA_get_frozen_gamma
 
      !------------------------------------------------------!
      ! This subroutine gets the frozen speed of sound and specific heat ratio.
      SUBROUTINE Nitrogen_NASA_get_frozen_gamma_sound_speed (rhoi, temp, gamma, c)
 
        USE mod_function_pointer_NASA,     ONLY: get_species_cv_NASA
 
        INTEGER :: is
        REAL(KIND=8) :: ki, sum1, sum2, rho, p_ov_rho, tmp, T
        REAL(KIND=8), DIMENSION(nb_ns) :: cv  
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), INTENT(OUT) :: gamma, c 
  
        ! Temperature
        T = temp(1)

        ! Frozen specific heats 
        CALL get_species_cv_NASA (T, cv)
   
        ! Inizialization
        sum1 = 0.d0   
        sum2 = 0.d0  
        rho  = 0.d0
        DO is = 1,nb_ns 
           tmp  = rhoi(is)
           rho  = rho + tmp
           sum1 = sum1 + tmp*Ri(is)
           sum2 = sum2 + tmp*cv(is)
        ENDDO
        p_ov_rho = sum1*T/rho
 
        ki    = sum1/sum2
        gamma = 1.d0 + ki
        c     = DSQRT(gamma*p_ov_rho)
  
      END SUBROUTINE Nitrogen_NASA_get_frozen_gamma_sound_speed
  
      !------------------------------------------------------!
      ! This subroutine computes the mixture static pressure (presence of free electrons is accounted for)
      SUBROUTINE Nitrogen_NASA_get_pressure (rhoi, temp, p)
   
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
  
      END SUBROUTINE Nitrogen_NASA_get_pressure
  
      !------------------------------------------------------!
      ! This subroutine computes the mixture density.
      SUBROUTINE Nitrogen_NASA_get_density (rhoi, rho)
  
        INTEGER :: is
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), INTENT(OUT) :: rho
  
        ! Mixture density
        rho = 0.d0
        DO is = 1,nb_ns 
           rho = rho + rhoi(is)
        ENDDO
  
      END SUBROUTINE Nitrogen_NASA_get_density
  
      !------------------------------------------------------!
      ! This subroutine computes the specific gas constant R.
      SUBROUTINE Nitrogen_NASA_get_Rgas (rhoi, R)
  
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
  
      END SUBROUTINE Nitrogen_NASA_get_Rgas
  
      !------------------------------------------------------!
      ! This subroutine gives the single species gas constants
      SUBROUTINE Nitrogen_NASA_get_Ri (out_Ri)
  
        INTEGER :: is
  
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: out_Ri
  
        DO is = 1,nb_ns 
           out_Ri(is) = Ri(is) 
        ENDDO
  
      END SUBROUTINE Nitrogen_NASA_get_Ri 
  
      !------------------------------------------------------!
      ! This subroutine gives the species molar masses
      SUBROUTINE Nitrogen_NASA_get_mi (out_mi)
  
        INTEGER :: is
  
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: out_mi
  
        DO is = 1,nb_ns 
           out_mi(is) = mass(is) 
        ENDDO
  
      END SUBROUTINE Nitrogen_NASA_get_mi

      !------------------------------------------------------!
      ! This subroutine computes the mixture number density 
      SUBROUTINE Nitrogen_NASA_get_mix_numb_dens (p, T, Te, xi, nb) 

        REAL(KIND=8), INTENT(IN) :: p, T, Te
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi
        REAL(KIND=8) :: nb

        nb = p/(kb*T)

      END SUBROUTINE Nitrogen_NASA_get_mix_numb_dens

      !------------------------------------------------------!
      ! This subroutine computes the species mass fraction
      SUBROUTINE Nitrogen_NASA_get_mass_fractions (rho, rhoi, yi)
  
        INTEGER :: is
  
        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi
  
        DO is = 1,nb_ns 
           yi(is) = rhoi(is)/rho
        ENDDO
  
      END SUBROUTINE Nitrogen_NASA_get_mass_fractions
 
      !------------------------------------------------------!
      ! This subroutine computes the species molar fractions
      SUBROUTINE Nitrogen_NASA_get_molar_fractions (rhoi, xi)

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

      END SUBROUTINE Nitrogen_NASA_get_molar_fractions

      !------------------------------------------------------!
      ! This subroutine computes the molar fractions from mass fractions. 
      ! The mixture molar mass is also computed and given in output. 
      SUBROUTINE Nitrogen_NASA_molar_frac_from_mass_frac (yi, mm, xi)

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

      END SUBROUTINE Nitrogen_NASA_molar_frac_from_mass_frac 
     
      !------------------------------------------------------!
      ! This subroutine computes the mass fractions from molar fractions. 
      ! The mixture molar mass is also computed and given in output.
      SUBROUTINE Nitrogen_NASA_mass_frac_from_molar_frac (xi, mm, yi)

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

      END SUBROUTINE Nitrogen_NASA_mass_frac_from_molar_frac 
 
      !------------------------------------------------------!
      ! This subroutine adds a small number to the composition respecting the mass constraint
      ! (sum of species molar fractions equal to one)
      SUBROUTINE Nitrogen_NASA_comp_tol(xi)

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

      END SUBROUTINE Nitrogen_NASA_comp_tol 

      !------------------------------------------------------!
      ! This subroutine computes the source term due to collisional and radiative processes.
      SUBROUTINE Nitrogen_NASA_get_source (rhoi, temp, s) 
  
        USE mod_function_pointer_NASA,     ONLY:  get_mass_prod_terms_NASA
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s
 
        CALL get_mass_prod_terms_NASA (rhoi, temp, s)        
  
      END SUBROUTINE Nitrogen_NASA_get_source 
  
      !------------------------------------------------------!
      ! This subroutine computes the source term due to collisional and radiative processes 
      ! and its Jacobian (with respect to primitive variables).
      SUBROUTINE Nitrogen_NASA_get_source_Jac (rhoi, temp, s, js) 
  
        USE mod_function_pointer_NASA,     ONLY:  get_mass_prod_terms_Jac_NASA
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s, js
       
        ! Source term and source term Jacobian (with respect to primitive variables)
        CALL get_mass_prod_terms_Jac_NASA (rhoi, temp, s, js)    

      END SUBROUTINE Nitrogen_NASA_get_source_Jac

      !------------------------------------------------------!
      ! This subroutine computes the species mass diffusion flux
      SUBROUTINE Nitrogen_NASA_compute_species_DiffFlux (p, T, Te, xi, diff_driv, Ji) 
 
        USE mod_nitrogen_NASA_CFD_transport,       ONLY: get_species_DiffFlux 
        
        REAL(KIND=8) :: nb

        REAL(KIND=8), INTENT(IN) :: p, T, Te
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, diff_driv        
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: Ji

        ! Compute mixture number density
        CALL Nitrogen_NASA_get_mix_numb_dens(p, T, Te, xi, nb) 
               
        ! Compute the species diffusion velocities
        CALL get_species_DiffFlux (T, nb, xi, diff_driv, Ji)

      END SUBROUTINE Nitrogen_NASA_compute_species_DiffFlux
  
      !------------------------------------------------------!
      ! This subrouotine computes the transport coefficients 
      SUBROUTINE Nitrogen_NASA_compute_transpCoeff (p, xi, temp, mu, kappa, lambda, Di, chi) 

        USE mod_nitrogen_NASA_CFD_transport,       ONLY: get_transport_coeff
 
        REAL(KIND=8) :: nb, mm
        REAL(KIND=8) :: T, lambda_tr, lambda_int 
        REAL(KIND=8), DIMENSION(nb_ns) :: xi_tol, yi_tol                       
 
        REAL(KIND=8), INTENT(IN) :: p
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, temp
        REAL(KIND=8), INTENT(OUT) :: mu, kappa
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: chi, Di, lambda

        ! Compute species molar fractions
        xi_tol = xi
        CALL Nitrogen_NASA_comp_tol(xi_tol)

        ! Compute mixture number density
        T = temp(1)
        CALL Nitrogen_NASA_get_mix_numb_dens(p, T, T, xi, nb) 

        ! Compute mixture transport coefficients
        CALL Nitrogen_NASA_mass_frac_from_molar_frac(xi_tol, mm, yi_tol)
        CALL get_transport_coeff (nb, xi_tol, yi_tol, temp, mu, kappa, lambda_tr, lambda_int, Di, chi)

        lambda(1) = lambda_tr + lambda_int

      END SUBROUTINE Nitrogen_NASA_compute_transpCoeff  
#endif 
  
  END MODULE mod_Nitrogen_NASA_library
!------------------------------------------------------------------------------!  
