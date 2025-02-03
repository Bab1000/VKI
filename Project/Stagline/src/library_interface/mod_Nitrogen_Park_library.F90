!------------------------------------------------------------------------------!
! This module implements an interface for the Nitrogen_Park thermodynamic library. It provides an
! implementation of functions and subroutines used to compute thermodynamic and transport properties and 
! source terms related to collisional and radiative processes.
! The interface is written such to have no data and variable dependency from the flow solver being used. 
! Therefore it can plugged to any code without big efforts. 
  MODULE mod_Nitrogen_Park_library 

#include"../config.h"

#ifdef NITROGEN_PARK
    IMPLICIT NONE
  
    INTEGER, SAVE :: nb_ns, nb_tvib, nb_te, nb_trot, nb_temp, nb_int_temp, nb_eq, nb_dim, posTr, posTv 
    REAL(KIND=8), PARAMETER :: kb = 1.380658d-23
    REAL(KIND=8), PARAMETER :: na = 6.0221367d23   
    REAL(KIND=8), PARAMETER :: ru = kb*na
    REAL(KIND=8), SAVE :: xi_tol
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: mass, Ri

    ! Subroutines of Nitrogen_Park thermodynamic library interface. In some cases they call 
    ! other subroutines for the model under use (even these are not dependent on flow solver data). 
    CONTAINS
  
      !------------------------------------------------------!
      ! This subroutine initializes the thermodynamic library in use.  
      SUBROUTINE Nitrogen_Park_initialize (input_ns, input_ntrot, input_ntvib, input_nte, input_ntemp, input_neq,  & 
                                         & input_ndim, input_solver, input_mixture, input_reaction, state_model,   &
                                         & input_transf, input_path, input_xi_tol) 
  
        USE mod_nitrogen_Park_initialize_CFD,            ONLY: initialize
        USE mod_nitrogen_Park_zeta,                      ONLY: set_pointer 
 
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
           WRITE(*,5)input_solver(1:length),'::Nitrogen_Park_library interface'   
           PRINT*
           WRITE(*,5)input_solver(1:length),'::Error in number of temperatures'     
           PRINT*
           STOP
        ENDIF

        neq = input_ns + 1 + input_ntrot + input_ntvib + input_nte + input_ndim
        IF (neq.NE.nb_eq) THEN
           PRINT*
           WRITE(*,5)input_solver(1:length),'::Nitrogen_Park_library interface'  
           PRINT* 
           WRITE(*,5)input_solver(1:length),'::Error in number of equations'     
           PRINT*
           STOP
        ENDIF
 
        ! Allocation of molecular mass and specific gas constant vectors
        ALLOCATE(mass(nb_ns), Ri(nb_ns))
        
        ! Nitrogen_Park library is initialized
        CALL initialize (input_solver, input_mixture, input_transf, input_reaction, input_path, & 
                       & nb_ns, nb_trot, nb_tvib, nb_te, nb_eq, nb_dim, mass)
  
        ! Tolerance on species molar fractions
        xi_tol = input_xi_tol

        ! Speficic gas constant of single species
        DO is = 1,nb_ns 
           Ri(is) = 1.d0/mass(is)
        ENDDO 
        Ri = ru*Ri

        ! Pointer initialization
        CALL set_pointer ()
  
        ! Useful indices 
        posTr = 1 + nb_trot
        posTv = posTr + nb_tvib

5     FORMAT(A,A)

      END SUBROUTINE Nitrogen_Park_initialize

      !------------------------------------------------------!
      ! This subroutine finalizes the Nitrogen Park thermodynamic library
      SUBROUTINE Nitrogen_Park_finalize ()

         USE mod_nitrogen_Park_initialize_CFD,            ONLY: finalize

         IF (ALLOCATED(mass))       DEALLOCATE(mass)
         IF (ALLOCATED(Ri))         DEALLOCATE(Ri)

         CALL finalize ()

      END SUBROUTINE Nitrogen_Park_finalize

      !------------------------------------------------------!
      ! This subroutine computes post-shock nonequilibrium conditions
      SUBROUTINE Nitrogen_Park_post_shock_neq (p1, u1, T1, p2, u2, T2, yi, xN)

        USE mod_nitrogen_Park_CFD_prop,                  ONLY: post_shock_neq

        REAL(KIND=8), INTENT(IN) :: p1, u1, T1
        REAL(KIND=8), INTENT(IN), OPTIONAL :: xN
        REAL(KIND=8), INTENT(OUT) :: p2, u2, T2
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi 

        CALL post_shock_neq (p1, u1, T1, p2, u2, T2, yi, xN)

      END SUBROUTINE Nitrogen_Park_post_shock_neq

      !------------------------------------------------------!
      ! This subroutine computes post-shock equilibrium conditions
      SUBROUTINE Nitrogen_Park_post_shock_eq (p1, u1, T1, p2, u2, T2, yi, rhoi, xN)

        USE mod_nitrogen_Park_CFD_eq,                  ONLY: post_shock_eq

        REAL(KIND=8), INTENT(IN) :: p1, u1, T1
        REAL(KIND=8), INTENT(IN), OPTIONAL :: xN
        REAL(KIND=8), INTENT(OUT) :: p2, u2, T2
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi, rhoi

        CALL post_shock_eq (p1, u1, T1, p2, u2, T2, yi, rhoi, xN)

      END SUBROUTINE Nitrogen_Park_post_shock_eq

      !------------------------------------------------------!
      ! This subroutine computes common factors for the ODE system in 
      ! case of thermla nonequilibrium
      SUBROUTINE Nitrogen_Park_compute_ODE_factors_tcneq (temp, omega, yi, fac_num, fac_den)

        USE mod_nitrogen_Park_initialize_CFD,      ONLY: pos_N2
        USE mod_nitrogen_Park_CFD_prop,            ONLY: ev_cvv_ho

        REAL(KIND=8) :: tmp1, tmp2
        REAL(KIND=8) :: Tr, Tv
        REAL(KIND=8) :: r, dexp0, den
        REAL(KIND=8) :: cv_vib, cv_rot, e_vib, e_rot

        REAL(KIND=8), INTENT(IN), DIMENSION(:) :: temp, omega, yi
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: fac_num, fac_den

        ! Temperatures 
        Tr = temp(posTr)
        Tv = temp(posTv)

        tmp1 = omega(pos_N2)
        tmp2 = yi(pos_N2)

        ! Rotation
        cv_rot = Ri(pos_N2) 
        e_rot  = cv_rot*Tr 

        ! Vibration
        CALL ev_cvv_ho (Tv, e_vib, cv_vib)

        ! Factors for rotational energy equation
        fac_num(MAX(1,posTr - 1)) = tmp1*e_rot
        fac_den(MAX(1,posTr - 1)) = tmp2*cv_rot

        ! Factors for vibrational energy equations
        fac_num(posTv - 1) = tmp1*e_vib
        fac_den(posTv - 1) = tmp2*cv_vib

      END SUBROUTINE Nitrogen_Park_compute_ODE_factors_tcneq

      !------------------------------------------------------!
      ! This subroutine write down the solution when solving the ODE system for
      ! studying the nonequilibrium flow behind a shock wave
      SUBROUTINE Nitrogen_Park_write_ODE_solution (flow_file, pop_file, k, xold, x, uold, m_dot, z)

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

      END SUBROUTINE Nitrogen_Park_write_ODE_solution 

      !------------------------------------------------------!
      ! This subroutine writes down the solution when solving the flow equations 
      ! by means of Finite Volume method
      SUBROUTINE Nitrogen_Park_write_fvmcc_solution (flow_file, pop_file, k, xold, x, uold, u, rhoi, temp)

        INTEGER :: is
        REAL(KIND=8) :: mm, fac
        REAL(KIND=8) :: p, rho
        REAL(KIND=8), SAVE :: time
        REAL(KIND=8), DIMENSION(2) :: xi
        REAL(KIND=8), DIMENSION(nb_ns) :: yi

        INTEGER, INTENT(IN) :: flow_file, pop_file, k
        REAL(KIND=8), INTENT(IN) :: xold, x, uold, u
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp

        ! Density 
        CALL Nitrogen_Park_get_density (rhoi, rho)

        ! Species mass fractions
        CALL Nitrogen_Park_get_mass_fractions (rho, rhoi, yi)
 
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

      END SUBROUTINE Nitrogen_Park_write_fvmcc_solution   

      !------------------------------------------------------!
      ! This subroutine computes the equilibrium composition given pressure and temperature
      SUBROUTINE Nitrogen_Park_compute_eq_composition (p, T, rhoi)

        USE mod_nitrogen_Park_CFD_eq,                 ONLY: eq_composition 

        INTEGER :: is
        REAL(KIND=8) :: xn, xn2
        REAL(KIND=8) :: mm, R, rho, ov_mm
        REAL(KIND=8), DIMENSION(nb_ns) :: yi, xi

        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: rhoi

        CALL eq_composition (p, T, xi)

        ! Mixture molecular mass
        mm = 0.d0
        DO is = 1,nb_ns 
           mm = mm + xi(is)*mass(is)
        ENDDO
        ov_mm = 1.d0/mm

        R = ru*ov_mm

        DO is = 1,nb_ns 
           yi(is) = xi(is)*mass(is)*ov_mm
        ENDDO

        ! Mixture density
        rho = p/(R*T)

        ! Species densities
        DO is = 1,nb_ns 
           rhoi(is) = rho*yi(is)
        ENDDO

      END SUBROUTINE Nitrogen_Park_compute_eq_composition

      !------------------------------------------------------!
      ! This subroutine computes the species equilibrium internal energy per unit mass 
      SUBROUTINE Nitrogen_Park_compute_eq_species_energy (T, e)

        USE mod_nitrogen_Park_CFD_eq,           ONLY: energy

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: e

         CALL energy (T, e)

      END SUBROUTINE Nitrogen_Park_compute_eq_species_energy 

      !------------------------------------------------------!
      ! This subroutine computes the species equilibrium enthalpy per unit mass 
      SUBROUTINE Nitrogen_Park_compute_eq_species_enthalpy (T, h)

        USE mod_nitrogen_Park_CFD_eq,           ONLY: enthalpy

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: h

        CALL enthalpy (T, h)

      END SUBROUTINE Nitrogen_Park_compute_eq_species_enthalpy 

      !------------------------------------------------------!
      ! This subroutine computes the species equilibrium entropy per unit mass. The
      ! entropy of mixing contributions is computed separately. 
      SUBROUTINE Nitrogen_Park_compute_eq_species_entropy (p, T, yi, xi, smix, s)

        USE mod_nitrogen_Park_CFD_eq,           ONLY: entropy

        REAL(KIND=8) :: tmp1
        REAL(KIND=8) :: QtN, QtN2

        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), INTENT(IN), DIMENSION(:) :: xi, yi
        REAL(KIND=8), INTENT(OUT) :: smix
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: s

        CALL entropy (p, T, yi, xi, smix, s)

      END SUBROUTINE Nitrogen_Park_compute_eq_species_entropy    

      !------------------------------------------------------!
      ! This subroutine computes the temperatures and thermodynamic data,
      ! additional data for Jacobians are given in output
      SUBROUTINE Nitrogen_Park_get_data (rho, rhoi, rho_eint, temp, c, gamma, p, alpha, beta, ei)
  
        USE mod_nitrogen_Park_initialize_CFD,  ONLY: pos_N2
        USE mod_function_pointer_Park,         ONLY: get_temp_Park, get_species_cv_Park
        USE mod_nitrogen_Park_CFD_prop,        ONLY: get_species_energy, ev_cvv_ho 
 
        IMPLICIT NONE
  
        INTEGER :: is
        INTEGER :: offset, pos 
        REAL(KIND=8) :: tmp
        REAL(KIND=8) :: sum1, sum2
        REAL(KIND=8) :: cvr, cvv, er, ev
        REAL(KIND=8) :: T, Tr, Tv, ki      
        REAL(KIND=8), DIMENSION(nb_ns) :: cv  
  
        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, rho_eint
        REAL(KIND=8), INTENT(OUT) :: c, gamma, p, alpha
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp, beta, ei
  
        ! Temperature computation
        CALL get_temp_Park (rhoi, rho_eint, temp)
        T = temp(1)
        CALL get_species_cv_Park (T, cv)
 
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
       
        ! Species energy 
        CALL get_species_energy (temp, ei(1:nb_ns))
        
        ! Additional data for thermal non equilibrium cases
        IF (nb_int_temp.GE.1) THEN

           ! Rotational and vibrational temperatures
           Tr = temp(posTr) 
           Tv = temp(posTv)  

           CALL ev_cvv_ho(Tv, ev, cvv)

           cvr = Ri(pos_N2)
           er  = cvr*Tr

           pos = MAX(2,posTr)
           IF (nb_trot.EQ.1) beta(pos)   = rhoi(pos_N2)*cvr
           IF (nb_tvib.EQ.1) beta(posTv) = rhoi(pos_N2)*cvv

           ! Nonequilibrium energy submatrix
           ! Nitrogen atom N 
           DO is = 1,nb_int_temp 
              ei(nb_ns + is) = 0.d0 
           ENDDO
 
           ! Nitrogen molecule N2
           offset = (pos_N2 - 1)*nb_int_temp
           ei(nb_ns + offset + pos   - 1) = er
           ei(nb_ns + offset + posTv - 1) = ev

        ENDIF
       
      END SUBROUTINE Nitrogen_Park_get_data 
  
      !------------------------------------------------------!
      ! This subroutine computes thermodynamic data needed to be stored
      ! for nonequilibrium gas flow computations. 
      SUBROUTINE Nitrogen_Park_get_thermodynamic_data (rho, rhoi, temp, c, gamma, p, alpha, beta, ei)
 
        USE mod_nitrogen_Park_initialize_CFD,  ONLY: pos_N2 
        USE mod_function_pointer_Park,         ONLY: get_species_cv_Park
        USE mod_nitrogen_Park_CFD_prop,        ONLY: get_species_energy, ev_cvv_ho 

        INTEGER :: is 
        INTEGER :: pos, offset
        REAL(KIND=8) :: tmp
        REAL(KIND=8) :: sum1, sum2
        REAL(KIND=8) :: cvr, cvv, er, ev
        REAL(KIND=8) :: T, Tr, Tv,  ki 
        REAL(KIND=8), DIMENSION(nb_ns) :: cv  
  
        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), INTENT(OUT) :: c, gamma, p, alpha
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: beta, ei

        ! Temperature   
        T = temp(1)

        ! Frozen specific heats
        CALL get_species_cv_Park (T, cv)  

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
       
        ! Species energy 
        CALL get_species_energy (temp, ei(1:nb_ns))

        ! Additional data for thermal non equilibrium cases
        IF (nb_int_temp.GE.1) THEN

           ! Rotational and vibrational temperatures
           Tr = temp(posTr) 
           Tv = temp(posTv)  

           CALL ev_cvv_ho(Tv, ev, cvv)

           cvr = Ri(pos_N2)
           er  = cvr*Tr

           pos = MAX(2,posTr)
           IF (nb_trot.EQ.1) beta(pos)   = rhoi(pos_N2)*cvr
           IF (nb_tvib.EQ.1) beta(posTv) = rhoi(pos_N2)*cvv

           ! Nonequilibrium energy submatrix
           ! Nitrogen atom N 
           DO is = 1,nb_int_temp 
              ei(nb_ns + is) = 0.d0 
           ENDDO
 
           ! Nitrogen molecule N2 
           offset = (pos_N2 - 1)*nb_int_temp
           ei(nb_ns + offset + pos   - 1) = er
           ei(nb_ns + offset + posTv - 1) = ev

        ENDIF
 
      END SUBROUTINE Nitrogen_Park_get_thermodynamic_data
 
      !------------------------------------------------------!
      ! This subroutine computes the energy density for total internal,
      ! rotational, vibrational and electronic components.
      SUBROUTINE Nitrogen_Park_get_energy_densities (rhoi, temp, rho_e)
  
        USE mod_nitrogen_Park_initialize_CFD,    ONLY: pos_N2
        USE mod_nitrogen_Park_CFD_prop,          ONLY: energy  
  
        INTEGER :: is
        REAL(KIND=8) :: er, ev, tmp
        REAL(KIND=8), DIMENSION(nb_ns) :: eint
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: rho_e
        
        CALL energy (temp, er, ev, eint)
        
        tmp = rhoi(pos_N2)
        rho_e(posTr) = tmp*er
        rho_e(posTv) = tmp*ev

        ! Internal energy per unit volume of the mixture
        tmp = 0.d0
        DO is = 1,nb_ns 
           tmp = tmp + rhoi(is)*eint(is)
        ENDDO
        rho_e(1) = tmp 

      END SUBROUTINE Nitrogen_Park_get_energy_densities
  
      !------------------------------------------------------!
      ! This subroutine computes the total internal, rotational, vibrational and electronic energy 
      ! per unit mass of the mixture.
      SUBROUTINE Nitrogen_Park_get_energy (rho, rhoi, temp, e)
  
        USE mod_nitrogen_Park_initialize_CFD,    ONLY: pos_N2
        USE mod_nitrogen_Park_CFD_prop,          ONLY: energy  
  
        INTEGER :: is
        REAL(KIND=8) :: er, ev, tmp
        REAL(KIND=8), DIMENSION(nb_ns) :: eint, yi
  
        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e
       
        ! Species mass fractions
        CALL Nitrogen_Park_get_mass_fractions (rho, rhoi, yi)
  
        CALL energy (temp, er, ev, eint)
  
        tmp = yi(pos_N2)
        e(posTr) = tmp*er
        e(posTv) = tmp*ev
  
        ! Static enthalpy per unit mass of the mixture
        tmp = 0.d0
        DO is = 1,nb_ns 
           tmp = tmp + yi(is)*eint(is)
        ENDDO
        e(1) = tmp  

      END SUBROUTINE Nitrogen_Park_get_energy
  
      !------------------------------------------------------!
      ! This subroutine computes the total internal, rotational, vibrational and electronic enthalpy 
      ! per unit mass of the mixture.
      SUBROUTINE Nitrogen_Park_get_enthalpy (rho, rhoi, temp, h)
 
        USE mod_nitrogen_Park_initialize_CFD,    ONLY: pos_N2 
        USE mod_nitrogen_Park_CFD_prop,          ONLY: enthalpy  
  
        INTEGER :: is
        REAL(KIND=8) :: er, ev, tmp
        REAL(KIND=8), DIMENSION(nb_ns) :: hint, yi
  
        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: h
       
        ! Species mass fractions
        CALL Nitrogen_Park_get_mass_fractions (rho, rhoi, yi)
  
        CALL enthalpy (temp, er, ev, hint)
  
        tmp = yi(pos_N2)
        h(posTr) = tmp*er 
        h(posTv) = tmp*ev
  
        ! Static enthalpy per unit mass of the mixture
        tmp = 0.d0
        DO is = 1,nb_ns 
           tmp = tmp + yi(is)*hint(is)
        ENDDO
        h(1) = tmp
  
      END SUBROUTINE Nitrogen_Park_get_enthalpy 
 
      !----------------------------------------------------!
      ! This subroutine computes the species internal energies and constant 
      ! volume specific heat per unit mass
      SUBROUTINE Nitrogen_Park_get_species_energy (temp, e) 
 
        USE mod_Nitrogen_Park_CFD_prop,    ONLY: get_species_energy

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e

        CALL get_species_energy (temp, e)

      END SUBROUTINE Nitrogen_Park_get_species_energy

      !----------------------------------------------------!
      ! This subroutine computes the species internal energies heat per unit mass
      SUBROUTINE Nitrogen_Park_get_species_energy_cv (temp, e, cv) 
 
        USE mod_Nitrogen_Park_CFD_prop,    ONLY: get_species_energy_cv

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e, cv

        CALL get_species_energy_cv (temp, e, cv)

      END SUBROUTINE Nitrogen_Park_get_species_energy_cv

      !----------------------------------------------------!
      ! This subroutine computes the species internal energies and constant 
      ! volume specific heat per unit mass
      SUBROUTINE Nitrogen_Park_get_species_enthalpy (temp, h) 
 
        USE mod_Nitrogen_Park_CFD_prop,    ONLY: get_species_energy

        INTEGER :: is
        REAL(KIND=8) :: T 

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: h

        CALL get_species_energy (temp, h)
 
        T = temp(1)
        DO is = 1,nb_ns 
           h(is) = h(is) + Ri(is)*T
        ENDDO        

      END SUBROUTINE Nitrogen_Park_get_species_enthalpy

      !----------------------------------------------------!
      ! This subroutine computes the species enthalpies and constant 
      ! pressure specific heat per unit mass
      SUBROUTINE Nitrogen_Park_get_species_enthalpy_cp (temp, h, cp) 
 
        USE mod_Nitrogen_Park_CFD_prop,    ONLY: get_species_energy_cv

        INTEGER :: is
        REAL(KIND=8) :: T, tmp1

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: h, cp

        CALL get_species_energy_cv (temp, h, cp)

        T = temp(1)
        DO is = 1,nb_ns 
           tmp1 = Ri(is)
           h(is)  = h(is) + tmp1*T
           cp(is) = cp(is) + tmp1
        ENDDO

      END SUBROUTINE Nitrogen_Park_get_species_enthalpy_cp

      !------------------------------------------------------!
      ! This subroutine computes the frozen specific heats at constant volume of each species
      SUBROUTINE Nitrogen_Park_get_species_frozen_cv (T, cv)  
 
        USE mod_function_pointer_Park,      ONLY: get_species_cv_Park   
 
        REAL(KIND=8), INTENT(IN) :: T 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv
  
        CALL get_species_cv_Park (T, cv)

      END SUBROUTINE Nitrogen_Park_get_species_frozen_cv 
  
      !------------------------------------------------------!
      ! This subroutine computes the frozen specific heats at constant pressure of each species
      SUBROUTINE Nitrogen_Park_get_species_frozen_cp (T, cp)  
  
        USE mod_function_pointer_Park,      ONLY: get_species_cv_Park   

        INTEGER :: i

        REAL(KIND=8), INTENT(IN) :: T 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cp
  
        CALL get_species_cv_Park (T, cp)

        DO i = 1,nb_ns 
           cp(i) = cp(i) + Ri(i)
        ENDDO 

      END SUBROUTINE Nitrogen_Park_get_species_frozen_cp

      !------------------------------------------------------!
      ! This subroutine computes the frozen specific heats at constant volume of each species and 
      ! the species energy components in thermal equilibrium with translation.
      SUBROUTINE Nitrogen_Park_get_species_frozen_energy_cv (temp, e, cv)  
 
        USE mod_function_pointer_Park,      ONLY: get_species_energy_cv_Park   
 
        REAL(KIND=8) :: T

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv, e
  
        T = temp(1)
        CALL get_species_energy_cv_Park (T, cv, e)

      END SUBROUTINE Nitrogen_Park_get_species_frozen_energy_cv 
  
      !------------------------------------------------------!
      ! This subroutine computes the frozen specific heats at constant pressure of each species and 
      ! the species enthalpy components in thermal equilibrium with translation.
      SUBROUTINE Nitrogen_Park_get_species_frozen_enthalpy_cp (temp, h, cp)  
  
        USE mod_function_pointer_Park,      ONLY: get_species_energy_cv_Park   

        INTEGER :: i

        REAL(KIND=8) :: T
        REAL(KIND=8) :: tmp

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cp, h
  
        T = temp(1)
        CALL get_species_energy_cv_Park (T, cp, h)

        DO i = 1,nb_ns 
           tmp   = Ri(i)
           h(i)  = h(i) + tmp*T
           cp(i) = cp(i) + tmp
        ENDDO 

      END SUBROUTINE Nitrogen_Park_get_species_frozen_enthalpy_cp

      !------------------------------------------------------!
      ! This subroutine computes the temperatures given the energy densities.
      SUBROUTINE Nitrogen_Park_get_temperatures (rhoi, rho_eint, temp)
  
        USE mod_function_pointer_Park,    ONLY: get_temp_Park
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rho_eint
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp
  
        CALL get_temp_Park (rhoi, rho_eint, temp) 
        
      END SUBROUTINE Nitrogen_Park_get_temperatures
 
      !------------------------------------------------------!
      ! This subroutine computes the equilibrium speed of sound
      SUBROUTINE Nitrogen_Park_get_eq_sound_speed (rhoi, p, T, c)
 
        USE mod_Nitrogen_Park_CFD_eq,    ONLY: energy

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
        CALL Nitrogen_Park_get_density (rhoi, rho)
        CALL Nitrogen_Park_get_mass_fractions (rho, rhoi, yi)

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
        
        CALL Nitrogen_Park_compute_eq_composition (p_pert, T, dens)
        CALL Nitrogen_Park_get_density (dens, rho_pert)

        drho_dp = (rho_pert - rho)/(eps*p)

        ! Equilibrium specific heat ratio
        T_pert = T*(1.d0 + eps)
        
        CALL Nitrogen_Park_compute_eq_composition (p, T_pert, dens)
        CALL Nitrogen_Park_get_density (dens, rho_pert)
        CALL Nitrogen_Park_get_mass_fractions (rho_pert, dens, yi)

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

      END SUBROUTINE Nitrogen_Park_get_eq_sound_speed 

      !------------------------------------------------------!
      ! This subroutine computes the equilibrium specific heat ratio
      SUBROUTINE Nitrogen_Park_get_eq_gamma (rhoi, p, T, gamma)
 
        USE mod_Nitrogen_Park_CFD_eq,    ONLY: energy

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
        CALL Nitrogen_Park_get_density (rhoi, rho)
        CALL Nitrogen_Park_get_mass_fractions (rho, rhoi, yi)

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
        
        CALL Nitrogen_Park_compute_eq_composition (p, T_pert, dens)
        CALL Nitrogen_Park_get_density (dens, rho_pert)
        CALL Nitrogen_Park_get_mass_fractions (rho_pert, dens, yi)

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

      END SUBROUTINE Nitrogen_Park_get_eq_gamma 

      !------------------------------------------------------!
      ! This subroutine computes the equilibrium specific heat ratio
      SUBROUTINE Nitrogen_Park_get_eq_gamma_sound_speed (rhoi, p, T, gamma, c)
 
         USE mod_Nitrogen_Park_CFD_eq,    ONLY: energy

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
        CALL Nitrogen_Park_get_density (rhoi, rho)
        CALL Nitrogen_Park_get_mass_fractions (rho, rhoi, yi)

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
        
        CALL Nitrogen_Park_compute_eq_composition (p_pert, T, dens)
        CALL Nitrogen_Park_get_density (dens, rho_pert)

        drho_dp = (rho_pert - rho)/(eps*p)

        ! Equilibrium specific heat ratio
        T_pert = T*(1.d0 + eps)
        
        CALL Nitrogen_Park_compute_eq_composition (p, T_pert, dens)
        CALL Nitrogen_Park_get_density (dens, rho_pert)
        CALL Nitrogen_Park_get_mass_fractions (rho_pert, dens, yi)

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

      END SUBROUTINE Nitrogen_Park_get_eq_gamma_sound_speed 

      !------------------------------------------------------!
      ! This subroutine computes the frozen speed of sound
      SUBROUTINE Nitrogen_Park_get_frozen_sound_speed (rhoi, temp, c)
 
        USE mod_function_pointer_Park,      ONLY: get_species_cv_Park 
 
        INTEGER :: is
        REAL(KIND=8) :: ki, sum1, sum2, R, gamma, rho, tmp, T
        REAL(KIND=8), DIMENSION(nb_ns) :: cv  
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp, rhoi
        REAL(KIND=8), INTENT(OUT) :: c 

        ! Temperature   
        T = temp(1)

        ! Frozen specific heats
        CALL Nitrogen_Park_get_density (rhoi, rho)
        CALL get_species_cv_Park (T, cv)  

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
  
      END SUBROUTINE Nitrogen_Park_get_frozen_sound_speed
 
      !------------------------------------------------------!
      ! This subroutine gets the frozen speed specific heat ratio.
      SUBROUTINE Nitrogen_Park_get_frozen_gamma (rhoi, temp, gamma)
 
        USE mod_function_pointer_Park,          ONLY: get_species_cv_Park
 
        INTEGER :: is
        REAL(KIND=8) :: ki, sum1, sum2, tmp, T
        REAL(KIND=8), DIMENSION(nb_ns) :: cv  
 
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), INTENT(OUT) :: gamma 
  
        ! Temperature 
        T = temp(1)

        ! Frozen specific heats 
        CALL get_species_cv_Park (T, cv)
   
        ! Inizialization
        sum1 = 0.d0   
        sum2 = 0.d0  
        DO is = 1,nb_ns 
           tmp  = rhoi(is)
           sum1 = sum1 + tmp*Ri(is)
           sum2 = sum2 + tmp*cv(is)
        ENDDO
  
        ki    = sum1/sum2
        gamma = 1.d0 + ki
  
      END SUBROUTINE Nitrogen_Park_get_frozen_gamma
 
      !------------------------------------------------------!
      ! This subroutine gets the frozen speed of sound and specific heat ratio.
      SUBROUTINE Nitrogen_Park_get_frozen_gamma_sound_speed (rhoi, temp, gamma, c)
  
        USE mod_function_pointer_Park,          ONLY: get_species_cv_Park

        INTEGER :: is
        REAL(KIND=8) :: ki, sum1, sum2, R, rho, tmp, T
        REAL(KIND=8), DIMENSION(nb_ns) :: cv  
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), INTENT(OUT) :: gamma, c 
  
        ! Temperature 
        T = temp(1)

        ! Frozen specific heats (constant volume)
        CALL get_species_cv_Park (T, cv)
        CALL Nitrogen_Park_get_density (rhoi, rho)
   
        ! Inizialization
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
  
      END SUBROUTINE Nitrogen_Park_get_frozen_gamma_sound_speed
  
      !------------------------------------------------------!
      ! This subroutine computes the mixture static pressure (presence of free electrons is accounted for)
      SUBROUTINE Nitrogen_Park_get_pressure (rhoi, temp, p)
   
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
  
      END SUBROUTINE Nitrogen_Park_get_pressure
  
      !------------------------------------------------------!
      ! This subroutine computes the mixture density.
      SUBROUTINE Nitrogen_Park_get_density (rhoi, rho)
  
        INTEGER :: is
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), INTENT(OUT) :: rho
  
        ! Mixture density
        rho = 0.d0
        DO is = 1,nb_ns 
           rho = rho + rhoi(is)
        ENDDO
  
      END SUBROUTINE Nitrogen_Park_get_density
  
      !------------------------------------------------------!
      ! This subroutine computes the specific gas constant R.
      SUBROUTINE Nitrogen_Park_get_Rgas (rhoi, R)
  
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
  
      END SUBROUTINE Nitrogen_Park_get_Rgas
  
      !------------------------------------------------------!
      ! This subroutine gives the species gas constants
      SUBROUTINE Nitrogen_Park_get_Ri (out_Ri)
  
        INTEGER :: is
  
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: out_Ri
  
        DO is = 1,nb_ns 
           out_Ri(is) = Ri(is) 
        ENDDO
  
      END SUBROUTINE Nitrogen_Park_get_Ri 
 
      !------------------------------------------------------!
      ! This subroutine gives the species molar masses
      SUBROUTINE Nitrogen_Park_get_mi (out_mi)
  
        INTEGER :: is
  
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: out_mi
  
        DO is = 1,nb_ns 
           out_mi(is) = mass(is) 
        ENDDO
  
      END SUBROUTINE Nitrogen_Park_get_mi

      !------------------------------------------------------!
      ! This subroutine computes the mixture number density 
      SUBROUTINE Nitrogen_Park_get_mix_numb_dens (p, T, Te, xi, nb) 

        REAL(KIND=8), INTENT(IN) :: p, T, Te
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi
        REAL(KIND=8) :: nb

        nb = p/(kb*T)

      END SUBROUTINE Nitrogen_Park_get_mix_numb_dens 
 
      !------------------------------------------------------!
      ! This subroutine computes the species mass fraction
      SUBROUTINE Nitrogen_Park_get_mass_fractions (rho, rhoi, yi)
  
        INTEGER :: is
  
        REAL(KIND=8), INTENT(IN) :: rho
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi
  
        DO is = 1,nb_ns 
           yi(is) = rhoi(is)/rho
        ENDDO
  
      END SUBROUTINE Nitrogen_Park_get_mass_fractions
 
      !------------------------------------------------------!
      ! This subroutine computes the species molar fractions
      SUBROUTINE Nitrogen_Park_get_molar_fractions (rhoi, xi)

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

      END SUBROUTINE Nitrogen_Park_get_molar_fractions

      !------------------------------------------------------!
      ! This subroutine computes the molar fractions from mass fractions. 
      ! The mixture molar mass is also computed and given in output. 
      SUBROUTINE Nitrogen_Park_molar_frac_from_mass_frac (yi, mm, xi)

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

      END SUBROUTINE Nitrogen_Park_molar_frac_from_mass_frac 
     
      !------------------------------------------------------!
      ! This subroutine computes the mass fractions from molar fractions. 
      ! The mixture molar mass is also computed and given in output.
      SUBROUTINE Nitrogen_Park_mass_frac_from_molar_frac (xi, mm, yi)

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

      END SUBROUTINE Nitrogen_Park_mass_frac_from_molar_frac 
 
      !------------------------------------------------------!
      ! This subroutine adds a small number to the composition respecting the mass constraint
      ! (sum of species molar fractions equal to one)
      SUBROUTINE Nitrogen_Park_comp_tol(xi)

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

      END SUBROUTINE Nitrogen_Park_comp_tol

      !------------------------------------------------------!
      ! This subroutine computes the source term due to collisional and radiative processes.
      SUBROUTINE Nitrogen_Park_get_source (rhoi, temp, s) 
  
        USE mod_function_pointer_Park,         ONLY: get_source_term
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s
  
        CALL get_source_term (rhoi, temp, s)        
  
      END SUBROUTINE Nitrogen_Park_get_source 
  
      !------------------------------------------------------!
      ! This subroutine computes the source term due to collisional and radiative processes 
      ! and its Jacobian (with respect to primitive variables).
      SUBROUTINE Nitrogen_Park_get_source_Jac (rhoi, temp, s, js)    

        USE mod_function_pointer_Park,         ONLY: get_source_term_Jac
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s, js

        ! Source term and source term Jacobian (with respect to primitive variables)
        CALL get_source_term_Jac (rhoi, temp, s, js)

      END SUBROUTINE Nitrogen_Park_get_source_Jac
 
      !------------------------------------------------------!
      ! This subroutine computes the species mass diffusion flux 
      SUBROUTINE Nitrogen_Park_compute_species_DiffFlux (p, T, Te, xi, diff_driv, Ji) 
 
        USE mod_nitrogen_Park_CFD_transport,       ONLY: get_species_DiffFlux 
        
        REAL(KIND=8) :: nb

        REAL(KIND=8), INTENT(IN) :: p, T, Te
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, diff_driv        
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: Ji 

        ! Compute mixture number density
        CALL Nitrogen_Park_get_mix_numb_dens(p, T, Te, xi, nb) 
               
        ! Compute the species mass diffusion flux
        CALL get_species_DiffFlux (T, nb, xi, diff_driv, Ji)

      END SUBROUTINE Nitrogen_Park_compute_species_DiffFlux
  
      !------------------------------------------------------!
      ! This subrouotine computes the transport coefficients 
      SUBROUTINE Nitrogen_Park_compute_transpCoeff (p, xi, temp, mu, kappa, lambda, Di, chi) 
                         
        USE mod_nitrogen_Park_initialize_CFD,      ONLY: posTr, posTv
        USE mod_nitrogen_Park_CFD_transport,       ONLY: get_transport_coeff, get_species_DiffFlux
  
        REAL(KIND=8) :: nb, mm
        REAL(KIND=8) :: lambda_tr, lambda_rot, lambda_vib
        REAL(KIND=8), DIMENSION(nb_ns) :: xi_tol, yi_tol

        REAL(KIND=8), INTENT(IN) :: p
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, temp
        REAL(KIND=8), INTENT(OUT) :: mu, kappa
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: chi, Di, lambda

        ! Compute species molar fractions
        xi_tol = xi
        CALL Nitrogen_Park_comp_tol(xi_tol)
  
        ! Compute mixture number density
        CALL Nitrogen_Park_get_mix_numb_dens(p, temp(1), temp(1), xi, nb)  

        ! Compute the trasport coefficients
        CALL Nitrogen_Park_mass_frac_from_molar_frac(xi_tol, mm, yi_tol)
        CALL get_transport_coeff (nb, xi_tol, yi_tol, temp, mu, kappa, lambda_tr, lambda_rot, lambda_vib, Di, chi) 

        ! Initialization
        lambda = 0.d0

        ! Fill thermal conductivity vector
        lambda(1)     = lambda_tr       
        lambda(posTr) = lambda(posTr) + lambda_rot
        lambda(posTv) = lambda(posTv) + lambda_vib  
        
      END SUBROUTINE Nitrogen_Park_compute_transpCoeff   
#endif
  
  END MODULE mod_Nitrogen_Park_library
!------------------------------------------------------------------------------!  
