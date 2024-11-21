!------------------------------------------------------------------------------!
! This module provides implementation of subroutines for dealing with thermodynamics of the N-N2 system
! when using the DSMC multi-temperature model.
  MODULE mod_nitrogen_DSMC_CFD_prop

  USE mod_nitrogen_DSMC_initialize_CFD,         ONLY: posTr, posTv   

  IMPLICIT NONE

  ! Subroutine for thermodynamics
  CONTAINS

    !------------------------------------------------------!
    ! This subroutine computes the translation partition function.
    SUBROUTINE Q_trans (T, mm, Q)

      USE mod_nitrogen_DSMC_initialize_CFD,     ONLY: fac_Q

      REAL(KIND=8), INTENT(IN)  :: T, mm
      REAL(KIND=8), INTENT(OUT) :: Q

      Q = fac_Q*(mm*T)**1.5
     
    END SUBROUTINE Q_trans

    !------------------------------------------------------!
    ! This subroutine computes the translation partition function and its
    ! derivative with respect to temperature.
    SUBROUTINE Q_trans_der (T, mm, Q, dQ_dT)

      USE mod_nitrogen_DSMC_initialize_CFD,     ONLY: fac_Q

      REAL(KIND=8), INTENT(IN)  :: T, mm
      REAL(KIND=8), INTENT(OUT) :: Q, dQ_dT

      Q = fac_Q*(mm*T)**1.5
      dQ_dT = 1.5d0*Q/T

    END SUBROUTINE Q_trans_der

    !------------------------------------------------------!
    ! This subroutine computes the internal partition function for the N2 molecule
    ! by means of the rigid-rotator harmonic oscillator approximation.
    SUBROUTINE Q_int_rr_ho (T, Q)

      USE mod_nitrogen_DSMC_initialize_CFD,     ONLY: theta_vib, theta_rot

      REAL(KIND=8) :: dexp_1       

      REAL(KIND=8), INTENT(IN)  :: T 
      REAL(KIND=8), INTENT(OUT) :: Q

      dexp_1 = DEXP(-theta_vib/T)
      Q = (T/(theta_rot*2.d0))*(1.d0/(1.d0 - dexp_1))

    END SUBROUTINE Q_int_rr_ho

    !------------------------------------------------------!
    ! This subroutine computes the internal partition function for the N2 molecule and its derivative 
    ! with respect to temperature by means of the rigid-rotator harmonic oscillator approximation.
    SUBROUTINE Q_int_rr_ho_der (T, Q, dQint_dT)

      USE mod_nitrogen_DSMC_initialize_CFD,     ONLY: theta_vib, theta_rot

      REAL(KIND=8) :: fac_1, fac_2, dexp_1       

      REAL(KIND=8), INTENT(IN)  :: T 
      REAL(KIND=8), INTENT(OUT) :: Q, dQint_dT

      dexp_1 = DEXP(-theta_vib/T)
      Q = (T/theta_rot/2.d0)*(1.d0/(1.d0 - dexp_1))

      fac_1 = theta_vib/T
      fac_2 = 1.d0 - dexp_1
      dQint_dT = 0.5d0/theta_rot/fac_2*(dexp_1*fac_1/fac_2 + 1.d0)

    END SUBROUTINE Q_int_rr_ho_der

    !------------------------------------------------------!
    ! This subroutine computes the frozen specific heats at constant volume (thermal equilibrium case).
    SUBROUTINE compute_frozen_cv_T (T, cv)

      USE mod_nitrogen_DSMC_initialize_CFD,     ONLY: pos_N2, theta_vib, cv_tr, cv_rot, Rn2

      REAL(KIND=8) :: cv_vib, den, dexp0, r 

      REAL(KIND=8), INTENT(IN) :: T
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv

      ! Translational component
      cv = cv_tr
      
      r = theta_vib/T
      dexp0 = DEXP(r)
      den   = dexp0 - 1.d0

      cv_vib = Rn2*dexp0*(r/den)**2 

      ! Adding rotational and vibrational contribution for the N2 molecule 
      cv(pos_N2) = cv(pos_N2) + cv_rot + cv_vib

    END SUBROUTINE compute_frozen_cv_T 

    !------------------------------------------------------!
    ! This subroutine computes the frozen specific heats at constant volume (thermal nonequilibrium case).
    SUBROUTINE compute_frozen_cv_T_Tv (T, cv)

      USE mod_nitrogen_DSMC_initialize_CFD,     ONLY: pos_N2, theta_vib, cv_tr, cv_rot

      REAL(KIND=8), INTENT(IN) :: T
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv

      ! Translational component
      cv = cv_tr
      
      ! Adding rotational contribution for the N2 molecule 
      cv(pos_N2) = cv(pos_N2) + cv_rot 

    END SUBROUTINE compute_frozen_cv_T_Tv 

    !------------------------------------------------------!
    ! This subroutine computes the frozen specific heats at constant volume (thermal nonequilibrium case).
    SUBROUTINE compute_frozen_cv_T_Tr (T, cv)

      USE mod_nitrogen_DSMC_initialize_CFD,     ONLY: theta_vib, cv_tr, cv_rot

      REAL(KIND=8), INTENT(IN) :: T
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv

      ! Translational component
      cv = cv_tr
      
    END SUBROUTINE compute_frozen_cv_T_Tr

    !------------------------------------------------------!
    ! This subroutine computes the frozen specific heats at constant volume (thermal nonequilibrium case).
    SUBROUTINE compute_frozen_cv_T_Tr_Tv (T, cv)

      USE mod_nitrogen_DSMC_initialize_CFD,     ONLY: theta_vib, cv_tr, cv_rot

      REAL(KIND=8), INTENT(IN) :: T
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv

      ! Translational component
      cv = cv_tr
      
    END SUBROUTINE compute_frozen_cv_T_Tr_Tv 

    !------------------------------------------------------!
    ! This subroutine computes the frozen specific heats at constant volume (thermal equilibrium case)
    SUBROUTINE compute_frozen_energy_cv_T (T, cv, e)

      USE mod_nitrogen_DSMC_initialize_CFD,     ONLY: pos_N, pos_N2, theta_vib, cv_tr, cv_rot, Rn2, hf_N

      REAL(KIND=8) :: cv_vib, den, dexp0, r 

      REAL(KIND=8), INTENT(IN) :: T
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv, e

      ! Translational component
      cv = cv_tr
      
      r = theta_vib/T
      dexp0 = DEXP(r)
      den   = dexp0 - 1.d0

      cv_vib = Rn2*dexp0*(r/den)**2 

      ! Adding rotational and vibrational contribution for the N2 molecule 
      cv(pos_N2) = cv(pos_N2) + cv_rot + cv_vib

      ! Species energy components in thermal equilibrium with translation
      e(pos_N)  = cv(pos_N)*T + hf_N
      e(pos_N2) = (cv_tr(pos_N2) + cv_rot)*T + Rn2*theta_vib/den

    END SUBROUTINE compute_frozen_energy_cv_T 

    !------------------------------------------------------!
    ! This subroutine computes the frozen specific heats at constant volume (thermal nonequilibrium case)
    SUBROUTINE compute_frozen_energy_cv_T_Tv (T, cv, e)

      USE mod_nitrogen_DSMC_initialize_CFD,     ONLY: pos_N, pos_N2, theta_vib, cv_tr, cv_rot, hf_N

      REAL(KIND=8), INTENT(IN) :: T
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv, e

      ! Translational component
      cv = cv_tr
      
      ! Adding rotational contribution for the N2 molecule 
      cv(pos_N2) = cv(pos_N2) + cv_rot 

      ! Species energy components in thermal equilibrium with translation
      e(pos_N)  = cv(pos_N)*T + hf_N
      e(pos_N2) = cv(pos_N2)*T

    END SUBROUTINE compute_frozen_energy_cv_T_Tv 

    !------------------------------------------------------!
    ! This subroutine computes the frozen specific heats at constant volume (thermal nonequilibrium case)
    SUBROUTINE compute_frozen_energy_cv_T_Tr (T, cv, e)

      USE mod_nitrogen_DSMC_initialize_CFD,     ONLY: pos_N, pos_N2, theta_vib, cv_tr, cv_rot, hf_N

      REAL(KIND=8), INTENT(IN) :: T
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv, e

      ! Translational component
      cv = cv_tr
      
      ! Species energy components in thermal equilibrium with translation
      e(pos_N)  = cv(pos_N)*T + hf_N
      e(pos_N2) = cv(pos_N2)*T

    END SUBROUTINE compute_frozen_energy_cv_T_Tr

    !------------------------------------------------------!
    ! This subroutine computes the frozen specific heats at constant volume (thermal nonequilibrium case)
    SUBROUTINE compute_frozen_energy_cv_T_Tr_Tv (T, cv, e)

      USE mod_nitrogen_DSMC_initialize_CFD,     ONLY: pos_N, pos_N2, theta_vib, cv_tr, cv_rot, hf_N

      REAL(KIND=8), INTENT(IN) :: T
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv, e

      ! Translational component
      cv = cv_tr
      
      ! Species energy components in thermal equilibrium with translation
      e(pos_N)  = cv(pos_N)*T + hf_N
      e(pos_N2) = cv(pos_N2)*T

    END SUBROUTINE compute_frozen_energy_cv_T_Tr_Tv 

    !------------------------------------------------------!
    ! This subroutine computes the frozen specific heats at constant pressure.
    SUBROUTINE frozen_cp (T, cp)

      USE mod_nitrogen_DSMC_initialize_CFD,          ONLY: nb_ns, Ri
      USE mod_function_pointer_DSMC,                 ONLY: get_species_cv_DSMC 

      INTEGER :: is

      REAL(KIND=8), INTENT(IN) :: T
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cp 

      CALL get_species_cv_DSMC (T, cp)

      DO is = 1,nb_ns 
         cp(is) = cp(is) + Ri(is)
      ENDDO

    END SUBROUTINE frozen_cp

    !------------------------------------------------------!
    ! This subroutine computes the vibrational energy and specific heat for the N2 molecule
    ! when using the harmonic oscillator model.
    SUBROUTINE ev_cvv_ho (Tv, ev, cvv)

      USE mod_nitrogen_DSMC_initialize_CFD,     ONLY: Rn2, theta_vib

      REAL(KIND=8) :: dexp0, Tr, den

      REAL(KIND=8), INTENT(IN) :: Tv
      REAL(KIND=8), INTENT(OUT) :: ev, cvv

      Tr = theta_vib/Tv
      dexp0 = DEXP(Tr)
      den   = dexp0 - 1.d0

      ! Energy and specific heat
      ev  = Rn2*theta_vib/den
      cvv = Rn2*dexp0*(Tr/den)**2 

    END SUBROUTINE ev_cvv_ho

    !------------------------------------------------------!
    ! This subroutine provides the species internal energy per unit mass 
    SUBROUTINE get_species_energy (temp, e)

      USE mod_nitrogen_DSMC_initialize_CFD,    ONLY: pos_N, pos_N2, nb_trot, nb_tvib, theta_vib, hf_n, Rn2, cv_tr, cv_rot 

      REAL(KIND=8) :: T, Tr, Tv, dexp_Tv, er, ev

      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp
      REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: e

      ! Temperatures 
      T  = temp(1)
      Tr = temp(posTr)
      Tv = temp(posTv)

      ! Rotational energy     
      er = cv_rot*Tr

      ! Vibrational energy
      IF ((nb_trot.EQ.1).AND.(nb_tvib.EQ.0)) THEN

         ev = 0.d0

      ELSE 

         dexp_Tv = DEXP(theta_vib/Tv)
         ev      = Rn2*theta_vib/(dexp_Tv - 1.d0)

      ENDIF      

      ! Internal energies
      e(pos_N)  = cv_tr(pos_N)*T  + hf_n
      e(pos_N2) = cv_tr(pos_N2)*T + er + ev 

    END SUBROUTINE get_species_energy 

    !------------------------------------------------------!
    ! This subroutine provides the species internal energy and constant volume
    ! specific heat per unit mass 
    SUBROUTINE get_species_energy_cv (temp, e, cv)

      USE mod_nitrogen_DSMC_initialize_CFD,    ONLY: pos_N, pos_N2, nb_temp, nb_trot, nb_tvib, theta_vib, &
                                                  &  hf_n, Rn2, cv_tr, cv_rot 

      REAL(KIND=8) :: den, T, Tr, Tv, dexp_Tv, er, ev

      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp
      REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: e, cv

      ! Temperatures 
      T  = temp(1)
      Tr = temp(posTr)
      Tv = temp(posTv)
    
      er      = cv_rot*Tr
      dexp_Tv = DEXP(theta_vib/Tv)
      ev      = Rn2*theta_vib/(dexp_Tv - 1.d0)
     
      IF ((nb_trot.EQ.1).AND.(nb_tvib.EQ.0)) ev = 0.d0
 
      ! Internal energies
      e(pos_N)  = cv_tr(pos_N)*T  + hf_n
      e(pos_N2) = cv_tr(pos_N2)*T + er + ev 

      ! Constant volume specific heat
      cv = cv_tr 
    
      ! Thermal equilibrium case
      IF (nb_temp.EQ.1) THEN

         den        = dexp_Tv - 1.d0
         cv(pos_N2) = cv(pos_N2) + cv_rot + Rn2*dexp_Tv*(theta_vib/(Tv*den))**2 

      ELSE
 
         ! Thermal nonequilibrium case (Tv ne T when using a separate temperature for rotation)
         IF ((nb_trot.EQ.0).AND.(nb_tvib.EQ.1)) THEN 

            cv(pos_N2) = cv(pos_N2) + cv_rot  

         ENDIF      
         
      ENDIF

    END SUBROUTINE get_species_energy_cv 

    !------------------------------------------------------!
    ! This subroutine provides the internal energy per unit mass of the N-N2 system
    SUBROUTINE energy (temp, er, ev, e)

      USE mod_nitrogen_DSMC_initialize_CFD,    ONLY: pos_N, pos_N2, nb_tvib, nb_trot, theta_vib, & 
                                                  &  hf_n, Rn2, cv_tr, cv_rot 

      REAL(KIND=8) :: T, Tr, Tv, dexp_Tv, erot

      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp
      REAL(KIND=8), INTENT(OUT) :: er, ev
      REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: e

      ! Temperatures 
      T  = temp(1)
      Tr = temp(posTr)
      Tv = temp(posTv)
    
      ! Rotational energy 
      er = cv_rot*Tr

      ! Vibrational energy (not considering in case of two temperature T and Tr)
      IF ((nb_trot.EQ.1).AND.(nb_tvib.EQ.0)) THEN

         ev = 0.d0

      ELSE 

         dexp_Tv = DEXP(theta_vib/Tv)
         ev      = Rn2*theta_vib/(dexp_Tv - 1.d0)

      ENDIF      

      ! Internal energies
      e(pos_N)  = cv_tr(pos_N)*T  + hf_n
      e(pos_N2) = cv_tr(pos_N2)*T + er + ev 

    END SUBROUTINE energy 

    !------------------------------------------------------!
    ! This subroutine provides the enthalpy per unit mass of the N-N2 system
    SUBROUTINE enthalpy (temp, er, ev, h)

      USE mod_nitrogen_DSMC_initialize_CFD,    ONLY: nb_ns, Ri 

      INTEGER :: is
      REAL(KIND=8) :: T

      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp
      REAL(KIND=8), INTENT(OUT) :: er, ev
      REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: h

      T = temp(1)
      CALL energy (temp, er, ev, h)

      DO is = 1,nb_ns 
         h(is) = h(is) + Ri(is)*T
      ENDDO

    END SUBROUTINE enthalpy  

    !------------------------------------------------------!
    ! This subroutine computes the temperature T for a given internal energy density (thermal equilibrium case).
    ! In this case non linear equation must be solved. The latter is solved by means of Newton-Raphson method.
    SUBROUTINE compute_T (rhoi, rho_eint, temp) 

      USE mod_nitrogen_DSMC_initialize_CFD,    ONLY: nb_ns, pos_N, pos_N2, hf_n, Rn, Rn2, cv_tr, theta_vib 

      INTEGER :: is, it
      INTEGER, PARAMETER :: it_max = 30
      REAL(KIND=8), PARAMETER :: tol = 1.d-6
      REAL(KIND=8) :: appr, exact
      REAL(KIND=8) :: den, dexp_T, ev, cv, tmp
      REAL(KIND=8) :: f, fp, g, r
      REAL(KIND=8) :: rho_N_hfn, rho_N2_Rn2
      REAL(KIND=8) :: T, Told
      REAL(KIND=8) :: res

      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, rho_eint
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp

      ! Useful data 
      rho_N_hfn  = rhoi(pos_N)*hf_n
      rho_N2_Rn2 = rhoi(pos_N2)*Rn2
      
      ! Initial temperature guess (vibration is considered half-excited)
      exact = rho_eint(1)
      g  = exact - rho_N_hfn
      fp  = 0.d0
      DO is = 1,nb_ns 
         fp = fp + rhoi(is)*cv_tr(is)
      ENDDO
      fp = fp + 1.5d0*rho_N2_Rn2 
      cv = fp - 0.5d0*rho_N2_Rn2

      T = g/fp
     
      ! Initialization
      it = 0
      res = 1.d0
      DO WHILE (res.GT.tol)

          Told = T 

          ! Vibrational energy (per unit volume) of N2 molecule
          r      = theta_vib/T
          dexp_T = DEXP(r)
          den    = dexp_T - 1.d0
          ev     = rho_N2_Rn2*theta_vib/den
 
          fp = 0.d0
          DO is = 1,nb_ns 
             fp = fp + rhoi(is)*cv_tr(is)
          ENDDO
          fp = fp + rho_N2_Rn2 + rho_N2_Rn2*dexp_T*(theta_vib/(T*den))**2

          tmp = cv*T + ev
          f   = tmp - g 
          T   = T - f/fp

          appr = tmp + rho_N_hfn 

          it = it + 1

          ! Residual computed in the internal energy density
          res = ABS((appr - exact)/exact)         

          IF (it.GT.it_max) THEN 
             PRINT*
             WRITE(*,'(A)')'solver_fvmcc_f90:excessive number of Newton iterators, in mod_nitrogen_DSMC_CFD_prop...'
             PRINT*
             PRINT*,'T ',T,'Told',Told,'res ',res
             PRINT*
             PAUSE
          ENDIF

      ENDDO
      
      ! Temperature value in output
      temp(1) = T

    END SUBROUTINE compute_T 
   
    !------------------------------------------------------!
    ! This subroutine computes the temperatures T and Tr for a given internal energy densities (thermal nonequilibrium case, 
    ! vibration is not taken into account). 
    SUBROUTINE compute_T_Tr (rhoi, rho_eint, temp) 

      USE mod_nitrogen_DSMC_initialize_CFD,    ONLY: pos_N, pos_N2, Rn2, Rn, hf_n 

      REAL(KIND=8) :: a, b, tmp
      REAL(KIND=8) :: T, Tr, rho_etr, rho_er
      REAL(KIND=8) :: rho_n, rho_n2 

      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, rho_eint
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp

      ! Initialization
      rho_etr = rho_eint(1)
      rho_er  = rho_eint(2)
      rho_n   = rhoi(pos_N)
      rho_n2  = rhoi(pos_N2) 
    
      IF (pos_N2.EQ.1) rho_n = 0.d0
 
      ! Rotational temperature  
      tmp = rho_n2*Rn2 
      Tr  = rho_er/tmp

      ! Translational temperature
      a = rho_etr - (rho_er + rho_n*hf_n)
      b = 1.5d0*(rho_n*Rn + tmp)

      T = a/b

      ! Temperature vector
      temp(1)     = T
      temp(posTr) = Tr

    END SUBROUTINE compute_T_Tr

    !------------------------------------------------------!
    ! This subroutine computes the temperatures T and Tv for a given internal energy densities (thermal nonequilibrium case). 
    SUBROUTINE compute_T_Tv (rhoi, rho_eint, temp) 

      USE mod_nitrogen_DSMC_initialize_CFD,    ONLY: pos_N, pos_N2, theta_vib, Rn2, Rn, hf_n 

      REAL(KIND=8) :: a, b
      REAL(KIND=8) :: T, Tv, rho_etr, rho_ev
      REAL(KIND=8) :: rho_n, rho_n2 

      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, rho_eint
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp

      ! Energy densoties
      rho_etr = rho_eint(1)
      rho_ev  = rho_eint(2)

      ! Species densities
      rho_n   = rhoi(pos_N)
      rho_n2  = rhoi(pos_N2) 

      IF (pos_N2.EQ.1) rho_n = 0.d0

      ! Vibrational temperature 
      Tv = theta_vib/DLOG(1.d0 + rho_n2*Rn2*theta_vib/rho_ev)

      ! Translational temperature
      a = rho_etr - (rho_ev + rho_n*hf_n)
      b = 0.5d0*(3.d0*rho_n*Rn + 5.d0*rho_n2*Rn2)

      T = a/b

      ! Temperature vector
      temp(1)     = T
      temp(posTv) = Tv

    END SUBROUTINE compute_T_Tv 

    !------------------------------------------------------!
    ! This subroutine computes the temperatures T, Tr and Tv for a given internal energy densities (thermal nonequilibrium case).
    SUBROUTINE compute_T_Tr_Tv (rhoi, rho_eint, temp) 

      USE mod_nitrogen_DSMC_initialize_CFD,    ONLY: pos_N, pos_N2, theta_vib, Rn2, Rn, hf_n 

      REAL(KIND=8) :: a, b
      REAL(KIND=8) :: T, Tr, Tv, rho_etr, rho_er, rho_ev
      REAL(KIND=8) :: rho_n, rho_n2 
      REAL(KIND=8) :: tmp

      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, rho_eint
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp

      ! Energy densities
      rho_etr = rho_eint(1)
      rho_er  = rho_eint(2)
      rho_ev  = rho_eint(3)

      ! Species densities
      rho_n   = rhoi(pos_N2)
      rho_n2  = rhoi(pos_N2) 

      IF (pos_N2.EQ.1) rho_n = 0.d0

      ! Rotational temperature
      tmp = Rn2*rho_n2
      Tr  = rho_er/tmp
      
      ! Vibrational temperature 
      Tv = theta_vib/DLOG(1.d0 + tmp*theta_vib/rho_ev)
      
      ! Translational temperature
      a = rho_etr - (rho_ev + rho_er + rho_n*hf_n)
      b = 1.5d0*(rho_n*Rn + tmp)

      T = a/b

      ! Temperature vector
      temp(1)     = T
      temp(posTr) = Tr
      temp(posTv) = Tv

    END SUBROUTINE compute_T_Tr_Tv 

    !------------------------------------------------------!
    ! This subroutine computes the nonequilibrium post-shock conditions
    SUBROUTINE post_shock_neq (p1, u1, T1, p2, u2, T2, yi, xN)

      USE mod_nitrogen_DSMC_initialize_CFD,         ONLY: pos_N2, pos_N, nb_ns, nb_temp, nb_trot, nb_tvib, urg, solver, & 
                                                       &  Ri, mm
      USE mod_Nitrogen_DSMC_CFD_eq,                 ONLY: eq_composition

      INTEGER :: is, length
      REAL(KIND=8), PARAMETER :: tol = 1.d-8
      REAL(KIND=8) :: f, fp, resR, resT, ratio, rhs, ratio_old, T_old
      REAL(KIND=8) :: mass, R
      REAL(KIND=8) :: tmp1, tmp2
      REAL(KIND=8) :: h1, h2, rho1, rho2, c1, g, gp1, M1, M1s, m_dot
      REAL(KIND=8), DIMENSION(nb_ns) :: xi, cvi, ei
      REAL(KIND=8), DIMENSION(nb_temp) :: temp
      REAL(KIND=8), DIMENSION(3) :: left, right, res

      REAL(KIND=8), INTENT(IN) :: p1, u1, T1
      REAL(KIND=8), INTENT(IN), OPTIONAL :: xN
      REAL(KIND=8), INTENT(OUT) :: p2, u2, T2
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi

      ! Useful data
      length = LEN_TRIM(solver)

      ! Molar fractions in the free-stream
      IF (PRESENT(xN)) THEN 

         xi(1) = xN
         xi(2) = 1.d0 - xN

      ELSE

         CALL eq_composition (p1, T1, xi)

      ENDIF
     
      ! Mixture molar mass 
      mass = 0.d0
      DO is = 1,nb_ns 
         mass = mass + xi(is)*mm(is)
      ENDDO
      R = urg/mass

      ! Species mass fractions
      DO is = 1,nb_ns 
         yi(is) = xi(is)*mm(is)/mass
      ENDDO

      ! Computation of post-shock conditions

      ! Thermal equilibrium case (a non linear system has to be solver)
      IF (nb_temp.EQ.1) THEN

         ! Mass, momentum and energy fluxes (pre-shock conditions)
         rho1  = p1/(R*T1)
         m_dot = rho1*u1

         h1   = 0.d0
         temp = T1
         CALL get_species_energy (temp, ei)
         DO is = 1,nb_ns 
            h1 = h1 + yi(is)*(ei(is) + Ri(is)*T1)
         ENDDO

         left(1) = m_dot 
         left(2) = m_dot*u1 + p1 
         left(3) = h1 + 0.5d0*u1**2
       
         ! Guess value for the density ratio 
         ratio = 0.1d0 
         resR  = 1.d0
         T2    = T1 
         DO WHILE (resR.GT.tol) 

            p2   = p1 + m_dot*u1*(1.d0 - ratio)
            rho2 = rho1/ratio
            
            rhs  = h1 + 0.5d0*u1**2*(1.d0 - ratio**2)
             
            ! Newton loop for post-shock temperature
            resT = 1.d0
            DO WHILE (resT.GT.tol) 

               f  = 0.d0
               fp = 0.d0

               temp = T2
               CALL get_species_energy_cv (temp, ei, cvi)
               DO is = 1,nb_ns
                  tmp1 = yi(is)
                  tmp2 = Ri(is) 
                  f    = f  + tmp1*(ei(is) + tmp2*T2)
                  fp   = fp + tmp1*(cvi(is) + tmp2)
               ENDDO

               T_old = T2 

               ! Post-shock temperature residual
               f  = f - rhs
               T2 = T2 - f/fp
               resT = ABS(T2 - T_old)/T_old

            ENDDO

            ratio_old = ratio
            
            rho2 = p2/(R*T2)
            
            ! Density ratio update and residual
            ratio = rho1/rho2
            resR = ABS(ratio - ratio_old)/ratio

         ENDDO

         u2 = u1*ratio

         ! Residuals
 
         ! Mass momentum and energy flux (post-shock)
         m_dot = rho2*u2
       
         h2   = 0.d0
         temp = T2
         CALL get_species_energy (temp, ei)
         DO is = 1,nb_ns 
            h2 = h2 + yi(is)*(ei(is) + Ri(is)*T2)
         ENDDO   
   
         right(1) = m_dot
         right(2) = m_dot*u2 + p2 
         right(3) = h2 + 0.5d0*u2**2

         DO is = 1,3 
            res(is) = ABS(right(is) - left(is))/left(is)*100.d0
         ENDDO
         
         WRITE(*,5)solver(1:length),':: Nitrogen DSMC -> thermal equilibrium post-shock conditions'
         PRINT*
         WRITE(*,10)'Residual on mass, momemtum and energy fluxes:'
         PRINT*
         WRITE(*,15)'Mass    ',res(1),' [%]'
         PRINT*
         WRITE(*,15)'Momentum',res(2),' [%]'
         PRINT*
         WRITE(*,15)'Energy  ',res(3),' [%]'
         PRINT*

      ELSE 

         ! Thermal nonequilibrium of rotation and vibration  
         IF ((nb_tvib.EQ.1).AND.(nb_trot.EQ.1)) THEN 

            g = 5.d0/3.d0 

         ! Thermal nonequilibrium of rotation only 
         ELSEIF ((nb_tvib.EQ.0).AND.(nb_trot.EQ.1)) THEN

            g = 5.d0/3.d0

         ! Thermal nonequilibrium of vibration only 
         ELSE 

            IF (pos_N2.EQ.1) THEN
               g = 7.d0/5.d0
            ELSE
               tmp1 = yi(pos_N)
               tmp2 = yi(pos_N2)
               g    = (5.d0*tmp1 + 7.d0*tmp2)/(3.d0*tmp1 + 5.d0*tmp2)
            ENDIF

         ENDIF

         gp1 = g + 1.d0
      
         c1  = SQRT(g*R*T1)
         M1  = u1/c1
         M1s = M1**2

         ! Pressure, velocity and temperature after the shock
         p2 = p1*(2.d0*g*M1s - g + 1.d0)/gp1
         u2 = u1 - c1*2.d0/gp1*(M1 - 1.d0/M1)
         T2 = T1*(2.d0*g*M1s - g + 1.d0)*(g - 1.d0 + 2.d0/M1s)/(gp1**2)

      ENDIF

5   FORMAT(A,A)
10  FORMAT(A)
15  FORMAT(A,E14.6,A)

    END SUBROUTINE post_shock_neq
    
  END MODULE mod_nitrogen_DSMC_CFD_prop
!------------------------------------------------------------------------------!

