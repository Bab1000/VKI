!------------------------------------------------------------------------------!
! This module provides subroutines for the computation of source of conservation equations for the N-N2 system
! when using Park multi-temperature model.
  MODULE mod_nitrogen_Park_CFD_source

  USE mod_nitrogen_Park_initialize_CFD,    ONLY: upi, una, urg, ukb, nb_ns, nb_tvib, nb_trot, nb_eq, nb_dim, nb_temp, & 
                                                 posTr, posTv, diss_imp_N2, vt_imp_mol, Tmin_kin

  IMPLICIT NONE

  REAL(KIND=8), PARAMETER :: rhoi_tol = 1.d-25
  REAL(KIND=8) :: Qn, Q_int

  ! Subroutine for the computations of source terms
  CONTAINS

    !------------------------------------------------------!
    ! This suroutine computes the source term due to collisional processes for the N-N2 system.
    SUBROUTINE source (rhoi, tvec, omega)

      USE mod_nitrogen_Park_initialize_CFD,           ONLY: mm_N, mm_N2, Rn, Rn2, fac_keq, c1, c2, c3, c4, theta_vib, & 
                                                            mw_a, mw_b, fac_Q, mu
      USE mod_nitrogen_Park_CFD_prop,                 ONLY: Q_trans, Q_int_rr_ho

      INTEGER :: is
      REAL(KIND=8), PARAMETER :: ratio = 2.5d0
      REAL(KIND=8), PARAMETER :: Zr = 23.3d0
      REAL(KIND=8), PARAMETER :: C  = 111.d0
      REAL(KIND=8), PARAMETER :: T0 = 300.55d0
      REAL(KIND=8), PARAMETER :: Tstar = 91.5d0
      REAL(KIND=8), PARAMETER :: mu0   = 1.781d-5
      REAL(KIND=8) :: sum1, sum2, tmp
      REAL(KIND=8) :: T, Tv, Tr, TTv, ln_T, ln_Tv, ln_TTv
      REAL(KIND=8) :: omega_VT, omega_CV, omega_RT, omega_CR
      REAL(KIND=8) :: tauc, visc, denom
      REAL(KIND=8) :: r, rv, dexp_r, dexp_rv, e, ev
      REAL(KIND=8) :: sigma, RT, T_m13, p, millikan, park, tau_VT_am, tau_VT_mm, tau_VT, tau_RT
      REAL(KIND=8) :: kfd, kfd_eq, keq, exp_1
      REAL(KIND=8) :: xn, xn2
      REAL(KIND=8), DIMENSION(nb_ns) :: rhoit

      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: tvec, rhoi
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: omega

      ! Temperatures and other useful data
      T  = tvec(1)
      Tr = tvec(posTr)
      Tv = tvec(posTv)

      ! Application of a fix in order to avoid numerical problems   
      IF (T.LT.Tmin_kin)  T = Tmin_kin
      IF (Tr.LT.Tmin_kin) Tr = Tmin_kin
      IF (Tv.LT.Tmin_kin) Tv = Tmin_kin

      TTv    = DSQRT(T*Tv)
      ln_T   = DLOG(T)
      ln_Tv  = DLOG(Tv)
      ln_TTv = DLOG(TTv) 

      ! Translational partition function of N
      Qn = fac_Q*(mm_N*T)**1.5
 
      ! Species concentrations (a fix is applied in order to avoid numerical problems)
      xn       = rhoi(1)/mm_N
      rhoit(1) = MAX(xn*mm_N,rhoi_tol)  

      xn2      = rhoi(2)/mm_N2
      rhoit(2) = MAX(xn2*mm_N2,rhoi_tol) 
     
      kfd    = c1*DEXP(c2*ln_TTv - c3/TTv)
      kfd_eq = c1*DEXP(c2*ln_T   - c3/T)

      ! Internal partition function of N2 (rigid rotator - harmonic oscillator model)
      CALL Q_int_rr_ho (T, Q_int)

      ! Equilibrium constant
      keq = fac_keq*Qn/Q_int*DEXP(-c3/T)

      ! Mass production terms (dissociation by N impact)
      exp_1    = xn*(kfd*xn2 - xn*xn*kfd_eq/keq)
      omega(1) = 2.d0*exp_1*mm_N
      omega(2) = - exp_1*mm_N2

      ! Dissociation by N2 impact (de-activated for comparison against the NASA Ames database for the N-N2 system)
      IF (diss_imp_N2.EQV..TRUE.) THEN

         tmp    = c4/c1
         kfd    = kfd*tmp
         kfd_eq = kfd_eq*tmp  

         exp_1  = xn2*(kfd*xn2 - xn*xn*kfd_eq/keq)
 
         omega(1) = omega(1) + 2.d0*exp_1*mm_N
         omega(2) = omega(2) - exp_1*mm_N2

      ENDIF

      omega(nb_ns + 1:nb_ns + nb_dim + 1) = 0.d0

      IF (nb_tvib.EQ.1) THEN

         ! Energy transfer source terms
         r  = theta_vib/T
         rv = theta_vib/Tv
         dexp_r  = DEXP(r)
         dexp_rv = DEXP(rv)
         e   = 1.d0/(dexp_r  - 1.d0)
         ev  = 1.d0/(dexp_rv - 1.d0)
      
         ! Limiting cross section for Park's correction [m^2]
         IF (T.GT.20000.d0) THEN
            sigma = 1.d-21*ratio*ratio
         ELSE
            tmp = 50000.d0/T
            sigma = 1.d-21*tmp*tmp
         ENDIF

         ! Static pressure
         p = (rhoit(1)*Rn + rhoit(2)*Rn2)*T

         ! Millikan and White relaxation time (plus Park correction for high temperatures) 
         T_m13     = T**(-0.3333333333)
         millikan  = DEXP(mw_a(1)*(T_m13 - mw_b(1)) - 18.42d0)*101325.d0/p
         park      = DSQRT(upi*mu*ukb*T/(8.d0*una))/(sigma*p)
         tau_VT_am = millikan + park

         ! VT transfer in N-N2 collisions only
         IF (vt_imp_mol.EQV..FALSE.) THEN

            tau_VT = tau_VT_am

         ! VT transfer also in N2-N2 collisions
         ELSE 

            millikan  = DEXP(mw_a(2)*(T_m13 - mw_b(2)) - 18.42d0)*101325.d0/p
            park      = DSQRT(0.5d0*upi*mm_N2*ukb*T/(8.d0*una))/(sigma*p)
            tau_VT_mm = millikan + park

            ! Frequency average the computation ofr VT energy transfer relaxation time
            sum1 = rhoit(1)/mm_N + rhoit(2)/mm_N2
            sum2 = rhoit(1)/(mm_N*tau_VT_am) + rhoit(2)/(mm_N2*tau_VT_mm)       
 
            tau_VT = sum1/sum2

         ENDIF

         tmp = theta_vib*Rn2

         Omega_VT = rhoit(2)*tmp*(e - ev)/tau_VT

         Omega_CV = tmp*ev*omega(2) 

         omega(nb_ns + nb_dim + posTv) = Omega_CV + Omega_VT

      ENDIF

      ! Parker model ror rotational relaxation
      IF (nb_trot.EQ.1) THEN

         ! Static pressure
         p = (rhoit(1)*Rn + rhoit(2)*Rn2)*T

         ! Dinamic viscosity (sutherland's law is used) 
         visc = mu0*((T0 + C)/(T + C))*(T/T0)**1.5

         ! Collision time 
         tauc = upi*visc/(4.d0*p)

         ! Parker expression for rotational relaxation time
         tmp    = Tstar/T
         denom  = (1 + 0.5d0*(upi**1.5)*(tmp**0.5) + (0.25d0*upi**2 + upi)*tmp)
         tau_RT = Zr*tauc/denom        

         omega_CR = Rn2*Tr*omega(2) 
         omega_RT = rhoit(2)*Rn2*(T - Tr)/tau_RT

         omega(nb_ns + nb_dim + posTr) = Omega_CR + Omega_RT

      ENDIF

    END SUBROUTINE source

    !------------------------------------------------------!
    ! This suroutine computes the source term and its Jacobian due to collisional processes for the N-N2 system.
    ! For sake of simplicity the Jacobian is computed with respect to primitive variables. The transformaton in order
    ! to obtaind its expression in terms of conservative variables is performed outside.
    SUBROUTINE source_Jac (rhoi, tvec, omega, domega_dp)

      USE mod_nitrogen_Park_initialize_CFD,           ONLY: mm_N, mm_N2, fac_keq, m_ratio, ov_m_ratio, c1, c2, c3, c4
      USE mod_nitrogen_Park_CFD_prop,                 ONLY: Q_trans_der, Q_int_rr_ho_der

      INTEGER :: is, start
      REAL(KIND=8) :: T, Tr, Tv, TTv, ln_T, ln_Tv, ln_TTv, omega_VT, omega_CV
      REAL(KIND=8) :: kfd, kbd, kfd_eq, keq, exp_1, exp_2, exp_3, tmp1, tmp2
      REAL(KIND=8) :: xn, xn2, xns 
      REAL(KIND=8) :: df_dT, df_dTv, dkf_dT, dkf_dTv, dkeq_dT, fac_dkf, dkf_eq_dT, dQn_dT, dQint_dT, dexp_eq
      
      REAL(KIND=8), DIMENSION(nb_ns + nb_temp) :: dOmegaV_dP 
      REAL(KIND=8), DIMENSION(nb_ns) :: rhoit

      REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: tvec, rhoi
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: omega
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: domega_dp

      ! Initialization
      omega     = 0.d0
      domega_dp = 0.d0

      ! Temperatures and other useful data
      T  = tvec(1)
      Tr = tvec(posTr) 
      Tv = tvec(posTv)

      ! Application of a fix in order to avoid numerical problems   
      IF (T.LT.Tmin_kin)  T = Tmin_kin
      IF (Tr.LT.Tmin_kin) Tr = Tmin_kin
      IF (Tv.LT.Tmin_kin) Tv = Tmin_kin

      TTv    = DSQRT(T*Tv)
      ln_T   = DLOG(T)
      ln_Tv  = DLOG(Tv)
      ln_TTv = DLOG(TTv) 

      ! Translational partition function of N atom and derivative with respect to temperature 
      CALL Q_trans_der (T, mm_N, Qn, dQn_dT)  
 
      ! Species concentrations (a fix is applied in order to avoid numerical problems)
      xn       = rhoi(1)/mm_N
      rhoit(1) = MAX(xn*mm_N,rhoi_tol) 

      xn2      = rhoi(2)/mm_N2
      rhoit(2) = MAX(xn2*mm_N2,rhoi_tol)

      kfd    = c1*DEXP(c2*ln_TTv - c3/TTv)
      kfd_eq = c1*DEXP(c2*ln_T   - c3/T)

      ! Internal partition function of N2 molecule and derivative with respect to temperature (rigid rotator - harmonic oscillator model)
      CALL Q_int_rr_ho_der (T, Q_int, dQint_dT)

      ! Equilibrium constant
      dexp_eq = DEXP(-c3/T)
      keq = fac_keq*Qn/Q_int*dexp_eq

      ! Backward dissociation rate coefficient
      kbd = kfd_eq/keq

      xns = xn**2
      tmp1 = xns*kbd
      tmp2 = kfd*xn2
      exp_1 = xn*(tmp2 - tmp1)
      exp_2 = tmp2 - 3.d0*tmp1 
      exp_3 = kfd*xn         

      omega(1) = 2.d0*exp_1*mm_N
      omega(2) = - exp_1*mm_N2
     
      omega(nb_ns + 1:nb_ns + nb_dim + 1) = 0.d0

      ! Computation of derivatives with respect to temperatures of rate coefficients
      ! Thermal nonequilibrium case
      IF (nb_temp.GT.1) THEN

        fac_dkf   = c1*0.5d0*DEXP((c2 - 2.d0)*ln_TTv  - c3/TTv)*(c2 + c3/TTv)
        dkf_eq_dT = c1*dexp_eq*DEXP((c2 - 1.d0)*ln_T)*(c2 + c3/T) 

        dkf_dT  = fac_dkf*Tv
        dkf_dTv = fac_dkf*T

      ! Thermal equilibrium case
      ELSE

        dkf_dT  = kfd*(c2/T + c3/T**2)
        dkf_dTv = dkf_dT
        dkf_eq_dT = dkf_dT

      ENDIF

      dkeq_dT = fac_keq*dexp_eq/Q_int*(c3/T**2*Qn + (dQn_dT - Qn/Q_int*dQint_dT))
      df_dT   = xn*(xn2*dkf_dT - xns/keq*(dkf_eq_dT - dkeq_dT*kbd))
      df_dTv  = xn*xn2*dkf_dTv

      ! N species continuity equation
      domega_dp(1) = 2.d0*exp_2
      domega_dp(2) = 2.d0*exp_3*m_ratio
      domega_dp(3) = 2.d0*mm_N*df_dT 

      ! N2 species continuity equation
      start = 3 + nb_temp - 1
      domega_dp(start + 1) = - exp_2*ov_m_ratio
      domega_dp(start + 2) = - exp_3
      domega_dp(start + 3) = - mm_N2*df_dT      
 
      IF (nb_tvib.EQ.1) THEN

         start = nb_ns + posTv 
         domega_dp(start)   = 2.d0*mm_N*df_dTv 
         domega_dp(2*start) = - mm_N2*df_dTv  
         
      ENDIF
 
      ! Dissociation by N2 impact (de-activated for comparison against the NASA Ames database for the N-N2 system)
      IF (diss_imp_N2.EQV..TRUE.) THEN
       
         ! Re-computation of rate coefficients and derivatives needed for the Jacobian
         tmp1    = c4/c1
         kfd     = kfd*tmp1
         kfd_eq  = kfd_eq*tmp1  
         kbd     = kfd_eq/keq
         dkf_dT  = dkf_dT*tmp1
         dkf_dTv = dkf_dTv*tmp1
         dkf_eq_dT = dkf_eq_dT*tmp1

         df_dT  = xn2*(xn2*dkf_dT - xns/keq*(dkf_eq_dT - dkeq_dT*kbd))
         df_dTv = xn2**2*dkf_dTv

         exp_1 = xn2*(kfd*xn2 - xns*kbd)
         exp_2 = - 2.d0*kbd*xn*xn2 
         exp_3 = 2.d0*kfd*xn2 - xns*kbd 
 
         omega(1) = omega(1) + 2.d0*exp_1*mm_N
         omega(2) = omega(2) - exp_1*mm_N2

         ! N species continuity equation
         domega_dp(1) = domega_dp(1) + 2.d0*exp_2 
         domega_dp(2) = domega_dp(2) + 2.d0*exp_3*m_ratio
         domega_dp(3) = domega_dp(3) + 2.d0*mm_N*df_dT 

         ! N2 species continuity equation
         start = 3 + nb_temp - 1
         domega_dp(start + 1) = domega_dp(start + 1) - exp_2*ov_m_ratio
         domega_dp(start + 2) = domega_dp(start + 2) - exp_3
         domega_dp(start + 3) = domega_dp(start + 3) - mm_N2*df_dT    

         IF (nb_tvib.EQ.1) THEN

            start = nb_ns + posTv 
            domega_dp(start)   = domega_dp(start) + 2.d0*mm_N*df_dTv 
            domega_dp(2*start) = domega_dp(2*start) - mm_N2*df_dTv  

         ENDIF

      ENDIF

      !! THE ROTATIONAL CONTRIBUTION MUST BE ADDED !!
      IF (nb_tvib.EQ.1) THEN

         CALL energy_transfer (rhoit, T, Tv, omega(2), Omega_VT, Omega_CV, domega_dp(5:8), dOmegaV_dP)

         omega(nb_ns + nb_dim + posTv) = Omega_CV + Omega_VT

         ! N2 vibrational energy conservation equation
         start = 2*(nb_ns + nb_temp)
         DO is = 1,nb_ns + nb_temp 
            domega_dP(start + is) =  dOmegaV_dP(is)
         ENDDO

      ENDIF

      !! THE ROTATIONAL CONTRIBUTION MUST BE ADDED !!
      IF (nb_trot.EQ.1) THEN

         PRINT*
         WRITE(*,'(A)')'solver_fvmcc_f90:in mod_nitrogen_Park_CFD_source, not implemented yet'
         PRINT*
         STOP

      ENDIF

    END SUBROUTINE source_Jac

    !------------------------------------------------------!
    ! This subroutine computes the source term due to VT and CV transfer in collisions between N and N2 species
    SUBROUTINE energy_transfer (rhoi, T, Tv, omega_n2, Omega_VT, Omega_CV, domega_n2_dP, dOmegaV_dP)

      USE mod_nitrogen_Park_initialize_CFD,         ONLY: mm_N, mm_N2, mu, theta_vib, mw_a, mw_b

      REAL(KIND=8), PARAMETER :: ratio = 2.5d0
      REAL(KIND=8) :: e, ev, sigma, p, millikan, park, tau_VT, tau_VT_am, tau_VT_mm, dtauVT_dn, dtauVT_dn2, dtauVT_dT
      REAL(KIND=8) :: dtauVT_am_dn, dtauVT_am_dn2, dtauVT_mm_dn, dtauVT_mm_dn2, dtauVT_am_dT, dtauVT_mm_dT
      REAL(KIND=8) :: sum1, sum2, dsum2_dn, dsum2_dn2, dsum2_dT, sum2s
      REAL(KIND=8) :: tmp1, tmp2, tmp3, tmp4
      REAL(KIND=8) :: cp, cpv, OmegaVT_div_tauVT, fac_1, fac_2, Tr, Trv, dexp_Tr, dexp_Trv, RT, T_m13
      REAL(KIND=8) :: tauVT_div_p, dmillikan_dT, dpark_dT, dsigma_dT

      REAL(KIND=8), INTENT(IN) :: T, Tv, omega_n2
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
      REAL(KIND=8), INTENT(OUT) :: Omega_VT, Omega_CV
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: domega_n2_dP
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: dOmegaV_dP

      ! Vibrational energy
      Tr  = theta_vib/T
      Trv = theta_vib/Tv
      dexp_Tr  = DEXP(Tr)
      dexp_Trv = DEXP(Trv)
      e   = 1.d0/(dexp_Tr  - 1.d0)
      ev  = 1.d0/(dexp_Trv - 1.d0)
     
      ! Useful factors
      fac_1 = theta_vib*urg/mm_N2
      fac_2 = urg/mm_N2

      ! Vibrational specific heats
      cp    = fac_2*dexp_Tr*(Tr*e)**2
      cpv   = fac_2*dexp_Trv*(Trv*ev)**2
 
      ! Limiting cross section for Park's correction [m^2]
      IF (T.GT.20000.d0) THEN
         sigma = 1.d-21*ratio*ratio
         dsigma_dT = 0.d0
      ELSE
         sigma = 1.d-21*(50000.d0/T)**2
         dsigma_dT = - 2.d0*sigma/T
      ENDIF

      ! Static pressure
      RT = urg*T
      p = (rhoi(1)/mm_N + rhoi(2)/mm_N2)*RT

      ! VT transfer in N-N2 collisions only
      IF (vt_imp_mol.EQV..FALSE.) THEN

         ! Millikan and White relaxation time (plus Park correction for high temperatures) 
         T_m13     = T**(-0.3333333333)
         millikan  = DEXP(mw_a(1)*(T_m13 - mw_b(1)) - 18.42d0)*101325.d0/p
         park      = DSQRT(upi*mu*ukb*T/(8.d0*una))/(sigma*p)
         tau_VT    = millikan + park

         Omega_VT = rhoi(2)*fac_1*(e - ev)/tau_VT
         Omega_CV = fac_1*ev*omega_n2 

         ! Useful variables
         OmegaVT_div_tauVT = Omega_VT/tau_VT
         tauVT_div_p = tau_VT/p

         ! Derivative of tau_VT (N-N2 interaction) with respect to species densities (N and N2) and temperature
         dtauVT_dn   = - tauVT_div_p*RT/mm_N 
         dtauVT_dn2  = - tauVT_div_p*RT/mm_N2  

         dmillikan_dT = - millikan/T*(mw_a(1)*T_m13/3.d0 + 1.d0)
         dpark_dT     = - park*(1.d0/sigma*dsigma_dT + 0.5d0/T)
         dtauVT_dT    =   dmillikan_dT + dpark_dT

       ! VT transfer also in N2-N2 collisions
       ELSE

         ! Useful data
         T_m13 = T**(-0.3333333333)
         tmp1  = 101325.d0/p
         tmp2  = T_m13/3.d0
         tmp3  = 1.d0/sigma*dsigma_dT
         tmp4  = 0.5d0/T

         ! N-N2 VT interaction
         millikan  = DEXP(mw_a(1)*(T_m13 - mw_b(1)) - 18.42d0)*tmp1
         park      = DSQRT(upi*mu*ukb*T/(8.d0*una))/(sigma*p)
         tau_VT_am = millikan + park

         ! Derivative of tau_VT (N-N2 interaction) with respect to species densities (N and N2) and temperature
         tauVT_div_p    = tau_VT_am/p

         dmillikan_dT = - millikan/T*(mw_a(1)*tmp2 + 1.d0)
         dpark_dT     = - park*(tmp3 + tmp4)

         dtauVT_am_dn  = - tauVT_div_p*RT/mm_N 
         dtauVT_am_dn2 = - tauVT_div_p*RT/mm_N2 
         dtauVT_am_dT  = dmillikan_dT + dpark_dT

         ! N2-N2 VT interaction
         millikan  = DEXP(mw_a(2)*(T_m13 - mw_b(2)) - 18.42d0)*tmp1
         park      = park*DSQRT(0.5d0*mm_N2/mu)
         tau_VT_mm = millikan + park  

         ! Derivative of tau_VT (N2-N2 interaction) with respect to species densities (N and N2) and temperature
         tauVT_div_p = tau_VT_mm/p

         dmillikan_dT = - millikan/T*(mw_a(2)*tmp2 + 1.d0)
         dpark_dT     = - park*(tmp3 + tmp4)

         dtauVT_mm_dn  = - tauVT_div_p*RT/mm_N 
         dtauVT_mm_dn2 = - tauVT_div_p*RT/mm_N2 
         dtauVT_mm_dT  = dmillikan_dT + dpark_dT

         ! Frequency average the for the computation VT energy transfer relaxation time
         sum1  = rhoi(1)/mm_N + rhoi(2)/mm_N2
         sum2  = rhoi(1)/(mm_N*tau_VT_am) + rhoi(2)/(mm_N2*tau_VT_mm)   
         sum2s = 1.d0/sum2**2    
 
         tau_VT = sum1/sum2

         dsum2_dn  = 1.d0/(mm_N*tau_VT_am)*(1.d0 - rhoi(1)/tau_VT_am*dtauVT_am_dn) - rhoi(2)/(mm_N2*tau_VT_mm**2)*dtauVT_mm_dn
         dsum2_dn2 = 1.d0/(mm_N2*tau_VT_mm)*(1.d0 - rhoi(2)/tau_VT_mm*dtauVT_mm_dn2) - rhoi(1)/(mm_N*tau_VT_am**2)*dtauVT_am_dn2
         dsum2_dT  = - (rhoi(1)/(mm_N*tau_VT_am**2)*dtauVT_am_dT +  rhoi(2)/(mm_N2*tau_VT_mm**2)*dtauVT_mm_dT)

         ! Final derivatives of relaxation time
         dtauVT_dn  = (sum2/mm_N - dsum2_dn*sum1)*sum2s
         dtauVT_dn2 = (sum2/mm_N2 - dsum2_dn2*sum1)*sum2s 
         dtauVT_dT  = - sum1*dsum2_dT*sum2s

         Omega_VT = rhoi(2)*fac_1*(e - ev)/tau_VT
         Omega_CV = fac_1*ev*omega_n2 

         OmegaVT_div_tauVT = Omega_VT/tau_VT

      ENDIF

      ! Vibrational energy per unit mass of N2 molecule (at Tv)
      ev = fac_1*ev

      ! Jacobian matrix row (vibrational energy conservation equation)
      ! Derivatives with respect to N and N2 species densities
      dOmegaV_dP(1) = - OmegaVT_div_tauVT*dtauVT_dn  + domega_n2_dP(1)*ev
      dOmegaV_dP(2) = - OmegaVT_div_tauVT*dtauVT_dn2 + Omega_VT/rhoi(2) + domega_n2_dP(2)*ev

      ! Derivatives with respect to T, Tr and Tv
      tau_VT = 1.d0/tau_VT
      dOmegaV_dP(3) =   rhoi(2)*cp*tau_VT  - OmegaVT_div_tauVT*dtauVT_dT + domega_n2_dP(3)*ev
      dOmegaV_dP(4) = 0.d0
      dOmegaV_dP(3 + nb_tvib) = - rhoi(2)*cpv*tau_VT + domega_n2_dP(4)*ev + omega_n2*cpv

    END SUBROUTINE energy_transfer

  END MODULE mod_nitrogen_Park_CFD_source
!------------------------------------------------------------------------------!
