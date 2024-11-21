!------------------------------------------------------------------------------!
! This module provides subroutines for the computation of source of conservation equations when the mixture 
! is made on lony N2 molecules.
  MODULE mod_nitrogen_Park_CFD_source_VT_RT

  USE mod_nitrogen_Park_initialize_CFD,    ONLY: upi, una, urg, ukb, nb_ns, nb_tvib, nb_trot, nb_eq, nb_dim, nb_temp, & 
                                                 posTr, posTv, pos_N2, Tmin

  IMPLICIT NONE

  REAL(KIND=8), PARAMETER :: rhoi_tol = 1.d-25

  ! Subroutine for the computations of source terms
  CONTAINS

    !------------------------------------------------------!
    ! This suroutine computes the source term due to collisional processes between N2 molecules
    SUBROUTINE VT_RT_transf (rhoi, tvec, omega)

      USE mod_nitrogen_Park_initialize_CFD,           ONLY: mm_N2, Rn2, theta_vib, mw_a, mw_b

      REAL(KIND=8), PARAMETER :: ratio = 2.5d0
      REAL(KIND=8), PARAMETER :: Zr = 23.3d0
      REAL(KIND=8), PARAMETER :: C  = 111.d0
      REAL(KIND=8), PARAMETER :: T0 = 300.55d0
      REAL(KIND=8), PARAMETER :: Tstar = 91.5d0
      REAL(KIND=8), PARAMETER :: mu0   = 1.781d-5
      REAL(KIND=8) :: tmp
      REAL(KIND=8) :: T, Tv, Tr, rhoit
      REAL(KIND=8) :: omega_VT, omega_RT
      REAL(KIND=8) :: a, b, denom, dexp_r, dexp_rv, e, ev, r, rv, Trank, tauc, visc
      REAL(KIND=8) :: sigma, T_m13, p, millikan, park
      REAL(KIND=8) :: tau_VT, tau_RT

      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: tvec, rhoi
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: omega

      ! Initialization
      omega = 0.d0

      ! Temperatures and other useful data
      T  = tvec(1)
      Tr = tvec(posTr)
      Tv = tvec(posTv)

      ! Application of a fix in order to avoid numerical problems   
      IF (T.LT.Tmin)  T = Tmin
      IF (Tr.LT.Tmin) Tr = Tmin
      IF (Tv.LT.Tmin) Tv = Tmin 

      rhoit = MAX(rhoi(pos_N2),rhoi_tol) 

      ! Static pressure 
      p = rhoit*Rn2*T

      ! VT energy transfer term (Landau-Teller model is used) 
      IF (nb_tvib.EQ.1) THEN

         ! Commom factors
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

         ! Millikan and White relaxation time (plus Park correction for high temperatures) 
         T_m13    = T**(-0.3333333333)
         millikan = DEXP(mw_a(pos_N2)*(T_m13 - mw_b(pos_N2)) - 18.42d0)*101325.d0/p
         park     = DSQRT(0.5d0*upi*mm_N2*ukb*T/(8.d0*una))/(sigma*p)
         tau_VT   = millikan + park
         Omega_VT = rhoit*theta_vib*Rn2*(e - ev)/tau_VT

         omega(nb_ns + nb_dim + posTv) = Omega_VT

     ENDIF

     ! RT energy transfer term (Parker model is used)
     IF (nb_trot.EQ.1) THEN

        ! Static pressure
        p = rhoit*Rn2*T

        ! Dinamic viscosity (sutherland's law is used) 
        visc = mu0*((T0 + C)/(T + C))*(T/T0)**1.5

        ! Collision time 
        tauc = upi*visc/(4.d0*p)

        ! Parker expression for rotational relaxation time
        tmp    = Tstar/T
        denom  = (1 + 0.5d0*(upi**1.5)*(tmp**0.5) + (0.25d0*upi**2 + upi)*tmp)
        tau_RT = Zr*tauc/denom 
 
        omega_RT = rhoit*Rn2*(T - Tr)/tau_RT

        omega(nb_ns + nb_dim + posTr) = Omega_RT

      ENDIF

    END SUBROUTINE VT_RT_transf

    !------------------------------------------------------!
    ! This suroutine computes the source term due to collisional processes between N2 molecules
    SUBROUTINE VT_RT_transf_Jac (rhoi, tvec, omega, domega_dp)

      USE mod_nitrogen_Park_initialize_CFD,           ONLY: mm_N2, Rn2, theta_vib, mw_a, mw_b

      REAL(KIND=8), PARAMETER :: ratio = 2.5d0
      REAL(KIND=8), PARAMETER :: Zr = 235.d-1
      REAL(KIND=8), PARAMETER :: C  = 111.d0
      REAL(KIND=8), PARAMETER :: T0 = 540.99d0
      REAL(KIND=8), PARAMETER :: mu0 = 1781.d-8
      REAL(KIND=8) :: rhoit, p, T_m13, RT, theta_R
      REAL(KIND=8) :: T, Tv, Tr, sigma, dsigma_dT
      REAL(KIND=8) :: cp, cpv, dexp_r, dexp_rv, e, ev, r, rv 
      REAL(KIND=8) :: millikan, park, tau_VT, ov_tau_VT
      REAL(KIND=8) :: OmegaVT_div_tauVT, tauVT_div_p, dtauVT_rho, dmillikan_dT, dpark_dT, dtauVT_dT, & 
                    & dtauVT_drho 
      REAL(KIND=8) :: omega_VT, new, Told

      REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: tvec, rhoi
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: omega
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: domega_dp

      ! Initialization
      omega = 0.d0
      domega_dp = 0.d0

      ! Temperatures and other useful data
      T  = tvec(1)
      Tr = tvec(posTr)
      Tv = tvec(posTv)

      ! Application of a fix in order to avoid numerical problems   
      IF (T.LT.Tmin)  T = Tmin
      IF (Tr.LT.Tmin) Tr = Tmin
      IF (Tv.LT.Tmin) Tv = Tmin 

      rhoit = MAX(rhoi(pos_N2),rhoi_tol) 

      ! Static pressure
      RT = Rn2*T 
      p  = rhoit*RT

      ! VT energy transfer term and related derivatives (Landau-Teller model is used) 
      IF (nb_tvib.EQ.1) THEN

         ! Vibrational energy
         r  = theta_vib/T
         rv = theta_vib/Tv
         dexp_r  = DEXP(r)
         dexp_rv = DEXP(rv)
         e   = 1.d0/(dexp_r  - 1.d0)
         ev  = 1.d0/(dexp_rv - 1.d0)
     
         ! Useful factor
         theta_R = theta_vib*Rn2

         ! Vibrational specific heats
         cp  = Rn2*dexp_r*(r*e)**2
         cpv = Rn2*dexp_rv*(rv*ev)**2
 
         ! Limiting cross section for Park's correction [m^2]
         IF (T.GT.20000.d0) THEN
            sigma = 1.d-21*ratio*ratio
            dsigma_dT = 0.d0
         ELSE
            sigma = 1.d-21*(50000.d0/T)**2
            dsigma_dT = - 2.d0*sigma/T
         ENDIF

         ! Millikan and White relaxation time (plus Park correction for high temperatures) 
         T_m13     = T**(-0.3333333333)
         millikan  = DEXP(mw_a(pos_N2)*(T_m13 - mw_b(pos_N2)) - 18.42d0)*101325.d0/p
         park      = DSQRT(0.5d0*upi*mm_N2*ukb*T/(8.d0*una))/(sigma*p)
         tau_VT    = millikan + park
         ov_tau_VT = 1.d0/tau_VT

         ! VT energy transfer term
         Omega_VT = rhoit*theta_R*(e - ev)*ov_tau_VT

         ! Useful variables for derivative evaluation
         OmegaVT_div_tauVT = Omega_VT*ov_tau_VT
         tauVT_div_p = tau_VT/p

         ! Derivative of tau_VT with respect to the density of N2 and temperature
         dtauVT_rho   = - tauVT_div_p*RT/mm_N2  
         dmillikan_dT = - millikan/T*(mw_a(pos_N2)*T_m13/3.d0 + 1.d0)
         dpark_dT     = - park*(1.d0/sigma*dsigma_dT + 0.5d0/T)
         dtauVT_dT    = dmillikan_dT + dpark_dT

         ! Source term
         omega(nb_ns + nb_dim + posTv) = Omega_VT

         ! Source term Jacobian 
         ! Derivative with respect to N2 species density
         domega_dp(4) = - OmegaVT_div_tauVT*dtauVT_rho + Omega_VT/rhoi(pos_N2)   

         ! Derivatives with respect to T, Tr and Tv
         domega_dp(5) = rhoit*cp*ov_tau_VT  - OmegaVT_div_tauVT*dtauVT_dT
         domega_dp(5 + nb_tvib) = - rhoit*cpv*ov_tau_VT

      ENDIF 

      ! RT energy transfer term and related derivatives (Parker model is used)
      IF (nb_trot.EQ.1) THEN

         PRINT*,'in "VT_RT_transf_Jac", not implemented yet...'
         STOP

      ENDIF  

    END SUBROUTINE VT_RT_transf_Jac

  END MODULE mod_nitrogen_Park_CFD_source_VT_RT
!------------------------------------------------------------------------------!
