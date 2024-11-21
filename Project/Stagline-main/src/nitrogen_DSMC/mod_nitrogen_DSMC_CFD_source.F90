!------------------------------------------------------------------------------!
! This module provides subroutines for the computation of source of conservation equations for the N-N2 system
! when using DSMC multi-temperature model.
  MODULE mod_nitrogen_DSMC_CFD_source

    USE mod_nitrogen_DSMC_initialize_CFD,        ONLY: pos_N2, pos_N2, posTr, posTv, nb_ns, nb_trot, nb_tvib, & 
                                                     & nb_temp, nb_dim, una, ukb, upi, urg

    IMPLICIT NONE

    REAL(KIND=8), PARAMETER :: rhoi_tol = 1.d-25
    REAL(KIND=8), PARAMETER :: Tmin     = 250.d0

    ! Subroutine for the computations of source terms
    CONTAINS

      !------------------------------------------------------!
      ! This suroutine computes the source term due to collisional processes for the N-N2 system.
      SUBROUTINE source (rhoi, tvec, omega)

        USE mod_nitrogen_DSMC_initialize_CFD,     ONLY: eta_ref, T_ref, exp_visc, theta_vib, par1, par2, ov_alpha, & 
                                                      & Zr_inf, Zrv, T_star, Ri, mm

        INTEGER :: is
        REAL(KIND=8) :: tmp1, tmp2, tmp3, tmp4
        REAL(KIND=8) :: T, Tr, Tv, p, eta
        REAL(KIND=8) :: Zrc
        REAL(KIND=8) :: e, ev
        REAL(KIND=8) :: omega_RT, omega_CR, omega_VT, omega_CV
        REAL(KIND=8) :: tau_RT, tau_VT
        REAL(KIND=8) :: r, rv, dexp_r, dexp_rv 
        REAL(KIND=8), DIMENSION(nb_ns) :: rhoit

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: tvec, rhoi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: omega

        ! Temperatures (a fix is applied in order to avoid numerical problems)
        T  = MAX(tvec(1),Tmin)
        Tr = MAX(tvec(posTr),Tmin)
        Tv = MAX(tvec(posTv),Tmin)

        ! Species densities (a fix is applied in order to avoid numerical problems)
        DO is = 1,nb_ns 
           rhoit(is) = MAX(rhoi(is),rhoi_tol)
        ENDDO

        ! Initialization
        omega = 0.d0

        ! Mixture static pressure 
        p = 0.d0
        DO is = 1,nb_ns 
           p = p + rhoit(is)*Ri(is)
        ENDDO
        p = p*T
        
        ! Dynamic viscosity 
        eta = eta_ref*(T/T_ref)**exp_visc

        ! Source term due to rotational nonequilibrium 
        IF (nb_trot.EQ.1) THEN

           ! Common data
           tmp1 = Ri(pos_N2)
           tmp2 = rhoit(pos_N2)

           ! Parker formula for Zrc 
           tmp3 = T_star/T
           tmp4 = DSQRT(tmp3)
           Zrc  = Zr_inf/(1.d0 + par1*tmp3 + par2*tmp4)   
           Zrc = 10.d0

           ! Relaxation time (Zrc*tauc)
           tau_RT = Zrc*(6.d0 - ov_alpha)*(4.d0 - ov_alpha)*eta/(p*30.d0)          

           omega_RT = tmp1*tmp2*(T - Tr)/tau_RT
           omega_CR = omega(pos_N2)*tmp1*Tr

           omega(nb_ns + nb_dim + posTr) = omega_RT + omega_CR

        ENDIF
      
        ! Source term due to vibrational nonequilibrium 
        IF (nb_tvib.EQ.1) THEN

           ! Common data
           tmp1 = Ri(pos_N2)*theta_vib
           tmp2 = rhoit(pos_N2)

           r       = theta_vib/T
           rv      = theta_vib/Tv
           dexp_r  = DEXP(r)
           dexp_rv = DEXP(rv)
           e       = 1.d0/(dexp_r  - 1.d0)
           ev      = 1.d0/(dexp_rv - 1.d0)

           ! Relaxation time (a fixed Zrv is assumed)
           tau_VT = Zrv*(6.d0 - ov_alpha)*(4.d0 - ov_alpha)*eta/(p*30.d0)          

           omega_VT = tmp1*tmp2*(e - ev)/tau_VT
           omega_CV = omega(pos_N2)*tmp1*ev

           omega(nb_ns + nb_dim + posTv) = omega_VT + omega_CV

        ENDIF
     
      END SUBROUTINE source

      !------------------------------------------------------!
      ! This suroutine computes the source term and its Jacobian due to collisional processes for the N-N2 system.
      ! For sake of simplicity the Jacobian is computed with respect to primitive variables. The transformaton in order
      ! to obtaind its expression in terms of conservative variables is performed outside.
      SUBROUTINE source_Jac (rhoi, tvec, omega, domega_dp)
   
        INTEGER :: is
        REAL(KIND=8) :: T, Tr, Tv
        REAL(KIND=8), DIMENSION(nb_ns) :: rhoit

        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: tvec, rhoi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: omega
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: domega_dp
  
        ! Temperatures (a fix is applied in order to avoid numerical problems)
        T  = MAX(tvec(1),Tmin)
        Tr = MAX(tvec(posTr),Tmin)
        Tv = MAX(tvec(posTv),Tmin)

        ! Species densities (a fix is applied in order to avoid numerical problems)
        DO is = 1,nb_ns 
           rhoit(is) = MAX(rhoi(is),rhoi_tol)
        ENDDO

        omega     = 0.d0
        domega_dp = 0.d0  

        PRINT*
        WRITE(*,'(A)')'in mod_nitrogen_DSMC_CFD_source.F90, not implemented yet...'
        PRINT*
        STOP

     END SUBROUTINE source_Jac

  END MODULE mod_nitrogen_DSMC_CFD_source
!------------------------------------------------------------------------------!
