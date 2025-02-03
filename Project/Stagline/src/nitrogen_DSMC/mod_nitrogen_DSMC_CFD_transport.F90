!------------------------------------------------------------------------------!
! This modules provides subroutine for computing transport properties and fluxes for the N-N2 system 
! when using the DSMC multi-temperature model. It also provides subroutines to be interfaced with 
! CFD codes for computing mass diffusion flux of species and diffusive components of heat flux vector
! (in this way model details are hidden to the CFD code that can be therefore written in a general and 
! flexible manner). 
  MODULE mod_nitrogen_DSMC_CFD_transport

    USE mod_nitrogen_DSMC_initialize_CFD,        ONLY: pos_N2, pos_N2, posTr, posTv, nb_ns, nb_trot, nb_tvib, nb_dim, & 
                                                     & una, ukb, upi, urg
 
    IMPLICIT NONE

    REAL(KIND=8), PARAMETER :: Tmin = 100.d0
    REAL(KIND=8), PARAMETER :: fac  = 15.d0/4.d0

    ! Subroutines dealing with transport properties and species diffusion velocities
    CONTAINS 

      !----------------------------------------------------!
      ! This subroutine computes transport coefficients. 
      SUBROUTINE get_transport_coeff (nb, xi, temp, mu, kappa, lambda_tr, lambda_rot, lambda_vib, Di, chi)

        USE mod_nitrogen_DSMC_initialize_CFD,        ONLY: Sc_ref, theta_vib, cv_rot, Ri
       
        INTEGER :: is
        REAL(KIND=8) :: T, Tr, Tv
        REAL(KIND=8) :: tmp1, tmp2, tmp3
        REAL(KIND=8) :: cv_vib
        
        REAL(KIND=8), INTENT(IN) :: nb
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, temp 
        REAL(KIND=8), INTENT(OUT) :: kappa, mu, lambda_tr, lambda_rot, lambda_vib
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: chi, Di
      
        ! Temperatures (a fix is applied in order to avoid numerical problems)
        T  = MAX(temp(1),Tmin)
        Tr = MAX(temp(posTr),Tmin) 
        Tv = MAX(temp(posTv),Tmin)
        
        ! No thermal diffusion 
        chi = 0.d0

        ! No bulk viscosity 
        kappa = 0.d0

        ! N2 mixture 
        IF (nb_ns.EQ.1) THEN

           ! Dynamic viscosity 
           CALL dyn_visc (T, mu)

           ! Thermal conductivity (translational component)
           lambda_tr = fac*mu*Ri(pos_N2)

           ! Eucken correction for internal degrees of freedom (rotation and vibration)
           IF ((nb_trot.EQ.1).AND.(nb_tvib.EQ.0))THEN

              cv_vib = 0.d0

           ELSE 

              ! Vibrational specific heat
              tmp1  = theta_vib/Tv
              tmp2  = DEXP(tmp1)
              tmp3  = tmp2 - 1.d0
 
              cv_vib = Ri(pos_N2)*tmp2*(tmp1/tmp3)**2 

           ENDIF
 
           ! The hypothesis of constant Schmidt number is done (in this way the
           ! D11 self-diffusion coefficient can be computed)
           tmp1       = Sc_ref*mu
           lambda_rot = cv_rot*tmp1  
           lambda_vib = cv_vib*tmp1

           ! Diffusion coefficient of N2 set to zero
           Di = 0.d0

        ! N-N2 mixture
        ELSE 

           PRINT*
           WRITE(*,'(A)')'in mod_nitrogen_DSMC_CFD:transport.F90, not implemented yet...'
           PRINT*

        ENDIF  

      END SUBROUTINE get_transport_coeff

      !----------------------------------------------------!
      ! This subroutine computes the species mass diffusion flux.
      SUBROUTINE get_species_DiffFlux (T, nb, xi, diff_driv, Ji)  

        REAL(KIND=8), INTENT(IN) :: T, nb
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, diff_driv
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: Ji

        ! N2 mixture
        IF (nb_ns.EQ.1) THEN

           ! No mass diffusion 
           Ji  = 0.d0

        ! N-N2 mixture
        ELSE 

           PRINT*
           WRITE(*,'(A)')'in mod_nitrogen_DSMC_CFD:transport.F90, not implemented yet...'
           PRINT*

        ENDIF

      END SUBROUTINE get_species_DiffFlux

      !----------------------------------------------------!
      ! This function computes the dynamic viscosity (to be used for mixture of pure N2) 
      SUBROUTINE dyn_visc (T, eta)

        USE mod_nitrogen_DSMC_initialize_CFD,        ONLY: eta_ref, T_ref, exp_visc

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8) :: eta

        eta = eta_ref*(T/T_ref)**exp_visc

      END SUBROUTINE dyn_visc 

  END MODULE mod_nitrogen_DSMC_CFD_transport
!------------------------------------------------------------------------------!
