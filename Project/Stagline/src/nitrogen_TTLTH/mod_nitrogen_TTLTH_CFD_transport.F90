!------------------------------------------------------------------------------!
! This modules provides subroutine for computing transport properties and fluxes for the N-N2 system 
! when using the TTLTH multi-temperature model. It also provides subroutines to be interfaced with 
! CFD codes for computing mass diffusion flux of species and diffusive components of heat flux vector
! (in this way model details are hidden to the CFD code that can be therefore written in a general and 
! flexible manner). 
  MODULE mod_nitrogen_TTLTH_CFD_transport

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(3) :: sigma_11, sigma_12, sigma_22, A_star, B_star, C_star

    ! Subroutines dealing with transport properties and fluxes
    CONTAINS 

      !----------------------------------------------------!
      ! This subroutine computes transport coefficients, mass diffusion flux and diffusive heat flux vector components. 
      ! Pre-kached values of collision integrals are used.
      SUBROUTINE get_transport_coeff_diff_heat_flux (ndim, rhoi, temp, grad_rhoi, grad_temp, mu, lambda, rhoi_Vdi, diff_heat)

        USE mod_nitrogen_TTLTH_initialize_CFD,        ONLY: posTr, posTv, nb_ns, una, ukb, upi, urg, mm_N, mm_N2, theta_vib, & 
                                                           Rn, Rn2, cv_tr, cv_rot, hf_n

        INTEGER :: i, j, pos_T, pos_sp
        REAL(KIND=8) :: fac1, fac2, fac3, fac4, tmp1, tmp2, x2
        REAL(KIND=8) :: nb, T, Tv, Tr, m_sum, m_prod, x_prod, ratio
        REAL(KIND=8) :: den, det_mu, det_lam, dexp0, cv_vib, ev, pres, rho, er
        REAL(KIND=8) :: lambda_tr, lambda_rot, lambda_vib
        REAL(KIND=8), DIMENSION(nb_ns) :: chi, nbi, xi, mui, mass, coeff_mu, coeff_lam, h
        REAL(KIND=8), DIMENSION(3) :: g_mu, g_lam, D_ij

        INTEGER, INTENT(IN) :: ndim
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp, grad_rhoi, grad_temp
        REAL(KIND=8), INTENT(OUT) :: mu
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: lambda, rhoi_Vdi, diff_heat

        ! Temperatures 
        T  = temp(1)
        Tr = temp(posTr)
        Tv = temp(posTv)

        ! Species molar fractions and number densities
        mass(1) = mm_N/una
        mass(2) = mm_N2/una

        nb = 0.d0
        DO i = 1,nb_ns 
           tmp1  = rhoi(i)/mass(i)
           nb    = nb + tmp1
           nbi(i) = tmp1
        ENDDO

        DO i = 1,nb_ns 
           xi(i) = nbi(i)/nb
        ENDDO

        ! Useful common factors
        m_sum  = mass(1) + mass(2)
        m_prod = mass(1)*mass(2)
        x_prod = xi(1)*xi(2)

        ! Collision integrals 
        CALL collision (T, sigma_11, sigma_12, sigma_22, A_star, B_star, C_star)

        tmp1 = DSQRT(upi*ukb*T) 
        DO i = 1,nb_ns 
 
           tmp2    = DSQRT(mass(i))

           mui(i)  = tmp1*tmp2/sigma_22(i)
           D_ij(i) = 2.d0*tmp1/(tmp2*nb*sigma_11(i))

        ENDDO

        D_ij(3) = tmp1*DSQRT(2.d0*m_sum/m_prod)/(nb*sigma_11(3))

        mui  = 5.d0/16.d0*mui
        D_ij = 3.d0/16.d0*D_ij 

        ! Linear system for viscosity and translational component of thermal conductivity
        
        ! Computation of diagonal terms (common factors are initially stored in order to speed up computations)
        fac1   = 16.d0/5.d0*x_prod*DSQRT(2.d0*m_prod/(upi*ukb*T*m_sum**3))
        fac3   = 1.d0/25.d0/(ukb*nb*D_ij(3))*x_prod/m_sum**2

        DO i = 1,nb_ns 

           x2   = xi(i)*xi(i) 
           tmp1 = x2/mui(i)
           tmp2 = 4.d0/15.d0*x2*mass(i)/(ukb*mui(i))
 
           j = 2 - i + 1

           fac2 = sigma_22(3)*mass(j)/mass(i) + 5.d0/3.d0*sigma_11(3)
           fac4 = 30.d0*mass(i)**2 + 25.d0*mass(j)**2 - 12.d0*mass(j)**2*B_star(3) + 16.d0*m_prod*A_star(3)   

           g_mu(i)  = tmp1 + fac1*fac2
           g_lam(i) = tmp2 + fac3*fac4  

        ENDDO

        ! Off diagonal term (linear systems are both symmetric)
        fac2 = sigma_22(3) - 5.d0/3.d0*sigma_11(3)
        fac4 = 16.d0*A_star(3) + 12.d0*B_star(3) - 55.d0

        g_mu(3)  = fac1*fac2
        g_lam(3) = fac3*fac4*m_prod  

        ! Solving linear systems by means of Cramer's rule (simple 2 x 2 algebraic system)
        fac1 = g_mu(3)
        fac2 = g_lam(3)

        det_mu  = g_mu(1)*g_mu(2) - fac1**2
        det_lam = g_lam(1)*g_lam(2) - fac2**2

        DO i = 1,nb_ns 

           j = 2 - i + 1

           tmp1 = xi(i)
           tmp2 = xi(j)

           coeff_mu(i)  = tmp1*g_mu(j) - tmp2*fac1 
           coeff_lam(i) = tmp1*g_lam(j) - tmp2*fac2 

        ENDDO

        coeff_mu  = coeff_mu/det_mu 
        coeff_lam = coeff_lam/det_lam

        ! Viscosity and translational component of thermal conductivity
        mu = 0.d0
        lambda_tr = 0.d0
        DO i = 1,nb_ns 
           fac1 = xi(i)
           mu = mu + coeff_mu(i)*fac1
           lambda_tr = lambda_tr + coeff_lam(i)*fac1
        ENDDO        

        ! Species enthalpies and N2 internal specific heats
        ratio = theta_vib/Tv
        dexp0 = DEXP(ratio)
        den   = dexp0 - 1.d0
 
        er     = cv_rot*Tr
        ev     = Rn2*theta_vib/den
        cv_vib = Rn2*dexp0*(ratio/den)**2 

        ! Eucken correction for internal degrees of freedom (rotation and vibration)
        tmp1 = xi(2)/D_ij(2) + xi(1)/D_ij(3)
        tmp1 = 1.d0/tmp1  
        tmp2 = rhoi(2) 

        lambda_rot = cv_rot*tmp1*tmp2
        lambda_vib = cv_vib*tmp1*tmp2

        ! Thermal conductivity vector 
        lambda    = 0.d0
        lambda(1) = lambda_tr
        lambda(posTr) = lambda(posTr) + lambda_rot
        lambda(posTv) = lambda(posTv) + lambda_vib
  
        h(1) = cv_tr(1)*T + Rn*T + hf_n 
        h(2) = cv_tr(2)*T + er + ev + Rn2*T 

        ! Thermal diffusion ratios 
        chi = 0.d0

        ! Mixture static pressure and density
        rho  = rhoi(1) + rhoi(2)
        pres = (rhoi(1)*Rn + rhoi(2)*Rn2)*T

        ! Species diffusion flux and diffusive heat flux vector
        tmp1 = rhoi(1)*rhoi(2)
        tmp2 = D_ij(3)*tmp1/(x_prod*pres*rho**2)
        DO i = 1,ndim

           pos_T  = 2*(i - 1)
           pos_sp = nb_ns*(i - 1)

           fac1 = tmp1*(Rn  - Rn2)*grad_temp(pos_T + 1)

           fac2 = T*(rhoi(2)*grad_rhoi(pos_sp + 1)*Rn -  rhoi(1)*grad_rhoi(pos_sp + 2)*Rn2)

           fac3 = - tmp2*(fac1 + fac2)
             
           rhoi_Vdi(pos_sp + 1) =  fac3
           rhoi_Vdi(pos_sp + 2) = -fac3   

           ! Rotational and vibrational energy conservation equations
           diff_heat(pos_T + posTr) = - fac3*er
           diff_heat(pos_T + posTv) = - fac3*ev

           ! Global energy conservation equation
           diff_heat(pos_T + 1) =  fac3*(h(1) - h(2))

        ENDDO
        
      END SUBROUTINE get_transport_coeff_diff_heat_flux

      !----------------------------------------------------!
      ! This subroutine computes transport coefficients (Dynamic viscosity and all components of thermal conductivity). 
      ! Pre-kached values are not used here.
      SUBROUTINE get_transport_coeff (rhoi, temp, mu, lambda_tr, lambda_rot, lambda_vib, lambda_tot)

        USE mod_nitrogen_TTLTH_initialize_CFD,         ONLY: posTv, nb_ns, una, ukb, upi, mm_N, mm_N2, theta_vib, Rn2, cv_rot

        INTEGER :: i, j
        REAL(KIND=8) :: fac1, fac2, fac3, fac4, tmp1, tmp2, x2
        REAL(KIND=8) :: nb, T, Tv, Tr, m_sum, m_prod, x_prod
        REAL(KIND=8) :: den, det_mu, det_lam, dexp0, cv_vib, ev
        REAL(KIND=8), DIMENSION(nb_ns) :: nbi, xi, mui, mass, coeff_mu, coeff_lam, h
        REAL(KIND=8), DIMENSION(3) :: g_mu, g_lam, D_ij, sigma_11, sigma_12, sigma_22, A_star, B_star, C_star
        
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), INTENT(OUT) :: mu, lambda_tr, lambda_rot, lambda_vib, lambda_tot

        ! Temperatures 
        T  = temp(1)
        Tv = temp(posTv)

        ! Collision integrals 
        CALL collision (T, sigma_11, sigma_12, sigma_22, A_star, B_star, C_star)

        ! Species molar fractions and number densities
        mass(1) = mm_N/una
        mass(2) = mm_N2/una

        nb = 0.d0
        DO i = 1,nb_ns 
           tmp1  = rhoi(i)/mass(i)
           nb    = nb + tmp1
           nbi(i) = tmp1
        ENDDO

        DO i = 1,nb_ns 
           xi(i) = nbi(i)/nb
        ENDDO

        ! Useful common factors
        m_sum  = mass(1) + mass(2)
        m_prod = mass(1)*mass(2)
        x_prod = xi(1)*xi(2)

        ! Species viscosity and diffusion coefficients
        tmp1 = DSQRT(upi*ukb*T) 
        DO i = 1,nb_ns 
 
           tmp2    = DSQRT(mass(i))

           mui(i)  = tmp1*tmp2/sigma_22(i)
           D_ij(i) = 2.d0*tmp1/(tmp2*nb*sigma_11(i))

        ENDDO

        D_ij(3) = tmp1*DSQRT(2.d0*m_sum/m_prod)/(nb*sigma_11(3)) 

        mui  = 5.d0/16.d0*mui
        D_ij = 3.d0/16.d0*D_ij 

        ! Linear system for viscosity and translational component of thermal conductivity
        
        ! Computation of diagonal terms (common factors are initially stored in order to speed up computations)
        fac1   = 16.d0/5.d0*x_prod*DSQRT(2.d0*m_prod/(upi*ukb*T*m_sum**3))
        fac3   = 1.d0/25.d0/(ukb*nb*D_ij(3))*x_prod/m_sum**2

        DO i = 1,nb_ns 

           x2   = xi(i)*xi(i) 
           tmp1 = x2/mui(i)
           tmp2 = 4.d0/15.d0*x2*mass(i)/(ukb*mui(i))
 
           j = 2 - i + 1

           fac2 = sigma_22(3)*mass(j)/mass(i) + 5.d0/3.d0*sigma_11(3)
           fac4 = 30.d0*mass(i)**2 + 25.d0*mass(j)**2 - 12.d0*mass(j)**2*B_star(3) + 16.d0*m_prod*A_star(3)   

           g_mu(i)  = tmp1 + fac1*fac2
           g_lam(i) = tmp2 + fac3*fac4  

        ENDDO

        ! Off diagonal term (linear systems are both symmetric)
        fac2 = sigma_22(3) - 5.d0/3.d0*sigma_11(3)
        fac4 = 16.d0*A_star(3) + 12.d0*B_star(3) - 55.d0

        g_mu(3)  = fac1*fac2
        g_lam(3) = fac3*fac4*m_prod  

        ! Solving linear systems by means of Cramer's rule (simple 2 x 2 algebraic system)
        fac1 = g_mu(3)
        fac2 = g_lam(3)

        det_mu  = g_mu(1)*g_mu(2) - fac1**2
        det_lam = g_lam(1)*g_lam(2) - fac2**2

        DO i = 1,nb_ns 

           j = 2 - i + 1

           tmp1 = xi(i)
           tmp2 = xi(j)

           coeff_mu(i)  = tmp1*g_mu(j) - tmp2*fac1 
           coeff_lam(i) = tmp1*g_lam(j) - tmp2*fac2 

        ENDDO
      
        coeff_mu  = coeff_mu/det_mu 
        coeff_lam = coeff_lam/det_lam

        ! Dynamic viscosity and translational component of thermal conductivity
        tmp1 = 0.d0
        tmp2 = 0.d0
        DO i = 1,nb_ns 
           fac1 = xi(i)
           tmp1 = tmp1 + coeff_mu(i)*fac1
           tmp2 = tmp2 + coeff_lam(i)*fac1
        ENDDO        

        mu        = tmp1
        lambda_tr = tmp2

        ! N2 internal specific heats
        Tr    = theta_vib/Tv
        dexp0 = DEXP(Tr)
        den   = dexp0 - 1.d0
 
        cv_vib = Rn2*dexp0*(Tr/den)**2 

        ! Eucken correction for internal degrees of freedom (rotation and vibration)
        tmp1 = xi(2)/D_ij(2) + xi(1)/D_ij(3)
        tmp1 = 1.d0/tmp1  
        tmp2 = rhoi(2) 

        lambda_rot = cv_rot*tmp1*tmp2
        lambda_vib = cv_vib*tmp1*tmp2

        ! Total thermal conductivity
        lambda_tot = lambda_tr + lambda_rot + lambda_vib

      END SUBROUTINE get_transport_coeff

      !----------------------------------------------------!
      ! This subroutine computes the species diffusion fluxes.
      SUBROUTINE get_species_diff_flux (ndim, temp, rhoi, grad_rhoi, grad_temp, rhoi_Vdi) 

        USE mod_nitrogen_TTLTH_initialize_CFD,             ONLY: nb_ns, una, ukb, upi, mm_N, mm_N2, Rn, Rn2

        INTEGER :: i, j, pos_T, pos_sp
        REAL(KIND=8) :: fac1, fac2, fac3, fac4, tmp1, tmp2
        REAL(KIND=8) :: nb, T, m_sum, m_prod, x_prod, pres, rho
        REAL(KIND=8), DIMENSION(nb_ns) :: nbi, xi, mass, h, chi
        REAL(KIND=8), DIMENSION(3) :: D_ij, sigma_11, sigma_12, sigma_22, A_star, B_star, C_star

        INTEGER, INTENT(IN) :: ndim
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp, rhoi, grad_rhoi, grad_temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: rhoi_Vdi

        ! Temperature 
        T  = temp(1)

        ! Collision integrals 
        CALL collision (T, sigma_11, sigma_12, sigma_22, A_star, B_star, C_star)

        ! Species molar fractions and number densities
        mass(1) = mm_N/una
        mass(2) = mm_N2/una

        nb = 0.d0
        DO i = 1,nb_ns 
           tmp1  = rhoi(i)/mass(i)
           nb    = nb + tmp1
           nbi(i) = tmp1
        ENDDO

        DO i = 1,nb_ns 
           xi(i) = nbi(i)/nb
        ENDDO

        ! Useful common factors
        m_sum  = mass(1) + mass(2)
        m_prod = mass(1)*mass(2)
        x_prod = xi(1)*xi(2)

        ! Species viscosity and diffusion coefficients
        tmp1 = DSQRT(upi*ukb*T) 
        DO i = 1,nb_ns 
 
           tmp2    = DSQRT(mass(i))
           D_ij(i) = 2.d0*tmp1/(tmp2*nb*sigma_11(i))

        ENDDO

        D_ij(3) = tmp1*DSQRT(2.d0*m_sum/m_prod)/(nb*sigma_11(3)) 

        D_ij = 3.d0/16.d0*D_ij 

        ! Thermal diffusion ratios 
        chi = 0.d0

        ! Mixture static pressure and density
        rho  = rhoi(1) + rhoi(2)
        pres = (rhoi(1)*Rn + rhoi(2)*Rn2)*T

        ! Species diffusion flux and diffusive heat flux vector
        tmp1 = rhoi(1)*rhoi(2)
        tmp2 = pres*rho 
        DO i = 1,ndim

           pos_T  = 2*(i - 1)
           pos_sp = nb_ns*(i - 1)

           fac1 = (ukb*tmp1/tmp2*(1.d0/mass(1) - 1.d0/mass(2)))*grad_temp(pos_T + 1)

           fac2 = ukb*T/tmp2*(rhoi(2)/mass(1)*grad_rhoi(pos_sp + 1) -  rhoi(1)/mass(2)*grad_rhoi(pos_sp + 2))

           fac3 = - D_ij(3)*(fac1 + fac2)/x_prod

           fac4 = fac3*tmp1/rho

           rhoi_Vdi(pos_sp + 1) =  fac4
           rhoi_Vdi(pos_sp + 2) = -fac4   

        ENDDO

      END SUBROUTINE get_species_diff_flux 

      !----------------------------------------------------!
      ! This subrotine computes the collision integrals at a given temperature for N-N, N2-N2 and N-N2 interactions.
      SUBROUTINE collision (T, sigma_11, sigma_12, sigma_22, A_star, B_star, C_star)

        USE mod_nitrogen_TTLTH_initialize_CFD,       ONLY: upi, q11, q12, q22, bs, cs, nb_ns

        INTEGER :: i, nb_inter
        INTEGER :: pos1, pos2, pos3, pos4
        REAL(KIND=8), PARAMETER :: fac = 1.d-20
        REAL(KIND=8) :: tmp1, tmp2
        REAL(KIND=8) :: lnt, lnt2, lnt3

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: sigma_11, sigma_12, sigma_22, A_star, B_star, C_star

        ! Useful data
        lnt  = DLOG(T)
        lnt2 = lnt*lnt
        lnt3 = lnt2*lnt

        ! Number of interactions
        nb_inter = nb_ns + 1
        
        pos1 = 1 
        pos2 = 2
        pos3 = 3
        pos4 = 4

        ! Computation of collision integrals
        DO i = 1,nb_inter 

           tmp1 = fac*DEXP(q11(pos1)*lnt3 + q11(pos2)*lnt2 + q11(pos3)*lnt + q11(pos4))   
           tmp2 = fac*DEXP(q22(pos1)*lnt3 + q22(pos2)*lnt2 + q22(pos3)*lnt + q22(pos4)) 

           sigma_11(i) = tmp1
           sigma_22(i) = tmp2
           A_star(i)   = tmp2/tmp1

           B_star(i)   = DEXP(bs(pos1)*lnt3 + bs(pos2)*lnt2 + bs(pos3)*lnt + bs(pos4)) 
           tmp2        = DEXP(cs(pos1)*lnt3 + cs(pos2)*lnt2 + cs(pos3)*lnt + cs(pos4)) 

           C_star(i)   = tmp2
           sigma_12(i) = tmp2/tmp1

           pos1 = pos1 + 4
           pos2 = pos2 + 4
           pos3 = pos3 + 4
           pos4 = pos4 + 4

        ENDDO 

      END SUBROUTINE collision

  END MODULE mod_nitrogen_TTLTH_CFD_transport
!------------------------------------------------------------------------------!
