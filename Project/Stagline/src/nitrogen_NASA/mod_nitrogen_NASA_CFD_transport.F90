!------------------------------------------------------------------------------!
! This modules provides subroutine for computing transport properties and fluxes for the N-N2 system 
! when using the NASA Ames ab-initio database. It also provides subroutines to be interfaced with 
! CFD codes for computing mass diffusion flux of species and diffusive components of heat flux vector
! (in this way model details are hidden to the CFD code that can be therefore written in a general and 
! flexible manner). 
  MODULE mod_nitrogen_NASA_CFD_transport

    USE mod_nitrogen_NASA_initialize_CFD,         ONLY: nb_ns, nb_ground, nb_inter

    IMPLICIT NONE

    ! Subroutines dealing with transport properties and fluxes
    CONTAINS 

      !----------------------------------------------------!
      ! This subroutine adds a small number to the composition respecting the mass constraint
      SUBROUTINE comp_tol (tol, xi, xi_tol)

        INTEGER :: i
        REAL(KIND=8) :: tmp
        REAL(KIND=8) :: sum_xi

        REAL(KIND=8), INTENT(IN) :: tol
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: xi_tol

        ! Initialization
        sum_xi = 0.d0
        DO i = 1,nb_ns
           tmp = xi(i) 
           xi_tol(i) = tmp + tol
           sum_xi    = sum_xi + tmp + tol
        ENDDO

        sum_xi = 1.d0/sum_xi
        DO i = 1,nb_ns 
           xi_tol(i) = xi_tol(i)*sum_xi
        ENDDO

      END SUBROUTINE comp_tol

      !----------------------------------------------------!
      ! This subroutine computes transport coefficients (dynamic viscosity and all components of thermal conductivity). 
      SUBROUTINE get_transport_coeff (nb, xi, yi, temp, mu, kappa, lambda_tr, lambda_int, Di, chi)
     
        USE mod_nitrogen_NASA_initialize_CFD,         ONLY: una, ukb, upi, urg, fac_mu, fac_Diff, mi, & 
                                                          & Ri, sqrt_mi, sqrt_mij, Tmin_CI
        USE mod_function_pointer_NASA,                ONLY: get_lambda_int

        INTEGER :: i, j, ii, jj, ij
        REAL(KIND=8) :: tmp1, tmp2, fac1, fac2 
        REAL(KIND=8) :: x_i, x_prod, m_i, m_j, miij, mjij,  m_sum
        REAL(KIND=8) :: miij_miij, miij_mjij, mjij_mjij
        REAL(KIND=8) :: Aij, Bij, w
        REAL(KIND=8) :: T, sqr_T
        REAL(KIND=8), DIMENSION(nb_ground) :: xig, mui, lambdai, xis, coeff_mu, coeff_lam
        REAL(KIND=8), DIMENSION(nb_inter) :: xi_xj_ov_nDij,  xi_xj_ov_Dij, Dij
        REAL(KIND=8), DIMENSION(nb_inter) :: Gij_mu, Gij_lam
        REAL(KIND=8), DIMENSION(nb_inter) :: Omega_11, Omega_12, Omega_22, Astar, Bstar, Cstar

        REAL(KIND=8), INTENT(IN) :: nb
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, yi, temp 
        REAL(KIND=8), INTENT(OUT) :: kappa, mu, lambda_tr, lambda_int
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: chi, Di

        ! Temperatures (a fix is applied in order to avoid numerical problems)
        T  = MAX(temp(1),Tmin_CI)
        sqr_T = DSQRT(T)
      
        ! Molar fractions of N and N2
        xig(1) = xi(1)
        tmp1 = 0.d0
        DO i = 2,nb_ns 
           tmp1 = tmp1 + xi(i) 
        ENDDO
        xig(2) = tmp1

        DO i = 1,nb_ground
           xis(i) = xig(i)**2
        ENDDO

        ! No thermal diffusion 
        chi = 0.d0

        ! No bulk viscosity
        kappa = 0.d0

        ! Compute collision integrals  
        CALL collision_table (T, Omega_11, Omega_12, Omega_22, Astar, Bstar, Cstar) 

        ! Compute the binary diffusion coefficients
        fac1 = fac_Diff/sqr_T
        DO i = 1,nb_ground 

           m_i = mi(i)

           DO j = i,nb_ground

              m_j = mi(j)
              x_prod = xig(i)*xig(j)

              ! Index
              ij = ((i - 1)*(2*nb_ground - i) + 2*j)/2

              ! Common factor
              tmp1 = fac1*Omega_11(ij)*x_prod*sqrt_mij(ij) 

              xi_xj_ov_nDij(ij) = tmp1
              Dij(ij) = x_prod/(nb*tmp1)

           ENDDO

        ENDDO
       
        ! Compute single species viscosity and translational component of thermal conductivity
        fac1 = fac_mu*sqr_T
        fac2 = 15.d0/4.d0
        DO i = 1,nb_ground

           ii = ((i - 1)*(2*nb_ground - i) + 2*i)/2 
           tmp1   = fac1*sqrt_mi(i)/Omega_22(ii)
           mui(i) = tmp1 
           lambdai(i) = fac2*tmp1*Ri(i)

        ENDDO 
 
        ! Assembly linear transport systems for viscosity and translational component of thermal conductivity
        ! Initialization
        Gij_mu  = 0.d0
        Gij_lam = 0.d0

        DO i = 1,nb_ground
           ii = ((i - 1)*(2*nb_ground - i) + 2*i)/2
           Gij_mu(ii)  = xis(i)/mui(i)
           Gij_lam(ii) = xis(i)/lambdai(i)
        ENDDO

        DO i = 1,nb_ground 

           m_i = mi(i)

           DO j = i + 1,nb_ground 

              m_j = mi(j)

              ! Indices
              ij  = ((i - 1)*(2*nb_ground - i) + 2*j)/2 
              ii  = ((i - 1)*(2*nb_ground - i) + 2*i)/2 
              jj  = ((j - 1)*(2*nb_ground - j) + 2*j)/2

              ! Collision integrals 
              Aij = Astar(ij)
              Bij = Bstar(ij)

              ! Common factors
              m_sum = 1.d0/(m_i + m_j)
              miij = m_i*m_sum 
              mjij = m_j*m_sum
              miij_miij = miij*miij
              miij_mjij = miij*mjij 
              mjij_mjij = mjij*mjij   

              fac1 = 2.d0*una*xi_xj_ov_nDij(ij)*m_sum
              fac2 = (1.d0/(25.d0*ukb))*xi_xj_ov_nDij(ij)

              ! Dynamic viscosity matrix
              Gij_mu(ij) = fac1*(- 1.d0 + 0.6d0*Aij)
              Gij_mu(ii) = Gij_mu(ii) + fac1*(1.d0 + 0.6d0*m_j/m_i*Aij)
              Gij_mu(jj) = Gij_mu(jj) + fac1*(1.d0 + 0.6d0*m_i/m_j*Aij)

              ! Translational component of thermal conductivity matrix
              Gij_lam(ij) = - fac2*miij_mjij*(55.d0 - 16.d0*Aij - 12.d0*Bij)
              Gij_lam(ii) = Gij_lam(ii) + fac2*(30.d0*miij_miij + 16.d0*miij_mjij*Aij + (25.d0 - 12.d0*Bij)*mjij_mjij)
              Gij_lam(jj) = Gij_lam(jj) + fac2*(30.d0*mjij_mjij + 16.d0*miij_mjij*Aij + (25.d0 - 12.d0*Bij)*miij_miij) 

           ENDDO

        ENDDO

        ! Solve the linear transport sytems for dynamic viscosity and translational component of thermal conductivity
        CALL solve_transp_syst(xig, Gij_mu, coeff_mu)
        CALL solve_transp_syst(xig, Gij_lam, coeff_lam)

        ! Compute dynamic viscosity and translational component of thermal conductivity from the solution of the 
        ! linear transport systems
        mu = 0.d0
        lambda_tr = 0.d0
        DO i = 1,nb_ground
           x_i = xig(i)
           mu        = mu + x_i*coeff_mu(i)
           lambda_tr = lambda_tr + x_i*coeff_lam(i)
        ENDDO

        ! Compute internal component of thermal conductivity (Eucken's correction is applied)
        CALL get_lambda_int(nb, T, xi, xig, Dij, lambda_int)

        ! Average diffusion coefficients
        ! Nitrogen atom N 
        i = 1
        tmp1 = 1.d0 - yi(i)
        tmp2 = 0.d0
        DO j = i + 1,nb_ns
           ij = ((i - 1)*(2*nb_ground - i) + 2*nb_ground)/2
           tmp2 = tmp2 + xi(j)/Dij(ij)
        ENDDO
         
        DO j = i + 1,nb_ground
           ij = ((i - 1)*(2*nb_ground - i) + 2*nb_ground)/2
           tmp2 = tmp2 + xig(j)/Dij(ij)
        ENDDO
     
        Di(i) = tmp1/tmp2

        ! Nitrogen molecule N2
        i = 2
        tmp1 = 1.d0 - yi(i)

        tmp2 = 0.d0
        DO j = 1,i - 1
           ij = ((i - 1)*(2*nb_ground - i) + 2*j)/2
           tmp2 = tmp2 + xi(j)/Dij(ij)
        ENDDO
              
        DO j = i + 1,nb_ns
           ij = ((i - 1)*(2*nb_ground - i) + 2*nb_ground)/2
           tmp2 = tmp2 + xi(j)/Dij(ij)
        ENDDO
             
        tmp1 = tmp1/tmp2
        DO i = 2,nb_ns 
           Di(i) = tmp1
        ENDDO
        tmp1 = tmp1/tmp2
 
      END SUBROUTINE get_transport_coeff
  
      !----------------------------------------------------!
      ! This subroutine computes the species mass diffusion flux
      SUBROUTINE get_species_DiffFlux (T, nb, xi, diff_driv, Ji)  

        USE mod_nitrogen_NASA_initialize_CFD,       ONLY: nb_dim, upi, urg, una, mi, Tmin_CI, dT_CI, T_table, & 
                                                        & Omega_11_table, fac_Diff, sqrt_mij

        INTEGER :: i, j, ij, ii, jj, d
        INTEGER :: left, right
        REAL(KIND=8) :: fac, tmp1, tmp2, tmp3, tmp4
        REAL(KIND=8) :: xprod, diag, rho
        REAL(KIND=8) :: a, b, m, dT_norm
        REAL(KIND=8) :: m_i, m_j, det_V
        REAL(KIND=8), DIMENSION(nb_ns) :: rhoi
        REAL(KIND=8), DIMENSION(nb_inter) :: Omega_11, Gij_V, Dij
        REAL(KIND=8), DIMENSION(nb_ground) :: xig, rhoig
        REAL(KIND=8), DIMENSION(nb_ground*nb_dim) :: diff_drivg, Vdig

        REAL(KIND=8), INTENT(IN) :: T, nb
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, diff_driv
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: Ji

        ! Densities and molar fractions of N and N2
        fac = nb/una
        xig(1)   = xi(1)
        rhoi(1)  = xi(1)*fac*mi(1)
        rhoig(1) = rhoi(1) 
        tmp1 = 0.d0
        fac  = fac*mi(2)
        DO i = 2,nb_ns 
           tmp2 = xi(i)
           tmp1 = tmp1 + tmp2
           rhoi(i) = tmp2*fac 
        ENDDO

        xig(2)   = tmp1
        rhoig(2) = tmp1*fac

        ! Mixture density
        rho = rhoig(1) + rhoig(2)

        ! Diffusion driving forces for N and N2
        DO d = 1,nb_dim 

           left  = (d - 1)*nb_ground + 1 
           right = left + nb_ground - 1

           diff_drivg(left) = diff_driv((d - 1)*nb_ns + 1)
           tmp1 = 0.d0
           DO i = 2,nb_ns 
              tmp1 = tmp1 +  diff_driv((d - 1)*nb_ns + i)
           ENDDO 
           diff_drivg(right) = tmp1
           
        ENDDO

        ! Computate the Omega(1,1) collision integrals (look-up tables are used)
        ! Table search
        left  = INT((T - Tmin_CI)/dT_CI) + 1 
        right = left + 1 

        ! Normalized temperature difference
        dT_norm = (T  - T_table(left))/dT_CI

        ! Loop over interactions
        DO i = 1,nb_ground

           DO j = i,nb_ground 

              ! Index 
              ij = ((i - 1)*(2*nb_ground - i) + 2*j)/2 

              a = Omega_11_table(left,ij)
              b = Omega_11_table(right,ij) 
              m = b - a
              Omega_11(ij) = a + m*dT_norm

           ENDDO

         ENDDO 

         ! Assembly the linear transport system for diffusion velocities (Stefan-Maxwell equations)
         ! Initialization
         Dij   = 0.d0
         Gij_V = 0.d0

         ! Binary diffusion coefficients 
         fac = DSQRT(T)/(nb*fac_Diff)
         DO i = 1,nb_ground

            DO j = i,nb_ground 

               ! Indices
               ij = ((i - 1)*(2*nb_ground - i) + 2*j)/2

               ! Binary diffusion coefficient
               Dij(ij) = fac/(Omega_11(ij)*sqrt_mij(ij))

            ENDDO

         ENDDO

         DO i = 1,nb_ground

            DO j = i + 1,nb_ground 

               ! Indices
               ij   = ((i - 1)*(2*nb_ground - i) + 2*j)/2
               ii   = ((i - 1)*(2*nb_ground - i) + 2*i)/2
               jj   = ((j - 1)*(2*nb_ground - j) + 2*j)/2

               ! Common factor
               fac = xig(i)*xig(j)/Dij(ij)
                 
               ! Stefan-Maxwell matrix entries
               Gij_V(ij) = - fac
               Gij_V(ii) = Gij_V(ii) + fac 
               Gij_V(jj) = Gij_V(jj) + fac 
                 
            ENDDO

        ENDDO
          
        ! Enforce mass conservation
        diag = 1.d0/((rho**2)*MAXVAL(Dij))
        
        ! Incorporate the mass conservation constraint in the Stefan-Maxwell equations
        DO i = 1,nb_ground
           DO j = i,nb_ground 
              ij = ((i - 1)*(2*nb_ground - i) + 2*j)/2
              Gij_V(ij) = Gij_V(ij) + rhoig(i)*rhoig(j)*diag
           ENDDO
        ENDDO

        ! Solve the linear transport system for species diffusion velocities
        tmp1 = rhoig(1)/(Gij_V(2)*rhoig(1)/rhoig(2) - Gij_V(1))
        tmp2 = 1.d0/rhoig(2)
        DO d = 1,nb_dim

           left  = (d - 1)*nb_ground + 1 
            
           tmp3 = diff_drivg(left)*tmp1
           Ji((d - 1)*nb_ns + 1) = tmp3
           tmp4 = tmp2*tmp3
           DO i = 2,nb_ns
              Ji((d - 1)*nb_ns + i) = - rhoi(i)*tmp4
           ENDDO

        ENDDO

      END SUBROUTINE get_species_DiffFlux 

      !----------------------------------------------------!
      ! This subroutine solves a generic linear tranport system for the N-N2 system. 
      SUBROUTINE solve_transp_syst(xi, Gij, yi)

        REAL(KIND=8) :: det

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, Gij
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi

        det   = 1.d0/(Gij(1)*Gij(3) - Gij(2)*Gij(2))
        yi(1) = (xi(1)*Gij(3) - xi(2)*Gij(2))*det
        yi(2) = (xi(2)*Gij(1) - xi(1)*Gij(2))*det

      END SUBROUTINE solve_transp_syst

      !----------------------------------------------------!
      ! This subroutine computes the component of thermal conductivity associated to internal  
      ! degrees of freedom by means of Eucken's correction in case of use of the VC model
      SUBROUTINE lambda_int_VC_rr(nb, T, xi, xig, Dij, lambda_int)

        USE mod_nitrogen_NASA_initialize_CFD,         ONLY: una, Rn2, theta_vib, mm_N2, ukb, EkJ

        INTEGER :: i, j, ij
        REAL(KIND=8) :: den1, den2, Tratio, dexp0
        REAL(KIND=8) :: cv_rot, cv_vib, Tv

        REAL(KIND=8), INTENT(IN) :: nb, T
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, xig, Dij
        REAL(KIND=8), INTENT(OUT) :: lambda_int

        i = 2
        den1 = 0.d0
        DO j = 1,i - 1
           ij = ((j - 1)*(2*nb_ground - j) + 2*i)/2  
           den1 = den1 + xi(j)/Dij(ij)
        ENDDO
        
        DO j = i,nb_ns 
           ij = ((i - 1)*(2*nb_ground - i) + 2*nb_ground)/2
           den1 = den1 + xi(j)/Dij(ij)
        ENDDO
        den1 = 1.d0/den1

        ! Rotational and vibrational specific heats
        cv_rot = Rn2

        ! Vibrational specific. This quantity is estimated by extracting from the 
        ! level population the vibrational temperature and using the harmonic 
        ! oscillator model in order to estimate the specific heat
        Tv = (EkJ(2) - EkJ(1))/(ukb*DLOG(ABS(xi(2)/xi(3))))
        Tratio = theta_vib/Tv
        dexp0  = DEXP(Tratio)
        den2   = dexp0 - 1.d0 
        cv_vib = Rn2*dexp0*(Tratio/den2)**2 

        ! Internal component of thermal conductivity
        lambda_int = (den1*nb*xig(2)*mm_N2/una)*(cv_rot + cv_vib)
        
      END SUBROUTINE lambda_int_VC_rr

      !----------------------------------------------------!
      ! This subroutine computes the component of thermal conductivity associated to internal  
      ! degrees of freedom by means of Eucken's correction in case of use of the VC model
      SUBROUTINE lambda_int_VC(nb, T, xi, xig, Dij, lambda_int)

        REAL(KIND=8), INTENT(IN) :: nb, T
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, xig, Dij
        REAL(KIND=8), INTENT(OUT) :: lambda_int

        lambda_int = 0.d0

        PRINT*
        WRITE(*,10)'In lambda_int_VC, not implemented yet!'
        PRINT*
        STOP

10    FORMAT(A)
        
      END SUBROUTINE lambda_int_VC

      !----------------------------------------------------!
      ! This subroutine computes the component of thermal conductivity associated to internal  
      ! degrees of freedom by means of Eucken's correction in case of use of the BRVC model
      SUBROUTINE lambda_int_BRVC(nb, T, xi, xig, Dij, lambda_int)

        USE mod_nitrogen_NASA_initialize_CFD,         ONLY: nb_bins, una, mm_N2, ukb, EkJ, T_min, T_max, inv_step, & 
                                                          & T_store, qint, cvint

        INTEGER :: i, j, ij
        INTEGER :: left, right
        REAL(KIND=8) :: tmp1, tmp2, den
        REAL(KIND=8) :: Qk0, Qk1
        REAL(KIND=8) :: cv_int, Tint, Tclip

        REAL(KIND=8), INTENT(IN) :: nb, T
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, xig, Dij
        REAL(KIND=8), INTENT(OUT) :: lambda_int

        i = 2
        den = 0.d0
        DO j = 1,i - 1
           ij = ((j - 1)*(2*nb_ground - j) + 2*i)/2  
           den = den + xi(j)/Dij(ij)
        ENDDO
        
        DO j = i,nb_ns 
           ij = ((i - 1)*(2*nb_ground - i) + 2*nb_ground)/2
           den = den + xi(j)/Dij(ij)
        ENDDO
        den = 1.d0/den

        ! Internal specific heat. This quantity is computed by assuming that the population
        ! distribution of internal energy levels is Boltzmann at its own internal temperatutre. 
        ! The former is computed from the slope of the population distribution 
        ! between the ground state and the first excited level.
        ! Value search
        left = INT((T - T_min)*inv_step) + 1
        tmp1 = (T - T_store(left))*inv_step 
        
        left  = (left - 1)*nb_bins 
        right = left + nb_bins

        ! Bin internal partition functions
        tmp2 = qint(left + 1)
        Qk0  = tmp2 + tmp1*(qint(right + 1) - tmp2)
        tmp2 = qint(left + 2)
        Qk1  = tmp2 + tmp1*(qint(right + 2) - tmp2)

        ! Internal temperature 
        Tint = (EkJ(2) - EkJ(1))/(ukb*DLOG(ABS(Qk1*xi(2)/(Qk0*xi(3))))) 

        ! Out of bound check 
        Tclip = MAX(Tint,T_min)
        Tclip = MIN(Tclip,T_max)

        ! Value search 
        left = INT((Tclip - T_min)*inv_step) + 1
        tmp1 = (Tclip - T_store(left))*inv_step 
        
        left  = (left - 1)*nb_bins 
        right = left + nb_bins

        ! Internal specific heat
        cv_int = 0.d0
        DO i = 1,nb_bins 
           tmp2 = cvint(left + i)
           cv_int = cv_int + xi(i + 1)*(tmp2 + tmp1*(cvint(right + i) - tmp2))
        ENDDO
        cv_int = cv_int/xig(2)

        ! Internal component of thermal conductivity
        lambda_int = (den*nb*xig(2)*mm_N2/una)*cv_int

      END SUBROUTINE lambda_int_BRVC

      !----------------------------------------------------!
      ! This subroutine computes the component of thermal conductivity associated to internal  
      ! degrees of freedom by means of Eucken's correction in case of use of the RVC model
      SUBROUTINE lambda_int_RVC(nb, T, xi, xig, Dij, lambda_int)

        REAL(KIND=8), INTENT(IN) :: nb, T
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, xig, Dij
        REAL(KIND=8), INTENT(OUT) :: lambda_int

        lambda_int = 0.d0

        PRINT*
        WRITE(*,10)'In lambda_int_RVC, not implemented yet!'
        PRINT*
        STOP

10    FORMAT(A)
        
      END SUBROUTINE lambda_int_RVC

      !----------------------------------------------------!
      ! This subroutine computes the component of thermal conductivity associated to internal  
      ! degrees of freedom by means of Eucken's correction in case of use of the MT_TTint model
      SUBROUTINE lambda_int_MT_TTint(nb, T, xi, xig, Dij, lambda_int)

        REAL(KIND=8), INTENT(IN) :: nb, T
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, xig, Dij
        REAL(KIND=8), INTENT(OUT) :: lambda_int

        lambda_int = 0.d0

        PRINT*
        WRITE(*,10)'In lambda_int_MT_TTint, not implemented yet!'
        PRINT*
        STOP

10    FORMAT(A)
        
      END SUBROUTINE lambda_int_MT_TTint

      !----------------------------------------------------!
      ! This subroutine computes the collision integrals for N-N, N-N2 and N2-N2 interactions based on look-up tables.
      SUBROUTINE collision_table (T, Omega_11, Omega_12, Omega_22, Astar, Bstar, Cstar)
 
        USE mod_nitrogen_NASA_initialize_CFD,       ONLY: Tmin_CI, dT_CI, T_table, Omega_11_table, Omega_12_table, & 
                                                        & Omega_22_table, Astar_table, Bstar_table, Cstar_table

        INTEGER :: i, j, ij
        INTEGER :: left, right
        REAL(KIND=8) :: a, b, m
        REAL(KIND=8) :: dT_norm

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: Omega_11, Omega_12, Omega_22, Astar, Bstar, Cstar

        ! Table search 
        left  = INT((T - Tmin_CI)/dT_CI) + 1 
        right = left + 1 

        ! Normalized temperature difference
        dT_norm = (T  - T_table(left))/dT_CI

        ! Loop over interactions
        DO i = 1,nb_ground

           DO j = i,nb_ground

              ! Index 
              ij = ((i - 1)*(2*nb_ground - i) + 2*j)/2 

              ! Omega(1,1)
              a = Omega_11_table(left,ij)
              b = Omega_11_table(right,ij) 
              m = b - a
              Omega_11(ij) = a + m*dT_norm

              ! Omega(1,2)
              a = Omega_12_table(left,ij)
              b = Omega_12_table(right,ij) 
              m = b - a
              Omega_12(ij) = a + m*dT_norm

              ! Omega(2,2)
              a = Omega_22_table(left,ij)
              b = Omega_22_table(right,ij) 
              m = b - a
              Omega_22(ij) = a + m*dT_norm

              ! Astar
              a = Astar_table(left,ij)
              b = Astar_table(right,ij) 
              m = b - a
              Astar(ij) = a + m*dT_norm

              ! Bstar
              a = Bstar_table(left,ij)
              b = Bstar_table(right,ij) 
              m = b - a
              Bstar(ij) = a + m*dT_norm

              ! Cstar
              a = Cstar_table(left,ij)
              b = Cstar_table(right,ij) 
              m = b - a
              Cstar(ij) = a + m*dT_norm

           ENDDO

        ENDDO

      END SUBROUTINE collision_table

  END MODULE mod_nitrogen_NASA_CFD_transport
!------------------------------------------------------------------------------!
