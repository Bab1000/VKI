!------------------------------------------------------------------------------!
! This modules provides subroutine for computing transport properties and fluxes for the N-N2 system 
! when using the Park multi-temperature model. It also provides subroutines to be interfaced with 
! CFD codes for computing mass diffusion flux of species and diffusive components of heat flux vector
! (in this way model details are hidden to the CFD code that can be therefore written in a general and 
! flexible manner). 
  MODULE mod_nitrogen_Park_CFD_transport

    USE mod_nitrogen_Park_initialize_CFD,         ONLY: nb_ns, nb_inter

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
      SUBROUTINE get_transport_coeff (nb, xi, yi, temp, mu, kappa, lambda_tr, lambda_rot, lambda_vib, Di, chi)
     
        USE mod_nitrogen_Park_initialize_CFD,         ONLY: pos_N, pos_N2, posTr, posTv, una, ukb, upi, urg, Tmin, fac_mu, & 
                                                          & fac_Diff, mi, Ri, sqrt_mi, sqrt_mij

        INTEGER :: i, j, ii, jj, ij
        REAL(KIND=8) :: tmp1, tmp2, fac1, fac2 
        REAL(KIND=8) :: x_i, x_prod, m_i, m_j, miij, mjij,  m_sum
        REAL(KIND=8) :: miij_miij, miij_mjij, mjij_mjij
        REAL(KIND=8) :: Aij, Bij, w
        REAL(KIND=8) :: T, Tr, Tv, sqr_T
        REAL(KIND=8), DIMENSION(nb_ns) :: mui, lambdai, xis, coeff_mu, coeff_lam
        REAL(KIND=8), DIMENSION(nb_inter) :: xi_xj_ov_nDij,  xi_xj_ov_Dij, Dij
        REAL(KIND=8), DIMENSION(nb_inter) :: Gij_mu, Gij_lam
        REAL(KIND=8), DIMENSION(nb_inter) :: Omega_11, Omega_12, Omega_22, Astar, Bstar, Cstar

        REAL(KIND=8), INTENT(IN) :: nb
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, yi, temp 
        REAL(KIND=8), INTENT(OUT) :: kappa, mu, lambda_tr, lambda_rot, lambda_vib
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: chi, Di

        ! Temperatures (a fix is applied in order to avoid numerical problems)
        T  = MAX(temp(1),Tmin)
        Tr = MAX(temp(posTr),Tmin) 
        Tv = MAX(temp(posTv),Tmin)
        sqr_T = DSQRT(T)
        
        DO i = 1,nb_ns 
           xis(i) = xi(i)**2
        ENDDO

        ! No thermal diffusion 
        chi = 0.d0

        ! No bulk viscosity 
        kappa = 0.d0

        ! Compute collision integrals  
        CALL collision_table (T, Omega_11, Omega_12, Omega_22, Astar, Bstar, Cstar) 

        ! Compute the binary diffusion coefficients
        fac1 = fac_Diff/sqr_T
        DO i = 1,nb_ns 

           DO j = i,nb_ns 

              x_prod = xi(i)*xi(j)

              ! Index
              ij = ((i - 1)*(2*nb_ns - i) + 2*j)/2

              ! Common factor
              tmp1 = fac1*Omega_11(ij)*x_prod*sqrt_mij(ij) 

              xi_xj_ov_nDij(ij) = tmp1
              Dij(ij) = x_prod/(nb*tmp1)

           ENDDO

        ENDDO
       
        ! Compute single species viscosity and translational component of thermal conductivity
        fac1 = fac_mu*sqr_T
        fac2 = 15.d0/4.d0
        DO i = 1,nb_ns 

           ii = ((i - 1)*(2*nb_ns - i) + 2*i)/2 
           tmp1   = fac1*sqrt_mi(i)/Omega_22(ii)
           mui(i) = tmp1 
           lambdai(i) = fac2*tmp1*Ri(i)

        ENDDO 
 
        ! Compute rotational and vibrational components of thermal conductivity
        CALL get_lambda_int(nb, Tr, Tv, xi, Dij, lambda_rot, lambda_vib)

        ! N2 mixture
        IF (nb_ns.EQ.1) THEN

           Di = 0.d0
           mu = mui(pos_N2)
           lambda_tr = lambdai(pos_N2)

        ! N-N2 mixture
        ELSE

           ! Assembly linear transport systems for viscosity and translational component of thermal conductivity
           ! Initialization
           Gij_mu  = 0.d0
           Gij_lam = 0.d0

           DO i = 1,nb_ns 
              ii = ((i - 1)*(2*nb_ns - i) + 2*i)/2
              Gij_mu(ii)  = xis(i)/mui(i)
              Gij_lam(ii) = xis(i)/lambdai(i)
           ENDDO

           DO i = 1,nb_ns 

              m_i = mi(i)

              DO j = i + 1,nb_ns 

                 m_j = mi(j)

                 ! Indices
                 ij  = ((i - 1)*(2*nb_ns - i) + 2*j)/2 
                 ii  = ((i - 1)*(2*nb_ns - i) + 2*i)/2 
                 jj  = ((j - 1)*(2*nb_ns - j) + 2*j)/2

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
           CALL solve_transp_syst(xi, Gij_mu, coeff_mu)
           CALL solve_transp_syst(xi, Gij_lam, coeff_lam)

           ! Compute dynamic viscosity and translational component of thermal conductivity from the solution of the 
           ! linear transport systems
           mu = 0.d0
           lambda_tr = 0.d0
           DO i = 1,nb_ns 
              x_i = xi(i)
              mu        = mu + x_i*coeff_mu(i)
              lambda_tr = lambda_tr + x_i*coeff_lam(i)
           ENDDO

           ! Average diffusion coefficients
           DO i = 1,nb_ns 

              ! Initialization
              tmp1 = 1.d0 - yi(i)
              tmp2 = 0.d0

              DO j = 1,i - 1
                 ij = ((i - 1)*(2*nb_ns - i) + 2*j)/2
                 tmp2 = tmp2 + xi(j)/Dij(ij)
              ENDDO
              
              DO j = i + 1,nb_ns
                 ij = ((i - 1)*(2*nb_ns - i) + 2*j)/2
                 tmp2 = tmp2 + xi(j)/Dij(ij)
              ENDDO
              
              Di(i) = tmp1/tmp2
             
           ENDDO

        ENDIF

      END SUBROUTINE get_transport_coeff

      !----------------------------------------------------!
      ! This subroutine computes the species diffusion flux
      SUBROUTINE get_species_DiffFlux (T, nb, xi, diff_driv, Ji)  

        USE mod_nitrogen_Park_initialize_CFD,       ONLY: nb_dim, pos_N, pos_N2, upi, urg, una, mi, Tmin, dT, T_table, & 
                                                        & Omega_11_table, fac_Diff, sqrt_mij

        INTEGER :: i, j, ij, ii, jj, d
        INTEGER :: left, right
        REAL(KIND=8) :: fac
        REAL(KIND=8) :: xprod, diag, rho
        REAL(KIND=8) :: a, b, m, dT_norm
        REAL(KIND=8) :: det_V, mm
        REAL(KIND=8), DIMENSION(nb_ns) :: rhoi
        REAL(KIND=8), DIMENSION(nb_inter) :: Omega_11, Gij_V, Dij

        REAL(KIND=8), INTENT(IN) :: T, nb
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, diff_driv
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: Ji

        ! N2 mixture
        IF (nb_ns.EQ.1) THEN

           ! No species diffusion 
           Ji = 0.d0

        ! N-N2 mixture
        ELSE 

           ! Species densities
           rho = 0.d0 
           fac = nb/una
           DO i = 1,nb_ns 
              rhoi(i) = xi(i)*mi(i)*fac
              rho     = rho + rhoi(i)
           ENDDO
           
           ! Computate the Omega(1,1) collision integrals (look-up tables are used)
           ! Table search
           left  = INT((T - Tmin)/dT) + 1 
           right = left + 1 

           ! Normalized temperature difference
           dT_norm = (T  - T_table(left))/dT

           ! Loop over interactions
           DO i = 1,nb_ns 

              DO j = i,nb_ns 

                 ! Index 
                 ij = ((i - 1)*(2*nb_ns - i) + 2*j)/2 

                 ! Omega(1,1)
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
           DO i = 1,nb_ns 

              DO j = i,nb_ns 

                 ! Indices
                 ij = ((i - 1)*(2*nb_ns - i) + 2*j)/2
                 
                 ! Binary diffusion coefficient
                 Dij(ij) = fac/(Omega_11(ij)*sqrt_mij(ij))

              ENDDO

           ENDDO

           DO i = 1,nb_ns 

              DO j = i + 1,nb_ns 

                 ! Indices
                 ij   = ((i - 1)*(2*nb_ns - i) + 2*j)/2
                 ii   = ((i - 1)*(2*nb_ns - i) + 2*i)/2
                 jj   = ((j - 1)*(2*nb_ns - j) + 2*j)/2

                 ! Common factor
                 fac = xi(i)*xi(j)/Dij(ij)
                 
                 ! Stefan-Maxwell matrix entries
                 Gij_V(ij) = - fac
                 Gij_V(ii) = Gij_V(ii) + fac 
                 Gij_V(jj) = Gij_V(jj) + fac 
                 
              ENDDO

          ENDDO
          
          ! Enforce mass conservation
          diag = 1.d0/((rho**2)*MAXVAL(Dij))
          
          ! Incorporate the mass conservation constraint in the Stefan-Maxwell equations
          DO i = 1,nb_ns 
             DO j = i,nb_ns 
                ij = ((i - 1)*(2*nb_ns - i) + 2*j)/2
                Gij_V(ij) = Gij_V(ij) + rhoi(i)*rhoi(j)*diag
             ENDDO
          ENDDO

          ! Solve the linear transport system for species diffusion velocities
          DO d = 1,nb_dim

             left  = (d - 1)*nb_ns + 1 
             right = left + nb_ns - 1

             CALL solve_transp_syst(-diff_driv(left:right), Gij_V, Ji(left:right))

             ! Compute the species diffusion flux 
             DO i = 1,nb_ns 
                Ji(left + i - 1) = rhoi(i)*Ji(left + i - 1)
             ENDDO

          ENDDO

        ENDIF

      END SUBROUTINE get_species_DiffFlux 

      !----------------------------------------------------!
      ! This subroutine solves a generic linear tranport system for the N-N2 system. 
      SUBROUTINE solve_transp_syst(xi, Gij, yi)

        REAL(KIND=8) :: det

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, Gij
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi

        det   = 1.d0/(Gij(1)*Gij(3) - Gij(2)**2)
        yi(1) = (xi(1)*Gij(3) - xi(2)*Gij(2))*det
        yi(2) = (xi(2)*Gij(1) - xi(1)*Gij(2))*det

      END SUBROUTINE solve_transp_syst

      !----------------------------------------------------!
      ! This subroutine computes the component of thermal conductivity associated to rotational and vibrational 
      ! degrees of freedom by means of Eucken's correction
      SUBROUTINE get_lambda_int(nb, Tr, Tv, xi, Dij, lambda_rot, lambda_vib)

        USE mod_nitrogen_Park_initialize_CFD,         ONLY: pos_N2, una, mm_N2, Ri, theta_vib, cv_rot

        INTEGER :: i, j, ij
        REAL(KIND=8) :: cv_vib
        REAL(KIND=8) :: fac, Tratio, dexp0, den

        REAL(KIND=8), INTENT(IN) :: nb, Tr, Tv
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi, Dij
        REAL(KIND=8), INTENT(OUT) :: lambda_rot, lambda_vib

        Tratio = theta_vib/Tv
        dexp0  = DEXP(Tratio)
        den    = dexp0 - 1.d0 
        cv_vib = Ri(pos_N2)*dexp0*(Tratio/den)**2

        i   = pos_N2
        den = 0.d0
        DO j = 1,i - 1
           ij = ((j - 1)*(2*nb_ns - j) + 2*i)/2  
           den = den + xi(j)/Dij(ij)
        ENDDO
        
        DO j = i,nb_ns 
           ij = ((i - 1)*(2*nb_ns - i) + 2*j)/2
           den = den + xi(j)/Dij(ij)
        ENDDO
        den = 1.d0/den

        fac = den*nb*mm_N2*xi(pos_N2)/una
        lambda_rot = fac*cv_rot
        lambda_vib = fac*cv_vib

      END SUBROUTINE get_lambda_int

      !----------------------------------------------------!
      ! This subroutine computes the collision integrals for N-N, N-N2 and N2-N2 interactions based on look-up tables.
      SUBROUTINE collision_table (T, Omega_11, Omega_12, Omega_22, Astar, Bstar, Cstar)
 
        USE mod_nitrogen_Park_initialize_CFD,       ONLY: Tmin, dT, T_table, Omega_11_table, Omega_12_table, Omega_22_table, &
                                                        & Astar_table, Bstar_table, Cstar_table

        INTEGER :: i, j, ij
        INTEGER :: left, right
        REAL(KIND=8) :: a, b, m
        REAL(KIND=8) :: dT_norm

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: Omega_11, Omega_12, Omega_22, Astar, Bstar, Cstar

        ! Table search 
        left  = INT((T - Tmin)/dT) + 1 
        right = left + 1 

        ! Normalized temperature difference
        dT_norm = (T  - T_table(left))/dT

        ! Loop over interactions
        DO i = 1,nb_ns 

           DO j = i,nb_ns 

              ! Index 
              ij = ((i - 1)*(2*nb_ns - i) + 2*j)/2 

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

      !----------------------------------------------------!
      ! This subrotine computes the collision integrals at a given temperature for N-N, N2-N2 and N2-N2 interactions.
      SUBROUTINE collision (T, Omega_11, Omega_12, Omega_22, Astar, Bstar, Cstar)

        USE mod_nitrogen_Park_initialize_CFD,       ONLY: upi, q11, q12, q22, bs, cs

        INTEGER :: i, j, ij
        INTEGER :: pos1, pos2, pos3, pos4
        REAL(KIND=8) :: tmp1, tmp2
        REAL(KIND=8) :: lnt, lnt2, lnt3

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: Omega_11, Omega_12, Omega_22, Astar, Bstar, Cstar

        ! Useful data
        lnt  = DLOG(T)
        lnt2 = lnt*lnt
        lnt3 = lnt2*lnt

        ! Computation of collision integrals
        DO i = 1,nb_ns

           DO j = i,nb_ns  

              ! Index 
              ij = ((i - 1)*(2*nb_ns - i) + 2*j)/2 

              pos1 = 4*(ij - 1) + 1
              pos2 = 4*(ij - 1) + 2
              pos3 = 4*(ij - 1) + 3
              pos4 = 4*(ij - 1) + 4
              
              tmp1 = 1.d-20*DEXP(q11(pos1)*lnt3 + q11(pos2)*lnt2 + q11(pos3)*lnt + q11(pos4))   
              tmp2 = 1.d-20*DEXP(q22(pos1)*lnt3 + q22(pos2)*lnt2 + q22(pos3)*lnt + q22(pos4)) 

              Omega_11(ij) = tmp1
              Omega_22(ij) = tmp2
              Astar(ij)    = tmp2/tmp1

              Bstar(ij) = DEXP(bs(pos1)*lnt3 + bs(pos2)*lnt2 + bs(pos3)*lnt + bs(pos4)) 
              tmp2      = DEXP(cs(pos1)*lnt3 + cs(pos2)*lnt2 + cs(pos3)*lnt + cs(pos4)) 

              Cstar(ij)    = tmp2
              Omega_12(ij) = tmp2*tmp1

           ENDDO

        ENDDO 

      END SUBROUTINE collision

  END MODULE mod_nitrogen_Park_CFD_transport
!------------------------------------------------------------------------------!
