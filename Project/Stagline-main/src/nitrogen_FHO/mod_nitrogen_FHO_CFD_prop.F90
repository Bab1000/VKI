!------------------------------------------------------------------------------!
! This module provides subroutines/function for the thermodynamic properties of FHO database for the N-N2 system.
  MODULE mod_nitrogen_FHO_CFD_prop

    IMPLICIT NONE

    ! Subroutine for thermodynamics
    CONTAINS 

      !------------------------------------------------------!
      ! This subroutine computes the translation partition function.
      SUBROUTINE Q_trans (T, mm, Q)

        USE mod_nitrogen_FHO_initialize_CFD,     ONLY: fac_Q

        REAL(KIND=8), INTENT(IN)  :: T, mm
        REAL(KIND=8), INTENT(OUT) :: Q

        Q = fac_Q*(mm*T)**1.5
     
      END SUBROUTINE Q_trans

      !------------------------------------------------------!
      ! This subroutine computes the translation partition function and its
      ! derivative with respect to temperature
      SUBROUTINE Q_trans_der (T, mm, Q, dQ_dT)

        USE mod_nitrogen_FHO_initialize_CFD,     ONLY: fac_Q

        REAL(KIND=8), INTENT(IN)  :: T, mm
        REAL(KIND=8), INTENT(OUT) :: Q, dQ_dT

        Q = fac_Q*(mm*T)**1.5
        dQ_dT = 1.5d0*Q/T

      END SUBROUTINE Q_trans_der 

      !----------------------------------------------------!
      ! This subroutine computes the rotational partition function
      SUBROUTINE Q_rot (T, Q)

        USE mod_nitrogen_FHO_initialize_CFD,     ONLY: theta_rot, ov_sigma

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), INTENT(OUT) :: Q

        Q = T*ov_sigma/theta_rot

      END SUBROUTINE Q_rot

      !----------------------------------------------------!
      ! This subroutine computes the rotational partition function and its
      ! derivative with respect to temperature
      SUBROUTINE Q_rot_der (T, Q, dQ_dT)

        USE mod_nitrogen_FHO_initialize_CFD,     ONLY: theta_rot, ov_sigma

        REAL(KIND=8) :: fac

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), INTENT(OUT) :: Q, dQ_dT

        fac   = ov_sigma/theta_rot
        Q     = T*fac
        dQ_dT = fac
        
      END SUBROUTINE Q_rot_der  

      !----------------------------------------------------!
      ! This subroutine computes the temperature given the internal energy density
      SUBROUTINE compute_T (rhoi, rho_eint, T)

        USE mod_nitrogen_FHO_initialize_CFD,     ONLY: nb_bins, hf_N, cv_tr_N, cv_tr_N2, R_N2, ek

        INTEGER :: i 
        REAL(KIND=8) :: cv_N, cv_N2, rhoN, rhoN2
        REAL(KIND=8) :: tmp, esum, lhs, rhs

        REAL(KIND=8), INTENT(IN) :: rho_eint
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), INTENT(OUT) :: T

        ! Useful quantities
        rhoN  = rhoi(1)

        rhoN2 = 0.d0
        esum  = 0.d0
        DO i = 1,nb_bins 
           tmp   = rhoi(i + 1)
           rhoN2 = rhoN2 + tmp
           esum  = esum  + tmp*ek(i) 
        ENDDO

        lhs = rho_eint - (esum + rhoN*hf_N)
        rhs = cv_tr_N*rhoN + (cv_tr_N2 + R_N2)*rhoN2         

        ! Temperature value
        T = lhs/rhs

      END SUBROUTINE compute_T

      !----------------------------------------------------!
      ! This subroutine computes the species internal energy per unit mass
      SUBROUTINE get_species_energy (T, e)

        USE mod_nitrogen_FHO_initialize_CFD,     ONLY: nb_bins, hf_N, cv_tr_N, cv_tr_N2, R_N2, ek

        INTEGER :: i
        REAL(KIND=8) :: etr

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), DIMENSION(:) :: e

        e(1) = cv_tr_N*T + hf_N 

        etr = (cv_tr_N2 + R_N2)*T
        DO i = 1,nb_bins 
           e(i + 1) = etr + ek(i)
        ENDDO 

      END SUBROUTINE get_species_energy

      !----------------------------------------------------!
      ! This subroutine computes the species internal energy and constant volume
      ! specific heat per unit mass
      SUBROUTINE get_species_energy_cv (T, e, cv)

        USE mod_nitrogen_FHO_initialize_CFD,     ONLY: nb_bins, hf_N, cv_tr_N, cv_tr_N2, R_N2, ek

        INTEGER :: i
        REAL(KIND=8) :: etr, cvtr

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), DIMENSION(:) :: e, cv

        e(1)  = cv_tr_N*T + hf_N 
        cv(1) = cv_tr_N

        cvtr = cv_tr_N2 + R_N2
        etr  = cvtr*T
        DO i = 1,nb_bins 
           e(i + 1)  = etr + ek(i)
           cv(i + 1) = cvtr
        ENDDO 

      END SUBROUTINE get_species_energy_cv

      !----------------------------------------------------!
      ! This subroutine computes the species constant volume specific heat per unit mass
      SUBROUTINE get_species_cv (T, cv)

        USE mod_nitrogen_FHO_initialize_CFD,     ONLY: nb_bins, cv_tr_N, cv_tr_N2, R_N2

        INTEGER :: i
        REAL(KIND=8) :: cvtr

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv

        cv(1) = cv_tr_N

        cvtr = cv_tr_N2 + R_N2
        DO i = 1,nb_bins 
           cv(i + 1) = cvtr
        ENDDO 

      END SUBROUTINE get_species_cv

      !----------------------------------------------------!
      ! This subroutine computes the species constant pressure specific heat per unit mass
      SUBROUTINE get_species_cp (T, cp)

        USE mod_nitrogen_FHO_initialize_CFD,     ONLY: nb_bins, R_N, R_N2

        INTEGER :: i

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cp

        CALL get_species_cv (T, cp)

        cp(1) = cp(1) + R_N
        DO i = 1,nb_bins
           cp(i + 1) = cp(i + 1) + R_N2
        ENDDO

      END SUBROUTINE get_species_cp

      !----------------------------------------------------!
      ! This subroutine computes post-shock nonequilibrium conditions
      SUBROUTINE post_shock_neq (p1, u1, T1, p2, u2, T2, yi, x_N)

        USE mod_nitrogen_FHO_initialize_CFD,      ONLY: nb_ns, mm_N, mm_N2, urg, cv_tr_N, cv_tr_N2
        USE mod_nitrogen_FHO_CFD_eq,              ONLY: eq_composition_vib

        INTEGER :: i
        REAL(KIND=8) :: tmp1, tmp2, tmp3
        REAL(KIND=8) :: cp, cv, g, mass, R
        REAL(KIND=8) :: c1, gp1, rho1, M1, M1s
        REAL(KIND=8), DIMENSION(nb_ns) :: rhoi, Ri, mi

        REAL(KIND=8), INTENT(IN) :: p1, u1, T1
        REAL(KIND=8), OPTIONAL, INTENT(IN) :: x_N
        REAL(KIND=8), INTENT(OUT) :: p2, u2, T2
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi

        ! Useful data
        Ri(1) = urg/mm_N
        mi(1) = mm_N 

        DO i = 2,nb_ns 
           Ri(i) = urg/mm_N2
           mi(i) = mm_N2
        ENDDO

        CALL eq_composition_vib (p1, T1, rhoi, x_N)

        rho1 = 0.d0
        DO i = 1,nb_ns 
           rho1 = rho1 + rhoi(i)
        ENDDO

        ! Species mass fractions
        DO i = 1,nb_ns 
           yi(i) = rhoi(i)/rho1
        ENDDO

        ! Mixture molar mass
        mass = 0.d0
        DO i = 1,nb_ns 
           mass = mass + yi(i)/mi(i)
        ENDDO
        mass = 1.d0/mass
        R    = urg/mass

        ! Frozen specific heats (constant pressure and volume)
        tmp1 = yi(1)
        cv  = tmp1*cv_tr_N 
        cp  = tmp1*(cv_tr_N + Ri(1))

        DO i = 2,nb_ns 
           tmp1 = yi(i)
           tmp2 = Ri(i)
           tmp3 = cv_tr_N2 + tmp2
           cv   = cv + tmp1*tmp3
           cp   = cp + tmp1*(tmp2 + tmp3)  
        ENDDO
 
        ! Specific heat ratio
        g   = cp/cv
        gp1 = g + 1.d0
      
        c1  = SQRT(g*R*T1)
        M1  = u1/c1
        M1s = M1**2

        ! Pressure, velocity and temperature after the shock
        p2 = p1*(2.d0*g*M1s - g + 1.d0)/gp1
        u2 = u1 - c1*2.d0/gp1*(M1 - 1.d0/M1)
        T2 = T1*(2.d0*g*M1s - g + 1.d0)*(g - 1.d0 + 2.d0/M1s)/(gp1**2) 

      END SUBROUTINE post_shock_neq

      !----------------------------------------------------!
      ! This subroutine computes the vibrational temperature for the energy level population. 
      ! For sake of robustness the first few iterations are performed out by means of bi-section method.
      ! Then, after the residual has dropped some order of magnitude, the faster Newton-Raphson method is used 
      SUBROUTINE compute_Tint (ni, Tint)

        USE mod_nitrogen_FHO_initialize_CFD,     ONLY: nb_bins, EkJ, ukb

        INTEGER :: i
        REAL(KIND=8), PARAMETER :: tol1 = 1.d-1, tol2 = 1.d-8
        REAL(KIND=8) :: ov_kb_T, Tint_old
        REAL(KIND=8) :: sum1l, sum2l, sum1r, sum2r, ov_kb_Tl, ov_kb_Tr, Tl, Tr, fl, fr, & 
                        tmp2l, tmp3l, tmp2r, tmp3r
        REAL(KIND=8) :: sum1, sum2, sum3, tmp1, tmp2, tmp3
        REAL(KIND=8) :: res, rhs, f, fp

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: ni 
        REAL(KIND=8), INTENT(OUT) :: Tint

        ! Computation of rhs
        rhs  = 0.d0
        sum1 = 0.d0
        sum2 = 0.d0
        DO i = 1,nb_bins 
           tmp1 = ni(i)
           sum1 = sum1 + tmp1
           sum2 = sum2 + tmp1*EkJ(i)
        ENDDO
        rhs = sum2/sum1

        ! Bi-section method (left and right bounds)
        Tl   = 300.d0
        Tr   = 30000.d0
        Tint = 0.5d0*(Tl + Tr)

        res = 1.d0
        DO WHILE (res.GT.tol1)
 
           sum1  = 0.d0
           sum2  = 0.d0

           sum1l = 0.d0
           sum2l = 0.d0

           sum1r = 0.d0
           sum2r = 0.d0
           
           ov_kb_Tl = 1.d0/(ukb*Tl)
           ov_kb_Tr = 1.d0/(ukb*Tr)
           ov_kb_T  = 1.d0/(ukb*Tint)

           DO i = 1,nb_bins 

              tmp1 = EkJ(i)

              tmp2  = DEXP(-tmp1*ov_kb_T)
              tmp3  = tmp2*tmp1

              tmp2l = DEXP(-tmp1*ov_kb_Tl)
              tmp3l = tmp2l*tmp1

              tmp2r = DEXP(-tmp1*ov_kb_Tr)
              tmp3r = tmp2r*tmp1

              ! Values at Tint
              sum1  = sum1 + tmp2
              sum2  = sum2 + tmp3

              ! Values at Tl
              sum1l = sum1l + tmp2l 
              sum2l = sum2l + tmp3l

              ! Values at Tr
              sum1r = sum1r + tmp2r 
              sum2r = sum2r + tmp3r

           ENDDO

           f  = sum2/sum1   - rhs 
           fl = sum2l/sum1l - rhs
           fr = sum2r/sum1r - rhs    

           ! Solution update
           Tint_old = Tint

           IF ((f*fr).LT.0.d0) THEN 

              Tl = Tint 
              Tint = 0.5d0*(Tl + Tr)              

           ELSEIF ((f*fl).LT.0.d0) THEN

              Tr = Tint 
              Tint = 0.5d0*(Tl + Tr)   

           ENDIF

           ! Residual 
           res = ABS(Tint - Tint_old)/Tint

        ENDDO

        ! Newton-Raphson method
        res  = 1.d0
        DO WHILE (res.GT.tol2)

           sum1 = 0.d0
           sum2 = 0.d0
           sum3 = 0.d0

           ov_kb_T = 1.d0/(ukb*Tint)
           DO i = 1,nb_bins 

              tmp1 = EkJ(i)
              tmp2 = DEXP(-tmp1*ov_kb_T)
              tmp3 = tmp2*tmp1

              sum1 = sum1 + tmp2 
              sum2 = sum2 + tmp3
              sum3 = sum3 + tmp3*tmp1 

           ENDDO

           Tint_old = Tint

           tmp1 = sum2/sum1
           f    = tmp1 - rhs
           fp   = (sum3/sum1 - tmp1**2)*ov_kb_T/Tint
           
           ! Solution update
           Tint = Tint - f/fp

           ! Residual 
           res = ABS(Tint - Tint_old)/Tint

        ENDDO

      END SUBROUTINE compute_Tint 

  END MODULE mod_nitrogen_FHO_CFD_prop
!------------------------------------------------------------------------------!
