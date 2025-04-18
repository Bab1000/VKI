!------------------------------------------------------------------------------!
! This modules provides subroutines for the computation of equilibrium properties of the N-N2 system 
! when using the Park multi-temperature model.
  MODULE mod_Nitrogen_Park_CFD_eq

    USE mod_nitrogen_Park_initialize_CFD,         ONLY: mm_N, mm_N2

    IMPLICIT NONE

    ! Subroutines for equilibrium properties 
    CONTAINS 

      !----------------------------------------------------!
      ! This subriutine computes the eqilibrium composition in terms of molar
      ! fractions
      SUBROUTINE eq_composition (p, T, xi)

        USE mod_nitrogen_Park_initialize_CFD,    ONLY: nb_ns, pos_N, pos_N2, gn, ukb, ed_n2

        REAL(KIND=8) :: delta2, rhs, kbT
        REAL(KIND=8) :: Qn, Qn2, Qint 

        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: xi

        IF (nb_ns.GT.1) THEN 

           ! Translational and internal partition functions of N and N2 species
           CALL tr_part_function  (T, mm_N, Qn)
           CALL tr_part_function  (T, mm_N2, Qn2)
           CALL int_part_function (T, Qint)
         
           kbT = ukb*T  
           rhs = (1.d0/p)*kbT*DEXP(-ed_n2/kbT)*((gn*Qn)**2)/(Qn2*Qint)
          
           delta2 = rhs*(rhs + 4.d0) 
          
           ! Molar fractions of N and N2 species
           xi(pos_N)  = 0.5d0*(- rhs + DSQRT(delta2))
           xi(pos_N2) = 1.d0 - xi(pos_N)

        ELSE 

          xi = 1.d0

        ENDIF

      END SUBROUTINE eq_composition

      !----------------------------------------------------!
      ! This subroutine computes the internal partition function for the N2 molecule
      SUBROUTINE int_part_function (T, Qint)

        USE mod_nitrogen_Park_initialize_CFD,     ONLY: theta_rot, theta_vib 

        INTEGER :: i
        REAL(KIND=8) :: ov_kbT, dexp_1

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), INTENT(OUT) :: Qint 

        ! Internal partition function
        dexp_1 = DEXP(-theta_vib/T)
        Qint   = (T/(theta_rot*2.d0))*(1.d0/(1.d0 - dexp_1))

      END SUBROUTINE int_part_function

      !----------------------------------------------------!
      ! This subroutine computes the translational patition function
      SUBROUTINE tr_part_function (T, mm, Qtr)

        USE mod_nitrogen_Park_initialize_CFD,     ONLY: fac_Q

        REAL(KIND=8), INTENT(IN)  :: T, mm
        REAL(KIND=8), INTENT(OUT) :: Qtr

        Qtr = fac_Q*(mm*T)**1.5

      END SUBROUTINE tr_part_function  

      !------------------------------------------------------!
      ! This subroutine provides the internal energy per unit mass of the N-N2 system
      SUBROUTINE energy (T, e)

        USE mod_nitrogen_Park_initialize_CFD,    ONLY: pos_N, pos_N2, theta_vib, hf_n, Rn2, cv_tr, cv_rot 

        INTEGER :: k
        REAL(KIND=8) :: er, ev, dexp_Tv, erot

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: e

        er      = cv_rot*T
        dexp_Tv = DEXP(theta_vib/T)
        ev      = Rn2*theta_vib/(dexp_Tv - 1.d0)
      
        ! Internal energies
        e(pos_N)  = cv_tr(pos_N)*T + hf_n
        e(pos_N2) = cv_tr(pos_N2)*T + er + ev 

      END SUBROUTINE energy 

      !------------------------------------------------------!
      ! This subroutine provides the static enthapy per unit mass of the N-N2 system
      SUBROUTINE enthalpy (T, h)

        USE mod_nitrogen_Park_initialize_CFD,    ONLY: nb_ns, Ri 

        INTEGER :: is

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: h

        CALL energy(T, h)

        DO is = 1,nb_ns 
           h(is) = h(is) + Ri(is)*T
        ENDDO

      END SUBROUTINE enthalpy

      !----------------------------------------------------!
      ! This subroutine computes the entropy per unit mass 
      ! of N and N2 species in equilibrium conditions. The entropy of 
      ! mixing contribution is computed as well. 
      SUBROUTINE entropy (p, T, yi, xi, smix, s)

        USE mod_nitrogen_Park_initialize_CFD,    ONLY: nb_ns, pos_N, pos_N2, gn, ukb, hf_n, mm_N, mm_N2, &
                                                    &  Rn, Rn2, Ri 

        INTEGER :: is
        REAL(KIND=8) :: tmp1, tmp2
        REAL(KIND=8) :: QtN, QtN2, QintN2
        REAL(KIND=8), DIMENSION(nb_ns) :: h
 
        REAL(KIND=8), INTENT(IN)  :: T, p
        REAL(KIND=8), INTENT(IN), DIMENSION(:) :: xi, yi
        REAL(KIND=8), INTENT(OUT) :: smix
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: s

        ! Translational and internal partition functions of N and N2 species     
        CALL tr_part_function (T, mm_N, QtN)
        CALL tr_part_function (T, mm_N2, QtN2)
        CALL int_part_function (T, QintN2)

        ! N and N2 enthalpies per unit mass
        CALL enthalpy (T, h)  

        ! Nuclear and electronic spin contribution
        QtN = gN*QtN

        ! Subtracting the formation enthalpy
        h(pos_N) = h(pos_N) - hf_n

        tmp1 = 1.d0/T 
        tmp2 = ukb*T/p

        s(pos_N)  = h(pos_N)*tmp1 + Rn*DLOG(tmp2*QtN) 
        s(pos_N2) = h(pos_N2)*tmp1 + Rn2*DLOG(tmp2*QtN2*QintN2)

        ! Entropy of mixing
        smix = 0.d0
        DO is = 1,nb_ns 
           smix = smix - Ri(is)*yi(is)*DLOG(xi(is))
        ENDDO 

      END SUBROUTINE entropy

      !----------------------------------------------------!
      ! This subroutine computes the equilibrium post-shock conditions 
      SUBROUTINE post_shock_eq (p1, u1, T1, p2, u2, T2, yi, rhoi, xN)

        USE mod_nitrogen_Park_initialize_CFD,         ONLY: nb_ns, pos_N, pos_N2, urg, Rn, Rn2, solver, mi

        INTEGER :: i, length
        REAL(KIND=8), PARAMETER :: tol = 1.d-8
        REAL(KIND=8) :: f, fl, fr, mm, mml, mmr
        REAL(KIND=8) :: ov_mm, ov_mml, ov_mmr, tmp
        REAL(KIND=8) :: res_r, res_T, ratio, ratio_old, T, Tl, Tr, T_old
        REAL(KIND=8) :: m_dot, R, rho1, rho2, h1, h2
        REAL(KIND=8), DIMENSION(nb_ns) :: h, h_left, h_right
        REAL(KIND=8), DIMENSION(nb_ns) :: yi_left, yi_right
        REAL(KIND=8), DIMENSION(nb_ns) :: xi, xi_left, xi_right 
        REAL(KIND=8), DIMENSION(3) :: left, right, res

        REAL(KIND=8), INTENT(IN) :: p1, u1, T1
        REAL(KIND=8), INTENT(IN), OPTIONAL :: xN
        REAL(KIND=8), INTENT(OUT) :: p2, u2, T2
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi, rhoi

        length = LEN_TRIM(solver)

        IF (PRESENT(xN)) THEN

           ! Atomic nitrogen in the free-stream
           xi(pos_N)  = xN 
           xi(pos_N2) = 1.d0 - xN

        ELSE 

           CALL eq_composition (p1, T1, xi)

        ENDIF
 
        ! Molecular mass and N and N2 species mass fractions
        mm = 0.d0
        DO i = 1,nb_ns
           mm = mm + xi(i)*mi(i)
        ENDDO
        ov_mm = 1.d0/mm        

        DO i = 1,nb_ns 
           yi(i) = xi(i)*mi(i)*ov_mm
        ENDDO

        R    = urg*ov_mm
        rho1 = p1/(R*T1)

        CALL enthalpy (T1, h)

        h1 = 0.d0
        DO i = 1,nb_ns 
           h1 = h1 + yi(i)*h(i)
        ENDDO
        m_dot = rho1*u1

        ! Mass, momentum and enthalpy fluxes (upstream of the shock)
        left(1) = m_dot
        left(2) = p1 + m_dot*u1
        left(3) = h1 + 0.5d0*u1**2

        ! Guess value of the density ratio
        ratio = 0.1d0
        res_r = 1.d0
        DO WHILE (res_r.GT.tol)

           ! Update of p2 and h2 values 
           p2  = p1 + m_dot*u1*(1.d0 - ratio)
           h2  = h1 + 0.5d0*u1**2*(1.d0 - ratio**2)

           ! The equation h(T2) = tmp is solved for T2
           ! For sake of robusteness the bisection method is used
           Tl  = 300.d0
           Tr  = 80000.d0
           T   = 0.5d0*(Tl + Tr)

           res_T = 1.d0
           DO WHILE (res_T.GT.tol)
 
              ! Species molar fractions and enthalpies per unit mass at T, Tl and Tr
              CALL eq_composition (p2, T, xi)
              CALL eq_composition (p2, Tl, xi_left)
              CALL eq_composition (p2, Tr, xi_right)
              
              CALL enthalpy (T, h)
              CALL enthalpy (Tl, h_left)
              CALL enthalpy (Tr, h_right)

              ! Mixture molecular mass at T, Tl and Tr
              mm  = 0.d0
              mml = 0.d0 
              mmr = 0.d0
              DO i = 1,nb_ns 
                 tmp = mi(i)
                 mm  = mm  + xi(i)*tmp
                 mml = mml + xi_left(i)*tmp
                 mmr = mmr + xi_right(i)*tmp
              ENDDO

              ov_mm  = 1.d0/mm
              ov_mml = 1.d0/mml
              ov_mmr = 1.d0/mmr 

              ! Species mass fractions at T, Tl and Tr
              DO i = 1,nb_ns 
                 tmp    = mi(i)
                 yi(i)  = xi(i)*tmp*ov_mm
                 yi_left(i)  = xi_left(i)*tmp*ov_mml
                 yi_right(i) = xi_right(i)*tmp*ov_mmr
              ENDDO

              ! Function at T, Tl and Tr
              f  = - h2
              fl = - h2
              fr = - h2
              DO i = 1,nb_ns 
                 f  = f  + yi(i)*h(i)
                 fl = fl + yi_left(i)*h_left(i)
                 fr = fr + yi_right(i)*h_right(i)
              ENDDO 
              
              ! Temperature update 
              T_old = T
              IF ((f*fr).LT.0.d0) THEN 

                 Tl = T 
                 T  = 0.5d0*(Tl + Tr)

              ELSEIF ((f*fl).LT.0.d0) THEN

                 Tr = T 
                 T  = 0.5d0*(Tl + Tr)

              ENDIF

              res_T = ABS(T - T_old)/T_old

           ENDDO          

           ! Density update 
           CALL eq_composition (p2, T, xi)            
           
           mm = 0.d0
           DO i = 1,nb_ns
              mm = mm + xi(i)*mi(i)
           ENDDO

           R    = urg/mm
           rho2 = p2/(R*T)
           
           ! Density ratio update 
           ratio_old = ratio
           ratio     = rho1/rho2
           res_r     = ABS(ratio - ratio_old)/ratio

        ENDDO

        u2 = u1*ratio
        p2 = p1 + m_dot*u1*(1.d0 - ratio)
        T2 = T

        CALL eq_composition (p2, T2, xi)            
           
        mm = 0.d0
        DO i = 1,nb_ns
           mm = mm + xi(i)*mi(i)
        ENDDO
        ov_mm = 1.d0/mm  

        R    = urg*ov_mm
        rho2 = p2/(R*T2)

        DO i = 1,nb_ns 
           yi(i) = xi(i)*mi(i)*ov_mm
        ENDDO 

        ! Species densities 
        rhoi = rho2*yi

        CALL enthalpy (T2, h)

        h2 = 0.d0
        DO i = 1,nb_ns 
           h2 = h2 + yi(i)*h(i)
        ENDDO
        m_dot = rho2*u2

        ! Mass, momentum and enthalpy fluxes (downstream of the shock)
        right(1) = m_dot
        right(2) = p2 + m_dot*u2
        right(3) = h2 + 0.5d0*u2**2

        DO i = 1,3
           res(i) = ABS(left(i) - right(i))/left(i)*100.d0
        ENDDO

        WRITE(*,5)solver(1:length),':: Nitrogen Park -> equilibrium post-shock conditions'
        PRINT*
        WRITE(*,10)'Residuals on mass, momentum and energy fluxes:'
        PRINT*
        WRITE(*,15)'Mass      ',res(1),' [%]'
        PRINT*
        WRITE(*,15)'Momentum  ',res(2),' [%]'
        PRINT*
        WRITE(*,15)'Energy    ',res(3),' [%]'
        PRINT* 
        WRITE(*,15)'p  ',p2,' [Pa]'
        WRITE(*,15)'u  ',u2,' [m/s]'
        WRITE(*,15)'T  ',T2,' [K]'
        PRINT*

5     FORMAT(A,A)
10    FORMAT(A)
15    FORMAT(A,E14.6,A)

      END SUBROUTINE post_shock_eq 

  END MODULE mod_Nitrogen_Park_CFD_eq
!------------------------------------------------------------------------------!
