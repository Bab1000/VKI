!------------------------------------------------------------------------------!
! This module provides subroutines for equilibrium properties of the NASA Ames
! database for the N-N2 ystem 
  MODULE mod_nitrogen_NASA_CFD_eq

    USE mod_nitrogen_NASA_initialize_CFD,     ONLY: mm_N, mm_N2

    IMPLICIT NONE

    ! Subroutines for equilibrium properties 
    CONTAINS 

      !----------------------------------------------------!
      ! This subroutine computes the macroscopic equilibrium composition (molar
      ! factions of N and N2 species)
      SUBROUTINE eq_composition (p, T, xn, xn2)

        USE mod_nitrogen_NASA_initialize_CFD,    ONLY: gn, ukb, Edis

        REAL(KIND=8) :: delta2, rhs, kbT
        REAL(KIND=8) :: Qn, Qn2, Qint 

        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), INTENT(OUT) :: xn, xn2

        ! Translational and internal partition functions of N and N2 species
        CALL tr_part_function  (T, mm_N, Qn)
        CALL tr_part_function  (T, mm_N2, Qn2)
        CALL int_part_function (T, Qint)

        kbT = ukb*T  
        rhs = (1.d0/p)*kbT*DEXP(-Edis/kbT)*((gn*Qn)**2)/(Qn2*Qint)

        delta2 = rhs**2 + 4.d0*rhs 

        ! Molar fractions of N and N2 species
        xn  = 0.5d0*(- rhs + DSQRT(delta2))
        xn2 = 1.d0 - xn

      END SUBROUTINE eq_composition

      !----------------------------------------------------!
      ! This subroutine computes the internal partition function (ro-vibrational
      ! levels) for the N2 molecule
      SUBROUTINE int_part_function (T, Qint)

        USE mod_nitrogen_NASA_initialize_CFD,     ONLY: model, nb_bins, levels, ukb, EvJ, EkJ, gvJ, theta_rot 

        INTEGER :: i
        REAL(KIND=8) :: ov_kbT

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), INTENT(OUT) :: Qint 

        ! Useful quantity 
        ov_kbT = 1.d0/(ukb*T)

        ! Internal partition function
        Qint = 0.d0

        ! Rigid rotor for rotational energy
        IF (model.EQ.'VC_rr') THEN

           DO i = 1,nb_bins
              Qint = Qint + DEXP(-EkJ(i)*ov_kbT)
           ENDDO
           Qint = 0.5d0*Qint*T/theta_rot

        ELSE 

           DO i = 1,levels 
              Qint = Qint + gvJ(i)*DEXP(-EvJ(i)*ov_kbT)
           ENDDO

        ENDIF

      END SUBROUTINE int_part_function

      !----------------------------------------------------!
      ! This subroutine computes the translational patition function
      SUBROUTINE tr_part_function (T, mm, Qtr)

        USE mod_nitrogen_NASA_initialize_CFD,     ONLY: fac_Q

        REAL(KIND=8), INTENT(IN)  :: T, mm
        REAL(KIND=8), INTENT(OUT) :: Qtr

        Qtr = fac_Q*(mm*T)**1.5

      END SUBROUTINE tr_part_function 

      !----------------------------------------------------!
      ! This subroutine computes the equilibrium composition when using the 
      ! Vibrational, bin Boltzmann and Full collisional models. This subroutine should not 
      ! be used in case of bin Uniform levels. 
      SUBROUTINE eq_composition_bins (p, T, rhoi, x_N)

        USE mod_nitrogen_NASA_initialize_CFD,    ONLY: nb_bins, levels, level_bin, ukb, una, urg,  &
                                                    &  degen, ek, delta_ek, gvJ, model, theta_rot

        INTEGER :: i, bin
        REAL(KIND=8) :: xn, xn2
        REAL(KIND=8) :: dexp0, EkJ, Qint
        REAL(KIND=8) :: fac1, fac2, ov_kbT, tmp1, tmp2
        REAL(KIND=8) :: mass, yN, yN2, nb_N2, p_N2, rho
        REAL(KIND=8), DIMENSION(nb_bins) :: Qk

        REAL(KIND=8), INTENT(IN) :: p, T
        REAL(KIND=8), INTENT(IN), OPTIONAL :: x_N
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: rhoi

        ! Molar fractions of N and N2 species
        IF (PRESENT(x_N)) THEN 

           xN  = x_N
           xN2 = 1.d0 - x_N
           
        ELSE 

           CALL eq_composition (p, T, xN, xN2)

        ENDIF

        ! Useful data
        ov_kbT = 1.d0/(ukb*T)
        fac1   = mm_N2/una

        ! Mixture molar mass
        mass = xN*mm_N + xN2*mm_N2
        
        ! Mass fractions of N and N2 species
        yN  = xN*mm_N/mass 
        yN2 = xN2*mm_N2/mass

        rho   = p/(urg*T)*mass
        p_N2  = yN2*rho*urg/mm_N2*T
        nb_N2 = p_N2*ov_kbT
        
        ! Evaluation of bin internal partition function
        Qk   = 0.d0
        Qint = 0.d0

        ! Rigid rotor used for rotational energy
        IF (model.EQ.'VC_rr') THEN

           tmp1 = 0.5d0*T/theta_rot
           DO i = 1,nb_bins
              EkJ   = ek(i)*fac1
              tmp2  = DEXP(-EkJ*ov_kbT)
              Qk(i) = tmp1
              Qint  = Qint + tmp2
           ENDDO
           Qint = Qint*tmp1

        ELSE

           DO i = 1,levels
              bin  = level_bin(i)
              EkJ  = ek(bin)*fac1
              tmp1 = degen(i)*DEXP(-delta_ek(i)*ov_kbT)
              tmp2 = tmp1*DEXP(-EkJ*ov_kbT)
              Qk(bin) = Qk(bin) + tmp1
              Qint    = Qint    + tmp2
           ENDDO

        ENDIF
        
        ! Partial densities of N and N2 bins
        rhoi(1) = rho*yN

        IF ((model.NE.'MT_TTint').AND.(model.NE.'MT_TrTv')) THEN

          fac2 = nb_n2*fac1
          Qint = 1.d0/Qint
          DO i = 1,nb_bins 
             EkJ   = ek(i)*fac1
             dexp0 = DEXP(-EkJ*ov_kbT)
             rhoi(i + 1) = fac2*(Qk(i)*dexp0*Qint)
          ENDDO
       
        ! Macroscopic models 
        ELSE 

          rhoi(2) = rho*yN2

        ENDIF

      END SUBROUTINE eq_composition_bins
 
      !----------------------------------------------------!
      ! This subroutine computes the internal energy per unit mass 
      ! of N and N2 species in equilibrium conditions
      SUBROUTINE energy (T, e)

        USE mod_nitrogen_NASA_initialize_CFD,     ONLY: levels, ukb, una, cv_tr_n, cv_tr_n2, hf_n, mm_N2, gvJ, EvJ

        INTEGER :: i
        REAL(KIND=8) :: eint
        REAL(KIND=8) :: tmp1, tmp2, tmp3, tmp4, tmp5 

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: e

        ! Nitrogen atom N
        e(1) = cv_tr_n*T + hf_n

        tmp1 = 0.d0
        tmp2 = 0.d0
        tmp3 = 1.d0/(ukb*T)
        DO i = 1,levels
           tmp4 = EvJ(i)
           tmp5 = gvJ(i)*DEXP(-tmp4*tmp3)
           tmp1 = tmp1 + tmp5
           tmp2 = tmp2 + tmp4*tmp5
        ENDDO

        ! Internal energy of N2 molecule 
        eint = tmp2*una/(tmp1*mm_N2) 

        ! Nitrogen molecule (translation and internal contributions)
        e(2) = cv_tr_n2*T + eint

      END SUBROUTINE energy

      !----------------------------------------------------!
      ! This subroutine computes the enthalpy per unit mass 
      ! of N and N2 species in equilibrium conditions 
      SUBROUTINE enthalpy (T, h)

        USE mod_nitrogen_NASA_initialize_CFD,     ONLY: Rn, Rn2

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: h

        CALL energy (T, h) 

        h(1) = h(1) + Rn*T
        h(2) = h(2) + Rn2*T

      END SUBROUTINE enthalpy

      !----------------------------------------------------!
      ! This subroutine computes the entropy per unit mass 
      ! of N and N2 species in equilibrium conditions. The entropy of 
      ! mixing contribution is computed as well. 
      SUBROUTINE entropy (p, T, yi, xi, smix, s)

        USE mod_nitrogen_NASA_initialize_CFD,     ONLY: nb_ns, ukb, hf_n, gn, Rn, Rn2, mm_N, mm_N2

        REAL(KIND=8) :: tmp1, tmp2
        REAL(KIND=8) :: QtN, QtN2, QintN2
        REAL(KIND=8), DIMENSION(nb_ns) :: h

        REAL(KIND=8), INTENT(IN) :: p, T
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
        h(1) = h(1) - hf_n

        tmp1 = 1.d0/T 
        tmp2 = ukb*T/p

        s(1) = h(1)*tmp1 + Rn*DLOG(tmp2*QtN) 
        s(2) = h(2)*tmp1 + Rn2*DLOG(tmp2*QtN2*QintN2)

        ! Entropy of mixing 
        smix = - Rn*yi(1)*DLOG(xi(1)) - Rn2*yi(2)*DLOG(xi(2))

      END SUBROUTINE entropy

      !----------------------------------------------------!
      ! This subroutine computes the equilibrium post-shock conditions 
      SUBROUTINE post_shock_eq (p1, u1, T1, p2, u2, T2, yi, rhoi, xN)

        USE mod_Nitrogen_NASA_initialize_CFD,           ONLY: nb_ns, nb_temp, urg, solver
        USE mod_function_pointer_NASA,                  ONLY: get_species_energy_NASA  

        INTEGER :: i, length
        REAL(KIND=8), PARAMETER :: tol = 1.d-8
        REAL(KIND=8) :: yN, yN2, mass
        REAL(KIND=8) :: rho1, rho2, h1, h2, R, m_dot
        REAL(KIND=8) :: f, fl, fr, h_left, h_right, rho, rho_left, rho_right, Tl, Tr, T, T_old
        REAL(KIND=8) :: ratio, ratio_old, res_r, res_T
        REAL(KIND=8) :: tmp1
        REAL(KIND=8), DIMENSION(3) :: left, right, res
        REAL(KIND=8), DIMENSION(nb_ns) :: rhoi_left, rhoi_right, yi_left, yi_right, Ri, mm, ei, ei_left, ei_right
        REAL(KIND=8), DIMENSION(nb_temp) :: temp

        REAL(KIND=8), INTENT(IN) :: p1, u1, T1
        REAL(KIND=8), INTENT(IN), OPTIONAL :: xN
        REAL(KIND=8), INTENT(OUT) :: p2, u2, T2
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi, rhoi

        length = LEN_TRIM(solver)

        ! Useful data
        Ri(1) = urg/mm_N
        mm(1) = mm_N
        DO i = 2,nb_ns 
           Ri(i) = urg/mm_N2
           mm(i) = mm_N2
        ENDDO

        ! Pre-shock species mass fractions
        CALL eq_composition_bins (p1, T1, rhoi, xN)

        rho1 = 0.d0
        DO i = 1,nb_ns 
           rho1 = rho1 + rhoi(i)
        ENDDO

        DO i = 1,nb_ns 
           yi(i) = rhoi(i)/rho1
        ENDDO

        ! Mixture molar mass
        mass = 0.d0
        DO i = 1,nb_ns 
           mass = mass + yi(i)/mm(i)
        ENDDO
        mass = 1.d0/mass
         
        ! Pre-shock density
        R    = urg/mass
        rho1 = p1/(R*T1)
        temp = T1
        CALL get_species_energy_NASA (temp, ei)

        h1 = 0.d0
        DO i = 1,nb_ns 
           h1 = h1 + yi(i)*(ei(i) + Ri(i)*T1)
        ENDDO
        m_dot   = rho1*u1

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
           Tr  = 50000.d0
           T   = 0.5d0*(Tl + Tr)

           res_T = 1.d0
           DO WHILE (res_T.GT.tol)
 
              ! Species mass fractions and enthalpies per unit mass at T, Tl and Tr
              CALL eq_composition_bins (p2, T, rhoi)
              CALL eq_composition_bins (p2, Tl, rhoi_left)
              CALL eq_composition_bins (p2, Tr, rhoi_right) 

              rho = 0.d0
              rho_left  = 0.d0
              rho_right = 0.d0
              DO i = 1,nb_ns 
                 rho = rho + rhoi(i)
                 rho_left  = rho_left + rhoi_left(i)
                 rho_right = rho_right + rhoi_right(i)
              ENDDO              

              DO i = 1,nb_ns 
                 yi(i) = rhoi(i)/rho
                 yi_left(i)  = rhoi_left(i)/rho_left
                 yi_right(i) = rhoi_right(i)/rho_right
              ENDDO

              temp = T
              CALL get_species_energy_NASA (temp, ei)
              temp = Tl
              CALL get_species_energy_NASA (temp, ei_left) 
              temp = Tr
              CALL get_species_energy_NASA (temp, ei_right)

              ! Function at T, Tl and Tr
              f  = - h2
              fl = - h2
              fr = - h2
              DO i = 1,nb_ns 
                 tmp1 = Ri(i)
                 f  = f  + yi(i)*(ei(i) + tmp1*T)
                 fl = fl + yi_left(i)*(ei_left(i) + tmp1*Tl)
                 fr = fr + yi_right(i)*(ei_right(i) + tmp1*Tr)
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
           CALL eq_composition_bins (p2, T, rhoi)
           
           rho2 = 0.d0
           DO i = 1,nb_ns 
              rho2 = rho2 + rhoi(i)
           ENDDO              

           ! Density ratio update 
           ratio_old = ratio
           ratio     = rho1/rho2
           res_r     = ABS(ratio - ratio_old)/ratio

        ENDDO

        u2 = u1*ratio
        p2 = p1 + m_dot*u1*(1.d0 - ratio)
        T2 = T

        temp = T2
        CALL eq_composition_bins (p2, T, rhoi)
        CALL get_species_energy_NASA (temp, ei)

        rho2 = 0.d0
        DO i = 1,nb_ns 
          rho2 = rho2 + rhoi(i)
        ENDDO    
        m_dot = rho2*u2   

        DO i = 1,nb_ns 
           yi(i) = rhoi(i)/rho2
        ENDDO

        h2 = 0.d0
        DO i = 1,nb_ns 
           h2 = h2 + yi(i)*(ei(i) + Ri(i)*T2)
        ENDDO

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

5     FORMAT(A,A)
10    FORMAT(A)
15    FORMAT(A,E14.6,A) 

      END SUBROUTINE post_shock_eq

  END MODULE mod_nitrogen_NASA_CFD_eq
!------------------------------------------------------------------------------!
