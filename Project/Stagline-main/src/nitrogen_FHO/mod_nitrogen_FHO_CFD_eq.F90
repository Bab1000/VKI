!------------------------------------------------------------------------------!
! This module provides subroutines for equilibrium properties of the FHO database for the N-N2 ystem 
  MODULE mod_nitrogen_FHO_CFD_eq

    USE mod_nitrogen_FHO_initialize_CFD,     ONLY: mm_N, mm_N2

    IMPLICIT NONE

    ! Subroutine for equilibrium properties
    CONTAINS 

      !----------------------------------------------------!
      ! This subroutine computes the macroscopic equilibrium composition (molar
      ! factions of N and N2 species)
      SUBROUTINE eq_composition (p, T, xn, xn2)

        USE mod_nitrogen_FHO_initialize_CFD,    ONLY: gn, ukb, Edis

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
      ! This subroutine computes the equilibrium composition   
      SUBROUTINE eq_composition_vib (p, T, rhoi, x_N)

        USE mod_nitrogen_FHO_initialize_CFD,    ONLY: nb_bins, ukb, una, urg, EkJ
                                                            

        INTEGER :: i
        REAL(KIND=8) :: xn, xn2
        REAL(KIND=8) :: dexp0, Qint
        REAL(KIND=8) :: fac1, fac2, ov_kbT, tmp1, tmp2
        REAL(KIND=8) :: mass, yN, yN2, nb_N2, p_N2, rho

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

        ! Partial densities of N and N2 vibrational levels
        rhoi(1) = rho*yN
        fac2    = nb_n2*fac1
        Qint    = 0.d0
        DO i = 1,nb_bins 
           dexp0 = DEXP(-EkJ(i)*ov_kbT)
           Qint  = Qint + dexp0
           rhoi(i + 1) = fac2*dexp0
        ENDDO

        Qint = 1.d0/Qint 
        DO i = 1,nb_bins 
           rhoi(i + 1) = rhoi(i + 1)*Qint
        ENDDO
        
      END SUBROUTINE eq_composition_vib

      !----------------------------------------------------!
      ! This subroutine computes the translational patition function
      SUBROUTINE tr_part_function (T, mm, Qtr)

        USE mod_nitrogen_FHO_initialize_CFD,     ONLY: fac_Q

        REAL(KIND=8), INTENT(IN)  :: T, mm
        REAL(KIND=8), INTENT(OUT) :: Qtr

        Qtr = fac_Q*(mm*T)**1.5

      END SUBROUTINE tr_part_function 

      !----------------------------------------------------!
      ! This subtroutine compute the internal partition fuction
      SUBROUTINE int_part_function (T, Qint)

        USE mod_nitrogen_FHO_initialize_CFD,     ONLY: nb_bins, ov_sigma, theta_rot, ukb, EkJ

        INTEGER :: i
        REAL(KIND=8) :: ov_kb_T
        REAL(KIND=8) :: Qrot, Qvib

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), INTENT(OUT) :: Qint 

        ! Rotational partition function
        Qrot = T*ov_sigma/theta_rot

        ! Vibrational partition function
        ov_kb_T = 1.d0/(ukb*T)
        Qvib    = 0.d0
        DO i = 1,nb_bins 
           Qvib = Qvib + DEXP(-EkJ(i)*ov_kb_T)
        ENDDO

        ! Internal partition function
        Qint = Qrot*Qvib

      END SUBROUTINE int_part_function
 
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
      ! This subroutine computes the internal energy per unit mass 
      ! of N and N2 species in equilibrium conditions
      SUBROUTINE energy (T, e)

        USE mod_nitrogen_FHO_initialize_CFD,    ONLY: nb_bins, cv_tr_N, cv_tr_N2, R_N2, hf_N, & 
                                                         &  ukb, una, mm_N2, EkJ

        INTEGER :: i
        REAL(KIND=8) :: tmp1, tmp2, tmp3, tmp4, tmp5 
        REAL(KIND=8) :: erot, evib

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: e

        ! N atom
        e(1) = cv_tr_N*T + hf_N 
        
        ! N2 molecule
        erot = R_N2*T 
         
        tmp1 = 0.d0
        tmp2 = 0.d0
        tmp3 = 1.d0/(ukb*T)
        DO i = 1,nb_bins 
           tmp4 = EkJ(i)
           tmp5 = DEXP(-tmp4*tmp3)
           tmp1 = tmp1 + tmp5
           tmp2 = tmp2 + tmp4*tmp5
        ENDDO
        
        ! Conversion [J] -> [J/kg]
        evib = tmp2*una/(tmp1*mm_N2)

        e(2) = cv_tr_N2*T + erot + evib

      END SUBROUTINE energy

      !----------------------------------------------------!
      ! This subroutine computes the enthalpy per unit mass 
      ! of N and N2 species in equilibrium conditions 
      SUBROUTINE enthalpy (T, h)

        USE mod_nitrogen_FHO_initialize_CFD,       ONLY: R_N, R_N2

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), INTENT(OUT), DIMENSION(:) :: h

        CALL energy (T, h)

        h(1) = h(1) + R_N*T
        h(2) = h(2) + R_N2*T

      END SUBROUTINE enthalpy

      !----------------------------------------------------!
      ! This subroutine computes the entropy per unit mass 
      ! of N and N2 species in equilibrium conditions. The entropy of 
      ! mixing contribution is computed as well. 
      SUBROUTINE entropy (p, T, yi, xi, smix, s) 

        USE mod_nitrogen_FHO_initialize_CFD,       ONLY: nb_ns, gN, ukb, hf_N, mm_N, mm_N2, R_N, R_N2 

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
        h(1) = h(1) - hf_N

        tmp1 = 1.d0/T 
        tmp2 = ukb*T/p

        s(1) = h(1)*tmp1 + R_N*DLOG(tmp2*QtN) 
        s(2) = h(2)*tmp1 + R_N2*DLOG(tmp2*QtN2*QintN2)

        ! Entropy of mixing 
        smix = - R_N*yi(1)*DLOG(xi(1)) - R_N2*yi(2)*DLOG(xi(2))

      END SUBROUTINE entropy

      !----------------------------------------------------!
      SUBROUTINE post_shock_eq (p1, u1, T1, p2, u2, T2, yi, rhoi, xN)

        USE mod_Nitrogen_FHO_initialize_CFD,           ONLY: nb_ns, urg, solver

        INTEGER :: i, length
        REAL(KIND=8), PARAMETER :: tol = 1.d-8
        REAL(KIND=8) :: yN, yN2, mass
        REAL(KIND=8) :: rho1, rho2, h1, h2, R, m_dot
        REAL(KIND=8) :: f, fl, fr, h_left, h_right, rho, rho_left, rho_right, Tl, Tr, T, T_old
        REAL(KIND=8) :: ratio, ratio_old, res_r, res_T
        REAL(KIND=8) :: tmp1
        REAL(KIND=8), DIMENSION(3) :: left, right, res
        REAL(KIND=8), DIMENSION(nb_ns) :: rhoi_left, rhoi_right, yi_left, yi_right, Ri, mm, ei, ei_left, ei_right

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
        CALL eq_composition_vib (p1, T1, rhoi, xN)

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

        CALL get_species_energy (T1, ei)

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
           Tr  = 60000.d0
           T   = 0.5d0*(Tl + Tr)

           res_T = 1.d0
           DO WHILE (res_T.GT.tol)
 
              ! Species mass fractions and enthalpies per unit mass at T, Tl and Tr
              CALL eq_composition_vib (p2, T, rhoi)
              CALL eq_composition_vib (p2, Tl, rhoi_left)
              CALL eq_composition_vib (p2, Tr, rhoi_right) 

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

              CALL get_species_energy (T, ei)
              CALL get_species_energy (Tl, ei_left)
              CALL get_species_energy (Tr, ei_right)

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
           CALL eq_composition_vib (p2, T, rhoi)
           
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

        CALL eq_composition_vib (p2, T, rhoi)
        CALL get_species_energy  (T2, ei)

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

        WRITE(*,5)solver(1:length),':: Nitrogen FHO -> equilibrium post-shock conditions'
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

  END MODULE mod_nitrogen_FHO_CFD_eq
!------------------------------------------------------------------------------!
