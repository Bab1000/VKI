!------------------------------------------------------------------------------!
! This module provides implementation of subroutines for dealing with thermodynamics 
! of the em, Ar, Arp system.
  MODULE mod_argon_CR_CFD_prop 
 
    IMPLICIT NONE

    ! Subroutines for thermodynamics
    CONTAINS 

      !------------------------------------------------------!
      ! This subroutine computes the translation partition function.
      SUBROUTINE Q_trans (T, mm, Q)

        USE mod_argon_CR_initialize_CFD,     ONLY: fac_Q

        REAL(KIND=8), INTENT(IN)  :: T, mm
        REAL(KIND=8), INTENT(OUT) :: Q

        Q = fac_Q*(mm*T)**1.5
     
      END SUBROUTINE Q_trans

      !------------------------------------------------------!
      ! This subroutine computes the translation partition function and its
      ! derivative with respect to temperature.
      SUBROUTINE Q_trans_der (T, mm, Q, dQ_dT)

        USE mod_argon_CR_initialize_CFD,     ONLY: fac_Q

        REAL(KIND=8), INTENT(IN)  :: T, mm
        REAL(KIND=8), INTENT(OUT) :: Q, dQ_dT

        Q = fac_Q*(mm*T)**1.5
        dQ_dT = 1.5d0*Q/T

      END SUBROUTINE Q_trans_der 

      !----------------------------------------------------!
      ! This subroutine computes the internal partition function 
      SUBROUTINE Q_int (T, species, Q)

        USE mod_argon_CR_initialize_CFD,     ONLY: ukb, levels_Ar, levels_Arp, gi_Ar, gi_Arp, & 
                                                 & EiJ_Ar, EiJ_Arp 

        INTEGER :: i
        REAL(KIND=8) :: ov_kbT

        REAL(KIND=8), INTENT(IN) :: T
        CHARACTER*(*), INTENT(IN) :: species
        REAL(KIND=8), INTENT(OUT) :: Q

        ! Initialization
        Q = 0.d0

        ! Useful quantity
        ov_kbT = 1.d0/(ukb*T)

        SELECT CASE(species)
 
          ! Ar 
          CASE('Ar')
            DO i = 1,levels_Ar
               Q = Q + gi_Ar(i)*DEXP(-EiJ_Ar(i)*ov_kbT)
            ENDDO

          ! Arp
          CASE('Arp')
            DO i = 1,levels_Arp
               Q = Q + gi_Arp(i)*DEXP(-EiJ_Arp(i)*ov_kbT)
            ENDDO

          CASE DEFAULT
            PRINT*
            WRITE(*,10)'in "mod_argon_CR_CFD_eq"'
            WRITE(*,10)'error in species selection for internal partition function'
            STOP

        END SELECT 

10    FORMAT(A)

      END SUBROUTINE Q_int

      !----------------------------------------------------!
      ! This subroutine computes the temperatures from the energy density when
      ! using the CR model 
      SUBROUTINE get_temperatures_CR (rhoi, rho_e, temp)

        USE mod_argon_CR_initialize_CFD,     ONLY: nb_ns, pos_em, pos_Ar, pos_Arp, pos_Te, levels_Ar,  & 
                                                 & levels_Arp, cv_tr_em, cv_tr_Ar, cv_tr_Arp,          & 
                                                 & ei_Ar, ei_Arp, hf_Arp, R_em, gamma_e_m1

        INTEGER :: i
        REAL(KIND=8) :: a, b, tmp
        REAL(KIND=8) :: rho, rho_em, rho_eint, rho_se, T, Te

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, rho_e
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp

        ! Initialization
        rho_eint = rho_e(1)
        rho_se   = rho_e(2)
        rho_em   = rhoi(1)

        rho = 0.d0
        DO i = 1,nb_ns 
           rho = rho + rhoi(i) 
        ENDDO

        ! Free electron temperature 
        Te = rho_se*(rho**(gamma_e_m1))/(rho_em*R_em)

        ! Heavy particle temperature
        a = 0.d0
        b = 0.d0
        DO i = 1,levels_Ar
           tmp = rhoi(pos_Ar + i - 1)
           a   = a + tmp*cv_tr_Ar
           b   = b + tmp*ei_Ar(i) 
        ENDDO

        DO i = 1,levels_Arp
           tmp = rhoi(pos_Arp + i - 1)
           a   = a + tmp*cv_tr_Arp
           b   = b + tmp*(ei_Arp(i) + hf_Arp)
        ENDDO

        T = (rho_eint - rho_em*cv_tr_em*Te - b)/a

        temp(1)      = T
        temp(pos_Te) = Te
         
      END SUBROUTINE get_temperatures_CR

      !----------------------------------------------------!
      ! This subroutine computes the temperatures from the energy density when
      ! using the CR model 
      SUBROUTINE get_temperatures_MT (rhoi, rho_e, temp)

        USE mod_argon_CR_initialize_CFD,     ONLY: nb_ns, pos_em, pos_Ar, pos_Arp, pos_Te, levels_Ar,  & 
                                                 & levels_Arp, cv_tr_em, cv_tr_Ar, cv_tr_Arp,          & 
                                                 & ei_Ar, ei_Arp, hf_Arp, R_em, gamma_e_m1

        INTEGER :: i
        REAL(KIND=8) :: a, b, tmp
        REAL(KIND=8) :: rho, rho_em, rho_eint, rho_se, T, Te

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, rho_e
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp

        ! Initialization
        rho_eint = rho_e(1)
        rho_se   = rho_e(2)
        rho_em   = rhoi(1)

        rho = 0.d0
        DO i = 1,nb_ns 
           rho = rho + rhoi(i) 
        ENDDO

        ! Free electron temperature
        Te = rho_se*(rho**(gamma_e_m1))/(rho_em*R_em)

        ! Heavy particle temperature
        a = 0.d0
        b = 0.d0
        tmp = rhoi(pos_Ar)
        a   = a + tmp*cv_tr_Ar

        tmp = rhoi(pos_Arp)
        a   = a + tmp*cv_tr_Arp
        b   = tmp*hf_Arp

        T = (rho_eint - rho_em*cv_tr_em*Te - b)/a

        temp(1)      = T
        temp(pos_Te) = Te

      END SUBROUTINE get_temperatures_MT

      !----------------------------------------------------!
      ! This subroutine computes the specific species energies when using the CR model
      SUBROUTINE energy_CR (temp, e) 

        USE mod_argon_CR_initialize_CFD,     ONLY: pos_em, pos_Ar, pos_Arp, pos_Te, levels_Ar,  & 
                                                 & levels_Arp, cv_tr_em, cv_tr_Ar, cv_tr_Arp,   & 
                                                 & ei_Ar, ei_Arp, hf_Arp

        INTEGER :: i
        REAL(KIND=8) :: tmp
        REAL(KIND=8) :: T, Te

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e

        ! Heavy particle and free electron temperatures 
        T  = temp(1)
        Te = temp(pos_Te)

        ! em 
        e(pos_em)  = cv_tr_em*Te

        ! Ar - electronic levels
        tmp = cv_tr_Ar*T
        DO i = 1,levels_Ar
           e(pos_Ar + i - 1) = tmp + ei_Ar(i) 
        ENDDO
       
        ! Arp - electronic levels
        tmp = cv_tr_Arp*T
        DO i = 1,levels_Arp
           e(pos_Arp + i - 1) = tmp + ei_Arp(i) + hf_Arp
        ENDDO 
        
      END SUBROUTINE energy_CR

      !----------------------------------------------------!
      ! This subroutine computes the specific species energies when using the MT model
      SUBROUTINE energy_MT (temp, e)

        USE mod_argon_CR_initialize_CFD,     ONLY: pos_em, pos_Ar, pos_Arp, pos_Te, levels_Ar, levels_Arp,  & 
                                                 & ukb, una, mm_Ar, mm_Arp, cv_tr_em, cv_tr_Ar, cv_tr_Arp,  & 
                                                 & gi_Ar, gi_Arp, EiJ_Ar, EiJ_Arp, hf_Arp

        INTEGER :: i
        REAL(KIND=8) :: exp_fac, e_el, Ei, ov_kb_Te, sum1, sum2
        REAL(KIND=8) :: T, Te

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e

        ! Heavy particle and free electron temperatures 
        T  = temp(1)
        Te = temp(pos_Te)

        ! Useful quantity 
        ov_kb_Te = 1.d0/(ukb*Te)

        ! em 
        e(pos_em)  = cv_tr_em*Te

        ! Ar 
        ! Electronic energy 
        sum1 = 0.d0
        sum2 = 0.d0
        DO i = 1,levels_Ar
           Ei = EiJ_Ar(i)
           exp_fac = gi_Ar(i)*DEXP(-Ei*ov_kb_Te) 
           sum1    = sum1 + exp_fac
           sum2    = sum2 + Ei*exp_fac
        ENDDO
        e_el = (sum2/sum1)*una/mm_Ar

        e(pos_Ar)  = cv_tr_Ar*T + e_el

        ! Arp 
        ! Electronic energy 
        sum1 = 0.d0
        sum2 = 0.d0
        DO i = 1,levels_Arp
           Ei = EiJ_Arp(i)
           exp_fac = gi_Arp(i)*DEXP(-Ei*ov_kb_Te) 
           sum1    = sum1 + exp_fac
           sum2    = sum2 + Ei*exp_fac
        ENDDO
        e_el = (sum2/sum1)*una/mm_Arp

        e(pos_Arp)  = cv_tr_Arp*T + e_el + hf_Arp 

      END SUBROUTINE energy_MT

      !----------------------------------------------------!
      ! This subroutine computes the specific species energy and constant volume specific
      ! heats when using the CR model
      SUBROUTINE energy_cv_CR (temp, e, cv) 

        USE mod_argon_CR_initialize_CFD,     ONLY: pos_em, pos_Ar, pos_Arp, pos_Te, levels_Ar,  & 
                                                 & levels_Arp, cv_tr_em, cv_tr_Ar, cv_tr_Arp,   & 
                                                 & ei_Ar, ei_Arp, hf_Arp

        INTEGER :: i
        REAL(KIND=8) :: tmp
        REAL(KIND=8) :: T, Te

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e, cv

        ! Heavy particle and free electron temperatures 
        T  = temp(1)
        Te = temp(pos_Te)

        ! em 
        e(pos_em)  = cv_tr_em*Te
        cv(pos_em) = cv_tr_em

        ! Ar - electronic levels
        tmp = cv_tr_Ar*T
        DO i = 1,levels_Ar
           e(pos_Ar + i - 1)  = tmp + ei_Ar(i)
           cv(pos_Ar + i - 1) = cv_tr_Ar 
        ENDDO
       
        ! Arp - electronic levels
        tmp = cv_tr_Arp*T
        DO i = 1,levels_Arp
           e(pos_Arp + i - 1)  = tmp + ei_Arp(i) + hf_Arp
           cv(pos_Arp + i - 1) = cv_tr_Arp
        ENDDO 
        
      END SUBROUTINE energy_cv_CR

      !----------------------------------------------------!
      ! This subroutine computes the specific species energy and constant volume specific
      ! heats when using the MT model
      SUBROUTINE energy_cv_MT (temp, e, cv)

        USE mod_argon_CR_initialize_CFD,     ONLY: pos_em, pos_Ar, pos_Arp, pos_Te, levels_Ar, levels_Arp,  & 
                                                 & ukb, una, mm_Ar, mm_Arp, cv_tr_em, cv_tr_Ar, cv_tr_Arp,  & 
                                                 & gi_Ar, gi_Arp, EiJ_Ar, EiJ_Arp, hf_Arp

        INTEGER :: i
        REAL(KIND=8) :: exp_fac, e_el, Ei, ov_kb_Te, sum1, sum2
        REAL(KIND=8) :: T, Te

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e, cv

        ! Heavy particle and free electron temperatures 
        T  = temp(1)
        Te = temp(pos_Te)

        ! Useful quantity 
        ov_kb_Te = 1.d0/(ukb*Te)

        ! em 
        e(pos_em)  = cv_tr_em*Te
        cv(pos_em) = cv_tr_em

        ! Ar 
        ! Electronic energy 
        !sum1 = 0.d0
        !sum2 = 0.d0
        !DO i = 1,levels_Ar
        !   Ei = EiJ_Ar(i)
        !   exp_fac = gi_Ar(i)*DEXP(-Ei*ov_kb_Te) 
        !   sum1    = sum1 + exp_fac
        !   sum2    = sum2 + Ei*exp_fac
        !ENDDO
        !e_el = (sum2/sum1)*una/mm_Ar

        cv(pos_Ar) = cv_tr_Ar
        e(pos_Ar)  = cv_tr_Ar*T !+ e_el

        ! Arp 
        ! Electronic energy 
        !sum1 = 0.d0
        !sum2 = 0.d0
        !DO i = 1,levels_Arp
        !   Ei = EiJ_Arp(i)
        !   exp_fac = gi_Arp(i)*DEXP(-Ei*ov_kb_Te) 
        !   sum1    = sum1 + exp_fac
        !   sum2    = sum2 + Ei*exp_fac
        !ENDDO
        !e_el = (sum2/sum1)*una/mm_Arp

        cv(pos_Arp) = cv_tr_Arp
        e(pos_Arp)  = cv_tr_Arp*T + hf_Arp !+ e_el

      END SUBROUTINE energy_cv_MT

      !----------------------------------------------------!
      ! This subroutine computes the post-shock nonequilibrium conditions
      SUBROUTINE post_shock_neq (p1, u1, T1, p2, u2, T2, yi)

        USE mod_argon_CR_initialize_CFD,     ONLY: nb_ns, levels_Ar, levels_Arp, pos_em, pos_Ar, pos_Arp, & 
                                                 & mm_species, urg
        USE mod_argon_CR_CFD_eq,             ONLY: eq_composition_CR

        INTEGER :: i 
        REAL(KIND=8) :: c1, g, gp1, mm, M1, M1s, rho1, R
        REAL(KIND=8), DIMENSION(nb_ns) :: rhoi

        REAL(KIND=8), INTENT(IN) :: p1, u1, T1
        REAL(KIND=8), INTENT(OUT) :: p2, u2, T2
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi

        ! Pre-shock densities and mass fractions
        CALL eq_composition_CR (p1, T1, rhoi)

        rho1 = 0.d0
        DO i = 1,nb_ns 
           rho1 = rho1 + rhoi(i)
        ENDDO
        
        rho1 = 1.d0/rho1
        DO i = 1,nb_ns 
           yi(i) = rhoi(i)*rho1
        ENDDO

        ! Mixture molar mass (inverse of)
        mm = 0.d0
        DO i = 1,nb_ns
           mm = mm + yi(i)/mm_species(i)
        ENDDO
    
        ! Mixture gas constant
        R = urg*mm

        ! Post-shock conditions (the free electron temperature is assumed to be
        ! frozen across the shock)
        ! Specific heat ratio (monoatomic gas-like case)
        g   = 5.d0/3.d0
        gp1 = g + 1.d0
        
        c1  = SQRT(g*R*T1)
        M1  = u1/c1
        M1s = M1**2

        ! Pressure, velocity and temperature after the shock
        p2 = p1*(2.d0*g*M1s - g + 1.d0)/gp1
        u2 = u1 - c1*2.d0/gp1*(M1 - 1.d0/M1)
        T2 = T1*(2.d0*g*M1s - g + 1.d0)*(g - 1.d0 + 2.d0/M1s)/(gp1**2) 

      END SUBROUTINE post_shock_neq 

  END MODULE mod_argon_CR_CFD_prop
!------------------------------------------------------------------------------!
