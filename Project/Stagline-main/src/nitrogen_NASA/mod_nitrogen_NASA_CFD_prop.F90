!------------------------------------------------------------------------------!
! This module provides subroutines/function for the thermodynamic properties of the NASA Ames database for the N-N2 system.
  MODULE mod_nitrogen_NASA_CFD_prop

    IMPLICIT NONE

      ! Subroutine for thermodynamics
      CONTAINS 

      !------------------------------------------------------!
      ! This subroutine computes the translation partition function.
      SUBROUTINE Q_trans (T, mm, Q)

        USE mod_nitrogen_NASA_initialize_CFD,     ONLY: fac_Q

        REAL(KIND=8), INTENT(IN)  :: T, mm
        REAL(KIND=8), INTENT(OUT) :: Q

        Q = fac_Q*(mm*T)**1.5
     
      END SUBROUTINE Q_trans

      !------------------------------------------------------!
      ! This subroutine computes the translation partition function and its
      ! derivative with respect to temperature.
      SUBROUTINE Q_trans_der (T, mm, Q, dQ_dT)

        USE mod_nitrogen_NASA_initialize_CFD,     ONLY: fac_Q

        REAL(KIND=8), INTENT(IN)  :: T, mm
        REAL(KIND=8), INTENT(OUT) :: Q, dQ_dT

        Q = fac_Q*(mm*T)**1.5
        dQ_dT = 1.5d0*Q/T

      END SUBROUTINE Q_trans_der

      !----------------------------------------------------!
      ! This subroutine computes the internal partition functions when using the 
      ! various models. In case of RVC model it reduces to the degeneracy of
      ! each single level.
      SUBROUTINE Q_int_bins (T, Qk)

         USE mod_nitrogen_NASA_initialize_CFD,     ONLY: model, nb_bins, levels, ukb, level_bin, degen, delta_ek

         INTEGER :: i, pos
         REAL(KIND=8) :: dexp0, kbT

         REAL(KIND=8), INTENT(IN) :: T
         REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: Qk

         ! Initialization
         Qk = 0.d0

         SELECT CASE(model)

           CASE('RVC','BRVC','VC')

             kbT = 1.d0/(ukb*T)
             DO i = 1,levels
                pos   = level_bin(i)
                dexp0 = degen(i)*DEXP(-delta_ek(i)*kbT) 
                Qk(pos) = Qk(pos) + dexp0
             ENDDO
 
           CASE DEFAULT 

             Qk = 1.d0

         END SELECT

      END SUBROUTINE Q_int_bins

      !----------------------------------------------------!
      ! This subroutine computes the gas temperature in case of use of the RVC model.
      SUBROUTINE compute_T_RVC (rhoi, rho_eint, temp) 
 
        USE mod_nitrogen_NASA_initialize_CFD,         ONLY: nb_bins, hf_n, cv_tr_n, cv_tr_n2, ek 

        INTEGER :: i
        REAL(KIND=8) :: rho_n, rho_n2
        REAL(KIND=8) :: tmp
        REAL(KIND=8) :: e_sum, cv_sum 

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, rho_eint
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp

        ! N and N2 partial densities and other useful data
        rho_n  = rhoi(1)

        rho_n2 = 0.d0
        e_sum  = 0.d0
        DO i = 1,nb_bins 
           tmp = rhoi(i + 1)
           rho_n2 = rho_n2 + tmp
           e_sum  = e_sum + tmp*ek(i)
        ENDDO
        
        e_sum  = e_sum + rho_n*hf_n
        cv_sum = rho_n*cv_tr_n + rho_n2*cv_tr_n2

        temp = (rho_eint(1) - e_sum)/cv_sum

      END SUBROUTINE compute_T_RVC

      !----------------------------------------------------!
      ! This subroutine computes the gas temperature in case of use of the BRVC model.
      SUBROUTINE compute_T_BRVC (rhoi, rho_eint, temp) 

        USE mod_nitrogen_NASA_initialize_CFD,         ONLY: nb_bins, hf_n, Rn2, cv_tr_n, cv_tr_n2,  &
                                                          & ek, eint, cvint, T_store, T_min, T_max, & 
                                                          & inv_step
 
        INTEGER, PARAMETER :: it_max = 100
        INTEGER :: i, it, left, right
        REAL(KIND=8), PARAMETER :: tol = 1.d-6
        REAL(KIND=8) :: f, fp, g, de
        REAL(KIND=8) :: rho_n, rho_n2, cv_sum, cv, e_sum
        REAL(KIND=8) :: tmp, tmp1, tmp2, tmp3
        REAL(KIND=8) :: res, T, T_old

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, rho_eint
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp

        ! Useful data
        rho_n = rhoi(1)

        rho_n2 = 0.d0
        e_sum  = 0.d0 
        DO i = 1,nb_bins 
           tmp1 = rhoi(i + 1)
           tmp2 = tmp1*ek(i)
           e_sum  = e_sum + tmp2
           rho_n2 = rho_n2 + tmp1
        ENDDO
         
        e_sum  = e_sum + rho_n*hf_n
        g      = rho_eint(1) - e_sum

        ! Approximative derivarive
        cv_sum = rho_n*cv_tr_n + rho_n2*cv_tr_n2     
        fp = cv_sum + rho_n2*Rn2

        ! Temperature guess value (rotational energy given by the rigid-rotator approximation)
        T = g/fp

        it  = 0
        res = 1.d0
        DO WHILE (res.GE.tol)

           ! Value search
           left = INT((T - T_min)*inv_step) + 1
           tmp1 = (T - T_store(left))*inv_step

           left  = (left - 1)*nb_bins 
           right = left + nb_bins

           tmp = 0.d0
           fp  = 0.d0
           DO i = 1,nb_bins 
              tmp2 = eint(left + i)
              tmp3 = cvint(left + i)
              de   = tmp2 + tmp1*(eint(right + i) - tmp2) 
              cv   = tmp3 + tmp1*(cvint(right + i) - tmp3) 
              tmp  = tmp  + de*rhoi(i + 1)
              fp   = fp + cv*rhoi(i + 1)
           ENDDO 
           fp = fp + cv_sum

           f = cv_sum*T + tmp

           ! Temperature update
           T_old = T
           T     = T - (f - g)/fp
          
           ! Out of bound check
           T = MAX(T,T_min)
           T = MIN(T,T_max) 
 
           it = it + 1 
           IF (it.GT.it_max) THEN
              PRINT*
              WRITE(*,10)'solver_fvmcc_F90:: excessive number of Newton iterator, in compute_T_BRVC ..'
              PRINT*
              PRINT*,'T',T,'Told',T_old
              PRINT*
              PRINT*,'Residual',res
              PRINT*
              STOP
           ENDIF

           res = ABS((T - T_old)/T)
           
        ENDDO 

        temp = T
        
10    FORMAT(A)

      END SUBROUTINE compute_T_BRVC

      !----------------------------------------------------!
      ! This subroutine computes the gas temperature in case of use of the VC model. 
      ! In this case a non-linear equation has to be solved. The Newton-Raphson method is used.
      SUBROUTINE compute_T_VC (rhoi, rho_eint, temp) 

        USE mod_nitrogen_NASA_initialize_CFD,         ONLY: nb_bins, hf_n, Rn2, cv_tr_n, cv_tr_n2,  &
                                                          & ek, eint, cvint, T_store, T_min, T_max, & 
                                                          & inv_step
 
        INTEGER, PARAMETER :: it_max = 100
        INTEGER :: i, it, left, right
        REAL(KIND=8), PARAMETER :: tol = 1.d-6
        REAL(KIND=8) :: f, fp, g, de
        REAL(KIND=8) :: rho_n, rho_n2, cv_sum, cv, e_sum
        REAL(KIND=8) :: tmp, tmp1, tmp2, tmp3
        REAL(KIND=8) :: res, T, T_old

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, rho_eint
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp

        ! Useful data
        rho_n = rhoi(1)

        rho_n2 = 0.d0
        e_sum  = 0.d0 
        DO i = 1,nb_bins 
           tmp1 = rhoi(i + 1)
           tmp2 = tmp1*ek(i)
           e_sum  = e_sum + tmp2
           rho_n2 = rho_n2 + tmp1
        ENDDO
         
        e_sum  = e_sum + rho_n*hf_n
        g      = rho_eint(1) - e_sum

        ! Approximative derivarive
        cv_sum = rho_n*cv_tr_n + rho_n2*cv_tr_n2     
        fp = cv_sum + rho_n2*Rn2

        ! Temperature guess value (rotational energy given by the rigid-rotator approximation)
        T = g/fp

        it  = 0
        res = 1.d0
        DO WHILE (res.GE.tol)

           ! Value search
           left  = INT((T - T_min)*inv_step) + 1
           tmp1  = (T - T_store(left))*inv_step

           left  = (left - 1)*nb_bins 
           right = left + nb_bins

           tmp = 0.d0
           fp  = 0.d0
           DO i = 1,nb_bins 
              tmp2 = eint(left + i)
              tmp3 = cvint(left + i)
              de   = tmp2 + tmp1*(eint(right + i) - tmp2) 
              cv   = tmp3 + tmp1*(cvint(right + i) - tmp3) 
              tmp  = tmp  + de*rhoi(i + 1)
              fp   = fp + cv*rhoi(i + 1)
           ENDDO 
           fp = fp + cv_sum

           f = cv_sum*T + tmp

           ! Temperature update
           T_old = T
           T     = T - (f - g)/fp
          
           ! Out of bound check
           T = MAX(T,T_min)
           T = MIN(T,T_max) 
 
           it = it + 1 
           IF (it.GT.it_max) THEN
              PRINT*
              WRITE(*,10)'solver_fvmcc_F90:: excessive number of Newton iterator, in compute_T_BRVC ..'
              PRINT*
              PRINT*,'T',T,'Told',T_old
              PRINT*
              PRINT*,'Residual',res
              PRINT*
              STOP
           ENDIF

           res = ABS((T - T_old)/T)
           
        ENDDO 

        temp = T

10    FORMAT(A)

      END SUBROUTINE compute_T_VC

      !----------------------------------------------------!
      ! This subroutine computes the gas temperature in case of use of the VC_rr model. 
      SUBROUTINE compute_T_VC_rr (rhoi, rho_eint, temp) 

        USE mod_nitrogen_NASA_initialize_CFD,         ONLY: nb_bins, hf_n, Rn2, cv_tr_n, cv_tr_n2, ek 
 
        INTEGER :: i
        REAL(KIND=8) :: rho_n, rho_n2, cv_sum, e_sum
        REAL(KIND=8) :: g, tmp, tmp1, tmp2

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, rho_eint
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp

        ! Useful data
        rho_n = rhoi(1)

        rho_n2 = 0.d0
        e_sum  = 0.d0 
        DO i = 1,nb_bins 
           tmp1 = rhoi(i + 1)
           tmp2 = tmp1*ek(i)
           e_sum  = e_sum + tmp2
           rho_n2 = rho_n2 + tmp1
        ENDDO
         
        e_sum  = e_sum + rho_n*hf_n
        tmp1    = rho_eint(1) - e_sum

        cv_sum = rho_n*cv_tr_n + rho_n2*cv_tr_n2     
        tmp2   = cv_sum + rho_n2*Rn2

        ! Temperature value
        temp = tmp1/tmp2
 
      END SUBROUTINE compute_T_VC_rr

      !----------------------------------------------------!
      ! This subroutine computes the gas temperatures in case of use of the MT_TTint model.
      SUBROUTINE compute_temp_MT_TTint (rhoi, rho_eint, temp) 
 
        USE mod_nitrogen_NASA_initialize_CFD,         ONLY: nb_bins, hf_n, cv_tr_n, cv_tr_n2 

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, rho_eint
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp

        ! N and N2 partial densities and other useful data
        temp = 0.d0
        
        PRINT*,'not implemented yet...'
        PRINT*,'in '          
        STOP

      END SUBROUTINE compute_temp_MT_TTint

      !----------------------------------------------------!
      ! This subroutine computes the species frozen specific heat in case of use of the RVC model.
      SUBROUTINE compute_frozen_cv_RVC (T, cv)

        USE mod_nitrogen_NASA_initialize_CFD,         ONLY: nb_bins, cv_tr_n, cv_tr_n2

        INTEGER :: i

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv

        ! Nitrogen atom N 
        cv(1) = cv_tr_n 

        DO i = 1,nb_bins 
           cv(i + 1) = cv_tr_n2
        ENDDO

      END SUBROUTINE compute_frozen_cv_RVC

      !----------------------------------------------------!
      ! This subroutine computes the species frozen specific heat in case of use of the BRVC.
      SUBROUTINE compute_frozen_cv_BRVC (T, cv)

        USE mod_nitrogen_NASA_initialize_CFD,         ONLY: nb_bins, cv_tr_n, cv_tr_n2, T_min, T_max, T_store,  & 
                                                          & cvint, inv_step 

        INTEGER :: i, left, right
        REAL(KIND=8) :: cvi
        REAL(KIND=8) :: tmp1, tmp2, Tclip

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv

        ! Out of bound check
        Tclip = MAX(T_min,T)
        Tclip = MIN(Tclip,T_max)

        ! Nitrogen atom N 
        cv(1) = cv_tr_n 
        
        ! Nitrogen molecule N2 (value search)
        ! Linear interpolation
        left = INT((Tclip - T_min)*inv_step) + 1
        tmp1 = (Tclip - T_store(left))*inv_step

        left  = (left - 1)*nb_bins 
        right = left + nb_bins

        DO i = 1,nb_bins 
           tmp2 = cvint(left + i)
           cvi  = tmp2 + tmp1*(cvint(right + i) - tmp2) 
           cv(i + 1) = cv_tr_n2 + cvi  
        ENDDO
            
      END SUBROUTINE compute_frozen_cv_BRVC

      !----------------------------------------------------!
      ! This subroutine computes the species frozen specific heat in case of use of the VC model.
      SUBROUTINE compute_frozen_cv_VC (T, cv)

        USE mod_nitrogen_NASA_initialize_CFD,         ONLY: nb_bins, cv_tr_n, cv_tr_n2, T_min, T_max, T_store,  & 
                                                          & cvint, inv_step 

        INTEGER :: i, left, right
        REAL(KIND=8) :: cvi
        REAL(KIND=8) :: tmp1, tmp2, Tclip

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv

        ! Out of bound check
        Tclip = MAX(T_min,T)
        Tclip = MIN(Tclip,T_max)

        ! Nitrogen atom N 
        cv(1) = cv_tr_n 
        
        ! Nitrogen molecule N2 (value search)
        ! Linear interpolation
        left = INT((Tclip - T_min)*inv_step) + 1
        tmp1 = (Tclip - T_store(left))*inv_step

        left  = (left - 1)*nb_bins 
        right = left + nb_bins

        DO i = 1,nb_bins 
           tmp2 = cvint(left + i)
           cvi  = tmp2 + tmp1*(cvint(right + i) - tmp2) 
           cv(i + 1) = cv_tr_n2 + cvi  
        ENDDO 

      END SUBROUTINE compute_frozen_cv_VC

      !----------------------------------------------------!
      ! This subroutine computes the species frozen specific heat in case of use of the VC_rr model.
      SUBROUTINE compute_frozen_cv_VC_rr (T, cv)

        USE mod_nitrogen_NASA_initialize_CFD,         ONLY: nb_bins, cv_tr_n, cv_tr_n2, Rn2  

        INTEGER :: i
        REAL(KIND=8) :: tmp

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv

        ! Nitrogen atom N 
        cv(1) = cv_tr_n 

        ! Nitrogen molecule N2 vibrational levels
        tmp = cv_tr_n2 + Rn2 
        
        DO i = 1,nb_bins        
           cv(i + 1) = tmp
        ENDDO

      END SUBROUTINE compute_frozen_cv_VC_rr

      !----------------------------------------------------!
      ! This subroutine computes the species frozen specific heat in case of use of the MT_TTint model.
      SUBROUTINE compute_frozen_cv_MT_TTint (T, cv)

        USE mod_nitrogen_NASA_initialize_CFD,         ONLY: cv_tr_n, cv_tr_n2

        INTEGER :: i

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cv

        ! Nitrogen atom N 
        cv(1) = cv_tr_n 

        ! Nitrogen molecule N2 
        cv(2) = cv_tr_n2

      END SUBROUTINE compute_frozen_cv_MT_TTint

      !----------------------------------------------------!
      ! This subroutine computes the species energy in case of use the RVC model.
      SUBROUTINE energy_RVC (temp, e)

        USE mod_nitrogen_NASA_initialize_CFD,         ONLY: nb_bins, hf_n, cv_tr_N, cv_tr_N2, ek

        INTEGER :: i
        REAL(KIND=8) :: T

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e

        ! Translational temperature 
        T = temp(1)

        ! Nitrogen atom N 
        e(1) = cv_tr_N*T + hf_n

         ! Nitroge molecule N2 (bin uniform state-to-state collisional model)
        DO i = 1,nb_bins 
           e(i + 1) = cv_tr_N2*T + ek(i)
        ENDDO

      END SUBROUTINE energy_RVC

      !----------------------------------------------------!
      ! This subroutine computes the nonequilibrium specific heat for the N2
      ! molecule when using the MT_TTint collisional model.
      SUBROUTINE get_neq_energy_cv_MT_TTint (T, e_int, cv_int)

        USE mod_nitrogen_NASA_initialize_CFD,         ONLY: nb_bins, ukb, una, mm_N2, EkJ, gk

        INTEGER :: i
        REAL(KIND=8) :: sum1, sum2, sum3, tmp1, tmp2
        REAL(KIND=8) :: ov_kbT, ov_kbT2, Elev, glev

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), INTENT(OUT) :: e_int, cv_int

        ! Useful quantities 
        ov_kbT  = 1.d0/(ukb*T)
        ov_kbT2 = ov_kbT/T

        ! Initialization
        sum1 = 0.d0
        sum2 = 0.d0
        sum3 = 0.d0

        DO i = 1,nb_bins 
           Elev = EkJ(i)
           glev = gk(i)
           tmp1 = glev*DEXP(-Elev*ov_kbT)
           tmp2 = tmp1*Elev
           sum1 = sum1 + tmp1
           sum2 = sum2 + tmp2
           sum3 = sum3 + Elev*tmp2
        ENDDO
        
        tmp1   = una/mm_N2
        tmp2   = sum2/sum1
        e_int  = tmp1*tmp2
        cv_int = tmp1*ov_kbT2*(sum3/sum1 - tmp2**2)

      END SUBROUTINE get_neq_energy_cv_MT_TTint

      !----------------------------------------------------!
      ! This subroutine computes the species energy in case of use of the BRVC model.
      SUBROUTINE energy_BRVC (temp, e)

        USE mod_nitrogen_NASA_initialize_CFD,         ONLY: nb_bins, hf_n, T_min, T_max, inv_step,   & 
                                                          & cv_tr_n, cv_tr_n2, T_store, ek, eint

        INTEGER :: i, left, right
        REAL(KIND=8) :: de, T
        REAL(KIND=8) :: tmp1, tmp2

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e

        ! Translational temperature 
        ! Out of bound check 
        T = MAX(temp(1),T_min)
        T = MIN(T,T_max)

        ! Nitrogen atom N 
        e(1) = cv_tr_n*T + hf_n

        ! Nitrogen molecule N2
        ! Linear interpolation (value search)
        left = INT((T - T_min)*inv_step) + 1
        tmp1 = (T - T_store(left))*inv_step

        left  = (left - 1)*nb_bins 
        right = left + nb_bins

        DO i = 1,nb_bins 
           tmp2 = eint(left + i)
           de   = tmp2 + tmp1*(eint(right + i) - tmp2) 
           e(i + 1) = cv_tr_n2*T + ek(i) + de 
        ENDDO

      END SUBROUTINE energy_BRVC

      !----------------------------------------------------!
      ! This subroutine computes the species energy in case of use of the VC model.
      SUBROUTINE energy_VC (temp, e)

        USE mod_nitrogen_NASA_initialize_CFD,         ONLY: nb_bins, hf_n, T_min, T_max, inv_step, & 
                                                          & cv_tr_n, cv_tr_n2, T_store, ek, eint

        INTEGER :: i, left, right
        REAL(KIND=8) :: de, T
        REAL(KIND=8) :: tmp1, tmp2

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e

        ! Translational temperature 
        ! Out of bound check 
        T = MAX(temp(1),T_min)
        T = MIN(T,T_max)

        ! Nitrogen atom N 
        e(1) = cv_tr_n*T + hf_n

        ! Nitrogen molecule N2
        ! Linear interpolation (value search)
        left = INT((T - T_min)*inv_step) + 1
        tmp1 = (T - T_store(left))*inv_step

        left  = (left - 1)*nb_bins 
        right = left + nb_bins

        DO i = 1,nb_bins 
           tmp2 = eint(left + i)
           de   = tmp2 + tmp1*(eint(right + i) - tmp2) 
           e(i + 1) = cv_tr_n2*T + ek(i) + de 
        ENDDO   

      END SUBROUTINE energy_VC

      !----------------------------------------------------!
      ! This subroutine computes the species energy in case of use of the VC_rr model.
      SUBROUTINE energy_VC_rr (temp, e)

        USE mod_nitrogen_NASA_initialize_CFD,         ONLY: nb_bins, hf_n, cv_tr_n, cv_tr_n2, Rn2, ek

        INTEGER :: i
        REAL(KIND=8) :: tmp, T

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e

        ! Translational temperature 
        T = temp(1)

        ! Nitrogen atom N 
        e(1) = cv_tr_n*T + hf_n

        ! Nitrogen molecule N2 vibrational levels 
        tmp = (cv_tr_n2 + Rn2)*T 
        DO i = 1,nb_bins 
           e(i + 1) = tmp + ek(i)
        ENDDO

      END SUBROUTINE energy_VC_rr

      !----------------------------------------------------!
      ! This subroutine computes the species energy in case of use the MT_TTint model.
      SUBROUTINE energy_MT_TTint (temp, e)

        USE mod_nitrogen_NASA_initialize_CFD,         ONLY: nb_bins, hf_n, cv_tr_N, cv_tr_N2, ek

        INTEGER :: i
        REAL(KIND=8) :: T

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e

        ! Translational temperature 
        T = temp(1)

        ! Nitrogen atom N 
        e(1) = cv_tr_N*T + hf_n

        ! Nitroge molecule N2
        e(2) = cv_tr_N2*T
       
        PRINT*,'Implemenntation to be finished!'
        STOP

      END SUBROUTINE energy_MT_TTint

      !----------------------------------------------------!
      ! This subroutine computes the species internal energies and specific heats at constant volume 
      ! for a given temperature value in case of use of the RVC collisional model
      SUBROUTINE energy_cv_RVC (T, e, cv)
 
         USE mod_nitrogen_NASA_initialize_CFD,         ONLY: nb_bins, hf_n, cv_tr_n, cv_tr_n2, ek

         INTEGER :: i
         REAL(KIND=8) :: e_tr_N2

         REAL(KIND=8), INTENT(IN) :: T
         REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e, cv

         ! Nitrogen atom N 
         e(1)  = cv_tr_N*T + hf_n
         cv(1) = cv_tr_N

         ! Nitrogen molecule N2 
         e_tr_N2 = cv_tr_N2*T
         DO i = 1,nb_bins 
            e(i + 1)  = e_tr_N2 + ek(i)
            cv(i + 1) = cv_tr_N2
         ENDDO

      END SUBROUTINE energy_cv_RVC

      !----------------------------------------------------!
      ! This subroutine computes the species internal energies and specific heats at constant volume 
      ! for a given temperature value in case of use of the bin BRVC model.
      SUBROUTINE energy_cv_BRVC (T, e, cv)
 
         USE mod_nitrogen_NASA_initialize_CFD,         ONLY: nb_bins, hf_n, T_min, T_max, inv_step,  & 
                                                           & cv_tr_n, cv_tr_n2, Rn2, T_store,        & 
                                                           & ek, eint, cvint

         INTEGER :: i, left, right
         REAL(KIND=8) :: cvi, de, e_tr_N2
         REAL(KIND=8) :: tmp1, tmp2, tmp3
         REAL(KIND=8) :: Tclip

         REAL(KIND=8), INTENT(IN) :: T
         REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e, cv

         ! Out of bound check
         Tclip = MAX(T_min,T)
         Tclip = MIN(Tclip,T_max)

         ! Nitrogen atom N 
         e(1)  = cv_tr_N*Tclip + hf_n
         cv(1) = cv_tr_N

         ! Nitrogen molecule N2
         ! Linear interpolation (value search)
         left = INT((Tclip - T_min)*inv_step) + 1
         tmp1 = (Tclip - T_store(left))*inv_step

         left  = (left - 1)*nb_bins 
         right = left + nb_bins

         e_tr_N2 = cv_tr_N2*Tclip
         DO i = 1,nb_bins 
            tmp2 = eint(left + i)
            tmp3 = cvint(left + i)
            de   = tmp2 + tmp1*(eint(right + i) - tmp2)
            cvi  = tmp3 + tmp1*(cvint(right + i) - tmp3) 
            e(i + 1)  = e_tr_N2 + ek(i) + de
            cv(i + 1) = cv_tr_N2 + cvi
         ENDDO

      END SUBROUTINE energy_cv_BRVC

      !----------------------------------------------------!
      ! This subroutine computes the species internal energies and specific heats at constant volume 
      ! for a given temperature value in case of use of the VC model.
      SUBROUTINE energy_cv_VC (T, e, cv)
 
         USE mod_nitrogen_NASA_initialize_CFD,         ONLY: nb_bins, hf_n, T_min, T_max, inv_step, & 
                                                           & cv_tr_n, cv_tr_n2, T_store, ek, eint, cvint

         INTEGER :: i, left, right
         REAL(KIND=8) :: cvi, de, e_tr_N2
         REAL(KIND=8) :: tmp1, tmp2, tmp3
         REAL(KIND=8) :: Tclip

         REAL(KIND=8), INTENT(IN) :: T
         REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e, cv

         ! Out of bound check
         Tclip = MAX(T_min,T)
         Tclip = MIN(Tclip,T_max)

         ! Nitrogen atom N 
         e(1)  = cv_tr_N*Tclip + hf_n
         cv(1) = cv_tr_N

         ! Nitrogen molecule N2
         ! Linear interpolation (value search)
         left = INT((Tclip - T_min)*inv_step) + 1
         tmp1 = (Tclip - T_store(left))*inv_step

         left  = (left - 1)*nb_bins 
         right = left + nb_bins

         e_tr_N2 = cv_tr_N2*Tclip
         DO i = 1,nb_bins 
            tmp2 = eint(left + i)
            tmp3 = cvint(left + i)
            de   = tmp2 + tmp1*(eint(right + i) - tmp2)
            cvi  = tmp3 + tmp1*(cvint(right + i) - tmp3) 
            e(i + 1)  = e_tr_N2 + ek(i) + de
            cv(i + 1) = cv_tr_N2 + cvi
         ENDDO

      END SUBROUTINE energy_cv_VC

      !----------------------------------------------------!
      ! This subroutine computes the species internal energies and specific heats at constant volume 
      ! for a given temperature value in case of use of the VC_rr.
      SUBROUTINE energy_cv_VC_rr (T, e, cv)
 
         USE mod_nitrogen_NASA_initialize_CFD,         ONLY: nb_bins, hf_n, cv_tr_n, cv_tr_n2, Rn2, ek

         INTEGER :: i, left, right
         REAL(KIND=8) :: a, b, de, cvint
         REAL(KIND=8) :: tmp1, tmp2

         REAL(KIND=8), INTENT(IN) :: T
         REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e, cv

         ! Nitrogen atom N 
         e(1)  = cv_tr_N*T + hf_n
         cv(1) = cv_tr_N

         ! Nitrogen molecule N2 vibrational levels 
         tmp1 = cv_tr_N2 + Rn2
         tmp2 = tmp1*T
         DO i = 1,nb_bins 
            e(i + 1)  = tmp2 + ek(i)
            cv(i + 1) = tmp1
         ENDDO

      END SUBROUTINE energy_cv_VC_rr

      !----------------------------------------------------!
      ! This subroutine computes the species internal energies and specific heats at constant volume 
      ! for a given temperature value in case of use of the MT_TTint collisional model
      SUBROUTINE energy_cv_MT_TTint (T, e, cv)
 
         USE mod_nitrogen_NASA_initialize_CFD,         ONLY: hf_n, cv_tr_n, cv_tr_n2

         INTEGER :: i

         REAL(KIND=8), INTENT(IN) :: T
         REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: e, cv

         ! Nitrogen atom N 
         e(1)  = cv_tr_N*T + hf_n
         cv(1) = cv_tr_N

         ! Nitrogen molecule N2 
         e(2)  = cv_tr_N2*T
         cv(2) = cv_tr_N2

      END SUBROUTINE energy_cv_MT_TTint

      !------------------------------------------------------!
      ! This subroutine computes the nonequilibrium post-shock conditions
      SUBROUTINE post_shock_neq (p1, u1, T1, p2, u2, T2, yi, xN)

        USE mod_nitrogen_NASA_initialize_CFD,         ONLY: nb_ns, nb_temp, urg, mm_N, mm_N2, solver, model
        USE mod_Nitrogen_NASA_CFD_eq,                 ONLY: eq_composition_bins
        USE mod_function_pointer_NASA,                ONLY: get_species_energy_cv_NASA, & 
                                                          & get_species_energy_NASA

        INTEGER :: is, length
        REAL(KIND=8), PARAMETER :: tol = 1.d-8
        REAL(KIND=8) :: f, fp, resR, resT, ratio, rhs, ratio_old, T_old
        REAL(KIND=8) :: mass, R
        REAL(KIND=8) :: tmp1, tmp2
        REAL(KIND=8) :: h1, h2, rho1, rho2, c1, g, gp1, M1, M1s, m_dot, xN2
        REAL(KIND=8), DIMENSION(nb_ns) :: mm, Ri, cvi, ei, rhoi
        REAL(KIND=8), DIMENSION(nb_temp) :: temp
        REAL(KIND=8), DIMENSION(3) :: left, right, res

        REAL(KIND=8), INTENT(IN) :: p1, u1, T1
        REAL(KIND=8), INTENT(IN), OPTIONAL :: xN
        REAL(KIND=8), INTENT(OUT) :: p2, u2, T2
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: yi

        ! Useful data
        length = LEN_TRIM(solver)

        Ri(1) = urg/mm_N
        mm(1) = mm_N

        DO is = 2,nb_ns
           Ri(is) = urg/mm_N2
           mm(is) = mm_N2
        ENDDO

        ! Species mass fractions
        CALL eq_composition_bins (p1, T1, rhoi, xN)

        rho1 = 0.d0
        DO is = 1,nb_ns 
           rho1 = rho1 + rhoi(is)
        ENDDO

        DO is = 1,nb_ns 
           yi(is) = rhoi(is)/rho1
        ENDDO

        ! Mixture molar mass 
        mass = 0.d0
        DO is = 1,nb_ns 
           mass = mass + yi(is)/mm(is)
        ENDDO
        mass = 1.d0/mass

        R = urg/mass

        ! Computation of post-chock conditions
        SELECT CASE(model)

          ! Application of simple Rankine-Hugoniot relations
          CASE ('RVC','MT_TTint','MT_TrTv')

            ! Specific heat ratio (monoatomic gas-like case)
            g = 5.d0/3.d0
            gp1 = g + 1.d0
      
            c1  = SQRT(g*R*T1)
            M1  = u1/c1
            M1s = M1**2

            ! Pressure, velocity and temperature after the shock
            p2 = p1*(2.d0*g*M1s - g + 1.d0)/gp1
            u2 = u1 - c1*2.d0/gp1*(M1 - 1.d0/M1)
            T2 = T1*(2.d0*g*M1s - g + 1.d0)*(g - 1.d0 + 2.d0/M1s)/(gp1**2) 
          
          ! Application of simple Rankine-Hugoniot relations
          CASE ('VC_rr')

            ! Specific heat ratio (diatomic gas-like case)
            IF (PRESENT(xN)) THEN

               xN2 =  1.d0 - xN
               g   = (7.d0*xN2 + 5.d0*xN)/(5.d0*xN2 + 3.d0*xN)  

            ELSE 

               g = 7.d0/5.d0

            ENDIF
            
            gp1 = g + 1.d0
      
            c1  = SQRT(g*R*T1)
            M1  = u1/c1
            M1s = M1**2

            ! Pressure, velocity and temperature after the shock
            p2 = p1*(2.d0*g*M1s - g + 1.d0)/gp1
            u2 = u1 - c1*2.d0/gp1*(M1 - 1.d0/M1)
            T2 = T1*(2.d0*g*M1s - g + 1.d0)*(g - 1.d0 + 2.d0/M1s)/(gp1**2) 

          ! A set of non linear equations must be solved in this case
          CASE('VC','BRVC')

             m_dot = rho1*u1

             h1   = 0.d0
             temp = T1
             CALL get_species_energy_NASA (temp, ei)
             DO is = 1,nb_ns 
                h1 = h1 + yi(is)*(ei(is) + Ri(is)*T1)
             ENDDO

             left(1) = m_dot 
             left(2) = m_dot*u1 + p1 
             left(3) = h1 + 0.5d0*u1**2
       
             ! Guess value for the density ratio 
             ratio = 0.1d0 
             resR  = 1.d0
             T2    = T1 
             DO WHILE (resR.GT.tol) 

                p2   = p1 + m_dot*u1*(1.d0 - ratio)
                rho2 = rho1/ratio
            
                rhs  = h1 + 0.5d0*u1**2*(1.d0 - ratio**2)
             
                ! Newton loop for post-shock temperature
                resT = 1.d0
                DO WHILE (resT.GT.tol) 

                   f  = 0.d0
                   fp = 0.d0

                   CALL get_species_energy_cv_NASA (T2, ei, cvi)
                   DO is = 1,nb_ns
                      tmp1 = yi(is) 
                      tmp2 = Ri(is)
                      f    = f  + tmp1*(ei(is) + tmp2*T2)
                      fp   = fp + tmp1*(cvi(is) + tmp2)
                   ENDDO

                   T_old = T2 

                  ! Post-shock temperature residual
                  f  = f - rhs
                  T2 = T2 - f/fp
                  resT = ABS(T2 - T_old)/T_old

                ENDDO

                ratio_old = ratio
            
                rho2 = p2/(R*T2)
            
                ! Density ratio update and residual
                ratio = rho1/rho2
                resR = ABS(ratio - ratio_old)/ratio

             ENDDO

             u2 = u1*ratio

             ! Residuals
 
             ! Mass momentum and energy flux (post-shock)
             m_dot = rho2*u2
       
             h2   = 0.d0
             temp = T2
             CALL get_species_energy_NASA (temp, ei)
             DO is = 1,nb_ns 
                h2 = h2 + yi(is)*(ei(is) + Ri(is)*T2)
             ENDDO   
   
             right(1) = m_dot
             right(2) = m_dot*u2 + p2 
             right(3) = h2 + 0.5d0*u2**2

             DO is = 1,3 
                res(is) = ABS(right(is) - left(is))/left(is)*100.d0
             ENDDO
         
             WRITE(*,5)solver(1:length),':: Nitrogen NASA -> non equilibrium post-shock conditions'
             PRINT*
             WRITE(*,10)'Residual on mass, momemtum and energy fluxes:'
             PRINT*
             WRITE(*,15)'Mass    ',res(1),' [%]'
             PRINT*
             WRITE(*,15)'Momentum',res(2),' [%]'
             PRINT*
             WRITE(*,15)'Energy  ',res(3),' [%]'
             PRINT*

        END SELECT

5     FORMAT(A,A)
10    FORMAT(A)
15    FORMAT(A,E14.6,A)

      END SUBROUTINE post_shock_neq

      !----------------------------------------------------!
      ! This subroutine computes the internal temperature for the energy level population in case of use of 
      ! bin RVC, BRVC or VC models. For sake of robustness the first few iterations are performed by means 
      ! of the bi-section method. Then, after the residual has dropped some order of magnitudes, 
      ! the faster Newton-Raphson method is applied.  
      SUBROUTINE compute_Tint_NonUnif (ni, temp_in, temp_out)

        USE mod_nitrogen_NASA_initialize_CFD,       ONLY: levels, nb_ns, nb_bins, level_bin, ukb, una, mm_N2, degen, EkJ, delta_ek

        INTEGER :: i, pos
        REAL(KIND=8), PARAMETER :: tol1 = 1.d-3, tol2 = 1.d-8
        REAL(KIND=8) :: res, rhs, tmp1, tmp2, tmp3, tmp2l, tmp3l, tmp2r, tmp3r
        REAL(KIND=8) :: ov_kb_T, ov_kb_Tint, ov_kb_Tl, ov_kb_Tr
        REAL(KIND=8) :: sum1, sum2, sum3, sum1l, sum2l, sum1r, sum2r, sum3r
        REAL(KIND=8) :: f, fl, fr, fp, Qint
        REAL(KIND=8) :: T, Tint, Tl, Tr, Tint_old
        REAL(KIND=8), DIMENSION(nb_bins) :: Qk

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: ni, temp_in
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp_out  

        ! Translational temperature 
        T = temp_in(1) 

        ! Computation of rhs 
        sum1 = 0.d0
        sum2 = 0.d0
        rhs  = 0.d0
        DO i = 1,nb_bins 
           tmp1 = ni(i)
           sum1 = sum1 + tmp1 
           sum2 = sum2 + tmp1*EkJ(i) 
        ENDDO
        rhs = sum2/sum1

        ! Computation of bin internal partition functions
        Qk = 0.d0
        ov_kb_T = 1.d0/(ukb*T)
        DO i = 1,levels
           pos = level_bin(i)
           Qk(pos) = Qk(pos) + degen(i)*DEXP(-delta_ek(i)*ov_kb_T)
        ENDDO

        ! Bi-section method (left and right bounds)
        Tl   = 300.d0 
        Tr   = 30000.d0
        Tint = 0.5d0*(Tl + Tr)        
        
        res = 1.d0
        DO WHILE (res.GT.tol1)

           sum1 = 0.d0
           sum2 = 0.d0

           sum1l = 0.d0
           sum2l = 0.d0

           sum1r = 0.d0
           sum2r = 0.d0

           ov_kb_Tint = 1.d0/(ukb*Tint)
           ov_kb_Tl   = 1.d0/(ukb*Tl)
           ov_kb_Tr   = 1.d0/(ukb*Tr)
           DO i = 1,nb_bins 

              tmp1 = EkJ(i)
              Qint = Qk(i)

              tmp2  = DEXP(-tmp1*ov_kb_Tint)
              tmp2l = DEXP(-tmp1*ov_kb_Tl)
              tmp2r = DEXP(-tmp1*ov_kb_Tr)

              tmp3  = Qint*tmp2
              tmp3l = Qint*tmp2l
              tmp3r = Qint*tmp2r

              ! Function at Tint
              sum1  = sum1 + Qint*tmp2
              sum2  = sum2 + tmp3*tmp1

              ! Function at Tl
              sum1l = sum1l + Qint*tmp2l
              sum2l = sum2l + tmp3l*tmp1

              ! Function at Tr
              sum1r = sum1r + Qint*tmp2r
              sum2r = sum2r + tmp3r*tmp1
              
           ENDDO
           
           Tint_old = Tint 

           f  = sum2/sum1   - rhs 
           fl = sum2l/sum1l - rhs
           fr = sum2r/sum1r - rhs

           ! Solution update
           IF ((f*fl).LT.0.d0) THEN 

              Tr = Tint
              Tint = 0.5d0*(Tl + Tr)

           ELSEIF ((f*fr).LT.0.d0) THEN

              Tl = Tint
              Tint = 0.5d0*(Tl + Tr)

           ENDIF

           ! Residual
           res = ABS(Tint - Tint_old)/Tint
           
        ENDDO
        
        ! Newton-Raphson method 
        res = 1.d0
        DO WHILE (res.GT.tol2)

           sum1 = 0.d0
           sum2 = 0.d0
           sum3 = 0.d0

           ov_kb_Tint = 1.d0/(ukb*Tint)
           DO i = 1,nb_bins

              tmp1 = EkJ(i)
              Qint = Qk(i)

              tmp2 = Qint*DEXP(-tmp1*ov_kb_Tint)
              tmp3 = tmp2*tmp1

              sum1 = sum1 + tmp2
              sum2 = sum2 + tmp3
              sum3 = sum3 + tmp3*tmp1

           ENDDO

           Tint_old = Tint

           tmp1 = sum2/sum1
           f    = tmp1 - rhs 
           fp   = (sum3/sum1 - tmp1**2)*ov_kb_Tint/Tint 

           ! Solution update
           Tint = Tint - f/fp

           ! Residual
           res = ABS(Tint - Tint_old)/Tint

        ENDDO

        ! Temperature output vector
        temp_out(1) = T
        temp_out(2) = Tint

      END SUBROUTINE compute_Tint_NonUnif 

      !----------------------------------------------------!
      ! This subroutine computes the internal temperature for the energy level population in case of use of 
      ! the VC_rr model. 
      SUBROUTINE compute_Tint_VC_rr (ni, temp_in, temp_out)

        USE mod_nitrogen_NASA_initialize_CFD,       ONLY: nb_bins, ukb, EkJ

        INTEGER :: i
        REAL(KIND=8), PARAMETER :: tol1 = 1.d-1, tol2 = 1.d-8
        REAL(KIND=8) :: T, Tint
        REAL(KIND=8) :: ov_kb_T, Tint_old
        REAL(KIND=8) :: sum1l, sum2l, sum1r, sum2r, ov_kb_Tl, ov_kb_Tr, Tl, Tr, fl, fr, & 
                        tmp2l, tmp3l, tmp2r, tmp3r
        REAL(KIND=8) :: sum1, sum2, sum3, tmp1, tmp2, tmp3
        REAL(KIND=8) :: res, rhs, f, fp

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: ni, temp_in
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: temp_out 

        ! Translational temperature 
        T = temp_in(1)

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

        ! Temperature output vector
        temp_out(1) = T
        temp_out(2) = Tint

      END SUBROUTINE compute_Tint_VC_rr

  END MODULE mod_nitrogen_NASA_CFD_prop
!------------------------------------------------------------------------------!
