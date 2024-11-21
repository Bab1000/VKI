!------------------------------------------------------------------------------!
! This module provides subroutines for the computation of source of conservation equations for the N-N2 system
! when using the Bari database.
  MODULE mod_nitrogen_Bari_CFD_source

    IMPLICIT NONE

    REAL(KIND=8), PARAMETER :: Tmin = 300.d0
    REAL(KIND=8), PARAMETER :: Park_fac = 0.23d0
    REAL(KIND=8), PARAMETER :: f1_Shat  = 0.1694d0 
    REAL(KIND=8), PARAMETER :: f2_Shat  = 5.5249d-5
    REAL(KIND=8), PARAMETER :: f3_Shat  = - 5.2927d-9
    REAL(KIND=8), PARAMETER :: f4_Shat  = 2.2742d-13

    CONTAINS 

      !------------------------------------------------------!
      ! This subroutine computes the source term due to collisional processes for the N-N2 system.
      SUBROUTINE source (rhoi, T, omega)
 
        USE mod_nitrogen_Bari_initialize_CFD,    ONLY: nb_bins, nb_ns, nb_eq, gn, ukb, una, Edis, mm_N, mm_N2,    & 
                                                    &  fac_exp, EkJ, a1_Da, a2_Da, a3_Da, a4_Da, a5_Da,           & 
                                                    &  a1_vta, a2_vta, a3_vta, a4_vta, a5_vta,                    & 
                                                    &  a_vtm, b_vtm, c_vtm, vTa, a_vv, b_vv, c_vv,                & 
                                                    &  vTm, vv, da, dm
        USE mod_nitrogen_Bari_CFD_prop,          ONLY: Q_trans, Q_rot  

        INTEGER :: i, i1, im1, ip, ip1, j, di 
        INTEGER :: pos, pos1   
        REAL(KIND=8) :: kf, kb, ov_keq
        REAL(KIND=8) :: Tclip, ds_T, ln_T, ov_lnT, ov_T, ov_T2, ov_T3, ov_T4
        REAL(KIND=8) :: exp_1, exp_2
        REAL(KIND=8) :: xn, xns, xn2, Qn, Qn2, Qrotn2
        REAL(KIND=8) :: tmp1, tmp2, tmp3, tmp4, tmp5
        REAL(KIND=8) :: const1, const2
        REAL(KIND=8) :: Shat_fac
        REAL(KIND=8), DIMENSION(nb_bins) :: xi, ek_const

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: omega

        ! Useful data
        IF (T.LT.Tmin) THEN
           Tclip = Tmin
        ELSE 
           Tclip = T 
        ENDIF

        ds_T   = DSQRT(Tclip)
        ln_T   = DLOG(Tclip)
        ov_lnT = 1.d0/ln_T
        ov_T   = 1.d0/Tclip
        ov_T2  = ov_T**2
        ov_T3  = ov_T2*ov_T
        ov_T4  = ov_T3*ov_T

        ! Species concentrations 
        xn = rhoi(1)/mm_N
        xns = xn*xn
 
        xn2 = 0.d0
        DO i = 1,nb_bins 
           tmp1  = rhoi(i + 1)/mm_N2
           xi(i) = tmp1
           xn2   = xn2 + tmp1
        ENDDO

        ! Internal and translational partition functions
        CALL Q_trans (Tclip, mm_N, Qn)
        CALL Q_trans (Tclip, mm_N2, Qn2) 
        CALL Q_rot   (Tclip, Qrotn2)

        ! Useful constants for the computation of backward rate coefficients (excitation and dissociation)
        tmp1   = 1.d0/(ukb*Tclip)
        const1 = -Edis*tmp1

        DO i = 1,nb_bins
           ek_const(i) = DEXP(EkJ(i)*tmp1)
        ENDDO

        tmp1 = 0.d0
        tmp2 = (una*Qn2*Qrotn2)/((gn*Qn)**2*DEXP(const1))  

        ! Initialization
        omega = 0.d0

        ! Collisional dissociation (N impact)
        ! N2(v) + N = 2N + N,   Da 
        IF ((da.EQV..TRUE.).AND.(dm.EQV..FALSE.)) THEN

           DO i = 1,nb_bins 
 
               kf     = fac_exp*DEXP(a1_Da(i) + a2_Da(i)*ov_T + a3_Da(i)*ov_T2 + a4_Da(i)*ov_T3 + a5_Da(i)*ln_T)
               ov_keq = tmp2/ek_const(i)
               kb     = kf*ov_keq 

               exp_1 = xn*(kf*xi(i) - kb*xns)

               tmp1         = tmp1 + exp_1
               omega(i + 1) = - exp_1*mm_N2 

           ENDDO

           omega(1) = 2.d0*tmp1*mm_N 

        ! Collisional dissociation (N and N2 impact)
        ! N2(v) + N  = 2N + N,    Da
        ! N2(v) + N2 = 2N + N2,   Dm     
        ELSEIF ((da.EQV..TRUE.).AND.(dm.EQV..TRUE.)) THEN 

           ! Shatalov factor
           tmp4 = Tclip**2  
           tmp5 = tmp4*Tclip
           Shat_fac = f1_Shat + f2_Shat*Tclip + f3_Shat*tmp4 + f4_Shat*tmp5

           DO i = 1,nb_bins 
 
              kf     = fac_exp*DEXP(a1_Da(i) + a2_Da(i)*ov_T + a3_Da(i)*ov_T2 + a4_Da(i)*ov_T3 + a5_Da(i)*ln_T)
              ov_keq = tmp2/ek_const(i)
              kb     = kf*ov_keq 

              exp_1 = xn*(kf*xi(i) - kb*xns)
              exp_2 = Shat_fac*exp_1*xn2/xn

              tmp3 = exp_1 + exp_2    
              tmp1 = tmp1 + tmp3
              omega(i + 1) = - tmp3*mm_N2 

           ENDDO

           omega(1) = 2.d0*tmp1*mm_N   

        ENDIF

        ! Collisional de-excitation VTa transfer 
        ! N2(v) + N = N2(v - dv) + N
        IF (vTa.EQV..TRUE.) THEN

           pos1 = 0      
           DO i = 1,nb_bins - 1

              tmp1 = xi(i)
              tmp2 = xn*mm_N2
              tmp3 = ek_const(i)  

              i1 = i + 1

              DO ip = i1,nb_bins

                 di = ip - i
 
                 ip1 = ip + 1
                 IF (di.LE.50) THEN

                    IF (ip.LT.nb_bins) THEN

                       ! First fit (1 <= dv <= 30)
                       IF (di.LE.30) THEN 

                          pos1 = pos1 + 1
                          kf     = fac_exp*DEXP(a1_vta(pos1) + a2_vta(pos1)*ov_T + a3_vta(pos1)*ov_T2 + a4_vta(pos1)*ov_T3 & 
                               & +  a5_vta(pos1)*ln_T)

                       ! Second fit (31 <= dv <= 50)
                       ELSE 

                          pos1 = pos1 + 1                     
                          kf     = fac_exp*DEXP(a1_vta(pos1) + a2_vta(pos1)*ov_T + a3_vta(pos1)*ov_T4 + a4_vta(pos1)*ov_lnT)

                       ENDIF

                    ! De-excitation from the last vibrational level v = 67
                    ELSE 

                       pos1 = pos1 + 1
                       kf     = fac_exp*DEXP(a1_vta(pos1) + a2_vta(pos1)*ov_T + a3_vta(pos1)*ov_T2 + a4_vta(pos1)*ov_T3 & 
                            & + a5_vta(pos1)*ln_T)

                    ENDIF

                    ov_keq = tmp3/ek_const(ip)
                    kb     = kf*ov_keq

                    exp_1 = tmp2*(kf*xi(ip) - kb*tmp1) 
                    omega(ip1)  = omega(ip1) - exp_1 
                    omega(i1)   = omega(i1)  + exp_1   

                 ELSE 

                    pos1 = pos1 + 1 
 
                 ENDIF

              ENDDO

           ENDDO

        ENDIF

        ! Collisional de-excitation VTm transfer 
        ! N2(v) + N2 = N2(v - 1) + N2
        IF (vTm.EQV..TRUE.) THEN 

           tmp1 = xn2*mm_N2

           DO i = 2,nb_bins

              ! Forward and backward rate coefficients
              im1 = i - 1
              kf  = a_vtm(im1)*DEXP(b_vtm(im1)*ln_T - c_vtm(im1)*ov_T) 
              kb  = kf/ek_const(i)*ek_const(im1)

              ! Excitation terms
              exp_1 = tmp1*(kf*xi(i) - kb*xi(im1)) 

              omega(i + 1) = omega(i + 1) - exp_1
              omega(i)     = omega(i)     + exp_1

            ENDDO 

        ENDIF

        ! Collisional excitation VV transfer
        ! N2(v) + N2(w - 1) = N2(v - 1) + N2(w)
        IF (vv.EQV..TRUE.) THEN 

           pos = 0

           ! i (vibrational level v) 
           DO i = 2,nb_bins 

              tmp1 = ek_const(i - 1)/ek_const(i)
              tmp2 = xi(i) 
              tmp3 = xi(i - 1) 

              ! j (vibrational level w)
              DO j = i,nb_bins
         
                 pos = pos + 1

                 ! Forward and backward rate coefficients 
                 kf     = a_vv(pos)*DEXP(b_vv(pos)*ln_T - c_vv(pos)*ov_T) 
                 ov_keq = tmp1*ek_const(j)/ek_const(j - 1)
                 kb     = kf*ov_keq
 
                 ! Excitation terms
                 exp_1 = mm_N2*(kf*tmp2*xi(j - 1) - kb*tmp3*xi(j))
                  
                 omega(i + 1) = omega(i + 1) - exp_1
                 omega(i)     = omega(i)     + exp_1
                 omega(j + 1) = omega(j + 1) + exp_1
                 omega(j)     = omega(j)     - exp_1 

              ENDDO

           ENDDO 

        ENDIF 

      END SUBROUTINE source
  
      !------------------------------------------------------!
      ! This subroutine computes the source term due to collisional processes for the N-N2 system.
      ! The source term Jacobian is also provided in output.
      SUBROUTINE source_Jac (rhoi, T, omega, js_omega)

        USE mod_nitrogen_Bari_initialize_CFD,     ONLY: nb_bins, nb_ns, nb_eq, gn, ukb, una, Edis, mm_N, mm_N2,      & 
                                                          &  fac_exp, fac_keq, EkJ, a1_Da, a2_Da, a3_Da, a4_Da, a5_Da,    & 
                                                          &  a1_vta, a2_vta, a3_vta, a4_vta, a5_vta, a_vtm, b_vtm, c_vtm, &
                                                          &  a_vv, b_vv, c_vv, vTa, vTm, vv, da, dm
        USE mod_nitrogen_Bari_CFD_prop,           ONLY: Q_trans, Q_trans_der, Q_rot_der  

        INTEGER :: i, i1, im1, ip, ip1, j, di
        INTEGER :: pos, pos1, pos2, pos3, pos4, posT
        REAL(KIND=8) :: Ev, Evm1
        REAL(KIND=8) :: fac, fac_der, kf, kb, ov_keq, dkf_dT, dkb_dT, dKeq_dT
        REAL(KIND=8) :: Tclip, ds_T, ln_T, ov_lnT, ov_lnT2, ov_T, ov_T2, ov_T3, ov_T4, ov_T5, ov_kbT2, ov_mmN
        REAL(KIND=8) :: exp_1, exp_2
        REAL(KIND=8) :: xn, xns, xn2, Qn, Qn2, Qrotn2, dQn_dT, dQrotn2_dT
        REAL(KIND=8) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9
        REAL(KIND=8) :: const1, const2, en_const
        REAL(KIND=8) :: fac1, fac2, fac3, fac4, fac5, fac6
        REAL(KIND=8) :: f1, f2, f3
        REAL(KIND=8) :: c2, c3, c4, c5
        REAL(KIND=8) :: Shat_fac, dShat_fac_dT
        REAL(KIND=8), DIMENSION(nb_bins) :: xi, ek_const
        
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: omega, js_omega

        ! Useful data
        IF (T.LT.Tmin)THEN 
           Tclip = Tmin
        ELSE 
           Tclip = T
        ENDIF
        ds_T    = DSQRT(Tclip)
        ln_T    = DLOG(Tclip)
        ov_lnT  = 1.d0/ln_T
        ov_lnT2 = ov_lnT**2
        ov_T    = 1.d0/Tclip
        ov_T2   = ov_T*ov_T
        ov_T3   = ov_T2*ov_T
        ov_T4   = ov_T3*ov_T
        ov_T5   = ov_T4*ov_T
        ov_kbT2 = ov_T2/ukb
        ov_mmN  = 1.d0/mm_N

        posT = nb_ns + 1

        ! Species concentrations
        tmp1 = rhoi(1) 
        xn   = rhoi(1)/mm_N
        xns  = xn*xn

        xn2 = 0.d0 
        DO i = 1,nb_bins 
           tmp1  = rhoi(i + 1)/mm_N2
           xi(i) = tmp1
           xn2   = xn2 + tmp1
        ENDDO

        ! Internal and translational partition functions
        CALL Q_trans_der (Tclip, mm_N, Qn, dQn_dT)
        CALL Q_trans     (Tclip, mm_N2, Qn2) 
        CALL Q_rot_der   (Tclip, Qrotn2, dQrotn2_dT)

        ! Useful constants for the computation of backward rate coefficients (excitation and dissociation)
        tmp1   = 1.d0/(ukb*Tclip)
        const1 = -Edis*tmp1

        en_const = DEXP(const1) 
        DO i = 1,nb_bins
           ek_const(i) = DEXP(EkJ(i)*tmp1)
        ENDDO

        ! Initialization
        omega    = 0.d0
        js_omega = 0.d0

        ! Dissociation (N impact)
        ! N2(v) + N = 3N,    Da
        IF ((da.EQV..TRUE.).AND.(dm.EQV..FALSE.)) THEN

           pos1 = posT

           tmp1 = 0.d0
           tmp2 = (una*Qn2*Qrotn2)/((gn*Qn)**2*DEXP(const1))  
           tmp3 = 0.d0
           tmp4 = 0.d0
           tmp5 = 1.d0/Qrotn2
           tmp6 = dQn_dT - Qn*dQrotn2_dT*tmp5

           DO i = 1,nb_bins 

              ! Common factors
              c2 = a2_Da(i)
              c3 = a3_Da(i)
              c4 = a4_Da(i)           
              c5 = a5_Da(i)

              tmp7 = xi(i)

              ! Forward and backward rate coefficients
              kf     = fac_exp*DEXP(a1_Da(i) + c2*ov_T + c3*ov_T2 + c4*ov_T3 + c5*ln_T)
              ov_keq = tmp2/ek_const(i)
              kb     = kf*ov_keq 

              ! Derivatives of forward and backward reaction reate coefficients with respect to temperature 
              fac_der = - c2*ov_T2 - 2.d0*c3*ov_T3 - 3.d0*ov_T4 + c5*ov_T 
              dKeq_dT = fac_keq*ek_const(i)*en_const*tmp5*(Qn*(Edis - EkJ(i))*ov_kbT2 + tmp6)

              dkf_dT  = kf*fac_der
              dkb_dT  = ov_keq*(dkf_dT - ov_keq*kf*dKeq_dT)

              ! Production terms
              exp_1 = xn*(kf*tmp7 - kb*xns)

              tmp1         = tmp1 + exp_1
              omega(i + 1) = - exp_1*mm_N2 

              ! Jacobian common factors
              fac1 = (kf*tmp7 - 3.d0*kb*xns)/mm_N 
              fac2 = xn*kf/mm_N2 
              fac3 = xn*(dkf_dT*tmp7 - dkb_dT*xns)

              ! Jacobian components
              tmp3 = tmp3 + fac1
              js_omega(i + 1) = 2.d0*fac2*mm_N

              js_omega(pos1 + 1)     = - fac1*mm_N2
              js_omega(pos1 + i + 1) = - fac2*mm_N2

              tmp4 = tmp4 + fac3
              js_omega(pos1 + posT) = - fac3*mm_N2

              ! Index update
              pos1 = pos1 + posT

           ENDDO

           ! Components relative to the equation for N (multiplication by a 2mm_N factor)
           fac = 2.d0*mm_N 

           omega(1)       = tmp1*fac 
           js_omega(1)    = tmp3*fac
           js_omega(posT) = tmp4*fac  

        ! Dissociation (N and N2 impact)
        ! N2(v) + N  = 3N,       Da
        ! N2(v) + N2 = 2N + N,   Dm
        ELSEIF ((da.EQV..TRUE.).AND.(dm.EQV..TRUE.)) THEN

           pos1 = posT

           tmp1 = 0.d0
           tmp2 = (una*Qn2*Qrotn2)/((gn*Qn)**2*DEXP(const1))  
           tmp3 = 0.d0
           tmp4 = 0.d0
           tmp5 = 1.d0/Qrotn2
           tmp6 = dQn_dT - Qn*dQrotn2_dT*tmp5

           ! Shatalov factor and derivative with respect to temperature
           tmp7         = Tclip**2  
           tmp8         = tmp7*Tclip
           Shat_fac     = f1_Shat + f2_Shat*Tclip + f3_Shat*tmp7 + f4_Shat*tmp8 
           dShat_fac_dT = f2_Shat + 2.d0*f3_Shat*Tclip + 3.d0*f4_Shat*tmp7

           DO i = 1,nb_bins 

              ! Common factors
              c2 = a2_Da(i)
              c3 = a3_Da(i)
              c4 = a4_Da(i)           
              c5 = a5_Da(i)

              tmp7 = xi(i)

              ! Forward and backward rate coefficients
              kf     = fac_exp*DEXP(a1_Da(i) + c2*ov_T + c3*ov_T2 + c4*ov_T3 + c5*ln_T)
              ov_keq = tmp2/ek_const(i)
              kb     = kf*ov_keq 

              ! Derivatives of forward and backward reaction reate coefficients with respect to temperature 
              fac_der = - c2*ov_T2 - 2.d0*c3*ov_T3 - 3.d0*ov_T4 + c5*ov_T 
              dKeq_dT = fac_keq*ek_const(i)*en_const*tmp5*(Qn*(Edis - EkJ(i))*ov_kbT2 + tmp6)

              dkf_dT  = kf*fac_der
              dkb_dT  = ov_keq*(dkf_dT - ov_keq*kf*dKeq_dT)

              ! Common factors for production terms and Jacobian
              tmp8 = kf*tmp7 - kb*xns
              tmp9 = dkf_dT*tmp7 - dkb_dT*xns

              ! Production terms
              exp_1 = tmp8*(xn + Shat_fac*xn2)

              tmp1 = tmp1 + exp_1
              omega(i + 1) = - exp_1*mm_N2 

              ! Jacobian common factors
              fac1 = (tmp8 - 2.d0*kb*xns)/mm_N
              fac2 = xn*kf/mm_N2
              fac3 = xn*tmp9         
              fac4 = - 2.d0*kb*xn*xn2*Shat_fac/mm_N
              fac5 = (tmp8 + xn2*kf)*Shat_fac/mm_N2
              fac6 = xn2*tmp9*Shat_fac + dShat_fac_dT*xn2*tmp8
              fac  = tmp8*Shat_fac/mm_N2
 
              ! Jacobian components
              ! Components relative to the equation for N
              f1  = fac1 + fac4
              f2  = fac2 + fac5
              f3  = fac3 + fac6   

              tmp3 = tmp3 + f1
              tmp4 = tmp4 + f3 

              DO j = 1,i - 1
                 js_omega(1 + j) = js_omega(1 + j) + fac
              ENDDO

              js_omega(1 + i) = js_omega(1 + i) + f2

              DO j = i + 1,nb_bins
                 js_omega(1 + j) = js_omega(1 + j) + fac
              ENDDO

              ! Components relative to the equation for N2(i) (vibrational level v = i - 1)
              js_omega(pos1 + 1) = - mm_N2*f1 

              fac = - mm_N2*fac           
              DO j = 1,i - 1
                 js_omega(pos1 + 1 + j) = fac
              ENDDO

              js_omega(pos1 + 1 + i) = - mm_N2*f2
    
              DO j = i + 1,nb_bins
                 js_omega(pos1 + 1 + j) = fac
              ENDDO

              js_omega(pos1 + posT) = - mm_N2*f3

              ! Index update
              pos1 = pos1 + posT

           ENDDO

           ! Components relative to the equation for N (multiplication by a 2mm_N factor)
           fac = 2.d0*mm_N

           omega(1)    = tmp1*fac  
           js_omega(1) = tmp3*fac
            
           DO i = 1,nb_bins
              js_omega(1 + i) = js_omega(1 + i)*fac
           ENDDO
           js_omega(posT) = tmp4*fac  

        ENDIF

        ! Collisional de-excitation VTa transfer 
        ! N2(v) + N = N2(v - dv) + N
        IF (vTa.EQV..TRUE.) THEN

           pos1 = 0     
           pos3 = posT 
           DO i = 1,nb_bins - 1

              ! Common factors
              tmp1 = xi(i)
              tmp2 = xn*mm_N2
              tmp3 = ek_const(i) 
              tmp4 = EkJ(i)*ov_kbT2 

              i1 = i + 1

              pos4 = i1*posT
              DO ip = i1,nb_bins

                 di = ip - i
 
                 ip1 = ip + 1
                 IF (di.LE.50) THEN

                    IF (ip.LT.nb_bins) THEN

                       ! First fit (1 <= dv <= 30)
                       IF (di.LE.30) THEN 

                          pos1 = pos1 + 1

                          ! Common factors
                          c2 = a2_vta(pos1)
                          c3 = a3_vta(pos1)
                          c4 = a4_vta(pos1)           
                          c5 = a5_vta(pos1)

                          kf      = fac_exp*DEXP(a1_vta(pos1) + c2*ov_T + c3*ov_T2 + c4*ov_T3 + c5*ln_T)
                          fac_der = - c2*ov_T2 - 2.d0*c3*ov_T3 - 3.d0*c4*ov_T4 + c5*ov_T 
                          dkf_dT  = fac_der*kf

                       ! Second fit (31 <= dv <= 50)
                       ELSE 

                          pos1 = pos1 + 1                     

                          ! Common factors
                          c2 = a2_vta(pos1)
                          c3 = a3_vta(pos1)
                          c4 = a4_vta(pos1)           

                          kf      = fac_exp*DEXP(a1_vta(pos1) + c2*ov_T + c3*ov_T4 + c4*ov_lnT)
                          fac_der = - c2*ov_T2 - 4.d0*c3*ov_T5 - 3.d0*c4*ov_lnT2*ov_T
                          dkf_dT  = fac_der*kf

                       ENDIF

                    ! De-excitation from the last vibrational level v = 67
                    ELSE 

                       pos1 = pos1 + 1

                       ! Common factors
                       c2 = a2_vta(pos1)
                       c3 = a3_vta(pos1)
                       c4 = a4_vta(pos1)           
                       c5 = a5_vta(pos1)

                       kf      = fac_exp*DEXP(a1_vta(pos1) + c2*ov_T + c3*ov_T2 + c4*ov_T3 + c5*ln_T)
                       fac_der = - c2*ov_T2 - 2.d0*c3*ov_T3 - 3.d0*ov_T4 + c5*ov_T  
                       dkf_dT  = fac_der*kf

                    ENDIF

                    ! Common factors
                    tmp5 = xi(ip)
                    tmp6 = ek_const(ip)

                    ! Excitation terms
                    ov_keq = tmp3/tmp6
                    kb     = kf*ov_keq

                    exp_2 = (kf*tmp5 - kb*tmp1)
                    exp_1 = tmp2*exp_2 
                    omega(ip1)  = omega(ip1) - exp_1 
                    omega(i1)   = omega(i1)  + exp_1

                    dKeq_dT = (tmp4 - EkJ(ip)*ov_kbT2)*tmp6/tmp3
                    dkb_dT  = ov_keq*(dkf_dT - ov_keq*kf*dKeq_dT)

                    ! Jacobian common factors
                    fac1 = xn*kb
                    fac2 = xn*kf
                    fac3 = tmp2*(dkf_dT*tmp5 - dkb_dT*tmp1) 
                    fac4 = mm_N2*exp_2*ov_mmN

                    ! Jacobian components
                    js_omega(pos3 + 1)    = js_omega(pos3 + 1)   + fac4
                    js_omega(pos3 + i1)   = js_omega(pos3 + i1)  - fac1
                    js_omega(pos3 + ip1)  = js_omega(pos3 + ip1) + fac2 
                    js_omega(pos3 + posT) = js_omega(pos3 + posT) + fac3

                    js_omega(pos4 + 1)    = js_omega(pos4 + 1)    - fac4 
                    js_omega(pos4 + i1)   = js_omega(pos4 + i1)   + fac1
                    js_omega(pos4 + ip1)  = js_omega(pos4 + ip1)  - fac2
                    js_omega(pos4 + posT) = js_omega(pos4 + posT) - fac3
           
                    ! Index update 
                    pos4 = pos4 + posT

                 ELSE 

                    ! Index update
                    pos1 = pos1 + 1 

                 ENDIF

            ENDDO

            ! Index update
            pos3 = pos3 + posT

          ENDDO

        ENDIF

        ! Collisional de-excitation VTm transfer 
        ! N2(v) + N2 = N2(v - 1) + N2
        IF (vTm.EQV..TRUE.) THEN 

           tmp1 = mm_N2*xn2
           tmp2 = xn2

           DO i = 2,nb_bins

              ! Forward and backward rate coefficients
              im1 = i - 1
              c2  = b_vtm(im1)
              c3  = c_vtm(im1)

              kf      = a_vtm(im1)*DEXP(c2*ln_T - c3*ov_T) 
              ov_keq  = ek_const(im1)/ek_const(i)
              kb      = kf*ov_keq

              ! Derivatives of forward and backward reaction rate coefficients with respect to temperature 
              fac_der = c2*ov_T + c3*ov_T2
              dkf_dT  = fac_der*kf
              dkb_dT  = ov_keq*(dkf_dT + kf*ov_kbT2*(EkJ(i) - EkJ(im1)))

              ! Common factors
              tmp3  = xi(i)
              tmp4  = xi(im1)
              tmp5  = kf*tmp3 - kb*tmp4   

              ! Excitation terms
              exp_1        = tmp1*tmp5 
              omega(i + 1) = omega(i + 1) - exp_1
              omega(i)     = omega(i)     + exp_1
              
              ! Common factors for Jacobian 
              fac1 = tmp5 + kf*tmp2
              fac2 = tmp5 - kb*tmp2
              fac3 = tmp5
              fac4 = tmp1*(dkf_dT*tmp3 - dkb_dT*tmp4) 

              ! Jacobian components relative to the equation for N2(i - 1) (vibrational level v - 1)
              pos1 = posT*(i - 1)
              DO j = 1,i - 2
                 js_omega(pos1 + 1 + j) = js_omega(pos1 + 1 + j) + fac3
              ENDDO

              js_omega(pos1 + 1 + i - 1) = js_omega(pos1 + 1 + i - 1) + fac2
              js_omega(pos1 + 1 + i)     = js_omega(pos1 + 1 + i)     + fac1

              DO j = i + 1,nb_bins
                 js_omega(pos1 + 1 + j) = js_omega(pos1 + 1 + j) + fac3
              ENDDO
             
              js_omega(pos1 + posT) = js_omega(pos1 + posT) + fac4

              ! Jacobian components relative to the equation for N2(i) (vibrational level v)
              pos1 = pos1 + posT
              DO j = 1,i - 2
                 js_omega(pos1 + 1 + j) = js_omega(pos1 + 1 + j) - fac3
              ENDDO

              js_omega(pos1 + 1 + i - 1) = js_omega(pos1 + 1 + i - 1) - fac2
              js_omega(pos1 + 1 + i)     = js_omega(pos1 + 1 + i)     - fac1

              DO j = i + 1,nb_bins
                 js_omega(pos1 + 1 + j) = js_omega(pos1 + 1 + j) - fac3
              ENDDO 

              js_omega(pos1 + posT) = js_omega(pos1 + posT) - fac4

            ENDDO 

        ENDIF
        
        ! Collisional excitation VV transfer
        ! N2(v) + N2(w - 1) = N2(v - 1) + N2(w)
        IF (vv.EQV..TRUE.) THEN 

           pos = 0

           ! i (vibrational level v) 
           DO i = 2,nb_bins 

              tmp1 = ek_const(i - 1)/ek_const(i)
              tmp2 = xi(i) 
              tmp3 = xi(i - 1) 

              Ev   = EkJ(i)
              Evm1 = EkJ(i - 1) 

              ! j (vibrational level w)
              DO j = i,nb_bins
               
                 pos  = pos + 1

                 tmp4 = xi(j - 1)
                 tmp5 = xi(j)

                 ! Forward and backward rate coefficients 
                 c2  = b_vv(pos)
                 c3  = c_vv(pos)

                 kf      = a_vv(pos)*DEXP(c2*ln_T - c3*ov_T) 
                 ov_keq  = tmp1*ek_const(j)/ek_const(j - 1)
                 kb      = kf*ov_keq
                 
                 ! Derivatives of forward and backward reaction rate coefficients with respect to temperature 
                 fac_der = c2*ov_T + c3*ov_T2
                 dkf_dT  = fac_der*kf
                 dkb_dT  = ov_keq*(dkf_dT + kf*ov_kbT2*(Ev + EkJ(j - 1) - Evm1 - EkJ(j))) 

                 ! Excitation terms
                 exp_1        = mm_N2*(kf*tmp2*tmp4 - kb*tmp3*tmp5)
                 omega(i + 1) = omega(i + 1) - exp_1    ! v
                 omega(i)     = omega(i)     + exp_1    ! v - 1
                 omega(j + 1) = omega(j + 1) + exp_1    ! w
                 omega(j)     = omega(j)     - exp_1    ! w - 1
                 
                 ! Jacobian common factors
                 fac1 = - kb*tmp5
                 fac2 = kf*tmp4
                 fac3 = kf*tmp2
                 fac4 = - kb*tmp3
                 fac5 = mm_N2*(dkf_dT*tmp2*tmp4 - dkb_dT*tmp3*tmp5)

                 ! Components relative to the equation for N2(i - 1) (vibrational level v - 1)
                 pos1 = posT*(i - 1)

                 js_omega(pos1 + i)     = js_omega(pos1 + i)     + fac1    ! v - 1
                 js_omega(pos1 + 1 + i) = js_omega(pos1 + 1 + i) + fac2    ! v 
                 js_omega(pos1 + j)     = js_omega(pos1 + j)     + fac3    ! w - 1
                 js_omega(pos1 + 1 + j) = js_omega(pos1 + 1 + j) + fac4    ! w
                 js_omega(pos1 + posT)  = js_omega(pos1 + posT)  + fac5    ! T

                 ! Components relative to the equation for N2(i) (vibrational level v)
                 pos1 = pos1 + posT

                 js_omega(pos1 + i)     = js_omega(pos1 + i)      - fac1   ! v - 1
                 js_omega(pos1 + 1 + i) = js_omega(pos1 + 1 + i)  - fac2   ! v 
                 js_omega(pos1 + j)     = js_omega(pos1 + j)      - fac3   ! w - 1
                 js_omega(pos1 + 1 + j) = js_omega(pos1 + 1 + j)  - fac4   ! w
                 js_omega(pos1 + posT)  = js_omega(pos1 + posT)   - fac5   ! T

                 ! Components relative to the equation for N2(j - 1) (vibrational level w - 1)
                 pos1 = posT*(j - 1)

                 js_omega(pos1 + i)     = js_omega(pos1 + i)     - fac1    ! v - 1
                 js_omega(pos1 + 1 + i) = js_omega(pos1 + 1 + i) - fac2    ! v 
                 js_omega(pos1 + j)     = js_omega(pos1 + j)     - fac3    ! w - 1
                 js_omega(pos1 + 1 + j) = js_omega(pos1 + 1 + j) - fac4    ! w
                 js_omega(pos1 + posT)  = js_omega(pos1 + posT)  - fac5    ! T

                 ! Components relative to the equation for N2(j) (vibrational level w)
                 pos1 = pos1 + posT

                 js_omega(pos1 + i)     = js_omega(pos1 + i)     + fac1    ! v - 1
                 js_omega(pos1 + 1 + i) = js_omega(pos1 + 1 + i) + fac2    ! v 
                 js_omega(pos1 + j)     = js_omega(pos1 + j)     + fac3    ! w - 1
                 js_omega(pos1 + 1 + j) = js_omega(pos1 + 1 + j) + fac4    ! w
                 js_omega(pos1 + posT)  = js_omega(pos1 + posT)  + fac5    ! T

              ENDDO

           ENDDO 

        ENDIF  

      END SUBROUTINE source_Jac 

  END MODULE mod_nitrogen_Bari_CFD_source
!------------------------------------------------------------------------------!
