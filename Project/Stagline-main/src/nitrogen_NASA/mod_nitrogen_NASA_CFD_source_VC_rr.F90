!------------------------------------------------------------------------------!
! This module provides subroutines for the computation of source of conservation equations for the N-N2 system
! when using the NASA Ames ab-initio database.
  MODULE mod_nitrogen_NASA_CFD_source_VC_rr

#define pre_KACHE

    IMPLICIT NONE

    REAL(KIND=8), PARAMETER :: Park_fac = 0.23d0
    REAL(KIND=8), PARAMETER :: f1_Shat  = 0.1694d0 
    REAL(KIND=8), PARAMETER :: f2_Shat  = 5.5249d-5
    REAL(KIND=8), PARAMETER :: f3_Shat  = - 5.2927d-9
    REAL(KIND=8), PARAMETER :: f4_Shat  = 2.2742d-13

    ! Subroutines for source term computation
    CONTAINS

      !------------------------------------------------------!
      ! This subroutine computes the source term due to collisional processes for the N-N2 system.
      ! This subroutine is called when using the VC_rr model. This subroutine ca be used also for computing
      ! the source term Jacobian by means of finite differences (in order to speed up its computation    
      ! pre-kached values of rate coeffiecients are used) 
      SUBROUTINE source_VC_rr (rhoi, temp, omega)

        USE mod_nitrogen_NASA_initialize_CFD,        ONLY: nb_ns, nb_bins, nb_eq, mm_N, mm_N2, Tmin_dis, Tmin_exc,  &      
                                                         & ad_kf, bd_kf, cd_kf, ae_kf, be_kf, ce_kf, a_vtm, b_vtm,  &
                                                         & c_vtm, a_vv, b_vv, c_vv, gn, EkJ, Edis, una, ukb, vTa,   & 
                                                         & vTm, vv, da, dm, theta_rot, nb_Da_proc, nb_vTa_proc,     & 
                                                         & nb_vTm_proc, nb_vv_proc
        USE mod_nitrogen_NASA_CFD_prop,              ONLY: Q_trans 

        INTEGER :: i, i1, im1, ip1, ip, j
        INTEGER :: pos, pos1, pos2
        REAL(KIND=8) :: kf, kb, ov_keq
        REAL(KIND=8) :: T, Tclip, ds_T, ln_T, ln_v, ov_T, exp_1, exp_2, exp_3
        REAL(KIND=8) :: xn, xns, xn2, Qn, Qn2, Qrot
        REAL(KIND=8) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
        REAL(KIND=8) :: const1, const2
        REAL(KIND=8) :: Shat_fac 
        REAL(KIND=8), DIMENSION(nb_bins) :: xi, ek_const

#ifdef pre_KACHE
        INTEGER, SAVE :: first_entry = 0
        REAL(KIND=8), SAVE :: T_old = 0.d0
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: ek_const_KACHE
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: kDaf_KACHE
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: kvTaf_KACHE
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: kvTmf_KACHE
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: kvvf_KACHE  
#endif

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: omega

        ! Temperature 
        T = temp(1)

        ! Fix applied in order to avoid numerical problems
        IF (T.LT.Tmin_dis) THEN 
           Tclip = Tmin_dis
        ELSE 
           Tclip = T
        ENDIF

        ds_T = DSQRT(Tclip)
        ln_T = DLOG(Tclip)
        ov_T = 1.d0/Tclip

#ifdef pre_KACHE
        IF (first_entry.EQ.0) THEN
           ALLOCATE(ek_const_KACHE(nb_bins)) 
           ALLOCATE(kDaf_KACHE(nb_Da_proc))
           ALLOCATE(kvTaf_KACHE(nb_vTa_proc))
           ALLOCATE(kvTmf_KACHE(nb_vTm_proc))
           ALLOCATE(kvvf_KACHE(nb_vv_proc))
        ENDIF
        first_entry = 1
#endif

        ! Computing pre-KACHE values 
#ifdef pre_KACHE
        IF (Tclip.NE.T_old) THEN
           tmp1 = 1.d0/(ukb*Tclip)
           DO i = 1,nb_bins
              ek_const_KACHE(i) = DEXP(EkJ(i)*tmp1)
           ENDDO
           
           ! Da process
           IF (Da.EQV..TRUE.) THEN
              DO i = 1,nb_Da_proc
                 kDaf_KACHE(i) = ad_kf(i)*DEXP(bd_kf(i)*ln_T - cd_kf(i)*ov_T)
              ENDDO
           ENDIF

           ! vTa process 
           IF (vTa.EQV..TRUE.) THEN
              DO i = 1,nb_vTa_proc
                 kvTaf_KACHE(i) = ae_kf(i)*DEXP(be_kf(i)*ln_T - ce_kf(i)*ov_T)
              ENDDO
           ENDIF

           ! vTm process
           IF (vTm.EQV..TRUE.) THEN
              DO i = 1,nb_vTm_proc
                 kvTmf_KACHE(i) = a_vtm(i)*DEXP(b_vtm(i)*ln_T - c_vtm(i)*ov_T)
              ENDDO
           ENDIF

           ! vv process
           IF (vv.EQV..TRUE.) THEN
              DO i = 1,nb_vv_proc
                 kvvf_KACHE(i) = a_vv(i)*DEXP(b_vv(i)*ln_T - c_vv(i)*ov_T)
              ENDDO
           ENDIF           

           T_old = Tclip
              
        ENDIF
#endif


        ! Species concentrations 
        xn = rhoi(1)/mm_N
        xns = xn*xn
 
        xn2 = 0.d0
        DO i = 1,nb_bins 
           tmp1  = rhoi(i + 1)/mm_N2
           xi(i) = tmp1 
           xn2   = xn2 + tmp1
        ENDDO

        ! Internal and translational partition function
        CALL Q_trans (Tclip, mm_N, Qn)
        CALL Q_trans (Tclip, mm_N2, Qn2)

        Qrot = 0.5d0*Tclip/theta_rot

        ! Useful constants for the computation of backward rate coefficients (excitation and dissociation)
        const1 = -Edis/(ukb*Tclip)
        const2 = 1.d0/(ukb*Tclip)

        DO i = 1,nb_bins
#ifdef pre_KACHE 
           ek_const(i) = ek_const_KACHE(i)
#else 
           ek_const(i) = DEXP(EkJ(i)*const2)
#endif
        ENDDO

        ! Initialization 
        omega = 0.d0

        ! Useful factors 
        tmp1 = 0.d0
        tmp2 = (una*Qn2*Qrot)/((gn*Qn)**2*DEXP(const1)) 

        ! Dissociation 
        IF ((da.EQV..TRUE.).AND.(dm.EQV..FALSE.)) THEN

           ! Collisional dissociation (N impact) 
           ! N2(v) + N  = 3N
           DO i = 1,nb_bins 

              ! Forward and backward rate coefficients
#ifdef pre_KACHE 
              kf = kDaf_KACHE(i)
#else 
              kf = ad_kf(i)*DEXP(bd_kf(i)*ln_T - cd_kf(i)*ov_T)
#endif
              ov_keq = tmp2/ek_const(i)
              kb     = kf*ov_keq 

              ! Production terms
              exp_1 = xn*(kf*xi(i) - kb*xns)

              tmp1         = tmp1 + exp_1
              omega(i + 1) = - exp_1*mm_N2

           ENDDO

           omega(1) = 2.d0*tmp1*mm_N 

        ELSEIF ((da.EQV..TRUE.).AND.(dm.EQV..TRUE.)) THEN

           ! Collisional dissociation (N and N2 impact) 
           ! N2(v) + N  = 3N
           ! N2(v) + N2 = 2N + N2 (Park correction factor is applied)

           ! Shatalov factor
           tmp4 = Tclip**2  
           tmp5 = tmp4*Tclip
           Shat_fac = f1_Shat + f2_Shat*Tclip + f3_Shat*tmp4 + f4_Shat*tmp5 

           DO i = 1,nb_bins 

              tmp3 = xi(i)

              ! Forward and backward rate coefficients
#ifdef pre_KACHE
              kf = kDaf_KACHE(i) 
#else 
              kf = ad_kf(i)*DEXP(bd_kf(i)*ln_T - cd_kf(i)*ov_T)
#endif
              ov_keq = tmp2/ek_const(i)
              kb     = kf*ov_keq 

              ! Production terms
              exp_1 = xn*(kf*tmp3 - kb*xns)
              exp_2 = Shat_fac*xn2*exp_1/xn

              tmp3 = (exp_1 + exp_2)
              tmp1 = tmp1 + tmp3
              omega(i + 1) = - tmp3*mm_N2

           ENDDO 
           
           omega(1) = 2.d0*tmp1*mm_N 

        ENDIF
        
        ! Collisional excitation (N impact)
        ! N2(v) + N = N(v + dv) + N
        IF (vTa.EQV..TRUE.) THEN  

           pos1 = 0
           DO i = 1,nb_bins - 1

              i1   = i + 1
              tmp1 = xi(i)
              tmp2 = xn*mm_N2
              tmp3 = ek_const(i)              

              DO ip = i1,nb_bins 
  
                 ip1  = ip + 1
                 pos2 = pos1 + ip - i
   
                 ! Forward and backward rate coefficients
#ifdef pre_KACHE
                 kf = kvTaf_KACHE(pos2)
#else 
                 kf = ae_kf(pos2)*DEXP(be_kf(pos2)*ln_T - ce_kf(pos2)*ov_T) 
#endif
                 ov_keq = ek_const(ip)/(tmp3)                 

                 kb  = kf*ov_keq

                 ! Excitation terms
                 exp_1 = tmp2*(kf*tmp1 - kb*xi(ip))
                
                 omega(i1)  = omega(i1)  - exp_1
                 omega(ip1) = omega(ip1) + exp_1

              ENDDO

              pos1 = pos1 + nb_bins - i

           ENDDO

        ENDIF

        ! Collisional excitation (translational-internal exchange) (N2 impact)
        ! N2(v) + N2 = N2(v - 1) + N2
        IF (vTm.EQV..TRUE.) THEN 

           tmp1 = xn2*mm_N2

           DO i = 2,nb_bins

              ! Forward and backward rate coefficients
              im1 = i - 1
#ifdef pre_KACHE
              kf = kvTmf_KACHE(im1)
#else
              kf = a_vtm(im1)*DEXP(b_vtm(im1)*ln_T - c_vtm(im1)*ov_T) 
#endif
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
#ifdef pre_KACHE
                 kf = kvvf_KACHE(pos)
#else 
                 kf = a_vv(pos)*DEXP(b_vv(pos)*ln_T - c_vv(pos)*ov_T) 
#endif
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

      END SUBROUTINE source_VC_rr

      !----------------------------------------------------!
      ! This subroutine computes the source term due to collisional processes for the N-N2 system.
      ! This subroutine is called when using the VC_rr model.
      ! The source term Jacobian is also provided in output.
      SUBROUTINE source_VC_rr_Jac (rhoi, temp, omega, js_omega)

        USE mod_nitrogen_NASA_initialize_CFD,        ONLY: nb_ns, nb_bins, nb_eq, mm_N, mm_N2, Tmin_dis, Tmin_exc,   & 
                                                         & ad_kf, bd_kf, cd_kf, ae_kf, be_kf, ce_kf, a_vtm, b_vtm,   & 
                                                         & c_vtm, a_vv, b_vv, c_vv, gn, EkJ, Edis, una, ukb,         & 
                                                         & fac_keq, theta_rot, vTa, vTm, vv, da, dm 
        USE mod_nitrogen_NASA_CFD_prop,              ONLY: Q_trans, Q_trans_der

        INTEGER :: i, im1, j, ip, i1, ip1
        INTEGER :: pos, pos1, pos2, pos3, pos4, posT
        REAL(KIND=8) :: Ev, Evm1
        REAL(KIND=8) :: kf, kb, ov_keq
        REAL(KIND=8) :: fac, fac_der, dkf_dT, dkb_dT, dKeq_dT, dkdT, Qrot, dQrot_dT
        REAL(KIND=8) :: T, Tclip, ds_T, ln_T, ov_T, ov_T2, ov_kbT2, ov_mmN, exp_1, exp_2
        REAL(KIND=8) :: xn, xns, xn2, Qn, Qn2, dQn_dT
        REAL(KIND=8) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9
        REAL(KIND=8) :: const1, const2, en_const
        REAL(KIND=8) :: fac1, fac2, fac3, fac4, fac5, fac6
        REAL(KIND=8) :: f1, f2, f3
        REAL(KIND=8) :: c2, c3
        REAL(KIND=8) :: Shat_fac, dShat_fac_dT
        REAL(KIND=8), DIMENSION(nb_bins) :: xi, ek_const
       
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: omega, js_omega

        ! Temperature 
        T = temp(1)

        ! Fix applied in order to avoid numerical problems
        IF (T.LT.Tmin_dis) THEN 
           Tclip = Tmin_dis
        ELSE 
           Tclip = T
        ENDIF

        ds_T    = DSQRT(Tclip)
        ln_T    = DLOG(Tclip)
        ov_T    = 1.d0/Tclip
        ov_T2   = ov_T**2
        ov_kbT2 = 1.d0/(ukb*Tclip**2)
        ov_mmN  = 1.d0/mm_N 

        posT = nb_ns + 1

        ! Species concentrations 
        xn   = rhoi(1)/mm_N
        xns  = xn*xn

        xn2 = 0.d0
        DO i = 1,nb_bins 
           tmp1  = rhoi(i + 1)/mm_N2
           xi(i) = tmp1
           xn2   = xn2 + tmp1
        ENDDO

        ! Internal and translational partition function (and derivatives with respect to temperature)
        CALL Q_trans_der (Tclip, mm_N, Qn, dQn_dT)
        CALL Q_trans (Tclip, mm_N2, Qn2)

        tmp1 = 0.5d0/theta_rot
        Qrot = tmp1*Tclip
        dQrot_dT = tmp1         

        ! Useful constants for the computation of backward rate coefficients (excitation and dissociation)
        const1 = -Edis/(ukb*Tclip)
        const2 = 1.d0/(ukb*Tclip)

        en_const = DEXP(const1) 
        DO i = 1,nb_bins
           ek_const(i) = DEXP(EkJ(i)*const2)
        ENDDO

        ! Initializaton
        omega    = 0.d0
        js_omega = 0.d0
        
        ! Dissociation (N impact)
        ! N2(v) + N = 3N
        IF ((da.EQV..TRUE.).AND.(dm.EQV..FALSE.)) THEN

           pos1 = posT

           tmp1 = 0.d0
           tmp2 = (una*Qn2*Qrot)/((gn*Qn)**2*en_const)
           tmp3 = 0.d0
           tmp4 = 0.d0
           tmp5 = 1.d0/Qrot
           tmp6 = dQn_dT - Qn*dQrot_dT*tmp5

           DO i = 1,nb_bins   

              ! Common factors
              c2   = bd_kf(i)
              c3   = cd_kf(i)
              tmp7 = xi(i)

              ! Forward and backward rate coefficients
              kf = ad_kf(i)*DEXP(c2*ln_T - c3*ov_T)

              ov_keq = tmp2/ek_const(i)
              kb     = kf*ov_keq 

              ! Derivatives of forward and backward reaction reate coefficients with respect to temperature
              fac_der = c2*ov_T + c3*ov_T2
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
          
           omega(1) = tmp1*fac
        
           js_omega(1)    = tmp3*fac
           js_omega(posT) = tmp4*fac         

        ! Dissociation (N and N2 impact)
        ! N2(v) + N  = 3N
        ! N2(v) + N2 = 2N + N        
        ELSEIF ((da.EQV..TRUE.).AND.(dm.EQV..TRUE.)) THEN

           pos1 = posT

           tmp1 = 0.d0
           tmp2 = (una*Qn2*Qrot)/((gn*Qn)**2*en_const)
           tmp3 = 0.d0
           tmp4 = 0.d0
           tmp5 = 1.d0/Qrot
           tmp6 = dQn_dT - Qn*dQrot_dT*tmp5

           ! Shatalov factor and derivative with respect to temperature
           tmp7         = Tclip**2  
           tmp8         = tmp7*Tclip
           Shat_fac     = f1_Shat + f2_Shat*Tclip + f3_Shat*tmp7 + f4_Shat*tmp8 
           dShat_fac_dT = f2_Shat + 2.d0*f3_Shat*Tclip + 3.d0*f4_Shat*tmp7 

           DO i = 1,nb_bins  

              ! Common factors
              c2   = bd_kf(i)
              c3   = cd_kf(i)
              tmp7 = xi(i)

              ! Forward and backward rate coefficients
              kf   = ad_kf(i)*DEXP(c2*ln_T - c3*ov_T)

              ov_keq = tmp2/ek_const(i)
              kb     = kf*ov_keq               

              ! Derivatives of forward and backward reaction reate coefficients with respect to temperature
              fac_der = c2*ov_T + c3*ov_T2
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

        ! vTa transfer 
        ! N2(v) + N = N2(v') + N
        IF (vTa.EQV..TRUE.) THEN

           ! Process not considered for T < Tmin_exc in order to avoid numerical problems 
           IF (T.GT.Tmin_exc) THEN 

              pos1 = 0
              pos3 = posT
              DO i = 1,nb_bins - 1

                 i1   = i + 1
                 tmp1 = xi(i)
                 tmp2 = xn*mm_N2
                 tmp3 = ek_const(i)   
                 tmp4 = EkJ(i)*ov_kbT2

                 pos4 = i1*posT
                 DO ip = i1,nb_bins 
  
                    ip1  = ip + 1
                    pos2 = pos1 + ip - i
  
                    ! Useful data
                    c2 = be_kf(pos2)
                    c3 = ce_kf(pos2)

                    ! Forward and backward rate coefficients
                    kf  = ae_kf(pos2)*DEXP(c2*ln_T - c3*ov_T)
                    ov_keq = ek_const(ip)/tmp3                 

                    kb  = kf*ov_keq

                    ! Derivatives of forward and backwards reaction rate coefficients
                    fac_der = c2*ov_T + c3*ov_T2

                    dKeq_dT = 1.d0/ov_keq*(EkJ(ip)*ov_kbT2 - tmp4)
                    dkf_dT  = kf*fac_der
                    dkb_dT  = ov_keq*(dkf_dT - ov_keq*kf*dKeq_dT)
                    
                    ! Excitation terms
                    exp_2 = (kf*tmp1 - kb*xi(ip))
                    exp_1 = tmp2*exp_2

                    omega(i1)  = omega(i1)  - exp_1
                    omega(ip1) = omega(ip1) + exp_1

                    ! Jacobian common factors
                    fac1 = kf*xn
                    fac2 = kb*xn
                    fac3 = tmp2*(dkf_dT*tmp1 - dkb_dT*xi(ip))
                    fac4 = mm_N2*exp_2*ov_mmN

                    ! Jacobian components 
                    js_omega(pos3 + 1)    = js_omega(pos3 + 1)   - fac4
                    js_omega(pos3 + i1)   = js_omega(pos3 + i1)  - fac1
                    js_omega(pos3 + ip1)  = js_omega(pos3 + ip1) + fac2 
                    js_omega(pos3 + posT) = js_omega(pos3 + posT) - fac3

                    js_omega(pos4 + 1)    = js_omega(pos4 + 1)   + fac4
                    js_omega(pos4 + i1)   = js_omega(pos4 + i1)  + fac1
                    js_omega(pos4 + ip1)  = js_omega(pos4 + ip1) - fac2 
                    js_omega(pos4 + posT) = js_omega(pos4 + posT) + fac3
           
                    ! Index update 
                    pos4 = pos4 + posT

                ENDDO

                ! Index update
                pos1 = pos1 + nb_bins - i
                pos3 = pos3 + posT

             ENDDO

           ENDIF

        ENDIF

        ! Collisional de-excitation vT transfer 
        ! N2(v) + N2 = N2(v - 1) + N2
        IF (vTm.EQV..TRUE.) THEN 

           tmp1 = mm_N2*xn2
           tmp2 = xn2

           DO i = 2,nb_bins

              ! Forward and backward rate coefficients
              im1 = i - 1
              c2  = b_vtm(im1)
              c3  = c_vtm(i - 1)

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
              
              ! Jacobian common factors 
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

        ! Collisional excitation vv transfer
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

      END SUBROUTINE source_VC_rr_Jac 

  END MODULE mod_nitrogen_NASA_CFD_source_VC_rr
!------------------------------------------------------------------------------!
