!------------------------------------------------------------------------------!
! This module provides subroutines for the computation of source of conservation equations for the N-N2 system
! when using the NASA Ames ab-initio database.
  MODULE mod_nitrogen_NASA_CFD_source_VC

#define pre_KACHE

    IMPLICIT NONE

    ! Subroutines for source term computation
    CONTAINS

      !------------------------------------------------------!
      ! This subroutine computes the source term due to collisional processes for the N-N2 system.
      ! This subroutine is called when using the VC model. This subroutine ca be used also for computing
      ! the source term Jacobian by means of finite differences (in order to speed up its computation    
      ! pre-kached values of rate coeffiecients are used) 
      SUBROUTINE source_VC (rhoi, temp, omega)

        USE mod_nitrogen_NASA_initialize_CFD,        ONLY: nb_ns, nb_bins, nb_eq, mm_N, mm_N2, ad_kf, bd_kf, cd_kf,   &  
                                                        &  ae_kf, be_kf, ce_kf, T_store, T_max, T_min, inv_step,      &
                                                        &  qint, gn, EkJ, Edis, una, ukb, fac_exp, Tmin_dis,          & 
                                                        &  Tmin_exc, nb_vTa_proc, nb_Da_proc, nb_exc_proc
        USE mod_nitrogen_NASA_CFD_prop,              ONLY: Q_trans

        INTEGER :: i, ip, i1, ip1, pos1, pos2
        INTEGER :: left, right
        REAL(KIND=8) :: kf, kfp, kb, ov_keq
        REAL(KIND=8) :: T, Tclip, ln_T, ln_v, ov_T, exp_1, exp_2, exp_3
        REAL(KIND=8) :: xn, xns, xn2, Qn, Qn2
        REAL(KIND=8) :: tmp1, tmp2, tmp3, const1, const2
        REAL(KIND=8), DIMENSION(nb_bins) :: xi, qk, ek_const

#ifdef pre_KACHE
        INTEGER, SAVE :: first_entry = 0
        REAL(KIND=8), SAVE :: T_old = 0.d0
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: ek_const_KACHE
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: kdf_KACHE
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: kef_KACHE  
#endif

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp, rhoi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: omega

        ! Temperature 
        T = temp(1)

        ! Fix applied in order to avoid numerical problems
        IF (T.LT.Tmin_dis) THEN 
           Tclip = Tmin_dis
        ELSE 
           Tclip = T
        ENDIF
        ln_T = DLOG(Tclip)
        ov_T = 1.d0/Tclip
         
        ! Allocating vectors for Kache values
#ifdef pre_KACHE
        IF (first_entry.EQ.0) THEN
           ALLOCATE(ek_const_KACHE(nb_bins)) 
           ALLOCATE(kdf_KACHE(nb_Da_proc))
           ALLOCATE(kef_KACHE(nb_exc_proc))
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
           
           DO i = 1,nb_Da_proc
             kdf_KACHE(i) = ad_kf(i)*DEXP(bd_kf(i)*ln_T - cd_kf(i)*ov_T)
           ENDDO

           DO i = 1,nb_exc_proc
             kef_KACHE(i) = ae_kf(i)*DEXP(be_kf(i)*ln_T - ce_kf(i)*ov_T)
           ENDDO
           
           T_old = Tclip
           
        ENDIF
#endif
        
        ! Species concentrations 
        xn = rhoi(1)/mm_N
        xns = xn*xn
 
        xn2 = 0.d0
        DO i = 1,nb_bins 
           tmp1  = rhoi(i + 1)/mm_N2 
           xn2   = xn2 + tmp1
           xi(i) = tmp1
        ENDDO

        ! Internal and translational partition function
        CALL Q_trans (Tclip, mm_N, Qn)
        CALL Q_trans (Tclip, mm_N2, Qn2)

        ! Value search
        left  = INT((Tclip - T_min)*inv_step) + 1
        tmp1  = (T - T_store(left))*inv_step 

        left  = (left - 1)*nb_bins 
        right = left + nb_bins

        ! Vibrational level rotational partition functions
        DO i = 1,nb_bins
           tmp2  = qint(left + i) 
           qk(i) = tmp2 + tmp1*(qint(right + i) - tmp2)
        ENDDO
         
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
        tmp2 = 1.d0/(una*Qn2)*(gn*Qn)**2*DEXP(const1) 

        ! Collisional dissociation (N impact) 
        ! N2(v) + N  = 3N
        DO i = 1,nb_bins 

           tmp3 = xi(i)

           ! Forward and backward rate coefficients
#ifdef pre_KACHE
           kf = kdf_KACHE(i)
#else
           kf = ad_kf(i)*DEXP(bd_kf(i)*ln_T - cd_kf(i)*ov_T)
#endif
           ov_keq = qk(i)/(tmp2*ek_const(i))
           kb     = kf*ov_keq 

           ! Production terms
           exp_1 = xn*(kf*tmp3 - kb*xns)

           tmp1         = tmp1 + exp_1
           omega(i + 1) = - exp_1*mm_N2

        ENDDO

        omega(1) = 2.d0*tmp1*mm_N 
       
        ! Collisional excitation (N impact)
        ! N2(v) + N = N2(v`) + N 
        IF (Tclip.GE.Tmin_exc) THEN

           pos1 = 0
           DO i = 1,nb_bins - 1

              i1   = i + 1
              tmp1 = xi(i)
              tmp2 = xn*mm_N2
              tmp3 = 1.d0/qk(i)*ek_const(i)              

              DO ip = i1,nb_bins 
  
                 ip1  = ip + 1
                 pos2 = pos1 + ip - i
   
                 ! Forward and backward rate coefficients
#ifdef pre_KACHE 
                 kf = kef_KACHE(pos2)
#else
                 kf = ae_kf(pos2)*DEXP(be_kf(pos2)*ln_T - ce_kf(pos2)*ov_T) 
#endif
                 ov_keq = ek_const(ip)/(tmp3*qk(ip))                 

                 kb  = kf*ov_keq

                 ! Excitation terms
                 exp_1 = tmp2*(kf*tmp1 - kb*xi(ip))
                
                 omega(i1)  = omega(i1)  - exp_1
                 omega(ip1) = omega(ip1) + exp_1
                                  
              ENDDO

              pos1 = pos1 + nb_bins - i

           ENDDO

        ENDIF
       
        ! Source term set to zero for momentum and global energy conservation equations 
        omega(nb_ns + 1:nb_eq) = 0.d0
        
      END SUBROUTINE source_VC

      !----------------------------------------------------!
      ! This subroutine computes the source term due to collisional processes for the N-N2 system.
      ! This subroutine is called when using the VC model.
      ! The source term Jacobian is also provided in output.
      SUBROUTINE source_VC_Jac (rhoi, temp, omega, js_omega)

        USE mod_nitrogen_NASA_initialize_CFD,        ONLY: nb_ns, nb_bins, nb_eq, mm_N, mm_N2, ad_kf, bd_kf, cd_kf, &  
                                                        &  ae_kf, be_kf, ce_kf, T_store, T_max, T_min, inv_step,    & 
                                                        &  qint, dqintdT, gn, EkJ, Edis, una, ukb, fac_keq,         & 
                                                        &  Tmin_dis, Tmin_exc
        USE mod_nitrogen_NASA_CFD_prop,              ONLY: Q_trans, Q_trans_der

        INTEGER :: i, ip, i1, ip1, pos1, pos2, pos3, pos4, posT
        INTEGER :: left, right
        REAL(KIND=8) ::  xn_tol, xn2_tol
        REAL(KIND=8) :: kf, kb, kfp, ov_keq
        REAL(KIND=8) :: fac_der, dkf_dT, dkb_dT, dKeq_dT, dkdT
        REAL(KIND=8) :: T, Tclip, ln_T, ov_T, ov_T2, ov_kbT2, ov_mmN, exp_1, exp_2
        REAL(KIND=8) :: xn, xns, Qn, Qn2, dQn_dT
        REAL(KIND=8) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
        REAL(KIND=8) :: c2, c3, fac1, fac2, fac3, fac4
        REAL(KIND=8) :: const1, const2, en_const
        REAL(KIND=8), DIMENSION(nb_bins) :: xi, qk, dqk_dT, ek_const
       
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

        ln_T    = DLOG(Tclip)
        ov_T    = 1.d0/Tclip
        ov_T2   = ov_T**2
        ov_kbT2 = 1.d0/(ukb*Tclip**2)
        ov_mmN  = 1.d0/mm_N 

        posT = nb_ns + 1

        ! Species concentrations 
        xn   = rhoi(1)/mm_N
        xns  = xn*xn

        DO i = 1,nb_bins 
           xi(i) = rhoi(i + 1)/mm_N2
        ENDDO

        ! Internal and translational partition function (and derivatives with respect to temperature)
        CALL Q_trans_der (Tclip, mm_N, Qn, dQn_dT)
        CALL Q_trans (Tclip, mm_N2, Qn2)

        ! Value search
        left = INT((Tclip - T_min)*inv_step) + 1
        tmp1 = (T - T_store(left))*inv_step

        left  = (left - 1)*nb_bins 
        right = left + nb_bins

        ! Bin or vibrational level internal partition functions
        DO i = 1,nb_bins
           tmp2 = qint(left + i)
           tmp3 = dqintdT(left + i)
           qk(i)     = tmp2 + tmp1*(qint(right + i) - tmp2)
           dqk_dT(i) = tmp3 + tmp1*(dqintdT(right + i) - tmp3)
        ENDDO
        
        ! Useful constants for the computation of backward rate coefficients (excitation and dissociation)
        const1 = -Edis/(ukb*Tclip)
        const2 = 1.d0/(ukb*Tclip)

        en_const = DEXP(const1) 
        DO i = 1,nb_bins
           ek_const(i) = DEXP(EkJ(i)*const2)
        ENDDO

        ! Collisional dissociation and pre-dissociation 
        ! N2(i) + N = 3N  , N2(i) = 2N
        pos1 = posT

        tmp1 = 0.d0
        tmp2 = 1.d0/(una*Qn2)*(gn*Qn)**2*en_const 
        tmp4 = 0.d0
        tmp5 = 0.d0

        ! Initializaton
        js_omega = 0.d0

        DO i = 1,nb_bins 

           ! Common factors
           c2   = bd_kf(i)
           c3   = cd_kf(i)
           tmp3 = xi(i)
           tmp6 = 1.d0/qk(i)

           ! Forward and backward rate coefficients
           kf   = ad_kf(i)*DEXP(c2*ln_T - c3*ov_T)
           kfp  = 0.d0

           ov_keq = qk(i)/(tmp2*ek_const(i))
           kb     = kf*ov_keq 

           ! Derivatives of forward and backward reaction reate coefficients with respect to temperature
           fac_der = c2*ov_T + c3*ov_T2
           dKeq_dT = fac_keq*ek_const(i)*en_const*tmp6*(Qn*(Edis - EkJ(i))*ov_kbT2 + dQn_dT - Qn*dqk_dT(i)*tmp6)

           dkf_dT  = kf*fac_der
           dkb_dT  = ov_keq*(dkf_dT - ov_keq*kf*dKeq_dT)
            
           ! Production terms
           exp_1 = xn*(kf*tmp3 - kb*xns)

           tmp1         = tmp1 + exp_1
           omega(i + 1) = - exp_1*mm_N2

           ! Jacobian common factors
           fac1 = (kf*tmp3 - 3.d0*kb*xns)/mm_N  
           fac2 = xn*kf/mm_N2 
           fac3 = xn*(dkf_dT*tmp3 - dkb_dT*xns) 

           ! Jacobian components
           tmp4 = tmp4 + fac1
           js_omega(i + 1) = 2.d0*fac2*mm_N

           js_omega(pos1 + 1)     = - fac1*mm_N2
           js_omega(pos1 + i + 1) = - fac2*mm_N2

           tmp5 = tmp5 + fac3
           js_omega(pos1 + posT) = - fac3*mm_N2

           ! Index update
           pos1 = pos1 + posT

        ENDDO

        omega(1) = 2.d0*tmp1*mm_N 
        
        js_omega(1)    = 2.d0*tmp4*mm_N
        js_omega(posT) = 2.d0*tmp5*mm_N         

        ! Collisional excitation N2(i) + N = N2(ip) + N ; ip > i (a fix is applied for T < Tmin_exc K)
        IF (Tclip.GE.Tmin_exc) THEN

           pos1 = 0
           pos3 = posT
           DO i = 1,nb_bins - 1

              i1   = i + 1
              tmp1 = xi(i)
              tmp2 = xn*mm_N2
              tmp3 = 1.d0/qk(i)*ek_const(i)   
              tmp4 = EkJ(i)*ov_kbT2
              tmp5 = dqk_dT(i)/qk(i)           

              pos4 = i1*posT
              DO ip = i1,nb_bins 
  
                 ip1  = ip + 1
                 pos2 = pos1 + ip - i
  
                 ! Useful data
                 c2 = be_kf(pos2)
                 c3 = ce_kf(pos2)

                 ! Forward and backward rate coefficients
                 kf  = ae_kf(pos2)*DEXP(c2*ln_T - c3*ov_T)
                 ov_keq = ek_const(ip)/(tmp3*qk(ip))                 

                 kb  = kf*ov_keq

                 ! Derivatives of forward and backwards reaction rate coefficients
                 fac_der = c2*ov_T + c3*ov_T2

                 dKeq_dT = 1.d0/ov_keq*(EkJ(ip)*ov_kbT2 - tmp4 + dqk_dT(ip)/qk(ip) - tmp5)
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

        ! Source term set to zero for momentum and global energy conservation equations 
        omega(nb_ns + 1:nb_eq) = 0.d0

      END SUBROUTINE source_VC_Jac 

  END MODULE mod_nitrogen_NASA_CFD_source_VC
!------------------------------------------------------------------------------!
