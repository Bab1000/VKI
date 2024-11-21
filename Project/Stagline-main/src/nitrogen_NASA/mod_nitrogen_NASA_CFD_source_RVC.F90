!------------------------------------------------------------------------------!
! This module provides subroutines for the computation of source of conservation equations for the N-N2 system
! when using the NASA Ames ab-initio database.
  MODULE mod_nitrogen_NASA_CFD_source_RVC

#define pre_KACHE

    IMPLICIT NONE

    ! Subroutines for source term computation
    CONTAINS

      !------------------------------------------------------!
      ! This subroutine computes the source term due to collisional processes for the N-N2 system.
      ! This subroutine is called when using the RVC collisional model. This subroutine ca be used also for computing
      ! the source term Jacobian by means of finite differences (in order to speed up its computation    
      ! pre-kached values of rate coeffiecients are used) 
      SUBROUTINE source_RVC (rhoi, temp, omega)

        USE mod_nitrogen_NASA_initialize_CFD,        ONLY: nb_ns, nb_bins, nb_eq, mm_N, mm_N2, ad_kf, bd_kf, cd_kf,   &  
                                                        &  ae_kf, be_kf, ce_kf, ap_kf, bp_kf, cp_kf, fac_keq,         & 
                                                        &  gn, gk, EkJ, Edis, una, ukb, Tmin_dis, Tmin_exc,           & 
                                                        &  nb_Da_proc, nb_exc_proc, exc_prod, exc_reac
        USE mod_nitrogen_NASA_CFD_prop,              ONLY: Q_trans 

        INTEGER :: i, reac, prod
        REAL(KIND=8) :: kf, kfp, kb, ov_keq
        REAL(KIND=8) :: T, Tclip, ln_T, ln_v, ov_T
        REAL(KIND=8) :: exp_1, exp_2, exp_3
        REAL(KIND=8) :: xn, xns, Qn, Qn2
        REAL(KIND=8) :: tmp1, tmp2, tmp3, const1, const2
        REAL(KIND=8), DIMENSION(nb_bins) :: xi, ek_const

#ifdef pre_KACHE
        INTEGER, SAVE :: first_entry = 0
        REAL(KIND=8), SAVE :: T_old = 0.d0
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: ek_const_KACHE
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: kdf_KACHE
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: kef_KACHE  
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
          
           ! Dissociation 
           DO i = 1,nb_Da_proc
             kdf_KACHE(i) = ad_kf(i)*DEXP(bd_kf(i)*ln_T - cd_kf(i)*ov_T)
           ENDDO

           ! Excitation
           DO i = 1,nb_exc_proc
             kef_KACHE(i) = ae_kf(i)*DEXP(be_kf(i)*ln_T - ce_kf(i)*ov_T)
           ENDDO
           
           T_old = Tclip
           
        ENDIF
#endif

        ! Species concentrations 
        xn  = rhoi(1)/mm_N
        xns = xn*xn
 
        tmp1 = 1.d0/mm_N2
        DO i = 1,nb_bins  
           xi(i) = rhoi(i + 1)*tmp1 
        ENDDO

        ! Internal and translational partition function
        CALL Q_trans (Tclip, mm_N, Qn)
        CALL Q_trans (Tclip, mm_N2, Qn2)

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
        ! N2(i) + N  = 3N
        DO i = 1,nb_bins 

           tmp3 = xi(i)

           ! Forward and backward rate coefficients
#ifdef pre_KACHE
           kf = kdf_KACHE(i)
#else
           kf = ad_kf(i)*DEXP(bd_kf(i)*ln_T - cd_kf(i)*ov_T)
#endif
           ov_keq = gk(i)/(tmp2*ek_const(i))
           kb     = kf*ov_keq 

           ! Production terms
           exp_1 = xn*(kf*tmp3 - kb*xns) 
           exp_2 = ap_kf(i)*(tmp3 - xns*ov_keq)
           exp_3 = exp_1 + exp_2           

           tmp1         = tmp1 + exp_3
           omega(i + 1) = - exp_3*mm_N2

        ENDDO

        omega(1) = 2.d0*tmp1*mm_N 

        ! Collisional excitation (N impact)
        ! N2(i) + N = N2(j) + N ; j > i
        tmp1 = mm_N2*xn
        DO i = 1,nb_exc_proc

           reac = exc_reac(i) 
           prod = exc_prod(i)

#ifdef pre_KACHE
           kf = kef_KACHE(i)
#else
           kf = ae_kf(i)*DEXP(be_kf(i)*ln_T - ce_kf(i)*ov_T)
#endif           
           kb     = kf*gk(reac)*ek_const(prod)/(gk(prod)*ek_const(reac))
           exp_1  = tmp1*(kf*xi(reac) - kb*xi(prod))

           omega(reac + 1) = omega(reac + 1) - exp_1
           omega(prod + 1) = omega(prod + 1) + exp_1

        ENDDO

      END SUBROUTINE source_RVC

      !----------------------------------------------------!
      ! This subroutine computes the source term due to collisional processes for the N-N2 system.
      ! This subroutine is called when using the RVC collisional model.
      SUBROUTINE source_RVC_Jac (rhoi, temp, omega, js_omega)

        USE mod_nitrogen_NASA_initialize_CFD,        ONLY: nb_ns, nb_bins, nb_eq, mm_N, mm_N2, ad_kf, bd_kf, cd_kf,     &  
                                                        &  ae_kf, be_kf, ce_kf, ap_kf, bp_kf, cp_kf, T_max, T_min, gn,  & 
                                                        &  EkJ, Edis, una, ukb, fac_keq, Tmin_dis, Tmin_exc, gk,        & 
                                                        &  nb_exc_proc, exc_reac, exc_prod
        USE mod_nitrogen_NASA_CFD_prop,              ONLY: Q_trans, Q_trans_der


        INTEGER :: i, reac, prod
        INTEGER :: pos1, pos2, posT
        REAL(KIND=8) :: kf, kb, kfp, ov_keq
        REAL(KIND=8) :: fac_der, dkf_dT, dkb_dT, dKeq_dT, dkdT
        REAL(KIND=8) :: T, Tclip, ln_T, ov_T, ov_T2, ov_kbT2, exp_1, exp_2
        REAL(KIND=8) :: xn, xns, Qn, Qn2, dQn_dT
        REAL(KIND=8) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
        REAL(KIND=8) :: c2, c3, fac1, fac2, fac3, fac4
        REAL(KIND=8) :: const1, const2, en_const
        REAL(KIND=8) :: xi_reac, xi_prod
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

        ln_T    = DLOG(Tclip)
        ov_T    = 1.d0/Tclip
        ov_T2   = ov_T**2
        ov_kbT2 = 1.d0/(ukb*Tclip**2)

        ! Useful quantity
        posT = nb_ns + 1

        ! Species concentrations 
        tmp1 = rhoi(1)
        xn   = rhoi(1)/mm_N
        xns  = xn*xn

        tmp1 = 1.d0/mm_N2
        DO i = 1,nb_bins 
           xi(i) = rhoi(i + 1)*tmp1
        ENDDO

        ! Translational partition function (and derivatives with respect to temperature)
        CALL Q_trans_der (Tclip, mm_N, Qn, dQn_dT)
        CALL Q_trans (Tclip, mm_N2, Qn2)

        ! Useful constants for the computation of backward rate coefficients (excitation and dissociation)
        const1 = -Edis/(ukb*Tclip)
        const2 = 1.d0/(ukb*Tclip)

        en_const = DEXP(const1) 
        DO i = 1,nb_bins
           ek_const(i) = DEXP(EkJ(i)*const2)
        ENDDO

        ! Collisional dissociation and pre-dissociation 
        ! N2(i) + N = 3N
        ! N2(i) = 2N
        pos1 = posT

        tmp1 = 0.d0
        tmp2 = 1.d0/(una*Qn2)*(gn*Qn)**2*en_const 
        tmp4 = 0.d0
        tmp5 = 0.d0

        ! Initialization
        omega    = 0.d0
        js_omega = 0.d0 

        DO i = 1,nb_bins 

           ! Common factors
           c2   = bd_kf(i)
           c3   = cd_kf(i)
           tmp3 = xi(i)
           tmp6 = 1.d0/gk(i)

           ! Forward and backward rate coefficients
           kf   = ad_kf(i)*DEXP(c2*ln_T - c3*ov_T)
           kfp  = 0.d0

           ov_keq = gk(i)/(tmp2*ek_const(i))
           kb     = kf*ov_keq 

           ! Derivatives of forward and backward reaction reate coefficients with respect to temperature
           fac_der = c2*ov_T + c3*ov_T2
           dKeq_dT = fac_keq*ek_const(i)*en_const*tmp6*(Qn*(Edis - EkJ(i))*ov_kbT2 + dQn_dT)

           dkf_dT  = kf*fac_der
           dkb_dT  = ov_keq*(dkf_dT - ov_keq*kf*dKeq_dT)
            
           ! Production terms
           exp_1 = xn*(kf*tmp3 - kb*xns)
           exp_2 = kfp*(tmp3 - xns*ov_keq)   

           tmp1         = tmp1 + (exp_1 + exp_2)
           omega(i + 1) = - (exp_1 + exp_2)*mm_N2

           ! Jacobian common factors
           fac1 = (kf*tmp3 - 3.d0*kb*xns)/mm_N - 2.d0*kfp*xn*ov_keq/mm_N 
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

        ! Collisional excitation (N impact)
        ! N2(i) + N = N2(j) + N ; j > i
        IF (T.GE.Tmin_exc) THEN

           tmp1 = mm_N2*xn
           tmp2 = mm_N2/mm_N
           DO i = 1,nb_exc_proc

              reac = exc_reac(i) 
              prod = exc_prod(i)

              ! Useful data
              c2 = be_kf(i)
              c3 = ce_kf(i)

              xi_reac = xi(reac)
              xi_prod = xi(prod)

              ! Forward and backward rate coefficients
              kf   = ae_kf(i)*DEXP(c2*ln_T - c3*ov_T)           
              tmp3 = gk(reac)*ek_const(prod)/(gk(prod)*ek_const(reac))
              kb   = kf*tmp3

              ! Derivatives of forward and backwards reaction rate coefficients
              fac_der = c2*ov_T + c3*ov_T2
              dkf_dT  = kf*fac_der
              dkb_dT  = dkf_dT*tmp3 + ov_kbT2*(EkJ(reac) - EkJ(prod))*kb

              exp_1 = (kf*xi_reac - kb*xi_prod)
              exp_2 = tmp1*exp_1

              ! Excitation terms
              omega(reac + 1) = omega(reac + 1) - exp_2
              omega(prod + 1) = omega(prod + 1) + exp_2

              ! Jacobian common factors
              fac1 = exp_1*tmp2
              fac2 = kf*xn
              fac3 = - kb*xn
              fac4 = tmp1*(dkf_dT*xi_reac - dkb_dT*xi_prod)

              ! Jacobian components
              ! Level i 
              pos1 = reac*posT
              js_omega(pos1 + 1)        = js_omega(pos1 + 1)        - fac1
              js_omega(pos1 + reac + 1) = js_omega(pos1 + reac + 1) - fac2
              js_omega(pos1 + prod + 1) = js_omega(pos1 + prod + 1) - fac3
              js_omega(pos1 + posT)     = js_omega(pos1 + posT)     - fac4

              ! Level j 
              pos2 = prod*posT
              js_omega(pos2 + 1)        = js_omega(pos2 + 1)        + fac1
              js_omega(pos2 + reac + 1) = js_omega(pos2 + reac + 1) + fac2
              js_omega(pos2 + prod + 1) = js_omega(pos2 + prod + 1) + fac3
              js_omega(pos2 + posT)     = js_omega(pos2 + posT)     + fac4

           ENDDO

        ENDIF

      END SUBROUTINE source_RVC_Jac 

  END MODULE mod_nitrogen_NASA_CFD_source_RVC
!------------------------------------------------------------------------------!
