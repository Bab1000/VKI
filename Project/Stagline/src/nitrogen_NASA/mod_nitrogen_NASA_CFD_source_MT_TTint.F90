!------------------------------------------------------------------------------!
! This module provides subroutines for the computation of source of conservation equations for the N-N2 system
! when using the NASA Ames ab-initio database.
  MODULE mod_nitrogen_NASA_CFD_source_MT_TTint

    IMPLICIT NONE

    CONTAINS 

    !------------------------------------------------------!
    ! This subroutine computes the source term due to collisional processes for the N-N2 system.
    ! This subroutine is called when using the MT_TTint collisional model. This subroutine ca be used also for computing
    ! the source term Jacobian by means of finite differences (in order to speed up its computation    
    ! pre-kached values of rate coeffiecients are used) 
    SUBROUTINE source_MT_TTint (rhoi, temp, omega)

       USE mod_nitrogen_NASA_initialize_CFD,         ONLY: nb_ns, nb_bins, nb_eq, Tmin_dis, Tmin_exc, mm_N, mm_N2,  & 
                                                         & ukb, una, EkJ, gk, gn, fac_Keq, ad_kf, bd_kf, cd_kf,     & 
                                                         & ae_kf, be_kf, ce_kf, exc_reac, exc_prod, nb_exc_proc,    & 
                                                         & Edis, ek, ap_kf
       USE mod_nitrogen_NASA_CFD_prop,               ONLY: Q_trans

       INTEGER :: i, reac, prod
       REAL(KIND=8) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
       REAL(KIND=8) :: T, Tint, QtN, Qint_T, Qint_Tint, ln_T, ov_kb_T, ov_kb_Tint, ov_T, ov_Tint
       REAL(KIND=8) :: kf, kf_eq, kb, kbp, kfp, kfp_eq, Keq
       REAL(KIND=8) :: xn, xn2, xns
       REAL(KIND=8) :: de, eL, eG, eL_p, eG_p
       REAL(KIND=8) :: Omega_CD, Omega_CE, Omega_CP
       REAL(KIND=8), DIMENSION(nb_bins) :: ek_const_T, ek_const_Tint

       REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
       REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: omega

       ! Translational and internal temperature 
       T    = temp(1)
       Tint = temp(2)

       ! Fix applied in order to avoid numerical problems
       IF (T.LT.Tmin_dis) THEN
         T = Tmin_dis
       ENDIF

       ! Fix applied in order to avoid numerical problems 
       IF (Tint.LT.Tmin_exc) THEN
         Tint = Tmin_exc
       ENDIF

       ! Translational partition function of N
       CALL Q_trans (T, mm_N, QtN)

       ! Concentrations of N and N2 species
       xn  = rhoi(1)/mm_N 
       xns = xn**2
       xn2 = rhoi(2)/mm_N2

       ! Useful quantities
       ln_T    = DLOG(T)
       ov_T    = 1.d0/T
       ov_Tint = 1.d0/Tint
       ov_kb_T    = ov_T/ukb
       ov_kb_Tint = ov_Tint/ukb
       Qint_T     = 0.d0
       Qint_Tint  = 0.d0
       DO i = 1,nb_bins
          tmp1 = -EkJ(i)
          tmp2 = gk(i)
          tmp3 = DEXP(ov_kb_T*tmp1)
          tmp4 = DEXP(ov_kb_Tint*tmp1)
          Qint_T    = Qint_T    + tmp2*tmp3
          Qint_Tint = Qint_Tint + tmp2*tmp4
          ek_const_T(i)    = tmp3 
          ek_const_Tint(i) = tmp4 
       ENDDO

       ! Initialization 
       omega = 0.d0

       ! Computation of forward and backward macroscopic reaction rate coefficients kf, kb
       ! and specific average energy lost/gained due to dissociation/recombination
       kf    = 0.d0
       kf_eq = 0.d0
       eL    = 0.d0
       eG    = 0.d0       

       DO i = 1,nb_bins 
          tmp1   = gk(i)*ad_kf(i)*DEXP(bd_kf(i)*ln_T - cd_kf(i)*ov_T)
          tmp2   = ek(i)
          tmp3   = tmp1*ek_const_Tint(i)
          tmp4   = tmp1*ek_const_T(i) 
          kf     = kf    + tmp3
          kf_eq  = kf_eq + tmp4
          eL     = eL + tmp2*tmp3
          eG     = eG + tmp2*tmp4
       ENDDO

       eL    = eL/kf
       eG    = eG/kf_eq
       kf    = kf/Qint_Tint
       kf_eq = kf_eq/Qint_T 

       ! Computation of forward and backward macroscopic reaction rate coefficients kf, kb
       ! and specific average energy lost/gained due to predissociation
       kfp    = 0.d0
       kfp_eq = 0.d0
       eL_p   = 0.d0
       eG_p   = 0.d0    
       DO i = 1,nb_bins
          tmp1 = gk(i)*ap_kf(i)
          tmp2 = tmp1*ek_const_Tint(i)
          tmp3 = tmp1*ek_const_T(i)
          tmp4 = ek(i)
          kfp    = kfp    + tmp2 
          kfp_eq = kfp_eq + tmp3
          eL_p   = eL_p + tmp2*tmp4
          eG_p   = eG_p + tmp3*tmp4
       ENDDO

       eL_p   = eL_p/kfp
       eG_p   = eG_p/kfp_eq
       kfp    = kfp/Qint_Tint
       kfp_eq = kfp_eq/Qint_T 

       ! Equilibrium constant for dissociation 
       Keq = fac_keq*QtN/Qint_T*DEXP(-Edis*ov_kb_T)
       kb  = kf_eq/Keq
       kbp = kfp_eq/Keq

       ! Mass production and energy transfer term due to dissociation and predissociation
       tmp1 = xn*(kf*xn2 - kb*xns) + kfp*xn2 - kbp*xns
       omega(1) = 2.d0*mm_N*tmp1
       omega(2) = - mm_N2*tmp1

       Omega_CD = - mm_N2*xn*(xn2*kf*eL - xns*kb*eG)
       Omega_CP = - mm_N2*(xn2*kfp*eL_p - xns*kbp*eG_p)

       ! Computation of energy transfer source term due to excitation/de-excitation
       Omega_CE = 0.d0
       DO i = 1,nb_exc_proc

          reac = exc_reac(i)
          prod = exc_prod(i)

          de = ek(prod) - ek(reac)
          kf = ae_kf(i)*DEXP(be_kf(i)*ln_T - ce_kf(i)*ov_T)
          kb = kf*ek_const_T(reac)/ek_const_T(prod)

          Omega_CE = Omega_CE + gk(reac)*(kf*ek_const_Tint(reac) - kb*ek_const_Tint(prod))*de

       ENDDO
       Omega_CE = Omega_CE*xn*xn2*mm_N2/Qint_Tint

       ! Energy transfer term for the internal energy equation
       omega(nb_eq) = Omega_CD + Omega_CE + Omega_CP
       
    END SUBROUTINE source_MT_TTint

    !----------------------------------------------------!
    ! This subroutine computes the source term due to collisional processes for the N-N2 system.
    ! This subroutine is called when using the RVC collisional model.
    SUBROUTINE source_MT_TTint_Jac (rhoi, temp, omega, js_omega)

      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi, temp
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: omega, js_omega

      omega    = 0.d0
      js_omega = 0.d0

      PRINT*,'not implemented yet'
      STOP

    END SUBROUTINE source_MT_TTint_Jac

  END MODULE mod_nitrogen_NASA_CFD_source_MT_TTint
!------------------------------------------------------------------------------!
