!------------------------------------------------------------------------------!
! This module provides subroutines for the computation of source of conservation equations 
! for the em, Ar, Arp system. The MT model is used.
  MODULE mod_argon_CR_CFD_source_MT

#define pre_KACHE

    IMPLICIT NONE

    REAL(KIND=8), PARAMETER :: Tfix = 300.d0

    CONTAINS 

      !----------------------------------------------------!
      ! This suroutine computes the source term due to collisional processes for the em, Ar, Arp system.
      SUBROUTINE source_MT (rhoi, temp, omega)

        USE mod_argon_CR_initialize_CFD,             ONLY: pos_em, pos_Ar, pos_Arp, pos_Te, nb_ns, nb_eq, una, ukb, uh, upi,   & 
                                                         & ue, ge, mm_Ar, mm_Arp, mm_em, deltaT, Tmin, Tmax, mm_species, tvec, & 
                                                         & Eion, ken, ahi_kf, bhi_kf, chi_kf, aei_kf, bei_kf, cei_kf,          & 
                                                         & flag_EI, flag_HI 
        USE mod_argon_CR_CFD_prop,                   ONLY: Q_trans, Q_int

        INTEGER :: i, ind1, ind1_p1
        REAL(KIND=8) :: a, b, c, m 
        REAL(KIND=8) :: T, Te, dTe
        REAL(KIND=8) :: xe, xAr, xArp, kf, kb, fac_Keq
        REAL(KIND=8) :: Qtr_em, Qtr_Ar_T, Qtr_Arp_T, Qtr_Ar_Te, Qtr_Arp_Te
        REAL(KIND=8) :: Qint_Ar_T, Qint_Arp_T, Qint_Ar_Te, Qint_Arp_Te 
        REAL(KIND=8) :: ov_mm_Ar, ov_mm_Arp, ov_mm_em
        REAL(KIND=8) :: prod_term
        REAL(KIND=8) :: nn, np, ne
	REAL(KIND=8) :: kei, sei, lam, ve
        REAL(KIND=8) :: Omega_EI, Omega_TE
        REAL(KIND=8), DIMENSION(nb_ns) :: xi
#ifdef pre_KACHE 
        INTEGER, SAVE :: first_entry = 0
        REAL(KIND=8), SAVE :: T_old  = 0.d0
        REAL(KIND=8), SAVE :: Te_old = 0.d0 
        REAL(KIND=8), SAVE :: Qtr_em_KACHE
        REAL(KIND=8), SAVE :: Qtr_Ar_T_KACHE, Qtr_Arp_T_KACHE
        REAL(KIND=8), SAVE :: Qtr_Ar_Te_KACHE, Qtr_Arp_Te_KACHE 
        REAL(KIND=8), SAVE :: Qint_Ar_T_KACHE, Qint_Arp_T_KACHE
        REAL(KIND=8), SAVE :: Qint_Ar_Te_KACHE, Qint_Arp_Te_KACHE
        REAL(KIND=8), SAVE :: kf_EI_KACHE, kb_EI_KACHE
        REAL(KIND=8), SAVE :: kf_HI_KACHE, kb_HI_KACHE
        REAL(KIND=8), SAVE :: kei_KACHE
#endif

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp, rhoi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: omega

#ifdef pre_KACHE 
        IF (first_entry.EQ.0) THEN
           first_entry = 1
        ENDIF
#endif

        ! Translational and free electron temperature 
        T  = temp(1)
        Te = temp(pos_Te)
        
        ! Temperature fix in order to avoid numerical problems
        IF (T.LT.Tfix)  T  = Tfix
        IF (Te.LT.Tfix) Te = Tfix
 
        ! Vector of species concentrations
        DO i = 1,nb_ns
           xi(i) = rhoi(i)/mm_species(i)
        ENDDO
         
        ! Molar fractions of em, Ar and Arp
        xe   = xi(pos_em)
        xAr  = xi(pos_Ar)
        xArp = xi(pos_Arp) 

        ! Inverted molar masses
        ov_mm_em  = 1.d0/mm_em  
        ov_mm_Ar  = 1.d0/mm_Ar
        ov_mm_Arp = 1.d0/mm_Arp

#ifdef pre_KACHE 
        IF (Te.NE.Te_old) THEN 

          CALL Q_trans (Te, mm_em, Qtr_em_KACHE)
          Qtr_em_KACHE = ge*Qtr_em_KACHE

          ne  = rhoi(pos_em)*una*ov_mm_em
          lam = 1.24*1.d7*DSQRT(Te**3/ne)
          sei = 5.85d-10*DLOG(lam)/(Te**2)
          ve  = DSQRT((8.d0*ukb*una*Te)/(upi*mm_em))
          kei_KACHE = ve*sei

          IF (flag_EI.EQV..TRUE.) THEN

             CALL Q_trans (Te, mm_Ar, Qtr_Ar_Te_KACHE)
             CALL Q_trans (Te, mm_Arp, Qtr_Arp_Te_KACHE)
             CALL Q_int (Te, 'Ar', Qint_Ar_Te_KACHE)
             CALL Q_int (Te, 'Arp', Qint_Arp_Te_KACHE)
             kf_EI_KACHE = aei_kf*(Te**bei_kf)*DEXP(-cei_kf/Te)
             kb_EI_KACHE = una*kf_EI_KACHE*Qtr_Ar_Te_KACHE*Qint_Ar_Te_KACHE/ &
                         & (Qtr_Arp_Te_KACHE*Qint_Arp_Te_KACHE*Qtr_em_KACHE)*DEXP(Eion/(ukb*Te))

          ENDIF

          Te_old = Te

        ENDIF

        IF (T.NE.T_old) THEN 

           IF (flag_HI.EQV..TRUE.) THEN

              CALL Q_trans (T, mm_Ar, Qtr_Ar_T_KACHE)
              CALL Q_trans (T, mm_Arp, Qtr_Arp_T_KACHE)
              CALL Q_int (T, 'Ar', Qint_Ar_T_KACHE)
              CALL Q_int (T, 'Arp', Qint_Arp_T_KACHE)
              kf_HI_KACHE = ahi_kf*(T**bhi_kf)*DEXP(-chi_kf/T)
              kb_HI_KACHE = una*kf*Qtr_Ar_T_KACHE*Qint_Ar_T_KACHE/ &
                          & (Qtr_Arp_T_KACHE*Qint_Arp_T_KACHE*Qtr_em_KACHE)*DEXP(Eion/(ukb*T))

           ENDIF

           T_old = T

        ENDIF

#else

        ! Translational partition functions of em
        CALL Q_trans (Te, mm_em, Qtr_em)

        ! Multiplication by ge = 2 (in order to account for the spin of the electron)
        Qtr_em = ge*Qtr_em

        ! Translational partition functions of Ar and Arp (at T and Te)
        CALL Q_trans (T, mm_Ar, Qtr_Ar_T)
        CALL Q_trans (T, mm_Arp, Qtr_Arp_T)
        CALL Q_trans (Te, mm_Ar, Qtr_Ar_Te)
        CALL Q_trans (Te, mm_Arp, Qtr_Arp_Te)

        ! Electronic partition functions of Ar and Arp (at T and Te)
        CALL Q_int (T, 'Ar', Qint_Ar_T)
        CALL Q_int (T, 'Arp', Qint_Arp_T)
        CALL Q_int (Te, 'Ar', Qint_Ar_Te)
        CALL Q_int (Te, 'Arp', Qint_Arp_Te)
#endif

        ! Initialization
        omega    = 0.d0
        Omega_EI = 0.d0
        
        ! Computation of source term vector due collisional-radiative processes starts here
        ! Electron impact ionization
        ! Ar + em <=> Arp + em + em
        IF (flag_EI.EQV..TRUE.) THEN

#ifdef pre_KACHE 
           kf = kf_EI_KACHE
           kb = kb_EI_KACHE
#else
           kf = aei_kf*(Te**bei_kf)*DEXP(-cei_kf/Te)
           kb = una*kf*Qtr_Ar_Te*Qint_Ar_Te/(Qtr_Arp_Te*Qint_Arp_Te*Qtr_em)*DEXP(Eion/(ukb*Te))
#endif

           ! Production term
           prod_term = xe*(kf*xAr - kb*xArp*xe) 

           ! em
           omega(pos_em) = prod_term 

           ! Ar
           omega(pos_Ar) = - prod_term

           ! Arp 
           omega(pos_Arp) = prod_term
        
           ! Energy lost by electrons due to electron-impact ionization
           Omega_EI = prod_term*Eion*una 

        ENDIF
         
        ! Heavy-particle impact ionization 
        ! Ar + Ar <=> Arp + em + Ar
        IF (flag_HI.EQV..TRUE.) THEN

#ifdef pre_KACHE
           kf = kf_HI_KACHE
           kb = kb_HI_KACHE
#else
           kf = ahi_kf*(T**bhi_kf)*DEXP(-chi_kf/T)
           kb = una*kf*Qtr_Ar_T*Qint_Ar_T/(Qtr_Arp_T*Qint_Arp_T*Qtr_em)*DEXP(Eion/(ukb*T))        
#endif     
   
           ! Production term
           prod_term = xAr*(kf*xAr - kb*xArp*xe) 
       
           ! em
           omega(pos_em) = omega(pos_em) + prod_term

           ! Ar 
           omega(pos_Ar) = omega(pos_Ar) - prod_term

           ! Arp
           omega(pos_Arp) = omega(pos_Arp) + prod_term

        ENDIF

        ! Multiplication of all species production terms by their molecular mass
        DO i = 1,nb_ns 
           omega(i) = omega(i)*mm_species(i)
        ENDDO

        ! Energy transfer term due to elastic collisions between free electrons and heavy particles
        ! Number density of Ar
        nn = rhoi(pos_Ar)*una*ov_mm_Ar

        ! Number density of em 
        ne = rhoi(pos_em)*una*ov_mm_em 
        
        ! Number density of Arp (neutral gas hypothesis is applied)
        np = ne

        ! Look-up table indices (for the computation of ken)
        ind1    = INT((Te - Tmin)/deltaT) + 1 
        ind1_p1 = ind1 + 1     
        
        ! Normalized temperature difference
        dTe = (Te - tvec(ind1))/deltaT 
 
#ifdef pre_KACHE
        kei = kei_KACHE
#else 
        lam = 1.24*1.d7*DSQRT(Te**3/ne)
        sei = 5.85d-10*DLOG(lam)/(Te**2)
        ve  = DSQRT((8.d0*ukb*una*Te)/(upi*mm_em))             
        kei = ve*sei
#endif
        a   = ken(ind1)
        b   = ken(ind1_p1) 
        m   = (b - a)
        c   = a + m*dTe
         
        Omega_TE = 3.d0*rhoi(pos_em)*ukb*(T - Te)*una*(c*nn*ov_mm_Ar + kei*np*ov_mm_Arp)       
        
        ! Source term for the free electron energy equation
        omega(nb_eq) = Omega_TE - Omega_EI
        
      END SUBROUTINE source_MT 

  END MODULE mod_argon_CR_CFD_source_MT
!------------------------------------------------------------------------------!
