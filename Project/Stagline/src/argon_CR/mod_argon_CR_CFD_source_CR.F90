!------------------------------------------------------------------------------!
! This module provides subroutines for the computation of source of conservation equations 
! for the em, Ar, Arp system. The CR model is used in order to treat electronic excited states of Ar and Arp atoms 
! as separate pseudo-species.
  MODULE mod_argon_CR_CFD_source_CR

#define pre_KACHE

    IMPLICIT NONE

    REAL(KIND=8), PARAMETER :: Tfix = 250.d0

    CONTAINS 

      !----------------------------------------------------!
      ! This suroutine computes the source term due to collisional processes for the em, Ar, Arp system.
      SUBROUTINE source_CR (rhoi, temp, omega)

        USE mod_argon_CR_initialize_CFD,             ONLY: pos_em, pos_Ar, pos_Te, levels_Ar, levels_Arp, nb_ns, nb_eq,    & 
                                                         & Arp_lev, urg, una, ukb, uh, upi, ue, ge, mm_Ar, mm_Arp, mm_em,  &
                                                         & deltaT, Tmin, Tmax, mm_species, tvec, gi_Ar, gi_arp, EiJ_Ar,    & 
                                                         & Im, Cij, Kij, Si, Vi, gi_Ar, ken, flag_EI, flag_HI, flag_EExc,  & 
                                                         & flag_HExc
        USE mod_argon_CR_CFD_prop,                   ONLY: Q_trans

        INTEGER :: i, j, pos
        INTEGER :: ind1, ind1_p1, ind2, ind2_p1
        REAL(KIND=8) :: xe, xi_Ar, xgs_Ar, T, Te, Te2, Te3, kf, kb, rho_em, Qtr_em
        REAL(KIND=8) :: a, b, c, m, dT, dTe
        REAL(KIND=8) :: ov_mm_Ar, ov_mm_Arp, ov_mm_em
        REAL(KIND=8) :: nn, np, ne, ntot
        REAL(KIND=8) :: Ei, gi_ov_exp_Ei, ov_kb_T, ov_kb_Te, fac_Keq, ov_keq
	REAL(KIND=8) :: kei, sei, lam, ve
        REAL(KIND=8) :: omega_em, omegai_Ar, Omega_TE, Omega_EI, Omega_EExc, prod_term 
        REAL(KIND=8), DIMENSION(nb_ns) :: xi
        REAL(KIND=8), DIMENSION(levels_Ar) :: ei_Ar_exp_T, ei_Ar_exp_Te, exp_Im_T, exp_Im_Te

#ifdef pre_KACHE
        INTEGER, SAVE :: first_entry = 0
        REAL(KIND=8), SAVE :: T_old  = 0.d0
        REAL(KIND=8), SAVE :: Te_old = 0.d0
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: ei_Ar_exp_T_KACHE, ei_Ar_exp_Te_KACHE 
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: exp_Im_T_KACHE, exp_Im_Te_KACHE
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: Si_KACHE, Vi_KACHE
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, SAVE :: Cij_KACHE, Kij_KACHE 
#endif

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: temp, rhoi
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: omega

        ! Translational and free electron temperature 
        T   = temp(1)
        Te  = temp(pos_Te)
        
        ! Application of a fix in order to avoid numerical problems
        IF (T.LT.Tfix)  T  = Tfix
        IF (Te.LT.Tfix) Te = Tfix
       
        ov_kb_T  = 1.d0/(ukb*T)
        ov_kb_Te = 1.d0/(ukb*Te)

        Te2 = Te*Te
        Te3 = Te2*Te

        ! Look-up table indices 
        ind1    = INT((T - Tmin)/deltaT) + 1 
        ind1_p1 = ind1 + 1     
        ind2    = INT((Te - Tmin)/deltaT) + 1 
        ind2_p1 = ind2 + 1     

        ! Normalized temperature differences (used for look-up table interpolation)
        dT  = (T  - tvec(ind1))/deltaT
        dTe = (Te - tvec(ind2))/deltaT

        ! Free electron density 
        rho_em = rhoi(pos_em)
        
        ! Vector of species concentrations
        DO i = 1,nb_ns
           xi(i) = rhoi(i)/mm_species(i)
        ENDDO
        
        ! Concentrations of em and Ar(1)
        xe     = xi(pos_em)
        xgs_Ar = xi(pos_Ar)

        ! Translational partition functions of em
        CALL Q_trans (Te, mm_em, Qtr_em)

        ! Multiplication by ge = 2 (in order to account for the spin of the electron)
        Qtr_em = ge*Qtr_em

        ! Inverted molar masses
        ov_mm_em  = 1.d0/mm_em  
        ov_mm_Ar  = 1.d0/mm_Ar
        ov_mm_Arp = 1.d0/mm_Arp

        ! Allocating vectors for Kache values
#ifdef pre_KACHE
        IF (first_entry.EQ.0) THEN
           ALLOCATE(ei_Ar_exp_T_KACHE(levels_Ar),ei_Ar_exp_Te_KACHE(levels_Ar))
           ALLOCATE(exp_Im_T_KACHE(levels_Ar),exp_Im_Te_KACHE(levels_Ar)) 
           ALLOCATE(Cij_KACHE(levels_Ar,levels_Ar),Kij_KACHE(levels_Ar,levels_Ar))
           ALLOCATE(Si_KACHE(levels_Ar),Vi_KACHE(levels_Ar))
           first_entry = 1  
        ENDIF
#endif
 
#ifdef pre_KACHE
        IF (T.NE.T_old) THEN
            
           ! Pre-Kached values of exponential factors at temperature T
           DO i = 1,levels_Ar
              exp_Im_T_KACHE(i)    = DEXP(Im(i)*ov_kb_T)
              ei_Ar_exp_T_KACHE(i) = DEXP(EiJ_Ar(i)*ov_kb_T)
           ENDDO

           ! Pre-kached values of forward rate coefficients for heavy-particle impact ionization 
           IF (flag_HI.EQV..TRUE.) THEN
              DO i = 1,levels_Ar
                 a  = Vi(i,ind1)
                 b  = Vi(i,ind1_p1)
                 m  = (b  - a)
                 Vi_KACHE(i) = a + m*dT
              ENDDO 
           ENDIF

           ! Pre-kached values of forward rate coefficients for heavy-particle impact excitation
           IF (flag_HExc.EQV..TRUE.) THEN 
              DO i = 1,levels_Ar - 1
                 DO j = i + 1,levels_Ar
                    a  = Kij(i,j,ind1)
                    b  = Kij(i,j,ind1_p1)
                    m  = (b  - a)
                    Kij_KACHE(j,i) = a + m*dT
                 ENDDO 
              ENDDO 
           ENDIF

           T_old = T

        ENDIF

        IF (Te.NE.Te_old) THEN
            
           ! Pre-Kached values of exponential factors at temperature Te
           DO i = 1,levels_Ar
              exp_Im_Te_KACHE(i)    = DEXP(Im(i)*ov_kb_Te)
              ei_Ar_exp_Te_KACHE(i) = DEXP(EiJ_Ar(i)*ov_kb_Te)
           ENDDO

           ! Pre-kached values of forward rate coefficients for electron impact ionization
           IF (flag_EI.EQV..TRUE.) THEN
              DO i = 1,levels_Ar
                 a  = Si(i,ind2)
                 b  = Si(i,ind2_p1)
                 m  = (b  - a)
                 Si_KACHE(i) = a + m*dTe
              ENDDO
           ENDIF

           ! Pre-kached values of forward rate coefficients for electron impact excitation
           IF (flag_EExc.EQV..TRUE.) THEN
              DO i = 1,levels_Ar - 1
                 DO j = i + 1,levels_Ar
                    a  = Cij(i,j,ind2)
                    b  = Cij(i,j,ind2_p1)
                    m  = (b  - a)
                    Cij_KACHE(j,i) = a + m*dTe
                 ENDDO 
              ENDDO
           ENDIF

           Te_old = Te

        ENDIF

#endif
        
        ! Exponential factors (for equilibrium constant computation) 
        DO i = 1,levels_Ar
#ifdef pre_KACHE
           exp_Im_T(i)     = exp_Im_T_KACHE(i)
           exp_Im_Te(i)    = exp_Im_Te_KACHE(i)
           ei_Ar_exp_T(i)  = ei_Ar_exp_T_KACHE(i)
           ei_Ar_exp_Te(i) = ei_Ar_exp_Te_KACHE(i)
#else
           exp_Im_T(i)     = DEXP(Im(i)*ov_kb_T)
           exp_Im_Te(i)    = DEXP(Im(i)*ov_kb_Te)
           ei_Ar_exp_T(i)  = DEXP(EiJ_Ar(i)*ov_kb_T)
           ei_Ar_exp_Te(i) = DEXP(EiJ_Ar(i)*ov_kb_Te)
#endif
        ENDDO
     
        ! Initialization
        omega      = 0.d0
        Omega_EExc = 0.d0
        Omega_EI   = 0.d0
        Omega_TE   = 0.d0

        ! Computation of the source term due collisional-radiative processes starts here
        ! Electron-impact excitation (Ar(i) + em <=> Ar(j) + em; j > i), kf = Cij
        IF (flag_EExc.EQV..TRUE.) THEN

           DO i = 1,levels_Ar - 1

              Ei    = EiJ_Ar(i)
              xi_Ar = xi(i + 1)  
              gi_ov_exp_Ei = gi_Ar(i)/ei_Ar_exp_Te(i)        
              omegai_Ar    = 0.d0
 
              DO j = i + 1,levels_Ar

                 ! Interpolation of the forward rate coefficient Cij
#ifdef pre_KACHE
                 kf = Cij_KACHE(j,i)
#else
                 a  = Cij(i,j,ind2)
                 b  = Cij(i,j,ind2_p1)
                 m  = (b  - a)
                 kf = a + m*dTe
#endif
                 ! Backward reaction rate coefficient Fji
                 kb = kf*ei_Ar_exp_Te(j)/gi_Ar(j)*gi_ov_exp_Ei 
             
                 ! Production term
                 prod_term = xe*(xi_Ar*kf - xi(j + 1)*kb)
              
                 omegai_Ar = omegai_Ar + prod_term

                 ! Level Ar(j)
                 omega(j + 1) = omega(j + 1) + prod_term

                 ! Thermal nonequilibrium
                 Omega_EExc = Omega_EExc - (EiJ_Ar(j) - Ei)*prod_term
                
              ENDDO
              
              omega(i + 1) = omega(i + 1) - omegai_Ar
          
           ENDDO

           Omega_EExc = Omega_EExc*una

        ENDIF       
        
        ! Heavy-particle impact excitation (Ar(i) + Ar(1) <=> Ar(j) + Ar(1); j > i), kf = Kij 
        IF (flag_HExc.EQV..TRUE.) THEN

           DO i = 1,levels_Ar - 1

              xi_Ar = xi(i + 1) 
              gi_ov_exp_Ei = gi_Ar(i)/ei_Ar_exp_T(i)         
              omegai_Ar    = 0.d0

              DO j = i + 1,levels_Ar

                 ! Interpolation of the forward reaction rate coefficient Kij
#ifdef pre_KACHE
                 kf = Kij_KACHE(j,i)
#else
                 a  = Kij(i,j,ind1)
                 b  = Kij(i,j,ind1_p1)
                 m  = (b - a)
                 kf = a + m*dT
#endif

                 ! Backward reaction rate coefficient Lji
                 kb = kf*ei_Ar_exp_T(j)/gi_Ar(j)*gi_ov_exp_Ei

                 ! Production term 
                 prod_term = xgs_Ar*(xi_Ar*kf - xi(j + 1)*kb)
              
                 ! Level Ar(i)
                 omegai_Ar = omegai_Ar + prod_term

                 ! Level Ar(j) 
                 omega(j + 1) = omega(j + 1) + prod_term

              ENDDO

              omega(i + 1) = omega(i + 1) - omegai_Ar

           ENDDO

        ENDIF
        
        ! Ionization with electron-impact (Ar(i) + em <=> Arp + em + em), kf = Si
        IF (flag_EI.EQV..TRUE.) THEN

           fac_Keq  = una*((mm_Ar/mm_Arp)**1.5)/Qtr_em
           omega_em = 0.d0
           DO i = 1,levels_Ar

              ! Level of Arp to be considered
              pos = 1 + levels_Ar + Arp_lev(i)

              ! Interpolation of the forward reaction coefficient Si
#ifdef pre_KACHE
              kf = Si_KACHE(i)
#else 
              a  = Si(i,ind2)
              b  = Si(i,ind2_p1)
              m  = (b - a)
              kf = a + m*dTe
#endif
              
              ! Backward reaction coefficient Oi
              kb = kf*(gi_Ar(i)/gi_Arp(Arp_lev(i)))*exp_Im_Te(i)*fac_Keq
           
              ! Production term
              prod_term = xe*(kf*xi(i + 1) - kb*xi(pos)*xe)  
            
              ! em
              omega_em = omega_em + prod_term 

              ! Level Ar(i)
              omega(i + 1) = omega(i + 1) - prod_term

              ! Level Arp (1 or 2)
              omega(pos) = omega(pos) + prod_term

              ! Thermal nonequilibrium
              Omega_EI = Omega_EI - Im(i)*prod_term 
                
           ENDDO

           ! em 
           omega(pos_em) = omega(pos_em) + omega_em
        
           Omega_EI = Omega_EI*una     

        ENDIF  
         
        ! Ionization with heavy-impact (Ar(i) + Ar(1) <=> Arp + em + Ar(1)),   kf = Vi
        IF (flag_HI.EQV..TRUE.) THEN

           fac_Keq  = una*((mm_Ar/mm_Arp)**1.5)/Qtr_em
           omega_em = 0.d0
           DO i = 1,levels_Ar

              ! Level of Arp to be considered
              pos = 1 + levels_Ar + Arp_lev(i)

              ! Interpolation of the forward reaction coefficient Vi
#ifdef pre_KACHE
              kf = Vi_KACHE(i) 
#else 
              a  = Vi(i,ind1)
              b  = Vi(i,ind1_p1)
              m  = (b - a)
              kf = a + m*dT
#endif

              ! Backward reaction coefficient Wi
              kb = kf*(gi_Ar(i)/gi_Arp(Arp_lev(i)))*exp_Im_T(i)*fac_Keq

              ! Production term
              prod_term = xgs_Ar*(kf*xi(i + 1) - kb*xi(pos)*xe) 

              ! em
              omega_em = omega_em + prod_term 

              ! Level Ar(i)
              omega(i + 1) = omega(i + 1) - prod_term

              ! Level Arp (1 or 2)
              omega(pos) = omega(pos) + prod_term
           
           ENDDO
        
           ! em 
           omega(pos_em) = omega(pos_em) + omega_em 

        ENDIF

        ! Multiplication of all species production terms by their molecular mass
        DO i = 1,nb_ns 
           omega(i) = omega(i)*mm_species(i)
        ENDDO

        ! Energy source term due to elastic collisions between electrons and heavy particles
        ! Number density of Ar
        nn = 0.d0
        DO i = 1,levels_Ar
           nn = nn + rhoi(pos_Ar + i - 1)*una*ov_mm_Ar
        ENDDO

        ! Number density of Arp (charge neutrality constraint is applied)
        ne  = rho_em*una*ov_mm_em
        np  = ne

        ! Electron-ion collision rate
        lam = 1.24d7*DSQRT(Te3/ne)
        sei = 5.85d-10*DLOG(lam)/Te2
        ve  = DSQRT((8.d0*urg*Te)/(upi*mm_em))                     
        kei = ve*sei

        ! Electron neutral collision rate (look-up table is used)
        a   = ken(ind2)
        b   = ken(ind2_p1) 
        m   = (b - a)
        c   = a + m*dTe
        
        Omega_TE = 3.d0*rho_em*(T - Te)*urg*(c*nn*ov_mm_Ar + kei*np*ov_mm_Arp)       
       
        ! Global source term for the free electron energy conservation equation
        omega(nb_eq) = Omega_EExc + Omega_TE + Omega_EI 

      END SUBROUTINE source_CR 

      !----------------------------------------------------!
      ! This subroutine provides an alternative formulation for the computation
      ! of the energy transfer term due to elastic collisions between free
      ! electrons and heavy particles (the physical model is the same as that
      ! implemented in the mutation library)
      SUBROUTINE get_Omega_TE (rhoi, T, Te, Omega_TE)   

        USE mod_argon_CR_initialize_CFD,        ONLY: nb_ns, pos_em, pos_Ar, pos_Arp, levels_Ar, levels_Arp, & 
                                                    & ue, ukb, upi, urg, ueps0, mm_Ar, mm_Arp, mm_species, Ri

        INTEGER :: is
        REAL(KIND=8), PARAMETER :: a1 = -1.67070199d-4
        REAL(KIND=8), PARAMETER :: a2 = 4.87662165d-3
        REAL(KIND=8), PARAMETER :: a3 = -6.34929831d-2
        REAL(KIND=8), PARAMETER :: a4 =  5.49816993d-1
        REAL(KIND=8), PARAMETER :: a5 = -7.91157696d-1
        REAL(KIND=8), PARAMETER :: b1 =  0.13742d0
        REAL(KIND=8), PARAMETER :: b2 = - 3.62064d0
        REAL(KIND=8), PARAMETER :: b3 =  3.14457d1
        REAL(KIND=8), PARAMETER :: b4 = - 8.69558d1
        REAL(KIND=8) :: Be, Bh, Ds, xe, mm, nu, nu_Ar, nu_Arp, nd, rho, ov_rho, tmp, p, ph, pe, rho_em
        REAL(KIND=8) :: efac, fac0, fac_Ar, fac_Arp, y_Ar, y_Arp
        REAL(KIND=8) :: q11_Ar, q11_Arp, sum_int
        REAL(KIND=8) :: lnte1, lnte2, lnte3, lnte4
        REAL(KIND=8) :: lntse1, lntse2, lntse3, lntse4, Tse
        REAL(KIND=8), DIMENSION(nb_ns) :: yi
        
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), INTENT(IN) :: T, Te
        REAL(KIND=8), INTENT(OUT) :: Omega_TE

        ! Free electron density 
        rho_em = rhoi(pos_em)

        ! Mixture density and pressure
        rho = rho_em
        pe  = rho_em*Ri(pos_em)*Te
 
        ph = 0.d0
        DO is = pos_em + 1,nb_ns 
           tmp = rhoi(is)
           rho = rho + tmp
           ph  = ph  + tmp*Ri(is)
        ENDDO
        ph = ph*T

        p = ph + pe

        ! Species mass fractions and mixture molar mass
        mm = 0.d0
        ov_rho = 1.d0/rho
        DO is = 1,nb_ns
           tmp    = rhoi(is)*ov_rho 
           yi(is) = tmp
           mm     = tmp/mm_species(is)
        ENDDO
        mm = 1.d0/mm
        
        ! Number density and free electron molar fraction
        xe = mm/mm_species(pos_em)*yi(pos_em)
        nd = p/(ukb*T*(1.d0 + xe*(Te/T -1.d0)))  

        ! Ar and Arp mass fractions
        y_Ar = 0.d0
        DO is = 1,levels_Ar 
           y_Ar = y_Ar + yi(pos_Ar + is - 1)
        ENDDO

        y_Arp = 0.d0
        DO is = 1,levels_Arp 
           y_Arp = y_Arp + yi(pos_Arp + is - 1)
        ENDDO

        ! Average closest impact parameter (e-Arp and e-e interaction)
        fac0 = ue*ue/(8.d0*upi*ueps0*ukb)
        Be   = fac0/Te

        ! Average closest impact parameter (Arp-Arp collision)
        Bh = fac0/T

        ! Debye shielding distance (e and Arp contribution)
        Ds = DSQRT(ueps0*ukb*Te/(2.d0*nd*xe*ue*ue))
        Ds = MIN(Ds,10000.d0*(Be + Bh))

        ! Non-dimensional temperature for charge-charge interactions
        Tse = MAX(Ds/(2.d0*Be),0.1d0)
        lntse1 = DLOG(Tse)
        lntse2 = lntse1*lntse1
        lntse3 = lntse2*lntse1
        lntse4 = lntse3*lntse1

        ! Collision integral
        efac    = upi*Ds*Ds/(Tse*Tse)
        q11_Arp = efac*DEXP(a1*lntse4 + a2*lntse3 + a3*lntse2 + a4*lntse1 + a5)
 
        ! Electron neutral interaction (em-Ar)
        lnte1 = DLOG(Te)
        lnte2 = lnte1*lnte1
        lnte3 = lnte2*lnte1 

        ! Collision integral
        q11_Ar = 1.d-20*DEXP(b1*lnte3 + b2*lnte2 + b3*lnte1 + b4)

        nu    = mm*y_Ar/mm_Ar*nd
        nu_Ar = SQRT(8.d0*urg*Te/(upi*mm_Ar))*nu*q11_Ar
 
        nu     = mm*y_Arp/mm_Arp*nd
        nu_Arp = SQRT(8.d0*urg*Te/(upi*mm_Arp))*nu*q11_Arp

        sum_int = nu_Ar/mm_Ar + nu_Arp/mm_Arp

        ! Energy transfer term
        Omega_TE = 3.d0*rho_em*urg*(T - Te)*sum_int

      END SUBROUTINE get_Omega_Te

      !----------------------------------------------------!
      SUBROUTINE source_CR_Jac

      END SUBROUTINE source_CR_Jac

  END MODULE mod_argon_CR_CFD_source_CR
!------------------------------------------------------------------------------!
