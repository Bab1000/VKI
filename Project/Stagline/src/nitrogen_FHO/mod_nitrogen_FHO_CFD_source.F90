!------------------------------------------------------------------------------!
! This module provides subroutines for the computation of source of conservation equations for the N-N2 system
! when using the FHO database.
  MODULE mod_nitrogen_FHO_CFD_source

    IMPLICIT NONE

    REAL(KIND=8), PARAMETER :: Tmin = 300.d0
    REAL(KIND=8), PARAMETER :: Park_fac = 0.23d0

    CONTAINS 

      !------------------------------------------------------!
      ! This subroutine computes the source term due to collisional processes for the N-N2 system.
      SUBROUTINE source (rhoi, T, omega)
!use mod_general_data, only : boite, use_boite
        USE mod_nitrogen_FHO_initialize_CFD,    ONLY: nb_bins, nb_ns, nb_eq, gn, ukb, una, Edis, mm_N, mm_N2, EkJ, &
                                                   &  ad_kf, bd_kf, cd_kf, a_n2_dis, b_n2_dis, c_n2_dis, &
                                                   &  ae_kf, be_kf, ce_kf, a_n2_exc, b_n2_exc, c_n2_exc, &
                                                   &  a_vv, b_vv, c_vv, vTa, vTm, vv, da, dm, v_Eff, ek, upi
        USE mod_nitrogen_FHO_CFD_prop,          ONLY: Q_trans, Q_rot  

        INTEGER :: i, j, i1, ip, ip1, pos1, pos
        REAL(KIND=8) :: kf, kf2, kb, kb2, ov_keq
        REAL(KIND=8) :: Tclip, ln_v, ln_T, ov_lnT, ov_T, ov_T2, ov_T3, ov_T4, exp_1, exp_2
        REAL(KIND=8) :: xn, xns, xn2, Qn, Qn2, Qrotn2
        REAL(KIND=8) :: tmp1, tmp2, tmp3, const1, const2, Z
        REAL(KIND=8), DIMENSION(nb_bins) :: xi, ek_const

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: omega
real(kind=8), dimension(size(omega)) :: omega_cp		
real(kind=8) :: OmegaCV, OmegaVT

integer	::	dv_max

        ! Useful data
        IF (T.LT.Tmin) THEN
           Tclip = Tmin
        ELSE 
           Tclip = T 
        ENDIF

        ln_T   = DLOG(Tclip)
        ov_lnT = 1.d0/ln_T
        ov_T   = 1.d0/Tclip
        ov_T2  = ov_T**2
        ov_T3  = ov_T2*ov_T
        ov_T4  = ov_T3*ov_T

        ! Species concentrations 
        xn  = rhoi(1)/mm_N
        xns = xn*xn
 
        xn2 = 0.d0
        DO i = 1,nb_bins 
           tmp1  = rhoi(i + 1)/mm_N2
           xi(i) = tmp1
           xn2   = xn2 + tmp1
        ENDDO

        ! Translational and rotational partition functions
        CALL Q_trans (Tclip, mm_N, Qn)
        CALL Q_trans (Tclip, mm_N2, Qn2) 
        CALL Q_rot   (Tclip, Qrotn2)

        ! Useful constants for the computation of backward rate coefficients (excitation and dissociation)
        tmp1   = 1.d0/(ukb*Tclip)
        const1 = -Edis*tmp1

        DO i = 1,nb_bins
           ek_const(i) = DEXP(EkJ(i)*tmp1)
        ENDDO

        ! Initialization
        omega = 0.d0
OmegaCV = 0.d0
OmegaVT = 0.d0

!print*, da, dm, vta, vtm, vv
!pause

        ! Collisional dissociation (N impact)
        ! N2(v) + N = 3N   
        IF ((da.EQV..TRUE.).AND.(dm.EQV..FALSE.)) THEN

           tmp1 = 0.d0
           tmp2 = (una*Qn2*Qrotn2)/((gn*Qn)**2*DEXP(const1))
 
           DO i = 1,nb_bins 

              kf     = ad_kf(i)*DEXP(bd_kf(i)*ln_T - cd_kf(i)*ov_T)
              ov_keq = tmp2/ek_const(i)
              kb     = kf*ov_keq 

              exp_1 = xn*(kf*xi(i) - kb*xns)

              tmp1 = tmp1 + exp_1
              omega(i + 1) = - exp_1*mm_N2 

           ENDDO

           omega(1) = 2.d0*tmp1*mm_N 
                       
        ! Collisional dissociation (N and N2 impact)
        ! N2(v) + N  = 3N   
        ! N2(v) + N2 = 2N + N2
        ELSEIF ((da.EQV..TRUE.).AND.(dm.EQV..TRUE.)) THEN

           tmp1 = 0.d0
           tmp2 = (una*Qn2*Qrotn2)/((gn*Qn)**2*DEXP(const1))
       
           DO i = 1,nb_bins

				 kf     = ad_kf(i)*DEXP(bd_kf(i)*ln_T - cd_kf(i)*ov_T)
				 kf2    = a_n2_dis(i)*DEXP(b_n2_dis(i)*ln_T - c_n2_dis(i)*ov_T)
                 ov_keq = tmp2/ek_const(i)
                 kb     = kf*ov_keq
                 kb2    = kf2*ov_keq

                 exp_1 = xn*(kf*xi(i) - kb*xns) + xn2*(kf2*xi(i) - kb2*xns)

                 tmp1 = tmp1 + exp_1
                 omega(i + 1) = - exp_1*mm_N2 
                 
				 OmegaCV = OmegaCV - mm_N2*exp_1*ek(i)
           ENDDO
 
           omega(1) = 2.d0*tmp1*mm_N
        ENDIF
 
		omega_cp = omega

dv_max = 2			          
        ! vTa transfer 
        ! N2(v) + N = N2(v') + N
        IF (vTa.EQV..TRUE.) THEN

           pos1 = 0  
           DO i = 1,nb_bins - 1

              tmp1 = xi(i)
              tmp2 = xn*mm_N2
              tmp3 = ek_const(i)  

              i1 = i + 1

              DO ip = i1,nb_bins

                 ip1 = ip + 1

                 pos1 = pos1 + 1
!if(ip-i <= dv_max) then
                 
				 kf   = ae_kf(pos1)*DEXP(be_kf(pos1)*ln_T - ce_kf(pos1)*ov_T)
                 ov_keq = tmp3/ek_const(ip)
                 kb     = kf*ov_keq

                 exp_1 = tmp2*(kf*xi(ip) - kb*tmp1) 
                 omega(ip1)  = omega(ip1) - exp_1 
                 omega(i1)   = omega(i1)  + exp_1   
!endif
!		OmegaVT = OmegaVT +(ek(i)-ek(ip))*exp_1
              ENDDO

           ENDDO 

        ENDIF

        ! vTm transfer 
        ! N2(v) + N2 = N2(v - 1) + N2
        IF (vTm.EQV..TRUE.) THEN

           pos1 = 0  
           DO i = 1,nb_bins - 1

              tmp1 = xi(i)
              tmp2 = xn2*mm_N2
              tmp3 = ek_const(i)  

              i1 = i + 1

              DO ip = i1,nb_bins

                 ip1 = ip + 1
                 pos1 = pos1 + 1

!if(ip-i <= dv_max) then 
                 kf   = a_n2_exc(pos1)*DEXP(b_n2_exc(pos1)*ln_T - c_n2_exc(pos1)*ov_T)
                 ov_keq = tmp3/ek_const(ip)
                 kb     = kf*ov_keq

                 exp_1 = tmp2*(kf*xi(ip) - kb*tmp1) 
                 omega(ip1)  = omega(ip1) - exp_1 
                 omega(i1)   = omega(i1)  + exp_1  
!endif
!OmegaVT = OmegaVT +(ek(i)-ek(ip))*exp_1
              ENDDO

           ENDDO 
        ENDIF
		!boite(1:nb_ns) = omega_cp
		!boite(nb_ns+1:2*nb_ns) = omega-omega_cp
		!boite(2*nb_ns+1) = OmegaCV
		!boite(2*nb_ns+2) = OmegaVT
        ! vv transfer 
        ! N2(v) + N2(w - 1) = N2(v - 1) + N2(w)
        IF (vv.EQV..TRUE.) THEN
		
		print*, 'Source : Attention, VV désactivé dans la bibliothèque FHO'
!			Z = upi*(4.2d-10)**2 * sqrt(8*ukb*T/(upi*mm_N2/2.d0))
!			
!			do i=1:nb_bins-1
!				
!				do j=i+1:nb_bins
!					
!					theta_b = 
!					
!					lambda = 2.d0**(-1.5) * sqrt(theta_b / T) * abs(DeltaE) / (ukb*theta_vib)
!					
!					S	= una * 9.37d-8*T* v / (1-xe*v) * (w+1) / (1-xe*(w+1)) * 0.5*(3-exp(-2/3*lambda))*exp(-2/3*lambda)
!					
!					kvv = Z * S * exp(-DeltaE/(2*ukb*T))
!					
!					
!				enddo
!				
!			enddo
!			
!           pos = 0
!
!           ! i (vibrational level v) 
!           DO i = 2,nb_bins 
!
!              tmp1 = ek_const(i - 1)/ek_const(i)
!              tmp2 = xi(i) 
!              tmp3 = xi(i - 1) 
!
!              ! j (vibrational level w)
!              DO j = i,nb_bins
!         
!                 pos = pos + 1
!
!                 ! Forward and backward rate coefficients 
!                 kf     = a_vv(pos)*DEXP(b_vv(pos)*ln_T - c_vv(pos)*ov_T) 
!                 ov_keq = tmp1*ek_const(j)/ek_const(j - 1)
!                 kb     = kf*ov_keq
! 
!                 ! Excitation terms
!                 exp_1 = mm_N2*(kf*tmp2*xi(j - 1) - kb*tmp3*xi(j))
!                  
!                 omega(i + 1) = omega(i + 1) - exp_1
!                 omega(i)     = omega(i)     + exp_1
!                 omega(j + 1) = omega(j + 1) + exp_1
!                 omega(j)     = omega(j)     - exp_1 
!
!              ENDDO
!
!           ENDDO 

        ENDIF

        ! Global momentum and energy equations
        omega(nb_ns + 1:nb_eq) = 0.d0 
		
		
      END SUBROUTINE source

      !----------------------------------------------------!
      ! This subroutine computes the source term due to collisional processes for the N-N2 system.
      ! The source term Jacobian is also provided in output.
      SUBROUTINE source_Jac (rhoi, T, omega, js_omega)

        USE mod_nitrogen_FHO_initialize_CFD,     ONLY: nb_bins, nb_ns, nb_eq, gn, ukb, una, Edis, mm_N, mm_N2, fac_keq, EkJ, &
		                                   &  ad_kf, bd_kf, cd_kf, a_n2_dis, b_n2_dis, c_n2_dis, &
                                                   &  ae_kf, be_kf, ce_kf, a_n2_exc, b_n2_exc, c_n2_exc, &
						   &  a_vv, b_vv, c_vv, vTa, vTm, vv, da, dm, v_Eff

        USE mod_nitrogen_FHO_CFD_prop,           ONLY: Q_trans, Q_trans_der, Q_rot_der         

        INTEGER :: i, j, ip, i1, ip1, pos, pos1, pos2, pos3, pos4, posT
        REAL(KIND=8) :: kf, kf_n2, kb, kb_n2, ov_keq, Ev, Evm1
        REAL(KIND=8) :: fac_der, fac_der_n2, dkf_dT, dkf_n2_dt, dkb_dT, dkb_n2_dt, dKeq_dT, dkdT, Qrot, dQrot_dT
        REAL(KIND=8) :: Tclip, ln_T, ov_T, ov_T2, ov_kbT2, ov_mmN, exp_1, exp_2
        REAL(KIND=8) :: xn, xns, Qn, Qn2, dQn_dT, xn2
        REAL(KIND=8) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7
        REAL(KIND=8) :: c2, c2_n2, c3, c3_n2, fac1, fac2, fac3, fac4, fac5
        REAL(KIND=8) :: const1, const2, en_const
        REAL(KIND=8), DIMENSION(nb_bins) :: xi, ek_const

        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: omega, js_omega
        
        real(kind=8), dimension(:), allocatable :: omg, omgp, rhoip
        real(kind=8) :: dXi, dT
        integer	::	szomg
        
        js_omega = 0.d0
        
        szomg = size(omega)
        
        allocate(omg(szomg))
        allocate(omgp(szomg))
        allocate(rhoip(szomg-2))
        !allocate(source(szomg))

        call source(rhoi,T,omg)
        
        do i=1,szomg-2
        	rhoip = rhoi
        	if(rhoi(i)<1.d-50) then
        		dXi = 1.d-50
        	else
        		dXi = 1.d-6*rhoi(i)
       		endif    
                rhoip    = rhoi
        	rhoip(i) = rhoi(i) + dXi	
        	call source(rhoip,T,omgp)
        	do j=1,szomg-2
                    js_omega((j-1)*(szomg-1)+i) = (omgp(j)-omg(j))/dXi
                enddo
        enddo
        
        dT = 1.d-6*T
        call source(rhoi,T+dT,omgp)
	do j=1,szomg-2
           js_omega((j-1)*(szomg-1)+szomg-1) = (omgp(j)-omg(j))/dT
        enddo

        omega = omg
        
        deallocate(omg)
        deallocate(omgp)
        deallocate(rhoip)
        !PRINT*,omega
        !PAUSE
!print*, " Jacobien"
		

        ! Useful data 
!        IF (T.LT.Tmin) THEN 
!           Tclip = Tmin
!        ELSE 
!           Tclip = T
!        ENDIF
!        ln_T    = DLOG(Tclip)
!        ov_T    = 1.d0/Tclip
!        ov_T2   = ov_T**2
!        ov_kbT2 = 1.d0/(ukb*Tclip**2)
!        ov_mmN  = 1.d0/mm_N 
!
!        posT = nb_ns + 1
!
!        ! Species concentrations 
!        xn   = rhoi(1)/mm_N
!        xns  = xn*xn
!		xn2  = 0.d0
!		
!        DO i = 1,nb_bins 
!           xi(i) = rhoi(i + 1)/mm_N2
!		   xn2 = xn2 + xi(i)
!        ENDDO
!
!        ! Internal and translational partition function (and derivatives with respect to temperature)
!        CALL Q_trans_der (Tclip, mm_N, Qn, dQn_dT)
!        CALL Q_trans (Tclip, mm_N2, Qn2)
!        CALL Q_rot_der (Tclip, Qrot, dQrot_dT)
!
!        ! Useful constants for the computation of backward rate coefficients (excitation and dissociation)
!        const1 = -Edis/(ukb*Tclip)
!        const2 = 1.d0/(ukb*Tclip)
!
!        en_const = DEXP(const1) 
!        DO i = 1,nb_bins
!           ek_const(i) = DEXP(EkJ(i)*const2)
!        ENDDO
!
!        ! Initialization
!        omega    = 0.d0
!        js_omega = 0.d0
!
!        ! Dissociation (N impact)
!        ! N2(v) + N = 3N
!        IF ((da.EQV..TRUE.).AND.(dm.EQV..FALSE.)) THEN
!
!           pos1 = posT
!
!           tmp1 = 0.d0
!           tmp2 = (una*Qn2*Qrot)/((gn*Qn)**2*en_const)  
!           tmp3 = 0.d0
!           tmp4 = 0.d0
!           tmp5 = 1.d0/Qrot
!           tmp6 = dQn_dT - Qn*dQrot_dT*tmp5
!
!           DO i = 1,nb_bins 
!
!              ! Common factors
!              c2   = bd_kf(i)
!              c3   = cd_kf(i)
!              tmp7 = xi(i)
!
!              ! Forward and backward rate coefficients
!              kf   = ad_kf(i)*DEXP(c2*ln_T - c3*ov_T)
!
!              ov_keq = tmp2/ek_const(i)
!              kb     = kf*ov_keq 
!
!              ! Derivatives of forward and backward reaction reate coefficients with respect to temperature
!              fac_der = c2*ov_T + c3*ov_T2
!              dKeq_dT = fac_keq*ek_const(i)*en_const*tmp5*(Qn*(Edis - EkJ(i))*ov_kbT2 + tmp6)
!
!              dkf_dT  = kf*fac_der
!              dkb_dT  = ov_keq*(dkf_dT - ov_keq*kf*dKeq_dT)
!            
!              ! Production terms
!              exp_1 = xn*(kf*tmp7 - kb*xns)
!
!              tmp1         = tmp1 + exp_1
!              omega(i + 1) = - exp_1*mm_N2
!
!              ! Jacobian common factors
!              fac1 = (kf*tmp7 - 3.d0*kb*xns)/mm_N  
!              fac2 = xn*kf/mm_N2 
!              fac3 = xn*(dkf_dT*tmp7 - dkb_dT*xns) 
!
!              ! Jacobian components
!              tmp3 = tmp3 + fac1
!              js_omega(i + 1) = 2.d0*fac2*mm_N
!
!              js_omega(pos1 + 1)     = - fac1*mm_N2
!              js_omega(pos1 + i + 1) = - fac2*mm_N2
!
!              tmp4 = tmp4 + fac3
!              js_omega(pos1 + posT) = - fac3*mm_N2
!
!              ! Index update
!              pos1 = pos1 + posT
!
!           ENDDO
!
!           ! Components relative to the equation for N (multiplication by a 2mm_N factor) 
!           tmp2 = 2.d0*mm_N
!
!           omega(1)       = tmp1*tmp2 
!
!           js_omega(1)    = tmp3*tmp2
!           js_omega(posT) = tmp4*tmp2  
!
!        ! Dissociation (N and N2 impact)
!        ! N2(v) + N  = 3N
!        ! N2(v) + N2 = 2N + N2
!        ELSEIF ((da.EQV..TRUE.).AND.(dm.EQV..TRUE.)) THEN
!
!           pos1 = posT
!
!           tmp1 = 0.d0
!           tmp2 = (una*Qn2*Qrot)/((gn*Qn)**2*en_const)  
!           tmp3 = 0.d0
!           tmp4 = 0.d0
!           tmp5 = 1.d0/Qrot
!           tmp6 = dQn_dT - Qn*dQrot_dT*tmp5
!
!           DO i = 1,nb_bins 
!
!              ! Common factors
!              c2   = bd_kf(i)
!              c3   = cd_kf(i)
!	      c2_n2 = b_n2_dis(i)
!	      c3_n2 = c_n2_dis(i)
!			  
!              tmp7 = xi(i)
!
!              ! Forward and backward rate coefficients
!              kf   = ad_kf(i)   *DEXP(c2*ln_T - c3*ov_T)
!			  kf_n2= a_n2_dis(i)*DEXP(c2_n2*ln_T - c3_n2*ov_T)
!
!              ov_keq = tmp2/ek_const(i)
!              kb     = kf*ov_keq
!			  kb_n2  = kf_n2*ov_keq 
!
!              ! Derivatives of forward and backward reaction reate coefficients with respect to temperature
!              fac_der = c2*ov_T + c3*ov_T2
!			  fac_der_n2 = c2_n2*ov_T + c3_n2*ov_T2
!              dKeq_dT = fac_keq*ek_const(i)*en_const*tmp5*(Qn*(Edis - EkJ(i))*ov_kbT2 + tmp6)
!
!              dkf_dT     = kf*fac_der
!			  dkf_n2_dT  = kf_n2*fac_der_n2
!              dkb_dT     = ov_keq*(dkf_dT -    ov_keq*kf*dKeq_dT)
!			  dkb_n2_dT  = ov_keq*(dkf_n2_dT - ov_keq*kf_n2*dKeq_dT)
!            
!              ! Production terms
!              exp_1 = xn*(kf*tmp7 - kb*xns) + xn2*(kf_n2*tmp7 - kb_n2*xns)
!
!              tmp1         = tmp1 + exp_1
!              omega(i + 1) = - exp_1*mm_N2
!
!              ! Jacobian common factors
!              fac1 = (kf*tmp7 - 3.d0*kb*xns - 2*kb_n2*xi(i)*xn)/mm_N  
!              fac2 = (kf*xn + kf_n2*xn2 + kf_n2*xi(i) - kb_n2*xns) / mm_N2
!              fac3 = xn*(dkf_dT*tmp7 - dkb_dT*xns) + xn2*(dkf_n2_dT*tmp7 - dkb_n2_dT*xns)
!			  
!              ! Jacobian components
!              tmp3 = tmp3 + fac1
!              js_omega(i + 1) = 2.d0*fac2*mm_N
!
!              js_omega(pos1 + 1)     = - fac1*mm_N2
!              js_omega(pos1 + i + 1) = - fac2*mm_N2
!
!              tmp4 = tmp4 + fac3
!              js_omega(pos1 + posT) = - fac3*mm_N2
!
!              ! Index update
!              pos1 = pos1 + posT
!
!           ENDDO
!
!           ! Components relative to the equation for N (multiplication by a 2mm_N factor) 
!           tmp2 = 2.d0*mm_N
!
!           omega(1)       = tmp1*tmp2 
!
!           js_omega(1)    = tmp3*tmp2
!           js_omega(posT) = tmp4*tmp2  
!
!        ENDIF
!
!        ! vTa transfer 
!        ! N2(v) + N = N2(v') + N
!        IF (vTa.EQV..TRUE. .AND. vTm.EQV..TRUE.) THEN
!
!           pos1 = 0
!           pos3 = posT
!           DO i = 1,nb_bins - 1
!
!              i1   = i + 1
!              tmp1 = xi(i)
!              tmp2 = xn*mm_N2
!              tmp3 = ek_const(i)   
!              tmp4 = EkJ(i)*ov_kbT2
!
!              pos4 = i1*posT
!              DO ip = i1,nb_bins 
!  
!                 ip1  = ip + 1
!                 pos2 = pos1 + ip - i
!  
!                 ! Useful data
!                 c2 = be_kf(pos2)
!                 c3 = ce_kf(pos2)
!				 c2_n2 = b_n2_exc(pos2)
!                 c3_n2 = c_n2_exc(pos2)
!
!                 ! Forward and backward rate coefficients
!                 kf     = ae_kf(pos2)*DEXP(c2*ln_T - c3*ov_T)
!                 kf_n2  = a_n2_exc(pos2)*DEXP(c2_n2*ln_T - c3_n2*ov_T)
!                 ov_keq = ek_const(ip)/tmp3                 
!
!                 kb     = kf*ov_keq
!				 kb_n2  = kf_n2*ov_keq
!
!                 ! Derivatives of forward and backwards reaction rate coefficients
!                 fac_der    = c2*ov_T + c3*ov_T2
!				 fac_der_n2 = c2_n2*ov_T + c3_n2*ov_T2
!
!                 dKeq_dT = 1.d0/ov_keq*(EkJ(ip)*ov_kbT2 - tmp4)
!                 dkf_dT     = kf*fac_der
!                 dkf_n2_dT  = kf_n2*fac_der_n2
!                 dkb_dT     = ov_keq*(dkf_dT - ov_keq*kf*dKeq_dT)
!				 dkb_n2_dT  = ov_keq*(dkf_n2_dT - ov_keq*kf_n2*dKeq_dT)
!
!                 ! Excitation terms
!                 exp_2 = (kf*tmp1 - kb*xi(ip))
!                 exp_1 = tmp2*exp_2
!                 
!                 omega(i1)  = omega(i1)  - exp_1
!                 omega(ip1) = omega(ip1) + exp_1
!
!                 ! Jacobian common factors
!                 fac1 = kf*xn + kf_n2*xn2 + (kf_n2*xi(i)-kb_n2*xi(ip)) 
!                 fac2 = kb*xn + kb_n2*xn2 - (kf_n2*xi(i)-kb_n2*xi(ip))
!                 fac3 = tmp2*(dkf_dT*tmp1 - dkb_dT*xi(ip) + dkf_n2_dT*tmp1 - dkb_n2_dT*xi(ip))
!                 fac4 = mm_N2*exp_2*ov_mmN
!
!                 ! Jacobian components 
!                 js_omega(pos3 + 1)    = js_omega(pos3 + 1)   - fac4
!                 js_omega(pos3 + i1)   = js_omega(pos3 + i1)  - fac1
!                 js_omega(pos3 + ip1)  = js_omega(pos3 + ip1) + fac2 
!                 js_omega(pos3 + posT) = js_omega(pos3 + posT) - fac3
!
!                 js_omega(pos4 + 1)    = js_omega(pos4 + 1)   + fac4
!                 js_omega(pos4 + i1)   = js_omega(pos4 + i1)  + fac1
!                 js_omega(pos4 + ip1)  = js_omega(pos4 + ip1) - fac2 
!                 js_omega(pos4 + posT) = js_omega(pos4 + posT) + fac3
!           
!                 ! Index update 
!                 pos4 = pos4 + posT
!
!              ENDDO
!
!              ! Index update
!              pos1 = pos1 + nb_bins - i
!              pos3 = pos3 + posT
!
!           ENDDO
!
!        ENDIF
!
!        ! vTm transfer 
!        ! N2(v) + N2 = N2(v - 1) + N2
!        IF (vTm.EQV..falsE.) THEN
!
!           PRINT*
!           WRITE(*,'(A)')'in mod_nitrogen_FHO_CFD_source.F90, not implemeted yet ...'
!           PRINT*
!           STOP
!
!        ENDIF
!
!        ! vv transfer 
!        ! N2(v) + N2(w - 1) = N2(v - 1) + N2(w)
!        IF (vv.EQV..TRUE.) THEN
!
!           pos = 0
!
!           ! i (vibrational level v) 
!           DO i = 2,nb_bins
!
!              tmp1 = ek_const(i - 1)/ek_const(i)
!              tmp2 = xi(i) 
!              tmp3 = xi(i - 1) 
!
!              Ev   = EkJ(i)
!              Evm1 = EkJ(i - 1) 
!
!              ! j (vibrational level w)
!              DO j = i,nb_bins 
!               
!                 pos  = pos + 1
!
!                 tmp4 = xi(j - 1)
!                 tmp5 = xi(j)
!
!                 ! Forward and backward rate coefficients 
!                 c2  = b_vv(pos)
!                 c3  = c_vv(pos)
!
!                 kf      = a_vv(pos)*DEXP(c2*ln_T - c3*ov_T) 
!                 ov_keq  = tmp1*ek_const(j)/ek_const(j - 1)
!                 kb      = kf*ov_keq
!                 
!                 ! Derivatives of forward and backward reaction rate coefficients with respect to temperature 
!                 fac_der = c2*ov_T + c3*ov_T2
!                 dkf_dT  = fac_der*kf
!                 dkb_dT  = ov_keq*(dkf_dT + kf*ov_kbT2*(Ev + EkJ(j - 1) - Evm1 - EkJ(j))) 
!                 
!                 ! Excitation terms
!                 exp_1        = mm_N2*(kf*tmp2*tmp4 - kb*tmp3*tmp5)
!                 omega(i + 1) = omega(i + 1) - exp_1    ! v
!                 omega(i)     = omega(i)     + exp_1    ! v - 1
!                 omega(j + 1) = omega(j + 1) + exp_1    ! w
!                 omega(j)     = omega(j)     - exp_1    ! w - 1
!                
!                 ! Jacobian common factors
!                 fac1 = - kb*tmp5
!                 fac2 = kf*tmp4
!                 fac3 = kf*tmp2
!                 fac4 = - kb*tmp3
!                 fac5 = mm_N2*(dkf_dT*tmp2*tmp4 - dkb_dT*tmp3*tmp5)
!                
!                 ! Components relative to the equation for N2(i - 1) (vibrational level v - 1)
!                 pos1 = posT*(i - 1)
!                 js_omega(pos1 + i)     = js_omega(pos1 + i)     + fac1    ! v - 1
!                 js_omega(pos1 + 1 + i) = js_omega(pos1 + 1 + i) + fac2    ! v 
!                 js_omega(pos1 + j)     = js_omega(pos1 + j)     + fac3    ! w - 1
!                 js_omega(pos1 + 1 + j) = js_omega(pos1 + 1 + j) + fac4    ! w
!                 js_omega(pos1 + posT)  = js_omega(pos1 + posT)  + fac5    ! T
!
!                 ! Components relative to the equation for N2(i) (vibrational level v)
!                 pos1 = pos1 + posT
!
!                 js_omega(pos1 + i)     = js_omega(pos1 + i)      - fac1   ! v - 1
!                 js_omega(pos1 + 1 + i) = js_omega(pos1 + 1 + i)  - fac2   ! v 
!                 js_omega(pos1 + j)     = js_omega(pos1 + j)      - fac3   ! w - 1
!                 js_omega(pos1 + 1 + j) = js_omega(pos1 + 1 + j)  - fac4   ! w
!                 js_omega(pos1 + posT)  = js_omega(pos1 + posT)   - fac5   ! T
!
!                 ! Components relative to the equation for N2(j - 1) (vibrational level w - 1)
!                 pos1 = posT*(j - 1)
!
!                 js_omega(pos1 + i)     = js_omega(pos1 + i)     - fac1    ! v - 1
!                 js_omega(pos1 + 1 + i) = js_omega(pos1 + 1 + i) - fac2    ! v 
!                 js_omega(pos1 + j)     = js_omega(pos1 + j)     - fac3    ! w - 1
!                 js_omega(pos1 + 1 + j) = js_omega(pos1 + 1 + j) - fac4    ! w
!                 js_omega(pos1 + posT)  = js_omega(pos1 + posT)  - fac5    ! T
!
!                 ! Components relative to the equation for N2(j) (vibrational level w)
!                 pos1 = pos1 + posT
!
!                 js_omega(pos1 + i)     = js_omega(pos1 + i)     + fac1    ! v - 1
!                 js_omega(pos1 + 1 + i) = js_omega(pos1 + 1 + i) + fac2    ! v 
!                 js_omega(pos1 + j)     = js_omega(pos1 + j)     + fac3    ! w - 1
!                 js_omega(pos1 + 1 + j) = js_omega(pos1 + 1 + j) + fac4    ! w
!                 js_omega(pos1 + posT)  = js_omega(pos1 + posT)  + fac5    ! T
!                 
!              ENDDO
!
!           ENDDO 
!           
!        ENDIF
!       write(9201,*) omega
      END SUBROUTINE source_Jac       
 
  END MODULE mod_nitrogen_FHO_CFD_source
!------------------------------------------------------------------------------!
