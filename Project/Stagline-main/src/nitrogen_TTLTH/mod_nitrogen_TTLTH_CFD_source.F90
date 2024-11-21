!------------------------------------------------------------------------------!
! This module provides subroutines for the computation of source of conservation equations for the N-N2 system
! when using TTLTH multi-temperature model.
  MODULE mod_nitrogen_TTLTH_CFD_source

  USE mod_nitrogen_TTLTH_initialize_CFD,    ONLY: upi, una, urg, ukb, nb_ns, nb_tvib, nb_trot, nb_eq, nb_dim, nb_temp, & 
                                                 posTr, posTv, diss_imp_N2, vt_imp_mol

  IMPLICIT NONE

  REAL(KIND=8), PARAMETER :: rhoi_tol = 1.d-25
  REAL(KIND=8), PARAMETER :: Tmin = 250.d0
  REAL(KIND=8) :: Qn, Q_int

  ! Subroutine for the computations of source terms
  CONTAINS

    !------------------------------------------------------!
    ! This suroutine computes the source term due to collisional processes for the N-N2 system.
    SUBROUTINE source (rhoi, tvec, omega)
!use mod_general_data, only: cell_number, boite, use_boite
      USE mod_nitrogen_TTLTH_initialize_CFD,   ONLY: mm_N, mm_N2, Rn, Rn2, fac_keq, c1, c2, c3, c4, theta_vib, theta_rot,& 
                                                   & mw_a, mw_b, fac_Q, mu, compteur_global, get_chemistry, get_VT,  & 
                                                   & get_VV, upi, ukb, uh, una, ekJ, n_T, tb_T, nb_tvib, &
                                                   & ad_kf, bd_kf, cd_kf, a_n2_dis,b_n2_dis,c_n2_dis, nb_bins, & 
                                                   & ed_n2, gn, gn2
      USE mod_nitrogen_TTLTH_CFD_prop,                 ONLY: Q_trans, Q_int_rr_ho, get_energy
      implicit none

      INTEGER :: is
      REAL(KIND=8), PARAMETER :: ratio = 2.5d0
      REAL(KIND=8), PARAMETER :: Zr = 235.d-1
      REAL(KIND=8), PARAMETER :: C  = 111.d0
      REAL(KIND=8), PARAMETER :: T0 = 540.99d0
      REAL(KIND=8), PARAMETER :: mu0 = 1781.d-8
      REAL(KIND=8) :: sum1, sum2, tmp
      REAL(KIND=8) :: T, Tr, TTv, ln_T, ln_Tv, ln_TTv, Tvl, Tvh, Tv
      REAL(KIND=8) :: omega_RT, omega_CR
      REAL(KIND=8) :: tauc, visc, a, b, Trank, denom
      REAL(KIND=8) :: r, rv, dexp_r, dexp_rv, e, ev, cv!, evl, cvl, evh, cvh, el, eh
      REAL(KIND=8) :: sigma, RT, T_m13, p, millikan, park, tau_VT_am, tau_VT_mm, tau_VT, tau_RT
      REAL(KIND=8) :: kfd, kfdinv, kfdh, kfd_eq, keq, exp_1, exp_2
      REAL(KIND=8) :: xn, xn2l, xn2h, xn2
      REAL(KIND=8), DIMENSION(nb_ns) :: rhoit

      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: tvec, rhoi
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: omega
	  
      real(kind=8), dimension(nb_tvib) :: aTv, xn2v, expr
	  real(kind=8), dimension(nb_tvib) :: kd_Ng, kd_N2g, kr_Ng, kr_N2g, VCd_N,	VCd_N2,	VCr_N,	VCr_N2
	  real(kind=8), dimension(nb_tvib, nb_tvib) :: KVT_N,	 KVT_N2,  VTPLUS_N, VTPLUS_N2, VTMOINS_N, VTMOINS_N2
	  real(kind=8), dimension(nb_tvib) :: omega_VT, Omega_CV,  VV_chi, VV_NRJ
!real(kind=8) :: kr_N, kr_N2, kd_N, kd_N2, OmCV_N_R, OmCV_N2_R, OmCV_N_D, OmCV_N2_D, OmVT_N, OmVT_N2!, omega_VT, Omega_CV
	integer :: i, j, szfv
	
!	      real(kind=8) :: kd_l_N, kd_l_N2, kr_l_N, kr_l_N2, kd_h_N, kd_h_N2, kr_h_N, kr_h_N2
!	  real(kind=8) :: VCd_l_N,VCd_l_N2,VCr_l_N,VCr_l_N2,VCd_h_N,VCd_h_N2,VCr_h_N,VCr_h_N2

	
	omega = 0.d0
	  
	  ! Temperatures and other useful data
      T      = tvec(1)
      Tr     = tvec(posTr)
	  aTv    = tvec(posTv:posTv+nb_tvib-1)
!print*, "source :", aTv
!pause
      ! Application of a fix in order to avoid numerical problems   
      IF (T.LT.Tmin)  T = Tmin
      IF (Tr.LT.Tmin) Tr = Tmin

!	  call get_chemistry_TTLTH(T,aTv(1),aTv(2), kd_l_N, kd_l_N2, kr_l_N, kr_l_N2, kd_h_N, kd_h_N2, kr_h_N, kr_h_N2, &
!	  VCd_l_N, VCd_l_N2, VCr_l_N, VCr_l_N2, VCd_h_N, VCd_h_N2, VCr_h_N, VCr_h_N2)

	  call get_chemistry(T, aTv, kd_Ng,	kd_N2g,	kr_Ng,	kr_N2g, &
								 VCd_N,	VCd_N2,	VCr_N,	VCr_N2)
									 
	  call get_VT(T, aTv,  KVT_N,	 KVT_N2,  VTPLUS_N, VTPLUS_N2, VTMOINS_N, VTMOINS_N2)
	  
	  call get_VV(T, aTv, rhoi,  VV_chi, VV_NRJ)

	  
    ! Chemical source terms
      ! Species concentrations (a fix is applied in order to avoid numerical problems)
      xn       = rhoi(1)/mm_N
      rhoit(1) = MAX(xn*mm_N,rhoi_tol) 

      xn2v      = rhoi(2:1+nb_tvib)/mm_N2
      xn2		= sum(xn2v)
	  rhoit(2) = MAX(xn2*mm_N2,rhoi_tol) 


      !!!!! Mass production terms
      ! Dissociation by N impact
	  do i=1,nb_tvib
		expr(i)    = xn*(kd_Ng(i)*xn2v(i) - xn*xn*kr_Ng(i))
		omega(1+i)  = - mm_N2*expr(i)
	  enddo
      omega(1) = 2.d0*mm_N*sum(expr)

      ! Dissociation by N2 impact (de-activated for comparison against the NASA Ames database for the N-N2 system)
      IF (diss_imp_N2.EQV..TRUE.) THEN
		  do i=1,nb_tvib
			expr(i)    = xn2*(kd_N2g(i)*xn2v(i) - xn*xn*kr_N2g(i))
			omega(1+i)  = omega(1+i) - mm_N2*expr(i)
		  enddo
		  omega(1) = omega(1) + 2.d0*mm_N*sum(expr)
      ENDIF

	  ! VT transfer
	  do i=1,nb_tvib
		  do j=1,nb_tvib
			omega(1+i)  = omega(1+i) + mm_N2*xn*(KVT_N(j,i)*xn2v(j) - KVT_N(i,j)*xn2v(i))
		  enddo
	  enddo
	  
	  if (vt_imp_mol.EQV..TRUE.) then
		  do i=1,nb_tvib
			  do j=1,nb_tvib
				omega(1+i)  = omega(1+i) + mm_N2*xn2*(KVT_N2(j,i)*xn2v(j) - KVT_N2(i,j)*xn2v(i))
			  enddo
		  enddo
	  endif

	  ! VV transfer
	  do i=1,nb_tvib
!		  omega(1+i)  = omega(1+i) + VV_chi(i)
	  enddo


	  omega(nb_ns + 1:nb_ns + nb_dim + 1) = 0.d0

      !!!! Energy transfer terms
      Omega_VT = 0.d0
	  Omega_CV = 0.d0
			
	  IF (nb_tvib.GE.1) THEN
         ! Static pressure
         !p = (rhoit(1)*Rn + rhoit(2)*Rn2)*T

		! VT transfer in N2-N only
		IF (vt_imp_mol.EQV..FALSE.) THEN
		 	do i=1, nb_tvib
			  Omega_VT(i) = Omega_VT(i) + xn*(	  sum(VTPLUS_N(1:nb_tvib,i) * xn2v(1:nb_tvib))  - sum(VTMOINS_N(i,1:nb_tvib))  * xn2v(i)  )
			  !Omega_VT(i) = Omega_VT(i) + xn*(	  VTPLUS_N(i,i) * xn2v(i)  - VTMOINS_N(i,i))  * xn2v(i)  
			  Omega_CV(i) = VCr_N(i) * xn**3 - VCd_N(i) * xn2v(i)*xn
			enddo
		! VT transfer also in N2-N2 collisions
		ELSE
		 	do i=1, nb_tvib
			  Omega_VT(i) = Omega_VT(i) + xn* (	  sum(VTPLUS_N(1:nb_tvib,i)  * xn2v(1:nb_tvib)) - sum(VTMOINS_N(i,1:nb_tvib))   * xn2v(i)  )
			  Omega_VT(i) = Omega_VT(i) + xn2*(	  sum(VTPLUS_N2(1:nb_tvib,i) * xn2v(1:nb_tvib)) - sum(VTMOINS_N2(i,1:nb_tvib))  * xn2v(i)  )
			  !Omega_VT(i) = Omega_VT(i) + xn* (	  VTPLUS_N(i,i) * xn2v(i)  - VTMOINS_N(i,i)  * xn2v(i)  )
			  !Omega_VT(i) = Omega_VT(i) + xn2*(	  VTPLUS_N2(i,i) * xn2v(i) - VTMOINS_N2(i,i)  * xn2v(i)  )
			  Omega_CV(i) = VCr_N(i) * xn**3 - VCd_N(i) * xn2v(i)*xn + VCr_N2(i) * xn**2*xn2 - VCd_N2(i) * xn2v(i)*xn2
			enddo
		ENDIF

        do i=1, nb_tvib
			omega(nb_ns + nb_dim + posTv + i - 1)     =   Omega_CV(i) +  Omega_VT(i)! + VV_NRJ(i)
			
		enddo
		
!		print*, "source :", VV_chi
!		print *, "source : ", nb_tvib, nb_trot, vt_imp_mol, Omega_VT, Omega_CV
!		pause
		!print *, "source : ", 			omega(nb_ns + nb_dim + posTv + 1 - 1:nb_ns + nb_dim + posTv + nb_tvib-1)/omega(1+1:1+nb_tvib) 
!print*, Omega_CV, Omega_VT
      ENDIF



      ! TTLTHer model for rotational relaxation
      IF (nb_trot.EQ.1) THEN

         ! Static pressure
         p = (rhoit(1)*Rn + rhoit(2)*Rn2)*T

         ! Dynamic viscosity 
         a     = 0.555d0*T0 + C
         Trank = T*5.d0/9.d0
         b     = 0.555d0*Trank + C
         visc  = mu0*(a/b)*(Trank/T0)**1.5

         ! Collision time 
         tauc = upi*visc/(4.d0*p)

         tmp   = 80.d0/T
         denom = (1 + 0.5d0*(upi**1.5)*(tmp**0.5) + (0.25d0*upi**2 + upi)*tmp)

         ! TTLTHer expression for rotational relaxation time
         tau_RT = Zr*tauc/denom        

         omega_CR = Rn2*Tr*omega(2) 
         omega_RT = rhoit(2)*Rn2*(T - Tr)/tau_RT

         omega(nb_ns + nb_dim + posTr) = Omega_CR + Omega_RT

      ENDIF
	  
    END SUBROUTINE source




    !------------------------------------------------------!
    ! This suroutine computes the source term and its Jacobian due to collisional processes for the N-N2 system.
    ! For sake of simplicity the Jacobian is computed with respect to primitive variables. The transformaton in order
    ! to obtaind its expression in terms of conservative variables is performed outside.
!    SUBROUTINE source_Jac (rhoi, tvec, omega, domega_dp)
!	  use mod_general_data,							   only: cell_number, nb_cells, use_boite
!      USE mod_nitrogen_TTLTH_initialize_CFD,           ONLY: mm_N, mm_N2, fac_keq, m_ratio, ov_m_ratio, c1, c2, c3, c4, &
!		compteur_global, EkJ, nb_bins, ukb, nlim
!      USE mod_nitrogen_TTLTH_CFD_prop,                 ONLY: Q_trans_der, Q_int_rr_ho_der
!
!      INTEGER :: is, start
!      REAL(KIND=8) :: T, Tr, aTvl, aTvh, TTv, ln_T, ln_Tv, ln_TTv, omega_VT, omega_CV
!      REAL(KIND=8) :: kfd, kbd, kfd_eq, keq, exp_1, exp_2, exp_3, tmp1, tmp2
!      REAL(KIND=8) :: xn, xn2, xns 
!      REAL(KIND=8) :: df_dT, df_dTv, dkf_dT, dkf_dTv, dkeq_dT, fac_dkf, dkf_eq_dT, dQn_dT, dQint_dT, & 
!                          dexp_eq
!      
!      REAL(KIND=8), DIMENSION(nb_ns + nb_temp) :: dOmegaV_dP 
!      REAL(KIND=8), DIMENSION(nb_ns) :: rhoit
!
!      REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: tvec, rhoi
!      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: omega
!      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: domega_dp
!
!	  real(KIND=8), DIMENSION(nb_temp) :: tvec_pr
!      real(kind=8), dimension(:), allocatable :: omg, omgp, rhoip
!      real(kind=8) :: dXi, dT, dTv, Qvl, Qvh, Tv
!      integer	::	szomg, i, j
!	  real(kind=8), dimension(3) :: rhoii, tvecc
!
!      ! Initialization
!      domega_dp = 0.d0
!	  szomg = size(omega)
!      	  
!      T      = tvec(1)
!      Tr     = tvec(posTr) 
!      aTvl     = tvec(posTv)
!	  aTvh     = tvec(posTv+1)
!
!	  allocate(omg(szomg))
!      allocate(omgp(szomg))
!      allocate(rhoip(nb_ns))
!	  omg     = 0.d0
!
!	  use_boite = .true.
!      call source (rhoi, tvec, omg)
!	  use_boite = .FALSE.
!	  
!	! Dérivées par rapport aux concentrations	  
!      do i=1,nb_ns ! boucle sur les espèces
!        	
!        	if(rhoi(i)<1.d-50) then
!        		dXi = 1.d-50
!        	else
!        		dXi = 1.d-6*rhoi(i)
!       		endif    
!            
!			rhoip    = rhoi
!        	rhoip(i) = rhoi(i) + dXi	
!        	
!			call source(rhoip,tvec,omgp)
!
!			do j=1,nb_ns ! boucle sur les sources à dériver: N, N2, Tv
!                    domega_dp((j-1)*(nb_ns+nb_temp)+i) = (omgp(j)-omg(j))/dXi		    ! N, N2 der, N N2
!            enddo
!			domega_dp(nb_ns    *(nb_ns+nb_temp)+i) = (omgp(nb_ns + nb_dim + posTv)    -omg(nb_ns + nb_dim + posTv))/dXi ! Tv    der, N N2
!			domega_dp((nb_ns+1)*(nb_ns+nb_temp)+i) = (omgp(nb_ns + nb_dim + posTv + 1)-omg(nb_ns + nb_dim + posTv + 1))/dXi ! Tv    der, N N2
!      enddo
!     
!	! Dérivées par rapport à T	   
!	   dT = 1.d-6*T
!		tvec_pr(1) = T+dT
!		tvec_pr(2) = aTvl
!		tvec_pr(3) = aTvh
!        
!		call source(rhoi,tvec_pr,omgp)
!	    
!		do j=1,nb_ns
!           domega_dp((j-1)*(nb_ns+nb_temp)+nb_ns+1) = (omgp(j)-omg(j))/dT
!        enddo
!		domega_dp(nb_ns*    (nb_ns+nb_temp)+nb_ns+1) = (omgp(nb_ns + nb_dim + posTv)-omg(nb_ns + nb_dim + posTv))/dT		! Tvl	
!		domega_dp((nb_ns+1)*(nb_ns+nb_temp)+nb_ns+1) = (omgp(nb_ns + nb_dim + posTv+1)-omg(nb_ns + nb_dim + posTv+1))/dT	! Tvh
!
!	! Dérivées par rapport à Tvl
!       if(aTvl==0.d0) then;	dTv = 1.d-9; else;	dTv = 1.d-6*aTvl; endif
!	   	tvec_pr(1) = T
!		tvec_pr(2) = aTvl+dTv
!        tvec_pr(3) = aTvh
!		
!		call source(rhoi,tvec_pr,omgp)
!	    
!		do j=1,nb_ns
!           domega_dp((j-1)*(nb_ns+nb_temp)+nb_ns+2) = (omgp(j)-omg(j))/dTv
!        enddo
!		domega_dp(nb_ns    *(nb_ns+nb_temp)+nb_ns+2) = (omgp(nb_ns + nb_dim + posTv)-omg(nb_ns + nb_dim + posTv))/dTv
!		domega_dp((nb_ns+1)*(nb_ns+nb_temp)+nb_ns+2) = (omgp(nb_ns + nb_dim + posTv+1)-omg(nb_ns + nb_dim + posTv+1))/dTv
!
!	! Dérivées par rapport à Tvh
!       if(aTvh==0.d0) then;	dTv = 1.d-9; else;	dTv = 1.d-6*aTvh; endif
!	   	tvec_pr(1) = T
!		tvec_pr(2) = aTvl
!        tvec_pr(3) = aTvh+dTv
!		
!		call source(rhoi,tvec_pr,omgp)
!	    
!		do j=1,nb_ns
!           domega_dp((j-1)* (nb_ns+nb_temp)+nb_ns+3) = (omgp(j)-omg(j))/dTv
!        enddo
!		domega_dp((nb_ns)*  (nb_ns+nb_temp)+nb_ns+3) = (omgp(nb_ns + nb_dim + posTv)  -omg(nb_ns + nb_dim + posTv))/dTv
!		domega_dp((nb_ns+1)*(nb_ns+nb_temp)+nb_ns+3) = (omgp(nb_ns + nb_dim + posTv+1)-omg(nb_ns + nb_dim + posTv+1))/dTv
!
!        omega = omg
!        
!        deallocate(omg)
!        deallocate(omgp)
!        deallocate(rhoip)
!
!!if(cell_number==1) then !.or.cell_number==nb_cells-1
!!print*, "sourceJac :: "
!!print*, domega_dp
!!print*,  omega(1), omega(2)+omega(3), omega(4:5), omega(6)+omega(7) 
!!print*,  rhoi(1), rhoi(2)+rhoi(3)
!!print*,  tvec
!!endif
!
!compteur_global = compteur_global + 1
!    END SUBROUTINE source_Jac


    !------------------------------------------------------!
    ! This suroutine computes the source term and its Jacobian due to collisional processes for the N-N2 system.
    ! For sake of simplicity the Jacobian is computed with respect to primitive variables. The transformaton in order
    ! to obtaind its expression in terms of conservative variables is performed outside.
    SUBROUTINE source_Jac (rhoi, tvec, omega, domega_dp)
!	  use mod_general_data,							   only: cell_number, nb_cells, use_boite
      USE mod_nitrogen_TTLTH_initialize_CFD,           ONLY: mm_N, mm_N2, fac_keq, m_ratio, ov_m_ratio, c1, c2, c3, c4, &
		compteur_global, EkJ, nb_bins, ukb, nb_tvib !, nlim
      USE mod_nitrogen_TTLTH_CFD_prop,                 ONLY: Q_trans_der, Q_int_rr_ho_der

      INTEGER :: is, start
      REAL(KIND=8) :: T, Tr, aTvl, aTvh, TTv, ln_T, ln_Tv, ln_TTv, omega_VT, omega_CV
      REAL(KIND=8) :: kfd, kbd, kfd_eq, keq, exp_1, exp_2, exp_3, tmp1, tmp2
      REAL(KIND=8) :: xn, xn2, xns 
      REAL(KIND=8) :: df_dT, df_dTv, dkf_dT, dkf_dTv, dkeq_dT, fac_dkf, dkf_eq_dT, dQn_dT, dQint_dT, & 
                          dexp_eq
      
      REAL(KIND=8), DIMENSION(nb_ns + nb_temp) :: dOmegaV_dP 
      REAL(KIND=8), DIMENSION(nb_ns) :: rhoit

      REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: tvec, rhoi
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: omega
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: domega_dp

	  real(KIND=8), DIMENSION(nb_temp) :: tvec_pr
      real(kind=8), dimension(:), allocatable :: omg, omgp, rhoip
      real(kind=8) :: dXi, dT, dTv, Qvl, Qvh, Tv
      integer	::	szomg, i, j
	  !real(kind=8), dimension(nb_ns) :: rhoii, tvecc
	  real(kind=8), dimension(nb_tvib) :: aTv

      ! Initialization
      domega_dp = 0.d0
	  szomg = size(omega)
      	  
      T      = tvec(1)
      Tr     = tvec(posTr) 
      aTv    = tvec(posTv:posTv+nb_tvib-1)
	  !aTvh     = tvec(posTv+1)

	  allocate(omg(szomg))
      allocate(omgp(szomg))
      allocate(rhoip(nb_ns))
	  omg     = 0.d0

!	  use_boite = .true.
      call source (rhoi, tvec, omg)
!	  use_boite = .FALSE.
	  
	! Dérivées par rapport aux concentrations	  
      do i=1,nb_ns ! boucle sur les espèces
        	
        	if(rhoi(i)<1.d-50) then
        		dXi = 1.d-50
        	else
        		dXi = 1.d-6*rhoi(i)
       		endif    
            
			rhoip    = rhoi
        	rhoip(i) = rhoi(i) + dXi	
        	
			call source(rhoip,tvec,omgp)

			do j=1,nb_ns ! boucle sur les sources à dériver: N, N2, Tv
				domega_dp((j-1)*(nb_ns+nb_temp)+i) = (omgp(j)-omg(j))/dXi		    ! N, N2 der, N N2
            enddo
			
			do j=1, nb_tvib
				domega_dp((nb_ns+j-1)    *(nb_ns+nb_temp)+i) = (omgp(nb_ns + nb_dim + posTv+j-1)    -omg(nb_ns + nb_dim + posTv+j-1))/dXi ! Tv    der, N N2
			enddo
			!domega_dp((nb_ns+1)*(nb_ns+nb_temp)+i) = (omgp(nb_ns + nb_dim + posTv + 1)-omg(nb_ns + nb_dim + posTv + 1))/dXi ! Tv    der, N N2
      enddo
     
	! Dérivées par rapport à T	   
	   dT = 1.d-6*T
		tvec_pr(1) = T+dT
		tvec_pr(2:1+nb_tvib) = aTv
		!tvec_pr(3) = aTvh
        
		call source(rhoi,tvec_pr,omgp)
	    
		do j=1,nb_ns
           domega_dp((j-1)*(nb_ns+nb_temp)+nb_ns+1) = (omgp(j)-omg(j))/dT
        enddo
		
		do j=1, nb_tvib
			domega_dp((nb_ns+j-1)*    (nb_ns+nb_temp)+nb_ns+1) = (omgp(nb_ns + nb_dim + posTv+j-1)-omg(nb_ns + nb_dim + posTv+j-1))/dT		! Tvl	
		enddo
		!domega_dp((nb_ns+1)*(nb_ns+nb_temp)+nb_ns+1) = (omgp(nb_ns + nb_dim + posTv+1)-omg(nb_ns + nb_dim + posTv+1))/dT	! Tvh

	! Dérivées par rapport aux temperatures Tv,i
	do i=1, nb_tvib
       if(aTv(i)==0.d0) then;	dTv = 1.d-9; else;	dTv = 1.d-6*aTv(i); endif
	   	tvec_pr(1) = T
		tvec_pr(2:1+nb_tvib) = aTv
        tvec_pr(1+i) = tvec_pr(1+i) + dTv
		
		call source(rhoi,tvec_pr,omgp)
	    
		do j=1,nb_ns
           domega_dp((j-1)*(nb_ns+nb_temp)+nb_ns+1+i) = (omgp(j)-omg(j))/dTv
        enddo
		
		do j=1, nb_tvib
			domega_dp((nb_ns+j-1)    *(nb_ns+nb_temp)+nb_ns+1+i) = (omgp(nb_ns + nb_dim + posTv+j-1)-omg(nb_ns + nb_dim + posTv+j-1))/dTv
		enddo
		!domega_dp((nb_ns+1)*(nb_ns+nb_temp)+nb_ns+2) = (omgp(nb_ns + nb_dim + posTv+1)-omg(nb_ns + nb_dim + posTv+1))/dTv
	enddo
	! Dérivées par rapport à Tvh
    !   if(aTvh==0.d0) then;	dTv = 1.d-9; else;	dTv = 1.d-6*aTvh; endif
	!   	tvec_pr(1) = T
!		tvec_pr(2) = aTvl
 !       tvec_pr(3) = aTvh+dTv
!		
!		call source(rhoi,tvec_pr,omgp)
!	    
!		do j=1,nb_ns
 !          domega_dp((j-1)* (nb_ns+nb_temp)+nb_ns+3) = (omgp(j)-omg(j))/dTv
  !      enddo
	!	domega_dp((nb_ns)*  (nb_ns+nb_temp)+nb_ns+3) = (omgp(nb_ns + nb_dim + posTv)  -omg(nb_ns + nb_dim + posTv))/dTv
	!	domega_dp((nb_ns+1)*(nb_ns+nb_temp)+nb_ns+3) = (omgp(nb_ns + nb_dim + posTv+1)-omg(nb_ns + nb_dim + posTv+1))/dTv

        omega = omg
        
        deallocate(omg)
        deallocate(omgp)
        deallocate(rhoip)

!if(cell_number==1) then !.or.cell_number==nb_cells-1
!print*, "sourceJac :: "
!print*, domega_dp
!print*,  omega(1), omega(2)+omega(3), omega(4:5), omega(6)+omega(7) 
!print*,  rhoi(1), rhoi(2)+rhoi(3)
!print*,  tvec
!endif

compteur_global = compteur_global + 1
    END SUBROUTINE source_Jac

    !------------------------------------------------------!
    ! This subroutine computes the source term due to VT and CV transfer in collisions between N and N2 species
    SUBROUTINE energy_transfer (rhoi, T, Tv, omega_n2, Omega_VT, Omega_CV, domega_n2_dP, dOmegaV_dP)

      USE mod_nitrogen_TTLTH_initialize_CFD,         ONLY: mm_N, mm_N2, mu, theta_vib, mw_a, mw_b
	 ! use mod_nitrogen_TTLTH_CFD_prop,				only: get_CV, get_VT

      REAL(KIND=8), PARAMETER :: ratio = 2.5d0
      REAL(KIND=8) :: e, ev, sigma, p, millikan, park, tau_VT, tau_VT_am, tau_VT_mm, dtauVT_dn, dtauVT_dn2, dtauVT_dT
      REAL(KIND=8) :: dtauVT_am_dn, dtauVT_am_dn2, dtauVT_mm_dn, dtauVT_mm_dn2, dtauVT_am_dT, dtauVT_mm_dT
      REAL(KIND=8) :: sum1, sum2, dsum2_dn, dsum2_dn2, dsum2_dT, sum2s
      REAL(KIND=8) :: tmp1, tmp2, tmp3, tmp4
      REAL(KIND=8) :: cp, cpv, OmegaVT_div_tauVT, fac_1, fac_2, Tr, Trv, dexp_Tr, dexp_Trv, RT, T_m13
      REAL(KIND=8) :: tauVT_div_p, dmillikan_dT, dpark_dT, dsigma_dT

      REAL(KIND=8), INTENT(IN) :: T, Tv, omega_n2
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi
      REAL(KIND=8), INTENT(OUT) :: Omega_VT, Omega_CV
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: domega_n2_dP
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: dOmegaV_dP
	  real(kind=8) :: Nn, Nn2
	  real(kind=8) :: OmVT_N, OMVT_N2, OMCV_N_R, OMCV_N2_R, OMCV_N_D, OMCV_N2_D

      Nn  = rhoi(1)/mm_N
	  Nn2 = rhoi(2)/mm_N2
	  
      ! Static pressure
      !RT = urg*T
      !p = (rhoi(1)/mm_N + rhoi(2)/mm_N2)*RT

      ! VT transfer in N-N2 collisions only
      !call get_VT(T,Tv,OmVT_N,OmVT_N2)
	  !call get_CV(T,Tv,OmCV_N_R,OmCV_N2_R,OmCV_N_D,OmCV_N2_D)
		 
	  IF (vt_imp_mol.EQV..FALSE.) THEN
		 
		! Omega_VT = OmVT_N * Nn2*Nn
        ! Omega_CV = OmCV_N_R * Nn**3 - OmCV_N_D * Nn2*Nn

       ! VT transfer also in N2-N2 collisions
	  ELSE

	!	 Omega_VT = OmVT_N * Nn2*Nn + OmVT_N2 * Nn2**2
       !  Omega_CV = OmCV_N_R * Nn**3 + OmCV_N2_R * Nn**2*Nn2 - OmCV_N_D * Nn2*Nn - OmCV_N2_D * Nn2**2

      ENDIF

      ! Jacobian matrix row (vibrational energy conservation equation)
      ! Derivatives with respect to N and N2 species densities
      !dOmegaV_dP(1) = - OmegaVT_div_tauVT*dtauVT_dn  + domega_n2_dP(1)*ev
      !dOmegaV_dP(2) = - OmegaVT_div_tauVT*dtauVT_dn2 + Omega_VT/rhoi(2) + domega_n2_dP(2)*ev

      ! Derivatives with respect to T, Tr and Tv
      !tau_VT = 1.d0/tau_VT
      !dOmegaV_dP(3) =   rhoi(2)*cp*tau_VT  - OmegaVT_div_tauVT*dtauVT_dT + domega_n2_dP(3)*ev
      !dOmegaV_dP(4) = 0.d0
      !dOmegaV_dP(3 + nb_tvib) = - rhoi(2)*cpv*tau_VT + domega_n2_dP(4)*ev + omega_n2*cpv

		 Omega_VT = 0.d0
         Omega_CV = 0.d0



       dOmegaV_dP = 0.d0
    END SUBROUTINE energy_transfer

  END MODULE mod_nitrogen_TTLTH_CFD_source
!------------------------------------------------------------------------------!
