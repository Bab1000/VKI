!------------------------------------------------------------------------------!
! This module provides functions and subroutines dealing with thermodynamic, transport properties 
! and source terms for the N-N2 system when using the TTLTH multi-temperature model.
! The initialization provided in this module should be use only when the library is interfaced with CFD codes
  MODULE mod_nitrogen_TTLTH_initialize_CFD

    IMPLICIT NONE

    INTEGER, SAVE :: posTr, posTv
    INTEGER, SAVE :: nb_ns, nb_tvib, nb_trot, nb_te, nb_dim, nb_eq, nb_temp, nb_bins
    REAL(KIND=8), SAVE :: ukb, una, urg, ue, uh, upi
    REAL(KIND=8), SAVE :: gn, gn2, mm_N, mm_N2, Rn, Rn2, mu, ed_n2, hf_n, theta_rot, theta_vib, fac_Q, fac_keq, fac_exp, & 
                       &  m_ratio, ov_m_ratio
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: ek, EkJ    
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: ad_kf, bd_kf, cd_kf, a_n2_dis, b_n2_dis, c_n2_dis, a_vv, b_vv, c_vv
    !!REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: ae_kf, be_kf, ce_kf, a_n2_exc, b_n2_exc, c_n2_exc
	REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, SAVE :: a_kvt_n, b_kvt_n, c_kvt_n, a_kvt_n2, b_kvt_n2, c_kvt_n2
	REAL(KIND=8), SAVE :: c1, c2, c3, c4
    REAL(KIND=8), SAVE :: cv_rot
    REAL(KIND=8), DIMENSION(2), SAVE :: mw_a, mw_b
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: cv_tr 
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: q11, q12, q22, bs, cs
    CHARACTER*80, SAVE :: solver
    LOGICAL, SAVE :: diss_imp_N2 = .FALSE., vt_imp_mol = .FALSE.
    integer, parameter :: in1=10
	integer :: compteur_global
	! Tables
	!REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: table_vt, table_cv, table_rates
	!REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE, SAVE :: table_rec
	!REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: table_diss, table_vt
	
	REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, SAVE :: table_energy
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: tb_T, tb_Tvl, tb_Tvh, tb_aTv
	
	integer :: n_T!, nlim
	integer, dimension(:), allocatable :: nivx	 

    ! Subroutine for the initialization of the TTLTH N-N2 library 
    CONTAINS

      !------------------------------------------------------!
      ! This subroutine initializes the TTLTH N-N2 library
      SUBROUTINE initialize (in_solver, in_mixture, in_transf, in_reaction, in_path, & 
                           & ns, ntrot, ntvib, nte, neq, ndim, mass)
 
         INTEGER :: i, comp, length, len_path

         INTEGER, INTENT(IN) :: ns, ntrot, ntvib, nte, ndim, neq
         CHARACTER*(*), INTENT(IN) :: in_solver, in_mixture, in_transf, in_reaction, in_path
         REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: mass
		 real(kind=8)	::	dum_double, dum2
		 integer		::	dum_int, longu, ios, j, comp2
		character(len=1)	::	dumchar ! Aurélien
         CHARACTER*200 :: path
         real(kind=8), dimension(:,:), allocatable :: a_temp, b_temp, c_temp 
		
		 
		 n_T = 50
		 
         ! Common data
         nb_ns   = ns
         nb_trot = ntrot
         nb_tvib = ntvib
         nb_te   = nte
         nb_dim  = ndim
         nb_eq   = neq 
         nb_temp = 1 + ntvib + ntrot

         posTr = 1 + nb_trot
         posTv = 1 + nb_trot + 1!nb_tvib

         ! Solver name
         solver = in_solver
         length = LEN_TRIM(solver)

         ! Path to the directory where rate coefficient data are stored
         path     = ADJUSTL(in_path)
         len_path = LEN_TRIM(path)

         WRITE(*,20)solver(1:length)//':: Nitrogen TTLTH library -> initialization'
         PRINT*

         ! Flags for chemical reactions and energy transfer terms
         SELECT CASE (in_mixture)

           CASE('complete')
             diss_imp_N2 = .TRUE.
             WRITE(*,20)solver(1:length),':: Nitrogen TTLTH -> N2 + N2 = 2N + N2 reaction also considered'
             PRINT*

           CASE('reduced','none')
             WRITE(*,20)solver(1:length),':: Nitrogen TTLTH _> only N2 + N = 3N reaction considered'
             PRINT*

           CASE('eq')
             WRITE(*,20)solver(1:length),':: Nitrogen TTLTH _> equilibrium conditions - no finite rate chemistry'
             PRINT*

           CASE DEFAULT 
             WRITE(*,20)solver(1:length),':: Nitrogen TTLTH -> in mod_initialize_CFD, error in mixture selection'
             PRINT* 
             STOP

         END SELECT 

         ! Flags for energy transfer mechanism
         SELECT CASE (in_transf)

           CASE('complete')
             vt_imp_mol = .TRUE.
             WRITE(*,20)solver(1:length),':: Nitrogen TTLTH -> VT transfer in N2-N2 collisions also considered'
             PRINT*

           CASE('reduced','none')
             WRITE(*,20)solver(1:length),':: Nitrogen TTLTH -> only VT transfer in N-N2 collisions considered'
             PRINT*

           CASE('eq')
             WRITE(*,20)solver(1:length),':: Nitrogen TTLTH _> equilibrium conditions - no finite rate energy trasnfer'
             PRINT*

           CASE DEFAULT 
             WRITE(*,20)solver(1:length),':: Nitrogen TTLTH -> in mod_initialize_CFD, error in transfer mechanism selection'
             PRINT* 
             STOP

         END SELECT 
         
         ! Physical constants
         ukb = 1.380658d-23   
         una = 6.0221367d23   
         urg = ukb*una
         ue  = 1.602191d-19
         uh  = 6.626075d-34
         upi = 3.14159265d0 
   
         ! Molecular masses of N and N2 species
         mm_N  = 14.0067d-3
         mm_N2 = 28.0134d-3 
  
         ! Reduced mass of the N-N2 system
         mu = mm_N*mm_N2/(mm_N + mm_N2)

         ! Useful factors 
         m_ratio = mm_N/mm_N2
         ov_m_ratio = 1.d0/m_ratio

         ! N and N2 specific gas constant
         Rn  = urg/mm_N
         Rn2 = urg/mm_N2

         ! Spin of N atom nucleous
         gn = 4.d0
		 gn2 = 1.d0

         ! Data of N2 molecule
         ! Characteristic vibrational and rotational temperatures  
         theta_rot = 2.88d0
         theta_vib = 3392.7d0

         ! Formation enthalpy (per unit mass) and dissociation energy (in J/mol)     
         ed_n2 = 9.7538d0
         ed_n2 = ed_n2*ue 
         !ed_n2 = 113272.d0*ukb
         !hf_n = (113272.d0*ukb/2.d0)*una/mm_N 
         !ed_n2 = 79849.d0 * 1.9863d-23
		 hf_n = (ed_n2/2.d0)*una/mm_N 
		
         fac_exp = una 
        ! Loading energy level file only for equilibrium computations
           OPEN(UNIT=in1,FILE=path(1:len_path)//'Ek.dat',STATUS='old',IOSTAT=ios)
 
           IF (ios.NE.0) THEN
              WRITE(*,20)solver(1:length),':: "Ek.dat" file not found, in mod_nitrogen_TTLTH_initialize_CFD ..'
              PRINT*
              STOP
           ENDIF

           ! Number of bins (or vibrational levels)

        READ(in1,'(A1,I2)')dumchar, nb_bins 

           ALLOCATE(ek(nb_bins), EkJ(nb_bins))

           DO i = 1,nb_bins
              READ(in1,*)dum_int, dum2, dum_double
              EkJ(i) = dum_double*ue
              ek(i)  = dum_double*ue*una/mm_N2
           ENDDO
		!print*, "alpha :: ", ekJ/ue
           CLOSE(in1)

		 allocate(nivx(1+nb_tvib))
		 nivx(1) = 0
		 nivx(2) = 13
		 nivx(3) = 28
		 nivx(4) = nb_bins
!		 nivx(5) = 35
!		 nivx(6) = nb_bins
!		 nivx(7) = nb_bins
		print*, "nb_bins", nb_bins 
		 if(size(nivx,1)-1.ne.nb_tvib .or. size(nivx,1)-1.ne.(nb_ns-1)) then
			print*, "Problème dans le nb de groupes de vibration", size(nivx,1),nb_tvib, nb_ns
			stop
		 endif


       ! Loading dissociation rate coefficients 
           ALLOCATE(ad_kf(nb_bins), bd_kf(nb_bins), cd_kf(nb_bins))
           ALLOCATE(a_n2_dis(nb_bins), b_n2_dis(nb_bins), c_n2_dis(nb_bins))
        
           OPEN(UNIT=in1,FILE=path(1:len_path)//'kd_n.dat',STATUS='old',IOSTAT=ios) 

           IF (ios.NE.0) THEN
              WRITE(*,20)solver(1:length),':: "kd_n.dat" file not found, in mod_nitrogen_TTLTH_initialize_CFD ..'
              PRINT*
              STOP 
           ENDIF

           DO i = 1,nb_bins
              READ(in1,*)dum_int,ad_kf(i),bd_kf(i),cd_kf(i)
           ENDDO

           CLOSE(in1)
  
           ! Conversion [m^3/s] -> [m^3/mol/s]
           ad_kf = ad_kf*fac_exp
 
           OPEN(UNIT=in1,FILE=path(1:len_path)//'kd_n2.dat',STATUS='old',IOSTAT=ios) 

           IF (ios.NE.0) THEN
              WRITE(*,20)solver(1:length),':: "kd_n2.dat" file not found, in mod_nitrogen_TTLTH_initialize_CFD ..'
              PRINT*
              STOP 
           ENDIF

           DO i = 1,nb_bins
              READ(in1,*)dum_int,a_n2_dis(i),b_n2_dis(i),c_n2_dis(i)
           ENDDO

           CLOSE(in1)

           ! Conversion [m^3/s] -> [m^3/mol/s]
           a_n2_dis = a_n2_dis*fac_exp

           WRITE(*,20)solver(1:length),':: Nitrogen FHO -> loaded dissociation rate coefficients'
           PRINT*

           ! Loading excitation rate coefficients (fit coefficients)
           allocate(a_kvt_n(nb_bins,nb_bins), b_kvt_n(nb_bins,nb_bins), c_kvt_n(nb_bins,nb_bins))
		   allocate(a_kvt_n2(nb_bins,nb_bins), b_kvt_n2(nb_bins,nb_bins), c_kvt_n2(nb_bins,nb_bins))
		   
		   open(UNIT=in1,FILE=path(1:len_path)//'ke_n.dat',STATUS='old',IOSTAT=ios) 
           do i = 1, nb_bins
              do j = 1, nb_bins
                 read(in1,*) dum_int, dum_int, a_kvt_n(i,j), b_kvt_n(i,j), c_kvt_n(i,j)
              enddo
           enddo
           close(in1)
		   a_kvt_n = a_kvt_n * fac_exp
           
           ! n2n2  
           OPEN(UNIT=in1,FILE=path(1:len_path)//'ke_n2.dat',STATUS='old',IOSTAT=ios) 

           do i = 1, nb_bins
              do j = 1, nb_bins
                 read(in1,*) dum_int, dum_int, a_kvt_n2(i,j), b_kvt_n2(i,j), c_kvt_n2(i,j)
              enddo
           enddo
           close(in1)
		   a_kvt_n2 = a_kvt_n2 * fac_exp


           WRITE(*,20)solver(1:length),':: Nitrogen FHO -> loaded de-excitation rate coefficients'
           PRINT*

           ! VV process 
          ! IF (vv.EQV..TRUE.) THEN

              comp = 0
              DO i = 2,nb_bins 
                 DO j = i,nb_bins
                    comp = comp + 1
                 ENDDO
              ENDDO 

              ALLOCATE(a_vv(comp), b_vv(comp), c_vv(comp))

              OPEN(UNIT=in1,FILE=path(1:len_path)//'kvv.dat',STATUS='old',IOSTAT=ios)

              IF (ios.NE.0) THEN
                 WRITE(*,20)solver(1:length),':: "kvv.dat" file not found, in mod_nitrogen_TTLTH_initialize_CFD ..'
                 PRINT*
                 STOP 
              ENDIF

              DO i = 1,comp
                 READ(in1,*)dum_int,dum_int,a_vv(i),b_vv(i),c_vv(i)
              ENDDO
 
              CLOSE(in1)

              WRITE(*,20)solver(1:length),':: Nitrogen Capitelli -> loaded VV process rate coefficients'
              PRINT*

              ! Conversion [cm^3/s] -> [m^3/mol/s]
              a_vv = a_vv*fac_exp*1.d-6 ! ATTENTION au facteur 1d6 !!!

          ! ENDIF

! Loading NASA's rates
!
!              ! Loading dissociation rate coefficients (fit coefficients)
!              ALLOCATE(ad_kf(nb_bins), bd_kf(nb_bins), cd_kf(nb_bins))       
!
!              OPEN(UNIT=in1,FILE='../sources/nitrogen_NASA/data/kdf.dat',STATUS='old',IOSTAT=ios) 
!
!              IF (ios.NE.0) THEN
!                 WRITE(*,20)solver(1:length),':: "kfd.dat" file not found, in mod_nitrogen_NASA_initialize_CFD ..'
!                 PRINT*
!                 STOP 
!              ENDIF
!
!              DO i = 1,nb_bins 
!                 READ(in1,*)dum_int,ad_kf(i),bd_kf(i),cd_kf(i)
!              ENDDO
!
!              CLOSE(in1)
!
!              WRITE(*,20)solver(1:length),':: Nitrogen NASA -> loaded Da process rate coefficients'
!              PRINT*
!              
!              ad_kf = ad_kf*una*1e-6
!              
!              ALLOCATE(a_n2_dis(nb_bins), b_n2_dis(nb_bins), c_n2_dis(nb_bins))
!              a_n2_dis = ad_kf
!              b_n2_dis = bd_kf
!              c_n2_dis = cd_kf
!              
!              ! Loading de-excitation rate coefficients (fit coefficients)
!              comp = 0
!              DO i = 1,nb_bins 
!                 comp = comp + (nb_bins - i)
!              ENDDO 
!
!              ALLOCATE(ae_kf(comp), be_kf(comp), ce_kf(comp))
!
!              OPEN(UNIT=in1,FILE='../sources/nitrogen_NASA/data/kef.dat',STATUS='old',IOSTAT=ios) 
!
!              IF (ios.NE.0) THEN
!                 WRITE(*,20)solver(1:length),':: "kef.dat" file not found, in mod_nitrogen_NASA_initialize_CFD ..'
!                 PRINT*
!                 STOP 
!              ENDIF
!
!              DO i = 1,comp
!                READ(in1,*)dum_int,dum_int,ae_kf(i),be_kf(i),ce_kf(i)
!              ENDDO
!
!              CLOSE(in1)
!
!              WRITE(*,20)solver(1:length),':: Nitrogen NASA -> loaded VTa process rate coefficients'
!              PRINT*
!
!              ! Conversion [cm^3/s] -> [m^3/mol/s]  
!              ae_kf = ae_kf*una*1e-6
!!molecules
!              ALLOCATE(a_n2_exc(comp), b_n2_exc(comp), c_n2_exc(comp))
!              
!              a_n2_exc = ae_kf
!              b_n2_exc = be_kf
!              c_n2_exc = ce_kf
		      
         ! Dissociation reaction N2 + N = N + N + N
         c1 = 0.049816d0
         c2 = -1.60d0
         c3 = 113200.d0 

         ! Dissociation reaction N2 + N2 = N2 + N + N (only the pre-exponential factor needs to be set)
         c4 = 0.0116237d0

         ! VT interactions (parameter for Millikan and White formula)
         ! N2-N  interaction
         mw_a(1) = 180.D0
         mw_b(1) = 0.0262D0

         ! N2-N2 interaction
         mw_a(2) = 221.D0 
         mw_b(2) = 0.0290D0

         ! Common factors
         fac_Q   = (2.d0*upi*ukb/uh**2.d0/una)**1.5
   !      fac_exp = 1.d-6*una
         fac_keq = 1.d0/una*gn**2*(mm_N/mm_N2)**1.5

         ! Unit conversion 
         c1 = c1*1.d-6*una
         c4 = c4*1.d-6*una

         ! Translational and rotational specific heats (constant volume)
         mass(1) = mm_N
         mass(2:nb_ns) = mm_N2

         ALLOCATE(cv_tr(nb_ns))

         DO i = 1,nb_ns 
            cv_tr(i) = 1.5d0/mass(i) 
         ENDDO
         cv_rot = 1.d0/mass(2)
   
         cv_tr  = cv_tr*urg
         cv_rot = cv_rot*urg

         ! Curve fits for collision integrals
         comp = 4*(nb_ns + 1)
         ALLOCATE(q11(comp), q12(comp), q22(comp), bs(comp), cs(comp))  
 
         ! N-N interaction
         q11(1) = -5.544d-3 
         q11(2) = 1.0789d-1   
         q11(3) = -9.011d-1 
         q11(4) = 5.8357d0

         q22(1) = -6.464d-3
         q22(2) = 1.42760d-1   
         q22(3) = -1.2781d0
         q22(4) = 7.2135d0

         bs(1) = 0.d0 
         bs(2) = 4.88d-3
         bs(3) = -5.95d-2 
         bs(4) = 3.21d-1

         cs = 0.d0

         ! N2-N2 interaction
         q11(5) = -5.04358d-3
         q11(6) = 1.08079d-1
         q11(7) = -9.98634d-1 
         q11(8) = 6.80089d0

         q22(5) = -6.99424d-3 
         q22(6) = 1.54912d-1 
         q22(7) = -1.33680d0  
         q22(8) = 7.66252d0

         bs(5) = -1.83418d-3
         bs(6) = 4.44907d-2
         bs(7) = -3.35033d-1
         bs(8) = 8.97329d-1

         cs = 0.d0 

         ! N-N2 interaction
         q11(9)  = -1.51789d-3
         q11(10) = 2.30672d-2
         q11(11) = -4.22508d-1 
         q11(12) = 5.46824d0

         q22(9)  = -3.51846d-3 
         q22(10) = 7.14265d-2 
         q22(11) = -7.70251d-1
         q22(12) = 6.34262d0

         bs(9)  = -1.61616d-3
         bs(10) = 4.59611d-2 
         bs(11) = -3.90469d-1 
         bs(12) = 1.13622d0

         cs = 0.d0

20  FORMAT(A,A)

compteur_global = 1
 
    END SUBROUTINE initialize
	

	subroutine interpol(x,y,tabx,taby)
			real(kind=8), intent(out)	::	x
			real(kind=8), intent(in)	::	y
			real(kind=8), dimension(:) :: tabx, taby
			
			real(kind=8) :: a
			integer :: sztbl, j
		
		sztbl = size(tabx)
		
		if (y < minval(taby) .or. y >= maxval(taby)) then
			write(*,*) "Interpol :: Variable hors du domaine"
			write(*,*) "Variable y : ", y, "min et max ", taby(1), taby(sztbl)
			write(*,*) "Variable x : ", x, "min et max ", tabx(1), tabx(sztbl)
!			write(*,*) "Using xmax instead od true x", tabx(sztbl)
!			x = tabx(sztbl)
			stop
		else
		
		j=1;
		do while(.not.( taby(j) <= y .and. y < taby(j+1)))
			j=j+1
		enddo
		
		a = (tabx(j+1)-tabx(j)) / (taby(j+1) - taby(j))
		
		x  = tabx(j) + a * (y-taby(j))
		
		endif

		!	print*, "Interpol :: y, x ", y, taby(1), taby(sztbl), x		
		
	end subroutine interpol

	subroutine interpol_2d(x,y,z,tabx,taby,matz,nx,ny)
			real(kind=8), intent(out)	::	z
			real(kind=8), intent(in)	::	x, y
			integer, intent(in)			:: nx, ny
			real(kind=8), intent(in), dimension(nx)		:: tabx
			real(kind=8), intent(in), dimension(ny)		:: taby
			real(kind=8), intent(in), dimension(nx,ny)		:: matz
			
			real(kind=8) :: a, b
			integer	::	i, j
		
		if (x < minval(tabx) .or. x >= maxval(tabx) .or. y < minval(taby) .or. y >= maxval(taby)) then
			write(*,*) "Interpol_2d :: Variable hors du domaine"
			write(*,*) "Variable x : ", x, "min et max ", tabx(1), tabx(nx)
			write(*,*) "Variable y : ", y, "min et max ", taby(1), taby(ny)
			stop
		endif
		
		i=1;	j=1;
		do while(.not.( tabx(i) <= x .and. x < tabx(i+1)))
			i=i+1
		enddo
		do while(.not.( taby(j) <= y .and. y < taby(j+1)))
			j=j+1
		enddo
		
		a = (matz(i+1,j)-matz(i,j)) / (tabx(i+1) - taby(i))
		b = (matz(i,j+1)-matz(i,j)) / (taby(j+1) - taby(j))
		
		z  = matz(i,j) + a * (x-tabx(i)) +  b * (y-taby(i))

		!	print*, "Interpol :: y, x ", y, taby(1), taby(sztbl), x		
		
	end subroutine interpol_2d

	
	subroutine create_tables(Tmin, Tmax, aTvmin, aTvmax, n)
			integer, intent(in)	::	n
			real(kind=8), intent(in)	::	Tmin, Tmax, aTvmin, aTvmax
			
			integer :: i, j, k, l, ct, szfv
			real(kind=8) :: a, av, T, aTv, aTv_p_dTv, ev, ev_p_dev
			!!real(kind=8) :: Q_int, fac_Q_loc, Qn, Qn2, Keq
			real(kind=8), dimension(:), allocatable :: fv, fv_p_dfv

		!n_T = n
		
		allocate(table_energy(n,nb_tvib*2), tb_T(n), tb_aTv(n))
		  a  = (Tmax/Tmin)**(1.d0/dble(n-1))
		  av = (aTvmax-aTvmin) / (n-1)
		
		do i=1, n
			tb_T(i)   = Tmin*a**(i-1)
			tb_aTv(i) = aTvmin + av*(i-1)
		enddo
		
		!do i=1,n
		!	T   = tb_T(i)
			do j=1,n
				aTv = tb_aTv(j)
				do i=1,nb_tvib
					szfv = nivx(i+1)-nivx(i)
					if(allocated(fv)) deallocate(fv); if(allocated(fv_p_dfv)) deallocate(fv_p_dfv);
					allocate(fv(szfv));
					allocate(fv_p_dfv(szfv));
					
					if(aTv==0.d0) then;	aTv_p_dTv = 1.d-9;		else;		aTv_p_dTv = 1.000001d0*aTv;		endif	
					fv			= dexp(aTv*		  EkJ(nivx(i)+1:nivx(i+1)) /ukb) / sum(dexp(aTv*		 EkJ(nivx(i)+1:nivx(i+1))	/ukb))
					fv_p_dfv	= dexp(aTv_p_dTv* EkJ(nivx(i)+1:nivx(i+1)) /ukb) / sum(dexp(aTv_p_dTv* EkJ(nivx(i)+1:nivx(i+1))	/ukb))					
					
					ev			= sum(fv*		Ek(nivx(i)+1:nivx(i+1)))
					ev_p_dev	= sum(fv_p_dfv*	Ek(nivx(i)+1:nivx(i+1)))
					
					! Ev
					table_energy(j,i) = ev
					! dEv/daTv
					table_energy(j,nb_tvib+i) = (ev_p_dev - ev) / (aTv_p_dTv - aTv)
				enddo					
			enddo
			
			print*, "table_energy : ", table_energy(1,1:nb_tvib), table_energy(n,1:nb_tvib)
			!pause
		!enddo

! Vérification
!	open(unit=1,file="tb_T.dat")
!		write(1,'(50(E14.6,1x))') tb_T
!	close(1)
!
!	open(unit=1,file="tb_Tv.dat")
!		write(1,'(50(E14.6,1x))') tb_aTv
!	close(1)

!	open(unit=1,file="tb_Tvh.dat")
!		write(1,'(50(E14.6,1x))') tb_aTvh
!	close(1)

!	open(unit=1,file="energyl.dat")
!		write(1,'(50(E14.6,1x))') table_energy(:,1)
!	close(1)
!
!	open(unit=1,file="energym.dat")
!		write(1,'(50(E14.6,1x))') table_energy(:,2)
!	close(1)
!!
!	open(unit=1,file="energyh.dat")
!		write(1,'(50(E14.6,1x))') table_energy(:,nb_tvib)
!	close(1)
!
  
print *, "Creation of Tables completed"
	end subroutine create_tables

	subroutine get_VT(T, aTv,  KVT_N,	 KVT_N2, &
							   VTPLUS_N, VTPLUS_N2, VTMOINS_N, VTMOINS_N2)
	  		real(kind=8), intent(in)	::	T
			real(kind=8), dimension(nb_tvib), intent(in)	::	aTv
			real(kind=8), dimension(nb_tvib,nb_tvib), intent(out)	::	KVT_N,	 KVT_N2, &
							   VTPLUS_N, VTPLUS_N2, VTMOINS_N, VTMOINS_N2	
			real(kind=8), dimension(nb_bins) :: fvk!, fvl
			integer	::	i, j, k, l, comp, szfv
			real(kind=8), dimension(nb_bins, nb_bins)	::	kup, kup_n2
	
		! Taux calculés par detailed balancing
		kup    = a_kvt_n    * T**b_kvt_n    * dexp(-c_kvt_n/T) 
		kup_n2 = a_kvt_n2	* T**b_kvt_n2	* dexp(-c_kvt_n2/T)

! Taux excitation calculés par DB
!		do i=1,nb_bins-1
!			do j=i+1,nb_bins
!				kup(i,j)    = kup(j,i)    * exp(-(ekJ(j)-ekj(i))/(ukB*T))
!				kup_n2(i,j) = kup_n2(j,i) * exp(-(ekJ(j)-ekj(i))/(ukB*T))
!			enddo
!		enddo
! Taux désexcitation calculés par DB
		do i=2,nb_bins
			do j=1,i-1
				kup(i,j)    = kup(j,i)    * exp(-(ekJ(j)-ekj(i))/(ukB*T))
				kup_n2(i,j) = kup_n2(j,i) * exp(-(ekJ(j)-ekj(i))/(ukB*T))
			enddo
		enddo
				
				! Taux NASA
!				kup=0
!				kup_n2=0
!				
!				comp=1
!				do i=1,nb_bins
!					do j=i+1,nb_bins
!						kup(i,j)    = ae_kf(comp)       * T**be_kf(comp)    * dexp(-ce_kf(comp)/T)
!						kup(j,i)    = kup(i,j) * dexp( -(ekJ(i)-ekj(j))/(ukB*T))
!						kup_n2(i,j) = a_n2_exc(comp)	* T**b_n2_exc(comp)	* dexp(-c_n2_exc(comp)/T)
!						kup_n2(j,i) = kup_n2(i,j) * dexp( -(ekJ(i)-ekj(j))/(ukB*T))
!						comp = comp+1
!					enddo
!				enddo

		KVT_N		= 0.d0
		KVT_N2		= 0.d0
		VTPLUS_N	= 0.d0
		VTPLUS_N2	= 0.d0
		VTMOINS_N	= 0.d0
		VTMOINS_N2	= 0.d0
		
		! Précalcul fonctions partition
		
		do k=1, nb_tvib
			szfv = nivx(k+1)-nivx(k)
			fvk(nivx(k)+1:nivx(k+1)) = dexp(aTv(k)*EkJ(nivx(k)+1:nivx(k+1))/ ukb) 
			fvk(nivx(k)+1:nivx(k+1)) = fvk(nivx(k)+1:nivx(k+1)) / sum( fvk(nivx(k)+1:nivx(k+1)) )
			!fvk(nivx(k)+1:nivx(k+1)) = dexp(aTv(k)*EkJ(nivx(k)+1:nivx(k+1))/ ukb) / sum(dexp(aTv(k)*EkJ(nivx(k)+1:nivx(k+1))/ ukb))
		enddo
						
		do k=1, nb_tvib

			do l=1, nb_tvib

				! Constantes de reaction
				do i=nivx(k)+1,nivx(k+1)
!					do j=nivx(l)+1,nivx(l+1)
!						KVT_N(k,l)  = KVT_N(k,l)  + kup(i,j)    * fvk(i)
!						KVT_N2(k,l) = KVT_N2(k,l) + kup_n2(i,j) * fvk(i)
!					enddo
						KVT_N(k,l)  = KVT_N(k,l)  + sum(kup	  (i,nivx(l)+1:nivx(l+1)) * fvk(i))
						KVT_N2(k,l) = KVT_N2(k,l) + sum(kup_n2(i,nivx(l)+1:nivx(l+1)) * fvk(i))
				enddo

				! Termes sources VT
				do i=nivx(k)+1,nivx(k+1)
!					do j=nivx(l)+1,nivx(l+1)
!						!VTPLUS_N(l,k)  = VTPLUS_N(l,k)   + kup(j,i)    * fvl(j)	* Ek(i) * mm_N2
!						!VTPLUS_N2(l,k) = VTPLUS_N2(l,k)  + kup_n2(j,i) * fvl(j)	* Ek(i) * mm_N2
!						VTPLUS_N(l,k)  = VTPLUS_N(l,k)   + kup(j,i)    * fvk(j)	* Ek(i) * mm_N2
!						VTPLUS_N2(l,k) = VTPLUS_N2(l,k)  + kup_n2(j,i) * fvk(j)	* Ek(i) * mm_N2
!						VTMOINS_N(k,l) = VTMOINS_N(k,l)  + kup(i,j)    * fvk(i)	* Ek(i) * mm_N2
!						VTMOINS_N2(k,l)= VTMOINS_N2(k,l) + kup_n2(i,j) * fvk(i)	* Ek(i) * mm_N2
!					enddo
						VTPLUS_N(l,k)  = VTPLUS_N(l,k)   + sum(kup	  (nivx(l)+1:nivx(l+1),i) * fvk(nivx(l)+1:nivx(l+1)))	* Ek(i) * mm_N2
						VTPLUS_N2(l,k) = VTPLUS_N2(l,k)  + sum(kup_n2 (nivx(l)+1:nivx(l+1),i) * fvk(nivx(l)+1:nivx(l+1)))	* Ek(i) * mm_N2
						VTMOINS_N(k,l) = VTMOINS_N(k,l)  + sum(kup	  (i,nivx(l)+1:nivx(l+1)))* fvk(i)	* Ek(i) * mm_N2
						VTMOINS_N2(k,l)= VTMOINS_N2(k,l) + sum(kup_n2 (i,nivx(l)+1:nivx(l+1)))* fvk(i)	* Ek(i) * mm_N2

				enddo
			enddo
		enddo
!	print *, "get_VT ",T, aTv, VTPLUS_N, VTMOINS_N, VTPLUS_N2, VTMOINS_N2			

	end subroutine get_VT

	subroutine get_VV(T, aTv, rhoi, VV_chi,	VV_NRJ)
	  		real(kind=8), intent(in)	::	T
			real(kind=8), dimension(nb_tvib), intent(in)	::	aTv
			real(kind=8), dimension(nb_tvib+1), intent(in)	::	rhoi
			real(kind=8), dimension(nb_tvib), intent(out)	::	VV_chi,	VV_NRJ
				
!			real(kind=8), dimension(nb_bins, nb_bins)	::	kup, kup_n2
			real(kind=8), dimension(nb_bins) :: fvk, xi, omega, ek_const
			real(kind=8)	::	tmp1, tmp2, tmp3, kf, ov_keq, kb, exp_1, ln_T, ov_T	
			integer	::	i, j, k, pos, comp, szfv
			
		ln_T = dlog(T)
		ov_T = 1.d0/T
		
		
		! Précalcul fonctions partition
		do k=1, nb_tvib
			szfv = nivx(k+1)-nivx(k)
			fvk(nivx(k)+1:nivx(k+1)) = dexp(aTv(k)*EkJ(nivx(k)+1:nivx(k+1))/ ukb) 
			fvk(nivx(k)+1:nivx(k+1)) = fvk(nivx(k)+1:nivx(k+1)) / sum( fvk(nivx(k)+1:nivx(k+1)) )
			
			xi(nivx(k)+1:nivx(k+1)) = rhoi(1+k)/mm_N2 * fvk(nivx(k)+1:nivx(k+1))
		enddo
		
		tmp1   = 1.d0/(ukb*T)
        
        DO i = 1,nb_bins
           ek_const(i) = DEXP(EkJ(i)*tmp1)
        ENDDO
		
		VV_chi = 0.d0
		VV_NRJ = 0.d0
		omega  = 0.d0
				
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
!print*, 'getVV', rhoi
			 ! Excitation terms
			 exp_1 = mm_N2*(kf*tmp2*xi(j - 1) - kb*tmp3*xi(j))
			  
			 omega(i)	= omega(i)		- exp_1
			 omega(i-1) = omega(i-1)    + exp_1
			 omega(j)	= omega(j)		+ exp_1
			 omega(j-1) = omega(j-1)    - exp_1 
			 
		  ENDDO

	   ENDDO 

		do k=1, nb_tvib
			VV_chi( k ) = sum(  omega( nivx(k)+1 : nivx(k+1) )  )
			VV_NRJ( k ) = sum(  omega( nivx(k)+1 : nivx(k+1) ) * Ek( nivx(k)+1 : nivx(k+1) )  )
		enddo
								
	end subroutine get_VV


	subroutine get_macro(T, Tv, ev, cv, kb_glo, kb_glo_n2, kd_glo, kd_glo_n2, cv_rec, cv_rec_n2, cv_dis, cv_dis_N2, vt, vt_n2)
			real(kind=8), intent(in)	::	T, Tv
			real(kind=8), intent(out)	::	kb_glo, kb_glo_n2, kd_glo, kd_glo_n2, cv_rec, cv_rec_n2, cv_dis, cv_dis_N2, vt, vt_n2
			
			integer :: i, j, k, l, ct
			real(kind=8) :: Tv_p_dTv, tmp1, tmp2, ev, ev_p_dev, cv
			real(kind=8) :: Q_int, fac_Q, Qn, Qn2, Keq
			real(kind=8), dimension(nb_bins) :: fv, fv_p_dfv, kd, kd_n2, kb, kb_n2
!			real(kind=8), dimension(nb_bins*(nb_bins-1)/2) :: kup, kdown, kup_n2, kdown_n2 
			real(kind=8), dimension(nb_bins, nb_bins)	::	kup, kup_n2
		
		! Tables de chimie
		! Tables de transferts CV VT
				fv = dexp(-EkJ/(ukb*Tv)) / sum(dexp(-EkJ/(ukb*Tv)))
				Tv_p_dTv = 1.000001d0*Tv
				fv_p_dfv = dexp(-EkJ/(ukb*Tv_p_dTv)) / sum(dexp(-EkJ/(ukb*Tv_p_dTv)))
				
				ev = sum(fv*Ek)
					ev_p_dev = sum(fv_p_dfv*Ek)
				cv = (ev_p_dev - ev) / (Tv_p_dTv - Tv)
				!if(i==1.AND.J==1) print *, ad_kf(1), ad_kf(nb_bins)
				kd    = ad_kf    * T**bd_kf    * dexp(-cd_kf/T)
				kd_n2 = a_n2_dis * T**b_n2_dis * dexp(-c_n2_dis/T)
					
					Q_int = (T/(theta_rot*2.d0))
					fac_Q = (2*upi*ukB/uh**2/uNa)**1.5;
					Qn    = fac_Q*(mm_N*T)**1.5;
					Qn2   = fac_Q*(mm_N2*T)**1.5;
					keq   = 1/uNa*(gn*Qn)**2/(gn2*Qn2*Q_int)*exp(-ed_n2/(ukB*T));
!print*, 'keq', keq, dexp(-EkJ(30)/(ukb*Tv))
				kb    = kd    * dexp(-EkJ/(ukb*T)) / keq
				kb_n2 = kd_n2 * dexp(-EkJ/(ukb*T)) / keq
				
				kb_glo		= sum(kb)
				kb_glo_n2	= sum(kb_n2)
				kd_glo		= sum(fv*kd)
				kd_glo_n2	= sum(fv*kd_n2)
				
				cv_rec		= mm_N2*sum(Ek*kb)
				cv_rec_N2	= mm_N2*sum(Ek*kb_n2)
				cv_dis		= mm_N2*sum(Ek*fv*kd)
				cv_dis_N2	= mm_N2*sum(Ek*fv*kd_n2)
				
!				kup    = ae_kf    * T**be_kf    * dexp(-ce_kf/T) 
!				kup_n2 = a_n2_exc * T**b_n2_exc * dexp(-c_n2_exc/T)
!if(i==50.AND.J==50) print *, "create_ta :: ", kD(50), table_cv(i,j,3)
!				ct = 1
!				tmp1 = 0.d0
!				tmp2 = 0.d0
!				do k=1, nb_bins-1
!					do l=k+1,nb_bins
!						kdown(ct)    = kup(ct)    * dexp((Ekj(l)-Ekj(k)) / (ukb*T))
!						kdown_n2(ct) = kup_n2(ct) * dexp((Ekj(l)-Ekj(k)) / (ukb*T))
!						tmp1 = tmp1 + mm_N2*(kup(ct)   * fv(k) - kdown(ct)    * fv(l)) * (Ek(l)-Ek(k))
!						tmp2 = tmp2 + mm_N2*(kup_n2(ct)* fv(k) - kdown_n2(ct) * fv(l)) * (Ek(l)-Ek(k))
!						!if(i==50.AND.J==50) print *, "create_ta :: ", k, l, kup(ct), kup_n2(ct)
!						ct = ct + 1
!					enddo
!				enddo
				tmp1 = 0.d0
				tmp2 = 0.d0

				kup    = a_kvt_n    * T**b_kvt_n    * dexp(-c_kvt_n/T) 
				kup_n2 = a_kvt_n2 * T**b_kvt_n2 * dexp(-c_kvt_n2/T)

				do i=1,nb_bins-1
					do j=i+1,nb_bins
						kup(i,j)    = kup(j,i)    * exp(-(ekJ(j)-ekj(i))/(ukB*T))
						kup_n2(i,j) = kup_n2(j,i) * exp(-(ekJ(j)-ekj(i))/(ukB*T))
					enddo
				enddo
				
				do i=1,nb_bins
					do j=1,nb_bins
						tmp1 = tmp1 + mm_N2*(kup(j,i)*fv(j)-kup(i,j)*fv(i))*Ek(i)
						tmp2 = tmp2 + mm_N2*(kup_n2(j,i)*fv(j)-kup_n2(i,j)*fv(i))*Ek(i)
					enddo
				enddo
				
				vt		= tmp1
				vt_n2	= tmp2


	end subroutine get_macro


	subroutine get_chemistry(T, aTv, kdg_N,	kdg_N2,	kr_N,	kr_N2, &
									 VCd_N,	VCd_N2,	VCr_N,	VCr_N2)
	  		real(kind=8), intent(in)	::	T
			real(kind=8), dimension(nb_tvib), intent(out)	::	kdg_N,	kdg_N2,	kr_N,	kr_N2, &
																VCd_N,	VCd_N2,	VCr_N,	VCr_N2, &
																aTv			
			real(kind=8) :: Q_int, fac_Q_loc, Qn, Qn2, Keq
			real(kind=8), dimension(nb_bins) :: kd, kd_n2, kb, kb_n2
			real(kind=8), dimension(nb_bins) :: fv
			
			integer	::	i, szfv
		


			kd    = ad_kf    * T**bd_kf    * dexp(-cd_kf/T)
			kd_n2 = a_n2_dis * T**b_n2_dis * dexp(-c_n2_dis/T)
	!print*, "chemistry :: ", kd			
				Q_int = (T/(theta_rot*2.d0))
				fac_Q_loc = (2*upi*ukB/uh**2/uNa)**1.5;
				Qn    = fac_Q_loc*(mm_N*T)**1.5;
				Qn2   = fac_Q_loc*(mm_N2*T)**1.5;
				keq   = 1/uNa*(gn*Qn)**2/(gn2*Qn2*Q_int)*exp(-ed_n2/(ukB*T));

!print*, 'keq', keq, dexp(-EkJ(30)/(ukb*Tv))
			kb    = kd    * dexp(-EkJ/(ukb*T)) / keq
			kb_n2 = kd_n2 * dexp(-EkJ/(ukb*T)) / keq

		! Constantes de groupes
			do i=1, nb_tvib
				szfv = nivx(i+1)-nivx(i)
				fv(nivx(i)+1:nivx(i+1)) = dexp(aTv(i)*EkJ(nivx(i)+1:nivx(i+1))	/ ukb)
				fv(nivx(i)+1:nivx(i+1)) = fv(nivx(i)+1:nivx(i+1)) / sum(fv(nivx(i)+1:nivx(i+1)))
				!fv(nivx(i)+1:nivx(i+1)) = dexp(aTv(i)*EkJ(nivx(i)+1:nivx(i+1))	/ ukb) / sum(dexp(aTv(i)*EkJ(nivx(i)+1:nivx(i+1))	/ ukb))

				!if(i==1.AND.J==1) print *, ad_kf(1), ad_kf(nb_bins)
				
				! Constantes de réactions
				kr_N(i)		    = sum(kb   (nivx(i)+1:nivx(i+1)))
				kr_N2(i)		= sum(kb_n2(nivx(i)+1:nivx(i+1)))
				
				kdg_N(i)		= sum(kd   (nivx(i)+1:nivx(i+1))	* fv(nivx(i)+1:nivx(i+1)))
				kdg_N2(i)		= sum(kd_n2(nivx(i)+1:nivx(i+1))	* fv(nivx(i)+1:nivx(i+1)))
			!print*, "chemistry :: ", kd(1:nlim)
			!print*, "chemistry :: ", fvl
				! Couplage chimie-vibration
				VCr_N(i)		= mm_N2 * sum(kb   (nivx(i)+1:nivx(i+1))		* ek(nivx(i)+1:nivx(i+1)))
				VCr_N2(i)		= mm_N2 * sum(kb_n2(nivx(i)+1:nivx(i+1))		* ek(nivx(i)+1:nivx(i+1)))
				
				VCd_N(i)		= mm_N2 * sum(kd(nivx(i)+1:nivx(i+1))			* fv(nivx(i)+1:nivx(i+1)) * ek(nivx(i)+1:nivx(i+1)))
				VCd_N2(i)		= mm_N2 * sum(kd_n2(nivx(i)+1:nivx(i+1))		* fv(nivx(i)+1:nivx(i+1)) * ek(nivx(i)+1:nivx(i+1)))
			enddo

	end subroutine get_chemistry

    !------------------------------------------------------!
    ! This subroutine de-allocates vectors and nullifies eventually used pointers to function subroutines. 
    SUBROUTINE finalize ()

      USE mod_function_pointer_TTLTH

      ! Data de-allocation
	  !IF (ALLOCATED(table_rec))      DEALLOCATE(table_rec)
	 ! IF (ALLOCATED(table_diss))     DEALLOCATE(table_diss)
	  !IF (ALLOCATED(table_vt))       DEALLOCATE(table_vt)
	  IF (ALLOCATED(tb_T))           DEALLOCATE(tb_T)
	  IF (ALLOCATED(tb_Tvl))         DEALLOCATE(tb_Tvl)
	  IF (ALLOCATED(tb_Tvh))         DEALLOCATE(tb_Tvh)
	  IF (ALLOCATED(tb_aTv))        DEALLOCATE(tb_aTv)
!	  IF (ALLOCATED(tb_aTvh))        DEALLOCATE(tb_aTvh)
	  IF (ALLOCATED(table_energy))   DEALLOCATE(table_energy)
	  
        IF (ALLOCATED(ek))          DEALLOCATE(ek)
        IF (ALLOCATED(EkJ))         DEALLOCATE(EkJ)
        IF (ALLOCATED(ad_kf))     DEALLOCATE(ad_kf)
        IF (ALLOCATED(bd_kf))     DEALLOCATE(bd_kf)
        IF (ALLOCATED(cd_kf))     DEALLOCATE(cd_kf)
        IF (ALLOCATED(a_n2_dis))  DEALLOCATE(a_n2_dis)
        IF (ALLOCATED(b_n2_dis))  DEALLOCATE(b_n2_dis)
        IF (ALLOCATED(c_n2_dis))  DEALLOCATE(c_n2_dis)
!        IF (ALLOCATED(ae_kf))     DEALLOCATE(ae_kf)
!        IF (ALLOCATED(be_kf))     DEALLOCATE(be_kf)
!        IF (ALLOCATED(ce_kf))     DEALLOCATE(ce_kf)
!        IF (ALLOCATED(a_n2_exc))  DEALLOCATE(a_n2_exc)
!        IF (ALLOCATED(b_n2_exc))  DEALLOCATE(b_n2_exc)
!!        IF (ALLOCATED(c_n2_exc))  DEALLOCATE(c_n2_exc)

        IF (ALLOCATED(a_kvt_n))   DEALLOCATE(a_kvt_n)
        IF (ALLOCATED(b_kvt_n))   DEALLOCATE(b_kvt_n)
        IF (ALLOCATED(c_kvt_n))   DEALLOCATE(c_kvt_n)
        IF (ALLOCATED(a_kvt_n2))  DEALLOCATE(a_kvt_n2)
        IF (ALLOCATED(b_kvt_n2))  DEALLOCATE(b_kvt_n2)
        IF (ALLOCATED(c_kvt_n2))  DEALLOCATE(c_kvt_n2)
		
		
      IF (ALLOCATED(q11))                             DEALLOCATE(q11)
      IF (ALLOCATED(q12))                             DEALLOCATE(q12)
      IF (ALLOCATED(q22))                             DEALLOCATE(q22)
      IF (ALLOCATED(bs))                              DEALLOCATE(bs)
      IF (ALLOCATED(cs))                              DEALLOCATE(cs)

      ! Nullify pointers
      IF (ASSOCIATED(get_temp_TTLTH))                  NULLIFY(get_temp_TTLTH)
      IF (ASSOCIATED(get_species_cv_TTLTH))            NULLIFY(get_species_cv_TTLTH)
      IF (ASSOCIATED(get_species_energy_cv_TTLTH))     NULLIFY(get_species_energy_cv_TTLTH)

    END SUBROUTINE finalize 

  END MODULE mod_nitrogen_TTLTH_initialize_CFD
!------------------------------------------------------------------------------!
