!------------------------------------------------------------------------------!
! This module provides the subroutine for initializing the nitrogen FHO library for the N-N2 system. 
  MODULE mod_nitrogen_FHO_initialize_CFD

    IMPLICIT NONE

    INTEGER, PARAMETER :: nb_ground = 2
    INTEGER, SAVE :: nb_ns, nb_dim, nb_eq, nb_bins
    INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: v_Eff   
    REAL(KIND=8), SAVE :: ukb, una, urg, ue, uh, upi
    REAL(KIND=8), SAVE :: gn, cv_tr_N, cv_tr_N2, mm_N, mm_N2, R_N, R_N2, hf_N
    REAL(KIND=8), SAVE :: Edis, sigma, ov_sigma, theta_rot
    REAL(KIND=8), SAVE :: fac_exp, fac_Q, fac_keq
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: ek, EkJ
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: ad_kf, bd_kf, cd_kf, a_n2_dis, b_n2_dis, c_n2_dis
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: ae_kf, be_kf, ce_kf, a_n2_exc, b_n2_exc, c_n2_exc
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: a_vv, b_vv, c_vv
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: q11, q12, q13, q22
    LOGICAL, SAVE :: vTa, vTm, vv, da, dm
    CHARACTER*80, SAVE :: solver

    CONTAINS

      !----------------------------------------------------!
      ! This subroutine initializes the nitrogen FHO library for the N-N2 system
      SUBROUTINE initialize (in_solver, in_reaction, in_path, ns, ndim, neq, mass)

        INTEGER :: i, j, comp, length, ios, comp2, longu
        INTEGER :: dum_int, len_path
        INTEGER, PARAMETER :: in1 = 10
        REAL(KIND=8) :: dum_double, dum2
        CHARACTER*20 :: line
        CHARACTER*200 :: path
        real(kind=8), dimension(:,:), allocatable :: a_temp, b_temp, c_temp 

        INTEGER, INTENT(IN) :: ns, ndim, neq
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: mass
        CHARACTER*(*), INTENT(IN) :: in_solver, in_reaction, in_path
		character(len=1)	::	dumchar ! Aurélien

        ! Solver name
        solver = in_solver
        length = LEN_TRIM(solver)

        ! Common data
        nb_ns   = ns
        nb_dim  = ndim
        nb_eq   = neq 

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
  
        ! Molecular mass vector
        mass(1) = mm_N
        DO i = 2,nb_ns
           mass(i) = mm_N2 
        ENDDO 

        ! N and N2 specific gas constant
        R_N  = urg/mm_N
        R_N2 = urg/mm_N2

        ! Translational components of species constant volume and pressure specific heats (N and N2)
        cv_tr_N  = 1.5d0*R_N 
        cv_tr_N2 = 1.5d0*R_N2 

        ! Spin of N atom nucleous
        gn = 4.d0 
 
        ! Symmetry factor and characteristic rotational temperature 
        sigma     = 2.d0
        ov_sigma  = 1.d0/sigma
        theta_rot = 2.88d0

        ! Dissociation energy [J] and formation enthalpy of N atom [J/Kg]
        Edis = 9.7538d0
        Edis = Edis*ue 
        !Edis = 79849.d0 * 1.9863d-23
		hf_N = (Edis/2.d0)*una/mm_N 

        ! Common factors
        fac_Q   = (2.d0*upi*ukb/uh**2/una)**1.5
        fac_exp = una
        fac_keq = 1.d0/una*gn**2*(mm_N/mm_N2)**1.5 

        ! Path to the directory where rate coefficient data are stored
        path     = ADJUSTL(in_path)
        len_path = LEN_TRIM(path)

        WRITE(*,20)solver(1:length),':: Nitrogen FHO library -> initialization'
        PRINT*
        
        ! Loading energy level file only for equilibrium computations

           OPEN(UNIT=in1,FILE=path(1:len_path)//'Ek.dat',STATUS='old',IOSTAT=ios)
 
           IF (ios.NE.0) THEN
              WRITE(*,20)solver(1:length),':: "Ek.dat" file not found, in mod_nitrogen_FHO_initialize_CFD ..'
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

           CLOSE(in1)

		! Loading reaction choice file
		
        IF (in_reaction.EQ.'eq') THEN

        ELSE

           ! Reading the reaction file (collisional processes)
           OPEN(UNIT=in1,FILE=path(1:len_path)//in_reaction(1:LEN_TRIM(in_reaction)),STATUS='old',IOSTAT=ios)

           IF (ios.NE.0) THEN
              WRITE(*,20)solver(1:length),':: "reaction" file not found, in mod_nitrogen_FHO_initialize_CFD ..'
              PRINT*
              STOP 
           ENDIF

           ! vT N2 + N 
           READ(in1,*)line 
           READ(in1,*)vTa

           ! vT N2 + N2 
           READ(in1,*)line 
           READ(in1,*)vTm

           ! vv  N2 + N2 
           READ(in1,*)line 
           READ(in1,*)vv

           ! Dissociation (N impact) 
           READ(in1,*)line 
           READ(in1,*)da

           ! Dissociation (N2 impact)
           READ(in1,*)line 
           READ(in1,*)dm  

           CLOSE(in1)

           ! Loading dissociation rate coefficients 
           ALLOCATE(ad_kf(nb_bins), bd_kf(nb_bins), cd_kf(nb_bins))
           ALLOCATE(a_n2_dis(nb_bins), b_n2_dis(nb_bins), c_n2_dis(nb_bins))
        
           OPEN(UNIT=in1,FILE=path(1:len_path)//'kd_n.dat',STATUS='old',IOSTAT=ios) 

           IF (ios.NE.0) THEN
              WRITE(*,20)solver(1:length),':: "kd_n.dat" file not found, in mod_nitrogen_FHO_initialize_CFD ..'
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
              WRITE(*,20)solver(1:length),':: "kd_n2.dat" file not found, in mod_nitrogen_FHO_initialize_CFD ..'
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
           allocate(a_temp(nb_bins,nb_bins), b_temp(nb_bins,nb_bins), c_temp(nb_bins,nb_bins))
           open(UNIT=in1,FILE=path(1:len_path)//'ke_n.dat',STATUS='old',IOSTAT=ios) 
           do i = 1, nb_bins
              do j = 1, nb_bins
                 read(in1,*) dum_int, dum_int, a_temp(i,j), b_temp(i,j), c_temp(i,j)
              enddo
           enddo
           close(in1)

           ! Conversion to Alessandro's format
           longu = nb_bins*(nb_bins-1)/2
           ALLOCATE(ae_kf(longu), be_kf(longu), ce_kf(longu))
           comp2 = 1
           do i = 1, nb_bins-1
              do j = i+1, nb_bins
                 ae_kf(comp2) = a_temp(j,i)
                 be_kf(comp2) = b_temp(j,i)
                 ce_kf(comp2) = c_temp(j,i)
                 comp2 = comp2 + 1
              enddo
           enddo
!comp2 = 1
!do i = 1, nb_bins-1
!	print*, i, " : ", ae_kf(comp2)
!	comp2 = comp2+1
!enddo
           ! Conversion [m^3/s] -> [m^3/mol/s]
           ae_kf = ae_kf*fac_exp

           ! n2n2  
           OPEN(UNIT=in1,FILE=path(1:len_path)//'ke_n2.dat',STATUS='old',IOSTAT=ios) 

           do i = 1, nb_bins
              do j = 1, nb_bins
                 read(in1,*) dum_int, dum_int, a_temp(i,j), b_temp(i,j), c_temp(i,j)
              enddo
           enddo
           close(in1)

           ! Conversion to Alessandro's format
           ALLOCATE(a_n2_exc(longu), b_n2_exc(longu), c_n2_exc(longu))
           comp2 = 1
           do i = 1, nb_bins-1
              do j = i+1, nb_bins
                 a_n2_exc(comp2) = a_temp(j,i)
                 b_n2_exc(comp2) = b_temp(j,i)
                 c_n2_exc(comp2) = c_temp(j,i)  !! BUG found
                 comp2 = comp2 + 1
               enddo
           enddo

           ! Conversion [m^3/s] -> [m^3/mol/s]
           a_n2_exc = a_n2_exc*fac_exp

           WRITE(*,20)solver(1:length),':: Nitrogen FHO -> loaded de-excitation rate coefficients'
           PRINT*

           ! VV process 
!           IF (vv.EQV..TRUE.) THEN
!
!              comp = 0
!              DO i = 2,nb_bins 
!                 DO j = i,nb_bins
!                    comp = comp + 1
!                 ENDDO
!              ENDDO 
!
!              ALLOCATE(a_vv(comp), b_vv(comp), c_vv(comp))
!
!              OPEN(UNIT=in1,FILE='../sources/nitrogen_FHO/data/kvv.dat',STATUS='old',IOSTAT=ios)
! 
!              IF (ios.NE.0) THEN
!                 WRITE(*,20)solver(1:length),':: "kvv.dat" file not found, in mod_nitrogen_Capitelli_initialize_CFD ..'
!                 PRINT*
!                 STOP 
!              ENDIF
!
!              DO i = 1,comp
!                 READ(in1,*)dum_int,dum_int,a_vv(i),b_vv(i),c_vv(i)
!              ENDDO
! 
!              CLOSE(in1)
!
!              WRITE(*,20)solver(1:length),':: Nitrogen FHO -> loaded VV process rate coefficients'
!              PRINT*
!
!              ! Conversion [m^3/s] -> [m^3/mol/s]
!              a_vv = a_vv*fac_exp 
!
!           ENDIF  

           ! VV process 
           IF (vv.EQV..TRUE.) THEN

              comp = 0
              DO i = 2,nb_bins 
                 DO j = i,nb_bins
                    comp = comp + 1
                 ENDDO
              ENDDO 

              ALLOCATE(a_vv(comp), b_vv(comp), c_vv(comp))

              OPEN(UNIT=in1,FILE=path(1:len_path)//'kvv.dat',STATUS='old',IOSTAT=ios)

              IF (ios.NE.0) THEN
                 WRITE(*,20)solver(1:length),':: "kvv.dat" file not found, in mod_nitrogen_FHO_initialize_CFD ..'
                 PRINT*
                 STOP 
              ENDIF

              DO i = 1,comp
                 READ(in1,*)dum_int,dum_int,a_vv(i),b_vv(i),c_vv(i)
              ENDDO
 
              CLOSE(in1)

              WRITE(*,20)solver(1:length),':: Nitrogen FHO -> loaded VV process rate coefficients'
              PRINT*

              ! Conversion [cm^3/s] -> [m^3/mol/s]
              a_vv = a_vv*fac_exp*1.d-6 ! ATTENTION au facteur 1d6 !!!

           ENDIF


        ENDIF
        
        ! Curve fit for collision integrals
        comp = 6*(nb_ground + 1)
        ALLOCATE(q11(comp), q12(comp), q13(comp), q22(comp))

        ! N-N interaction  
        q11(1) = 1.9939d1 
        q11(2) = 9.0652d0   
        q11(3) = -8.4875d-1 
        q11(4) = 2.3809d0
        q11(5) = 3.8782d-2
        q11(6) = 4.8674d-1

        q22(1) = 2.3736d1
        q22(2) = 9.1364d2   
        q22(3) = -7.7981d-1
        q22(4) = 2.7874d0
        q22(5) = 4.9002d-2
        q22(6) = 4.6224d-1
        
        q12(1) = 9.3402d0 
        q12(2) = 1.0845d3   
        q12(3) = -1.2761d0 
        q12(4) = 7.9736d-1
        q12(5) = 3.5324d-2
        q12(6) = 4.4711d-1

        q13(1) = 5.5929d0
        q13(2) = 2.1522d3   
        q13(3) = -1.6966d0
        q13(4) = 4.4782d-1
        q13(5) = 2.4128d-2
        q13(6) = 4.4567d-1

        ! N2-N2 interaction 
        q11(7)  = 4.2447d2  
        q11(8)  = 2.5868d4   
        q11(9)  = -1.0330d0 
        q11(10) = 2.7547d1
        q11(11) = 1.d0
        q11(12) = 4.2785d-1

        q22(7)  = 7.5032d2
        q22(8)  = 4.3769d4   
        q22(9)  = 9.6115d-1
        q22(10) = 5.5575d1
        q22(11) = 7.3450d-1
        q22(12) = 4.7747d-1
        
        q12(7)  = 2.7590d2 
        q12(8)  = 4.6513d4   
        q12(9)  = -1.4489d0 
        q12(10) = 1.5340d1
        q12(11) = 8.8021d-1
        q12(12) = 4.1708d-1

        q13(7)  = 1.3522d2 
        q13(8)  = 6.1370d4   
        q13(9)  = -1.8121d0
        q13(10) = 7.5032d0
        q13(11) = 4.3667d-1
        q13(12) = 4.2663d-1
 
        ! N-N2 interaction  
        q11(13) = 2.3258d1  
        q11(14) = 1.9547d3   
        q11(15) = -1.2673d0 
        q11(16) = 9.9578d-1
        q11(17) = 1.0382d-1
        q11(18) = 4.4037d-1

        q22(13) = 2.0822d1 
        q22(14) = 9.8291d2   
        q22(15) = -1.0120d0
        q22(16) = 1.1670d0
        q22(17) = 6.8362d-2
        q22(18) = 4.3962d-1
        
        q12(13) = 2.3258d1 
        q12(14) = 1.4173d4   
        q12(15) = -1.8460d0 
        q12(16) = 1.2542d0
        q12(17) = 7.0239d-2
        q12(18) = 4.9160d-1

        q13(13) = 1.4489d-1
        q13(14) = 4.8346d2   
        q13(15) = -2.3505d0
        q13(16) = 9.0340d-3
        q13(17) = 3.2446d-4
        q13(18) = 5.3192d-1
        
20    FORMAT(A,A)

      END SUBROUTINE initialize 

      !------------------------------------------------------!
      ! This subroutine de-allocates vectors and nullifies eventually used pointers to function subroutines. 
      SUBROUTINE finalize ()

        ! Vector de-allocation 
        IF (ALLOCATED(ek))          DEALLOCATE(ek)
        IF (ALLOCATED(EkJ))         DEALLOCATE(EkJ)
        IF (ALLOCATED(ad_kf))     DEALLOCATE(ad_kf)
        IF (ALLOCATED(bd_kf))     DEALLOCATE(bd_kf)
        IF (ALLOCATED(cd_kf))     DEALLOCATE(cd_kf)
        IF (ALLOCATED(a_n2_dis))  DEALLOCATE(a_n2_dis)
        IF (ALLOCATED(b_n2_dis))  DEALLOCATE(b_n2_dis)
        IF (ALLOCATED(c_n2_dis))  DEALLOCATE(c_n2_dis)
        IF (ALLOCATED(ae_kf))     DEALLOCATE(ae_kf)
        IF (ALLOCATED(be_kf))     DEALLOCATE(be_kf)
        IF (ALLOCATED(ce_kf))     DEALLOCATE(ce_kf)
        IF (ALLOCATED(a_n2_exc))  DEALLOCATE(a_n2_exc)
        IF (ALLOCATED(b_n2_exc))  DEALLOCATE(b_n2_exc)
        IF (ALLOCATED(c_n2_exc))  DEALLOCATE(c_n2_exc)
        IF (ALLOCATED(a_vv))      DEALLOCATE(a_vv)
        IF (ALLOCATED(b_vv))      DEALLOCATE(b_vv)
        IF (ALLOCATED(c_vv))      DEALLOCATE(c_vv)
        IF (ALLOCATED(v_Eff))       DEALLOCATE(v_Eff)
        IF (ALLOCATED(q11))         DEALLOCATE(q11)
        IF (ALLOCATED(q12))         DEALLOCATE(q12)
        IF (ALLOCATED(q22))         DEALLOCATE(q22)
        IF (ALLOCATED(q13))         DEALLOCATE(q13)

      END SUBROUTINE finalize 

  END MODULE mod_nitrogen_FHO_initialize_CFD
!------------------------------------------------------------------------------!
