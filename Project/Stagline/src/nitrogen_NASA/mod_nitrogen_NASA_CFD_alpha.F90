!------------------------------------------------------------------------------!
! This module provides the subroutine for initializing the nitrogen NASA library for the N-N2 system. 
! It is designed such that the vibrationally specific bin Boltzmann and full collisional models can be used.
  MODULE mod_nitrogen_NASA_initialize_CFD

    IMPLICIT NONE

    INTEGER, PARAMETER :: levels    = 9390
    INTEGER, PARAMETER :: nb_ground = 2
    INTEGER, SAVE :: nb_ns, nb_trot, nb_tvib, nb_temp, nb_dim, nb_eq, nb_bins, nb_inter, nb_points, table_points
    INTEGER, SAVE :: nb_exc_proc, nb_Da_proc, nb_vTa_proc, nb_Vtm_proc, nb_vv_proc
    INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: level_bin
    INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: exc_prod, exc_reac
    REAL(KIND=8), PARAMETER :: step  = 10.d0
    REAL(KIND=8), PARAMETER :: inv_step = 1.d0/step
    REAL(KIND=8), PARAMETER :: T_min = 20.d0
    REAL(KIND=8), PARAMETER :: T_max = 90000.d0
    REAL(KIND=8), PARAMETER :: Tmin_CI = 250.d0
    REAL(KIND=8), PARAMETER :: Tmax_CI = 90000.d0
    REAL(KIND=8), PARAMETER :: dT_CI   = 10.d0
    REAL(KIND=8), SAVE :: gn, hf_n, Rn, Rn2, mm_N, mm_N2, cv_tr_n, cv_tr_n2, cp_tr_n, cp_tr_n2
    REAL(KIND=8), SAVE :: ukb, una, urg, ue, uh, upi, Edis 
    REAL(KIND=8), SAVE :: fac_Q, fac_exp, fac_keq, fac_mu, fac_Diff
    REAL(KIND=8), SAVE :: Tmin_dis, Tmin_exc
    REAL(KIND=8), SAVE :: theta_rot, theta_vib
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: T_store
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: T_table
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: Ri, mi, sqrt_mi, sqrt_mij
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: ek, delta_ek, degen, gk, EkJ
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: EvJ, gvJ 
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: ad_kf, bd_kf, cd_kf
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: ap_kf, bp_kf, cp_kf 
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: ae_kf, be_kf, ce_kf 
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: a_vtm, b_vtm, c_vtm 
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: a_vv, b_vv, c_vv  
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: q11, q12, q22, bs, cs
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: dqintdT, qint, cvint, eint
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, SAVE :: Omega_11_table, Omega_12_table, Omega_22_table
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, SAVE :: Astar_table, Bstar_table, Cstar_table
    CHARACTER*80, SAVE :: solver, model
    LOGICAL, SAVE :: vTa, vTm, vv, da, dm, exa, pred

    CONTAINS

      !----------------------------------------------------!
      ! This subroutine initializes the nitrogen NASA library for the N-N2 system
      SUBROUTINE initialize (in_solver, in_mixture, in_reaction, in_path, ns, ntrot, ntvib, ndim, neq, mass)

        INTEGER, PARAMETER :: in1 = 10
        INTEGER :: i, j, ij, ios, np, p
        INTEGER :: len_path
        INTEGER :: comp, dum_int, length, pos
        INTEGER :: pos1, pos2, pos3, pos4
        REAL(KIND=8) :: dum_double 
        REAL(KIND=8) :: edisgr, T, lnt, lnt2, lnt3
        REAL(KIND=8) :: Omega_11, Omega_12, Omega_22, Astar, Bstar, Cstar
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dqintdT_bins, qint_bins, eint_bins, cvint_bins
        CHARACTER*20 :: line
        CHARACTER*200 :: path

        INTEGER, INTENT(IN) :: ns, ntrot, ntvib, ndim, neq
        CHARACTER*(*), INTENT(IN) :: in_solver, in_mixture, in_reaction, in_path
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: mass

        ! Solver name
        solver = in_solver
        length = LEN_TRIM(solver)

        ! Common data
        nb_ns   = ns
        nb_dim  = ndim
        nb_eq   = neq 
        nb_trot = ntrot
        nb_tvib = ntvib
        nb_temp = 1 + nb_trot + nb_tvib

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
 
        ! Characteristic rotational and vibrational temperatures (from the rovibrational levels)
        theta_rot = 2.84311389931011d0
        theta_vib = 3352.20535548991d0
 
        ! Molecular mass vector
        mass(1) = mm_N
        DO i = 2,nb_ns
           mass(i) = mm_N2 
        ENDDO 

        ! N and N2 specific gas constants
        Rn  = urg/mm_N
        Rn2 = urg/mm_N2

        ! Degeneracy of N atom (nuclear and spin contributions)
        gn = 12.d0
 
        IF ((in_mixture.EQ.'VC_rr').OR.(in_mixture.EQ.'vc_rr')) THEN
           gn = 4.d0
        ENDIF

        ! Formation enthalpy of N atom [J/Kg]
        edisgr = -9.75377d0
        Edis   = -edisgr*ue
        hf_n   = (Edis/2.d0)*una/mm_N

        ! Common factors
        fac_Q   = (2.d0*upi*ukb/uh**2/una)**1.5
        fac_exp = 1.d-6*una
        fac_keq = 1.d0/una*gn**2*(mm_N/mm_N2)**1.5 

        ! Translational components of species constant volume and pressure specific heats (N and N2)
        cv_tr_n = 1.5d0*Rn 
        cp_tr_n = cv_tr_n + Rn

        cv_tr_n2 = 1.5d0*Rn2 
        cp_tr_n2 = cv_tr_n2 + Rn2

        ! Data allocation
        ALLOCATE(EvJ(levels), gvJ(levels))

        WRITE(*,20)solver(1:length),':: Nitrogen NASA library -> initialization'
        PRINT*

        ! Path to the directory where rate coefficient data are stored
        path     = ADJUSTL(in_path)
        len_path = LEN_TRIM(path) 
         
        ! Reading the rovibrational energy levels (vJ order)
        OPEN(UNIT=in1,FILE=path(1:len_path)//'LevOri.dat',STATUS='old',IOSTAT=ios) 

        CALL check_file (ios, 'LevOri.dat')

        DO i = 1,levels
           READ(in1,*)dum_int,dum_int,j,dum_double
           EvJ(i) = (dum_double - edisgr)*ue
           IF (MOD(j,2).EQ.0) THEN
              gvJ(i) = (2.d0*j + 1.d0)*6.d0
           ELSE 
              gvJ(i) = (2.d0*j + 1.d0)*3.d0
           ENDIF      
        ENDDO

        CLOSE(in1)

        ! Minimum temperatures for evaluation of dissociation and excitation rate coefficients
        OPEN(UNIT=in1,FILE=path(1:len_path)//'T_fix.dat',STATUS='old',IOSTAT=ios) 

        CALL check_file (ios, 'T_fix.dat')

        READ(in1,*)Tmin_dis
        READ(in1,*)Tmin_exc

        CLOSE(in1)

        ! Physical model being used
        model = in_mixture

        ! Data allocation for the model selected 
        SELECT CASE(model)
 
          !------------------------------------------------!
          ! RVC (Rovibrational collisional model)
          CASE('RVC','rvc')
            WRITE(*,20)solver(1:length),':: Nitrogen NASA library -> RVC model is being used'
            PRINT*

            ! Model name 
            model = 'RVC'

            ! Loading energy data file(s). A checking on the number of species specified in input is also performed.
            ! Base energy file "Ek.dat"
            OPEN(UNIT=in1,FILE=path(1:len_path)//'Ek.dat',STATUS='old',IOSTAT=ios)

            CALL check_file (ios, 'Ek.dat')

            ! Number of bins 
            READ(in1,*)nb_bins 

            ! Checking on the number of bins 
            CALL check_nb_bins (nb_ns, nb_bins)

            ! Number of dissociation processes
            nb_Da_proc = nb_bins

            ! Memory allocation
            ALLOCATE(level_bin(nb_bins), degen(nb_bins), delta_ek(nb_bins))
            ALLOCATE(ek(nb_bins), gk(nb_bins), EkJ(nb_bins))
            ALLOCATE(ad_kf(nb_bins), bd_kf(nb_bins), cd_kf(nb_bins))  
            ALLOCATE(ap_kf(nb_bins), bp_kf(nb_bins), cp_kf(nb_bins))

            DO i = 1,nb_bins 
              READ(in1,*)dum_int,gk(i),ek(i)
            ENDDO
            degen = gk    

            DO i = 1,nb_bins 
               level_bin(i) = i
            ENDDO
 
            ! Conversion [ev] -> [J/kg]
            EkJ = ek*ue
            ek  = ek*ue*una/mm_N2
            delta_ek = 0.d0

            CLOSE(in1)

            ! Dissociation rate coefficients 
            OPEN(UNIT=in1,FILE=path(1:len_path)//'kdf.dat',STATUS='old',IOSTAT=ios) 

            CALL check_file (ios, 'kdf.dat')

            DO i = 1,nb_bins 
               READ(in1,*)dum_int,ad_kf(i),bd_kf(i),cd_kf(i)
            ENDDO

            CLOSE(in1)

            WRITE(*,20)solver(1:length),':: Nitrogen NASA -> loaded dissociation rate coefficients'
            PRINT*

            ! Excitation rate coefficients
            OPEN(UNIT=in1,FILE=path(1:len_path)//'kef.dat',STATUS='old',IOSTAT=ios) 

            CALL check_file (ios, 'kef.dat')
 
            ! Number of excitation processes
            READ(in1,*)nb_exc_proc

            ! Memory allocation
            ALLOCATE(ae_kf(nb_exc_proc), be_kf(nb_exc_proc), ce_kf(nb_exc_proc))
            ALLOCATE(exc_prod(nb_exc_proc), exc_reac(nb_exc_proc))

            DO i = 1,nb_exc_proc
               READ(in1,*)exc_reac(i),exc_prod(i),ae_kf(i),be_kf(i),ce_kf(i)
            ENDDO

            CLOSE(in1)

            WRITE(*,20)solver(1:length),':: Nitrogen NASA -> loaded excitation rate coefficients'
            PRINT*

            ! Predissociation rate coefficients
            OPEN(UNIT=in1,FILE=path(1:len_path)//'kpf.dat',STATUS='old',IOSTAT=ios) 

            CALL check_file (ios, 'kpf.dat')

            DO i = 1,nb_bins 
               READ(in1,*)dum_int,ap_kf(i)
               ap_kf(i) = ap_kf(i)/2.419d-17
            ENDDO

            ! Conversion [cm^3/s] -> [m^3/mol/s]  
            ad_kf = ad_kf*fac_exp 
            ae_kf = ae_kf*fac_exp 

            CLOSE(in1)

            WRITE(*,20)solver(1:length),':: Nitrogen NASA -> loaded predissociation rate coefficients'
            PRINT*

          !------------------------------------------------!
          ! BRVC (Boltzmann rovibrational collisional model)
          CASE('BRVC','brvc')
            WRITE(*,20)solver(1:length),':: Nitrogen NASA library -> BRVC model is being used'
            PRINT*

            ! Model name
            model = 'BRVC'
            
            ! Number of point for thermodynamic properties look-up tables 
            np = INT((T_max - T_min)/step) + 1 
            nb_points = np

            ! Read collisional processes taken into account
            OPEN(UNIT=in1,FILE=path(1:len_path)//in_reaction(1:LEN_TRIM(in_reaction)), &
               & STATUS='old',IOSTAT=ios)

            CALL check_file (ios, in_reaction)

            ! Dissociation (N impact)
            READ(in1,*)line 
            READ(in1,*)da

            ! Excitation (N impact)
            READ(in1,*)line 
            READ(in1,*)exa

            ! Predissociation
            READ(in1,*)line 
            READ(in1,*)pred
 
            CLOSE(in1)

            ! Loading energy data file(s). A checking on the number of species specified in input is also performed.
            ! Base energy file "Ek.dat"
            OPEN(UNIT=in1,FILE=path(1:len_path)//'Ek.dat',STATUS='old',IOSTAT=ios)

            CALL check_file (ios, 'Ek.dat')

            ! Number of bins 
            READ(in1,*)nb_bins 

            ! Checking on the number of bins 
            CALL check_nb_bins (nb_ns, nb_bins)

            ! Number of collisional excitation processes       
            nb_Da_proc   = nb_bins    
            nb_exc_proc  = 0
            DO i = 1,nb_bins 
              nb_exc_proc = nb_exc_proc + (nb_bins - i)
            ENDDO 

            ! Memory allocation
            ALLOCATE(ek(nb_bins), gk(nb_bins), EkJ(nb_bins))
            ALLOCATE(level_bin(levels), degen(levels), delta_ek(levels))
            ALLOCATE(dqintdT(nb_bins*np), qint(nb_bins*np), eint(nb_bins*np), cvint(nb_bins*np), T_store(np)) 
            ALLOCATE(dqintdT_bins(nb_bins), qint_bins(nb_bins), eint_bins(nb_bins), cvint_bins(nb_bins))
            IF (da.EQV..TRUE.) THEN 
               ALLOCATE(ad_kf(nb_bins), bd_kf(nb_bins), cd_kf(nb_bins))
               ALLOCATE(ap_kf(nb_bins), bp_kf(nb_bins), cp_kf(nb_bins))   
            ENDIF
            IF (exa.EQV..TRUE.) THEN 
               ALLOCATE(ae_kf(nb_exc_proc), be_kf(nb_exc_proc), ce_kf(nb_exc_proc))
            ENDIF

            DO i = 1,nb_bins 
              READ(in1,*)dum_int,gk(i),ek(i)
            ENDDO
     
            ! Conversion [ev] -> [J/kg]
            EkJ = ek*ue
            ek  = ek*ue*una/mm_N2

            CLOSE(in1)

            ! Bin internal energy file
            OPEN(UNIT=in1,FILE=path(1:len_path)//'deltaEki.dat',STATUS='old',IOSTAT=ios)

            CALL check_file (ios, 'deltaEki.dat')

            DO i = 1,levels 
               READ(in1,*)level_bin(i),degen(i),delta_ek(i)
            ENDDO
 
            CLOSE(in1)

            IF (MAXVAL(level_bin).GT.nb_bins) THEN
               WRITE(*,20)solver(1:length),':: error in "deltaEki.dat" file, nb_bins NE MAXVAL(level_bin)'
               PRINT*
               STOP    
            ENDIF

            ! Conversion [eV] -> [J]
            delta_ek = ue*delta_ek

            ! Computing needed thermodynamic properties
            T   = T_min
            i   = 1
            pos = 0

            DO WHILE (T.LE.T_max)

               T_store(i) = T
               CALL int_en_cp (T, dqintdT_bins, qint_bins, eint_bins, cvint_bins)

               DO j = 1,nb_bins 
                  dqintdT(pos + j) = dqintdT_bins(j)
                  qint(pos + j)    = qint_bins(j)
                  eint(pos + j)    = eint_bins(j)
                  cvint(pos + j)   = cvint_bins(j)
               ENDDO

               T   = T + step
               i   = i + 1
               pos = pos + nb_bins 

            ENDDO

            ! Dissociation rate coefficients 
            ! N2(i) + N = 2N + N 
            IF (da.EQV..TRUE.) THEN

              OPEN(UNIT=in1,FILE=path(1:len_path)//'kdf.dat',STATUS='old',IOSTAT=ios) 

              CALL check_file (ios, 'kdf.dat')

              DO i = 1,nb_bins 
                 READ(in1,*)dum_int,ad_kf(i),bd_kf(i),cd_kf(i)
              ENDDO

              CLOSE(in1)

              ! Conversion [cm^3/s] -> [m^3/mol/s]  
              ad_kf = ad_kf*fac_exp 

              WRITE(*,20)solver(1:length),':: Nitrogen NASA -> loaded dissociation rate coefficients'
              PRINT*

            ENDIF

            ! Predissociation rate coefficients
            ! N2(i) = 2N
            IF (pred.EQV..TRUE.) THEN

               OPEN(UNIT=in1,FILE=path(1:len_path)//'kpf.dat',STATUS='old',IOSTAT=ios) 

               CALL check_file (ios, 'kpf.dat')

               DO i = 1,nb_bins 
                  READ(in1,*)dum_int,ap_kf(i),bp_kf(i),cp_kf(i)
                  ap_kf(i) = ap_kf(i)/2.419d-17
               ENDDO

               WRITE(*,20)solver(1:length),':: Nitrogen NASA -> loaded predissociation rate coefficients'
               PRINT* 

            ELSE

               ap_kf = 0.d0
               bp_kf = 0.d0
               cp_kf = 0.d0

            ENDIF

            ! Excitation rate coefficients
            ! N2(i) + N = N2(ip) + N2,  ip > i
            IF (exa.EQV..TRUE.) THEN

              OPEN(UNIT=in1,FILE=path(1:len_path)//'kef.dat',STATUS='old',IOSTAT=ios) 

              CALL check_file (ios, 'kef.dat')

              DO i = 1,nb_exc_proc
                READ(in1,*)dum_int,dum_int,ae_kf(i),be_kf(i),ce_kf(i)
              ENDDO

              CLOSE(in1)

              WRITE(*,20)solver(1:length),':: Nitrogen NASA -> loaded excitation rate coefficients'
              PRINT*

              ! Conversion [cm^3/s] -> [m^3/mol/s]  
              ae_kf = ae_kf*fac_exp 

            ENDIF

          !------------------------------------------------!
          ! VC (Vibrational collisional) model
          CASE('VC','vc')
            WRITE(*,20)solver(1:length),':: Nitrogen NASA library -> VC model is being used'
            PRINT*

            ! Model name
            model = 'VC'
 
            ! Number of point for thermodynamic properties look-up tables 
            np = INT((T_max - T_min)/step) + 1 
            nb_points = np

            ! Loading energy data file(s). A checking on the number of species specified in input is also performed.
            ! Base energy file "Ek.dat"
            OPEN(UNIT=in1,FILE=path(1:len_path)//'Ek.dat',STATUS='old',IOSTAT=ios)

            CALL check_file (ios, 'Ek.dat')

            ! Number of bins 
            READ(in1,*)nb_bins 

            ! Checking on the number of bins 
            CALL check_nb_bins (nb_ns, nb_bins)

            ! Number of collisional excitation processes  
            nb_Da_proc  = nb_bins            
            nb_vTa_proc = 0
            DO i = 1,nb_bins 
              nb_vTa_proc = nb_vTa_proc + (nb_bins - i)
            ENDDO 
            nb_exc_proc = nb_vTa_proc

            ! Memory allocation
            ALLOCATE(ek(nb_bins), gk(nb_bins), EkJ(nb_bins))
            ALLOCATE(level_bin(levels), degen(levels), delta_ek(levels))
            ALLOCATE(dqintdT(nb_bins*np), qint(nb_bins*np), eint(nb_bins*np), cvint(nb_bins*np), T_store(np)) 
            ALLOCATE(dqintdT_bins(nb_bins), qint_bins(nb_bins), eint_bins(nb_bins), cvint_bins(nb_bins))
            ALLOCATE(ad_kf(nb_Da_proc), bd_kf(nb_Da_proc), cd_kf(nb_Da_proc))  
            ALLOCATE(ae_kf(nb_vTa_proc), be_kf(nb_vTa_proc), ce_kf(nb_vTa_proc))

            DO i = 1,nb_bins 
              READ(in1,*)dum_int,gk(i),ek(i)
            ENDDO
     
            ! Conversion [ev] -> [J/kg]
            EkJ = ek*ue
            ek  = ek*ue*una/mm_N2

            CLOSE(in1)

            ! Bin internal energy file
            OPEN(UNIT=in1,FILE=path(1:len_path)//'deltaEki.dat',STATUS='old',IOSTAT=ios)

            CALL check_file (ios, 'deltaEki.dat')

            DO i = 1,levels 
               READ(in1,*)level_bin(i),degen(i),delta_ek(i)
            ENDDO
 
            CLOSE(in1)

            IF (MAXVAL(level_bin).GT.nb_bins) THEN
               WRITE(*,20)solver(1:length),':: error in "deltaEki.dat" file, nb_bins NE MAXVAL(level_bin)'
               PRINT*
               STOP    
            ENDIF

            ! Conversion [eV] -> [J]
            delta_ek = ue*delta_ek

            ! Computing needed thermodynamic properties
            T   = T_min
            i   = 1
            pos = 0

            DO WHILE (T.LE.T_max)

               T_store(i) = T
               CALL int_en_cp (T, dqintdT_bins, qint_bins, eint_bins, cvint_bins)

               DO j = 1,nb_bins 
                  dqintdT(pos + j) = dqintdT_bins(j)
                  qint(pos + j)    = qint_bins(j)
                  eint(pos + j)    = eint_bins(j)
                  cvint(pos + j)   = cvint_bins(j)
               ENDDO

               T   = T + step
               i   = i + 1
               pos = pos + nb_bins 

            ENDDO

            ! Dissociation rate coefficients 
            OPEN(UNIT=in1,FILE=path(1:len_path)//'kdf.dat',STATUS='old',IOSTAT=ios) 

            CALL check_file (ios, 'kdf.dat')

            DO i = 1,nb_bins 
               READ(in1,*)dum_int,ad_kf(i),bd_kf(i),cd_kf(i)
            ENDDO

            CLOSE(in1)

            WRITE(*,20)solver(1:length),':: Nitrogen NASA -> loaded Da rate coefficients'
            PRINT*

            ! Excitation rate coefficients
            OPEN(UNIT=in1,FILE=path(1:len_path)//'kef.dat',STATUS='old',IOSTAT=ios) 

            CALL check_file (ios, 'kef.dat')

            DO i = 1,nb_exc_proc
              READ(in1,*)dum_int,dum_int,ae_kf(i),be_kf(i),ce_kf(i)
            ENDDO

            CLOSE(in1)

            WRITE(*,20)solver(1:length),':: Nitrogen NASA -> loaded VTa rate coefficients'
            PRINT*

            ! Conversion [cm^3/s] -> [m^3/mol/s]  
            ad_kf = ad_kf*fac_exp 
            ae_kf = ae_kf*fac_exp 

          !------------------------------------------------!
          ! VC_rr (Vibrational collisional) model (use is made of the rigid rotor assumption) 
          CASE('VC_rr','vc_rr')
            WRITE(*,20)solver(1:length),':: Nitrogen NASA library -> VC_rr model is being used'
            PRINT*

            ! Model name
            model = 'VC_rr'

            ! Read collisional processes taken into account
            OPEN(UNIT=in1,FILE=path(1:len_path)//in_reaction(1:LEN_TRIM(in_reaction)), & 
               & STATUS='old',IOSTAT=ios)

            CALL check_file (ios, in_reaction)

            ! vt N2 + N
            READ(in1,*)line 
            READ(in1,*)vTa

            ! vt N2 + N2
            READ(in1,*)line 
            READ(in1,*)vTm

            ! vv N2 + N2
            READ(in1,*)line 
            READ(in1,*)vv
 
            ! Dissociation (N impact)
            READ(in1,*)line 
            READ(in1,*)da
 
            ! Dissociation (N2 impact)
            READ(in1,*)line 
            READ(in1,*)dm

            CLOSE(in1)

            ! Number of Da, vTa, vTm and vv processes
            nb_Da_proc  = nb_ns - 1 
            nb_vTm_proc = nb_ns - 2
            nb_vTa_proc = 0
            nb_vv_proc  = 0

            nb_vTa_proc = 0
            DO i = 1,nb_ns - 1
               nb_vTa_proc = nb_vTa_proc + (nb_ns - 1 - i)
            ENDDO 

            DO i = 2,nb_ns - 1
               DO j = i,nb_ns - 1
                  nb_vv_proc = nb_vv_proc + 1
               ENDDO
            ENDDO 

            ! Loading energy data file(s). A checking on the number of species specified in input is also performed.
            ! Vibrational energy file "Ek.dat"
            OPEN(UNIT=in1,FILE=path(1:len_path)//'Ek.dat',STATUS='old',IOSTAT=ios)

            CALL check_file (ios,'Ek.dat')

            ! Number of bins (or vibrational levels)
            READ(in1,*)nb_bins 

            ! Checking on the number of bins 
            CALL check_nb_bins (nb_ns, nb_bins)
            
            ALLOCATE(ek(nb_bins), gk(nb_bins), EkJ(nb_bins))
            ALLOCATE(ad_kf(nb_Da_proc), bd_kf(nb_Da_proc), cd_kf(nb_Da_proc))   
            ALLOCATE(ae_kf(nb_vTa_proc), be_kf(nb_vTa_proc), ce_kf(nb_vTa_proc))
            ALLOCATE(a_vtm(nb_vTm_proc), b_vtm(nb_vTm_proc), c_vtm(nb_vTm_proc))
            ALLOCATE(a_vv(nb_vv_proc), b_vv(nb_vv_proc), c_vv(nb_vv_proc))

            DO i = 1,nb_bins 
               READ(in1,*)dum_int,gk(i),ek(i)
            ENDDO
     
            ! Conversion [ev] -> [J/kg]
            EkJ = ek*ue
            ek  = ek*ue*una/mm_N2

            ! Da processes (N and N2 impact)
            IF (Da.EQV..TRUE.) THEN

               OPEN(UNIT=in1,FILE=path(1:len_path)//'kdf.dat',STATUS='old',IOSTAT=ios) 

               CALL check_file (ios,'kdf.dat')

               DO i = 1,nb_Da_proc
                  READ(in1,*)dum_int,ad_kf(i),bd_kf(i),cd_kf(i)
               ENDDO

               CLOSE(in1)

               WRITE(*,20)solver(1:length),':: Nitrogen NASA -> loaded Da and Dm rate coefficients'
               PRINT*
  
            ENDIF

            ! VTa processes (N impact)
            IF (vTa.EQV..TRUE.) THEN

               OPEN(UNIT=in1,FILE=path(1:len_path)//'kef.dat',STATUS='old',IOSTAT=ios) 

               CALL check_file (ios,'kef.dat')

               DO i = 1,nb_vTa_proc 
                  READ(in1,*)dum_int,dum_int,ae_kf(i),be_kf(i),ce_kf(i)
               ENDDO

               CLOSE(in1)

               WRITE(*,20)solver(1:length),':: Nitrogen NASA -> loaded VTa rate coefficients'
               PRINT*
  
            ENDIF

            ! VTm processes (N2 impact)
            IF (vTm.EQV..TRUE.) THEN

               OPEN(UNIT=in1,FILE=path(1:len_path)//'kvtm.dat',STATUS='old',IOSTAT=ios) 

               CALL check_file (ios,'kvtm.dat')

               DO i = 1,nb_vTm_proc
                  READ(in1,*)dum_int,a_vtm(i),b_vtm(i),c_vtm(i)
               ENDDO

               CLOSE(in1)

               WRITE(*,20)solver(1:length),':: Nitrogen NASA -> loaded VTm rate coefficients'
               PRINT*
  
            ENDIF

            ! VV processes 
            IF (vv.EQV..TRUE.) THEN

               OPEN(UNIT=in1,FILE=path(1:len_path)//'kvv.dat',STATUS='old',IOSTAT=ios) 

               CALL check_file (ios,'kvv.dat')

               DO i = 1,nb_vv_proc
                  READ(in1,*)dum_int,dum_int,a_vv(i),b_vv(i),c_vv(i)
               ENDDO

               CLOSE(in1)

               WRITE(*,20)solver(1:length),':: Nitrogen NASA -> loaded VV rate coefficients'
               PRINT*
  
            ENDIF

            ! Conversion [cm^3/s] -> [m^3/mol/s]  
            ad_kf = ad_kf*fac_exp 
            ae_kf = ae_kf*fac_exp  
            a_vtm = a_vtm*fac_exp  
            a_vv  = a_vv*fac_exp

          !------------------------------------------------!
          ! Multi-temperature model (2 temperatures T-Tint)
          CASE('MT_TTint','mt_TTint')
            WRITE(*,20)solver(1:length),':: Nitrogen NASA library -> MT_TTint macroscopic model is being used'
            PRINT*

            ! Model name  
            model = 'MT_TTint'

            ! Number of bins (energy levels of N2 molecule)
            nb_bins = levels 

            ! Loading energy data file(s). A checking on the number of species specified in input is also performed.
            ! Base energy file "Ek.dat"
            OPEN(UNIT=in1,FILE=path(1:len_path)//'Ek.dat',STATUS='old',IOSTAT=ios)

            CALL check_file (ios, 'Ek.dat')

            ! Number of bins 
            READ(in1,*)nb_bins 

            ! Number of dissociation processes
            !nb_exc_proc = 17229850 
            nb_Da_proc  = nb_bins

            ! Memory allocation
            ALLOCATE(level_bin(nb_bins), degen(nb_bins), delta_ek(nb_bins))
            ALLOCATE(ek(nb_bins), gk(nb_bins), EkJ(nb_bins))
            ALLOCATE(ad_kf(nb_bins), bd_kf(nb_bins), cd_kf(nb_bins))  
            ALLOCATE(ap_kf(nb_bins), bp_kf(nb_bins), cp_kf(nb_bins))

            DO i = 1,nb_bins 
              READ(in1,*)dum_int,gk(i),ek(i)
            ENDDO
            degen = gk    

            DO i = 1,nb_bins 
               level_bin(i) = i
            ENDDO
 
            ! Conversion [ev] -> [J/kg]
            EkJ = ek*ue
            ek  = ek*ue*una/mm_N2
            delta_ek = 0.d0

            CLOSE(in1)

            ! Dissociation rate coefficients 
            OPEN(UNIT=in1,FILE=path(1:len_path)//'kdf.dat',STATUS='old',IOSTAT=ios) 

            CALL check_file (ios, 'kdf.dat')

            DO i = 1,nb_bins 
               READ(in1,*)dum_int,ad_kf(i),bd_kf(i),cd_kf(i)
            ENDDO

            CLOSE(in1)

            WRITE(*,20)solver(1:length),':: Nitrogen NASA -> loaded dissociation rate coefficients'
            PRINT*

            ! Excitation rate coefficients
            OPEN(UNIT=in1,FILE=path(1:len_path)//'kef.dat',STATUS='old',IOSTAT=ios) 

            CALL check_file (ios, 'kef.dat')

            ! Number of excitation processes 
            READ(in1,*)nb_exc_proc

            ! Memotry allocation 
            ALLOCATE(ae_kf(nb_exc_proc), be_kf(nb_exc_proc), ce_kf(nb_exc_proc))
            ALLOCATE(exc_prod(nb_exc_proc), exc_reac(nb_exc_proc))

            DO i = 1,nb_exc_proc
               READ(in1,*)exc_reac(i),exc_prod(i),ae_kf(i),be_kf(i),ce_kf(i)
            ENDDO

            CLOSE(in1)

            WRITE(*,20)solver(1:length),':: Nitrogen NASA -> loaded excitation rate coefficients'
            PRINT*

            ! Predissociation rate coefficients
            OPEN(UNIT=in1,FILE=path(1:len_path)//'kpf.dat',STATUS='old',IOSTAT=ios) 

            CALL check_file (ios, 'kpf.dat')

            DO i = 1,nb_bins 
               READ(in1,*)dum_int,ap_kf(i)
               ap_kf(i) = ap_kf(i)/2.419d-17
            ENDDO

            ! Conversion [cm^3/s] -> [m^3/mol/s]  
            ad_kf = ad_kf*fac_exp 
            ae_kf = ae_kf*fac_exp 

            WRITE(*,20)solver(1:length),':: Nitrogen NASA -> loaded predissociation rate coefficients'
            PRINT*

          !------------------------------------------------!
          ! Multi-temperature model (3 temperatures T-Tr-Tv)
          CASE('MT_TTrTv','mt_TTrTv')
            WRITE(*,20)solver(1:length),':: Nitrogen NASA library -> MT_TTrTv macroscopic model is being used'
            PRINT*

            ! Model name  
            model = 'MT_TTrTv'

            PRINT*
            WRITE(*,10)'not implemented yet...'
            WRITE(*,10)'in "mod_nitrogen_NASA_CFD_initialize"...'
            STOP

          !------------------------------------------------!
          ! Equilibrium model (no finite rate kinetics)
          CASE('EQ','eq')
            WRITE(*,20)solver(1:length),':: Nitrogen NASA library -> Full equilibrium model (vJ levels) is being used'
            PRINT*

          !------------------------------------------------!
          CASE DEFAULT
            WRITE(*,20)solver(1:length),':: error in model selection, in mod_nitrogen_NASA_initialize_CFD ..'
            PRINT*
            STOP
          
        END SELECT 

        ! Curve fit for collision integrals
        nb_inter = (nb_ground + 1)*nb_ground/2
        comp     = 4*nb_inter
        ALLOCATE(q11(comp), q12(comp), q22(comp), bs(comp), cs(comp), mi(nb_ground), Ri(nb_ground)) 

        mi(1) = mm_N
        mi(2) = mm_N2

        Ri(1) = Rn
        Ri(2) = Rn2

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

        ! Generate look-up tables for collision integrals
        table_points = (Tmax_CI - Tmin_CI)/(dT_CI) + 1
        ALLOCATE(T_table(table_points))
        ALLOCATE(Omega_11_table(table_points,nb_inter))
        ALLOCATE(Omega_12_table(table_points,nb_inter))
        ALLOCATE(Omega_22_table(table_points,nb_inter))
        ALLOCATE(Astar_table(table_points,nb_inter))          
        ALLOCATE(Bstar_table(table_points,nb_inter))
        ALLOCATE(Cstar_table(table_points,nb_inter))
        ALLOCATE(sqrt_mi(nb_ground),sqrt_mij(nb_inter))

        ! Common factors for viscosity and binary diffusion coefficients
        fac_mu   = (5.d0/16.d0)*SQRT(upi*ukb/una)
        fac_Diff = (16.d0/3.d0)/SQRT(2.d0*upi*urg)

        ! Square roots of molecular masses and reduced molecular masses
        DO i = 1,nb_ground 

           sqrt_mi(i) = DSQRT(mi(i))

           DO j = i,nb_ground 

              ij = ((i - 1)*(2*nb_ground - i) + 2*j)/2 

              sqrt_mij(ij) = SQRT(mi(i)*mi(j)/(mi(i) + mi(j))) 

           ENDDO

        ENDDO

        ! Temperature table
        T_table(1) = Tmin_CI
        DO i = 2,table_points
           T_table(i) = T_table(i - 1) + dT_CI
        ENDDO

        WRITE(*,20)solver(1:length),':: Nitrogen NASA -> generating tables for collision integrals'
        PRINT*

        ! Collision integral look-up tables
        DO p = 1,table_points 

           lnt  = DLOG(T_table(p))
           lnt2 = lnt*lnt
           lnt3 = lnt2*lnt

           DO i = 1,nb_ground 

              DO j = i,nb_ground

                 ! Index
                 ij = ((i - 1)*(2*nb_ground - i) + 2*j)/2 

                 pos1 = 4*(ij - 1) + 1
                 pos2 = 4*(ij - 1) + 2
                 pos3 = 4*(ij - 1) + 3
                 pos4 = 4*(ij - 1) + 4

                 ! Collision integrals
                 Omega_11 = 1.d-20*DEXP(q11(pos1)*lnt3 + q11(pos2)*lnt2 + q11(pos3)*lnt + q11(pos4))   
                 Omega_22 = 1.d-20*DEXP(q22(pos1)*lnt3 + q22(pos2)*lnt2 + q22(pos3)*lnt + q22(pos4)) 
                 Astar    = Omega_22/Omega_11
                 Bstar    = DEXP(bs(pos1)*lnt3 + bs(pos2)*lnt2 + bs(pos3)*lnt + bs(pos4)) 
                 Cstar    = DEXP(cs(pos1)*lnt3 + cs(pos2)*lnt2 + cs(pos3)*lnt + cs(pos4)) 
                 Omega_12 = Bstar*Cstar

                 ! Fill the collision integral look-up tables 
                 Omega_11_table(p,ij) = Omega_11
                 Omega_12_table(p,ij) = Omega_12
                 Omega_22_table(p,ij) = Omega_22
                 Astar_table(p,ij)    = Astar
                 Bstar_table(p,ij)    = Bstar
                 Cstar_table(p,ij)    = Cstar

              ENDDO

           ENDDO

        ENDDO
 
10    FORMAT(I4)
20    FORMAT(A,A)
25    FORMAT(A)

      END SUBROUTINE initialize 

      !----------------------------------------------------!
      ! This subroutine check if the number of levels of the N2 molecule
      ! is consistent with what provided in the solver input file
      SUBROUTINE check_nb_bins (ns, nbins) 

        INTEGER, INTENT(IN) :: ns, nbins

        IF (nbins.NE.(ns - 1)) THEN 
          WRITE(*,10)solver(1:LEN_TRIM(solver))//' :: number of bins in Ek.dat file set to uncorrect value'
          PRINT*
          WRITE(*,10)solver(1:LEN_TRIM(solver))//' :: in mod_nitrogen_NASA_initialize_CFD ..'
          PRINT*
          STOP
       ENDIF 
        
10    FORMAT(A)

      END SUBROUTINE check_nb_bins

      !----------------------------------------------------!
      ! This subroutine checks if a given file has been opened correctly 
      SUBROUTINE check_file (ios, filename)

        INTEGER, INTENT(IN) :: ios
        CHARACTER*(*), INTENT(IN) :: filename

        IF (ios.NE.0) THEN 
           WRITE(*,10)solver(1:LEN_TRIM(solver))//':: '//'"'//filename(1:LEN_TRIM(filename))//'"'//& 
                    & ' file not found, in mod_nitrogen_NASA_initialize_CFD ..'
           PRINT*
           STOP
        ENDIF

10    FORMAT(A)

      END SUBROUTINE check_file 

      !----------------------------------------------------!
      ! This subroutine computes the internal energy and internal specific heat
      ! contributions needed when employing the vibrational collisional or
      ! boltzmann bin models
      SUBROUTINE int_en_cp (T, dqintdT, qint, eint, cpint)

        INTEGER :: i, pos
        REAL(KIND=8) :: fac1, fac2
        REAL(KIND=8) :: kbT, dE, g, dexp0
        REAL(KIND=8), DIMENSION(nb_bins) :: sum1, sum2, sum3

        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: dqintdT, qint, eint, cpint
 
        ! Initialization
        sum1 = 0.d0
        sum2 = 0.d0
        sum3 = 0.d0        

        kbT = 1.d0/(ukb*T)
        DO i = 1,levels
           pos   = level_bin(i)
           dE    = delta_ek(i) 
           g     = degen(i) 
           dexp0 = g*DEXP(-dE*kbT) 
           fac1  = dexp0*dE
           fac2  = fac1*dE
           sum1(pos) = sum1(pos) + dexp0
           sum2(pos) = sum2(pos) + fac1
           sum3(pos) = sum3(pos) + fac2
        ENDDO
 
        fac1 = 1.d0/(ukb*T**2)
        fac2 = una/mm_N2 
        DO i = 1,nb_bins 
           dqintdT(i)  = sum2(i)
           qint(i)     = sum1(i)
           eint(i)     = sum2(i)/sum1(i)
           cpint(i)    = fac1*(sum3(i)/sum1(i) - (sum2(i)/sum1(i))**2)
        ENDDO

        ! Conversion from [J] -> [J/kg]
        dqintdT = fac1*dqintdT
        eint    = fac2*eint
        cpint   = fac2*cpint

      END SUBROUTINE int_en_cp 

      !------------------------------------------------------!
      ! This subroutine de-allocates vectors and nullifies eventually used pointers to function subroutines. 
      SUBROUTINE finalize ()
 
        USE mod_function_pointer_NASA

        ! Data de-allocation
        IF (ALLOCATED(q11))                             DEALLOCATE(q11)
        IF (ALLOCATED(q12))                             DEALLOCATE(q12)
        IF (ALLOCATED(q22))                             DEALLOCATE(q22)
        IF (ALLOCATED(bs))                              DEALLOCATE(bs)
        IF (ALLOCATED(cs))                              DEALLOCATE(cs)
        IF (ALLOCATED(mi))                              DEALLOCATE(mi)
        IF (ALLOCATED(Ri))                              DEALLOCATE(Ri)
        IF (ALLOCATED(sqrt_mi))                         DEALLOCATE(sqrt_mi)
        IF (ALLOCATED(sqrt_mij))                        DEALLOCATE(sqrt_mij)
        IF (ALLOCATED(T_table))                         DEALLOCATE(T_table)
        IF (ALLOCATED(Omega_11_table))                  DEALLOCATE(Omega_11_table)
        IF (ALLOCATED(Omega_12_table))                  DEALLOCATE(Omega_12_table)
        IF (ALLOCATED(Omega_22_table))                  DEALLOCATE(Omega_22_table) 
        IF (ALLOCATED(ad_kf))                           DEALLOCATE(ad_kf)
        IF (ALLOCATED(bd_kf))                           DEALLOCATE(bd_kf)
        IF (ALLOCATED(cd_kf))                           DEALLOCATE(cd_kf)
        IF (ALLOCATED(ap_kf))                           DEALLOCATE(ap_kf)
        IF (ALLOCATED(bp_kf))                           DEALLOCATE(bp_kf)
        IF (ALLOCATED(cp_kf))                           DEALLOCATE(cp_kf)
        IF (ALLOCATED(ae_kf))                           DEALLOCATE(ae_kf)
        IF (ALLOCATED(be_kf))                           DEALLOCATE(be_kf)
        IF (ALLOCATED(ce_kf))                           DEALLOCATE(ce_kf)
        IF (ALLOCATED(a_vtm))                           DEALLOCATE(a_vtm)
        IF (ALLOCATED(b_vtm))                           DEALLOCATE(b_vtm)
        IF (ALLOCATED(c_vtm))                           DEALLOCATE(c_vtm)
        IF (ALLOCATED(a_vv))                            DEALLOCATE(a_vv)
        IF (ALLOCATED(b_vv))                            DEALLOCATE(b_vv)
        IF (ALLOCATED(c_vv))                            DEALLOCATE(c_vv)
        IF (ALLOCATED(T_store))                         DEALLOCATE(T_store)
        IF (ALLOCATED(ek))                              DEALLOCATE(ek)
        IF (ALLOCATED(EkJ))                             DEALLOCATE(EkJ)
        IF (ALLOCATED(EvJ))                             DEALLOCATE(EvJ)
        IF (ALLOCATED(gvJ))                             DEALLOCATE(gvJ)
        IF (ALLOCATED(gk))                              DEALLOCATE(gk)
        IF (ALLOCATED(exc_reac))                        DEALLOCATE(exc_reac)
        IF (ALLOCATED(exc_prod))                        DEALLOCATE(exc_prod)

        ! Nullify pointers
        IF (ASSOCIATED(get_temp_NASA))                  NULLIFY(get_temp_NASA)
        IF (ASSOCIATED(get_species_cv_NASA))            NULLIFY(get_species_cv_NASA)
        IF (ASSOCIATED(get_species_energy_NASA))        NULLIFY(get_species_energy_NASA)  
        IF (ASSOCIATED(get_mass_prod_terms_NASA))       NULLIFY(get_mass_prod_terms_NASA)
        IF (ASSOCIATED(get_mass_prod_terms_Jac_NASA))   NULLIFY(get_mass_prod_terms_Jac_NASA) 

      END SUBROUTINE finalize 

  END MODULE mod_nitrogen_NASA_initialize_CFD
!------------------------------------------------------------------------------!
