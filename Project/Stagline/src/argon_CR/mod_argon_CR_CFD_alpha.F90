!------------------------------------------------------------------------------!
! This module provides functions and subroutines dealing with thermodynamic, transport properties 
! and source terms for the electronic specific CR model for the em,Ar(i),Arp system.
! The initialization provided in this module should be use only when the library is interfaced with CFD codes
  MODULE mod_argon_CR_initialize_CFD

    IMPLICIT NONE

    INTEGER, SAVE :: levels_Ar, levels_Arp, lev_Ar, lev_Arp
    INTEGER, SAVE :: nb_ns, nb_tvib, nb_trot, nb_te, nb_dim, nb_eq, nb_temp, nb_points
    INTEGER, SAVE :: pos_em, pos_Ar, pos_Arp, pos_Te
    INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: Arp_lev
    REAL(KIND=8), PARAMETER :: gamma_e    = 5.d0/3.d0
    REAL(KIND=8), PARAMETER :: gamma_e_m1 = gamma_e - 1.d0 
    REAL(KIND=8), PARAMETER :: ge = 2.d0
    REAL(KIND=8), PARAMETER :: jc_Arp1 = 1.5d0
    REAL(KIND=8), PARAMETER :: jc_Arp2 = 0.5d0
    REAL(KIND=8), SAVE :: ukb, una, urg, ue, uh, upi, ueps0
    REAL(KIND=8), SAVE :: fac_Q, Eion, hf_Arp
    REAL(KIND=8), SAVE :: ahi_kf, bhi_kf, chi_kf, aei_kf, bei_kf, cei_kf 
    REAL(KIND=8), SAVE :: mm_em, mm_Ar, mm_Arp, cv_tr_em, cv_tr_Ar, cv_tr_Arp, cp_tr_em, cp_tr_Ar, & 
                        & cp_tr_Arp, R_em, R_Ar, R_Arp 
    REAL(KIND=8), SAVE :: deltaT, Tmin, Tmax
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: jc, ken, kin, tvec, Im, mm_species, Ri
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: gi_Ar, gi_Arp, ei_Ar, EiJ_Ar, ei_Arp, EiJ_Arp
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, SAVE :: Si, Vi
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: Cij, Kij
    LOGICAL, SAVE :: flag_EI, flag_HI, flag_EExc, flag_HExc
    CHARACTER*80, SAVE :: solver, model

    ! Subroutine for the initialization of the Argon CR library 
    CONTAINS

      !------------------------------------------------------!
      ! This subroutine initializes the Argon CR library
      ! The species order is the following 
      ! 1 em
      ! 2 Ar(i)
      ! 3 Arp(i)
      SUBROUTINE initialize (in_solver, in_mixture, in_transf, in_reaction, in_path, ns, ntrot, ntvib, nte, neq, ndim, mass)
 
         INTEGER, PARAMETER :: in1 = 10
         INTEGER :: i, j, s, ios, dum_int
         INTEGER :: length, length_i, length_j, len_path, len_file
         REAL(KIND=8) :: dum_double, T
         CHARACTER*10 :: i_char, j_char
         CHARACTER*20 :: filename
         CHARACTER*80 :: dum_char
         CHARACTER*200 :: path

         INTEGER, INTENT(IN) :: ns, ntrot, ntvib, nte, ndim, neq
         CHARACTER*(*), INTENT(IN) :: in_solver, in_mixture, in_transf, in_reaction, in_path
         REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: mass
 
         ! Common data
         nb_ns   = ns
         nb_trot = ntrot
         nb_tvib = ntvib
         nb_te   = nte
         nb_dim  = ndim
         nb_eq   = neq 
         nb_temp = 1 + ntvib + ntrot

         ! Physical constants 
         ukb   = 1.380658d-23   
         una   = 6.0221367d23   
         urg   = ukb*una
         ue    = 1.602191d-19
         uh    = 6.626075d-34
         upi   = 3.14159265d0 
         ueps0 = 8.854188d-12

         ! Solver name
         solver = in_solver
         length = LEN_TRIM(solver)
          
         ! Molecular masses of em, Ar, and Arp species
         mm_em  = 0.00055d-3
         mm_Ar  = 39.948d-3
         mm_Arp = 39.9475d-3 

         ! Common factor for translational partition faction
         fac_Q = (2.d0*upi*ukb/uh**2.d0/una)**1.5

         ! Specific gas constants of em, Ar and Arp species
         R_em  = urg/mm_em
         R_Ar  = urg/mm_Ar
         R_Arp = urg/mm_Arp

         ! Constant volume and pressure specific heats
         cv_tr_em  = 1.5d0*R_em
         cv_tr_Ar  = 1.5d0*R_Ar
         cv_tr_Arp = 1.5d0*R_Arp

         cp_tr_em  = cv_tr_em + R_em
         cp_tr_Ar  = cv_tr_Ar + R_Ar
         cp_tr_Arp = cv_tr_Arp + R_Arp

         ! Physical model in use
         model = in_mixture

         ! Initialization for collisional processes flag
         flag_HI = .FALSE.
         flag_EI = .FALSE.
         flag_EExc = .FALSE.
         flag_HExc = .FALSE.

         SELECT CASE(model)

           !-----------------------------------------------!
           ! Use of electronic specific Collisional-Radiative (CR) model
           CASE('CR','cr')
             WRITE(*,20)solver(1:length),':: Argon CR library -> CR model is being used'
             PRINT*

             ! Number of electronic levels 
             levels_Arp = 2

             ! Number of electronic levels Ar and Arp atoms 
             levels_Ar = nb_ns - 1 - levels_Arp

             lev_Ar  = levels_Ar
             lev_Arp = levels_Arp 

             ! Model name
             model = 'CR'

             ! Path to data directory
             path     = in_path
             len_path = LEN_TRIM(path)  

             ! Check on the number of species 
             CALL check_nb_species (nb_ns, (1 + levels_Ar + levels_Arp)) 

             ! Posistions of em, Ar and Arp in the species list 
             pos_em  = 1
             pos_Ar  = pos_em + 1
             pos_Arp = pos_Ar + levels_Ar

             ! Position of free electron temperature in the temperature vector
             pos_Te = 2

             mass(pos_em) = mm_em
             DO i = 1,levels_Ar
                mass(pos_Ar + i - 1) = mm_Ar
             ENDDO 

             DO i = 1,levels_Arp
                mass(pos_Arp + i - 1) = mm_Arp
             ENDDO 
            
             ! Read the transfer file (processes to be accounted for)
             OPEN(UNIT=in1,FILE=path(1:len_path)//'kin_scheme',STATUS='unknown',IOSTAT=ios)

             ! Check file status 
             CALL check_file (ios, 'kin_scheme')

             ! Heavy-particle impact ionization
             READ(in1,*)dum_char
             READ(in1,*)flag_HI

             ! Electron impact ionization
             READ(in1,*)dum_char
             READ(in1,*)flag_EI

             ! Heavy-particle impact excitation
             READ(in1,*)dum_char
             READ(in1,*)flag_HExc

             ! Electron impact excitation
             READ(in1,*)dum_char
             READ(in1,*)flag_EExc

             CLOSE(in1)

             ! Memory allocation
             ALLOCATE(mm_species(nb_ns), Ri(nb_ns))
             ALLOCATE(gi_Ar(levels_Ar),gi_Arp(levels_Arp))
             ALLOCATE(ei_Ar(levels_Ar),ei_Arp(levels_Arp),EiJ_Ar(levels_Ar),EiJ_Arp(levels_Arp))
             ALLOCATE(jc(levels_Ar),Arp_lev(levels_Ar))
             ALLOCATE(Im(levels_Ar))

             ! Vector storing the species molecular masses
             mm_species = mass

             ! Vector storing the species gas constants
             DO i = 1,nb_ns 
                Ri(i) = urg/mm_species(i)
             ENDDO

             ! Open and read the energy level file for Ar
             OPEN(UNIT=in1,FILE=path(1:len_path)//'Ei_Ar',STATUS='unknown',IOSTAT=ios)     

             ! Check file status 
             CALL check_file (ios, 'Ei_Ar')

             READ(in1,*)dum_int 
             DO i = 1,levels_Ar
                READ(in1,*)gi_Ar(i),EiJ_Ar(i)
             ENDDO
             
             CLOSE(in1)

             ! Conversions [eV] -> [J], [J] -> [J/kg]
             EiJ_Ar = EiJ_Ar*ue
             ei_Ar  = EiJ_Ar*una/mm_Ar

             ! Open and read the energy level file for Arp
             OPEN(UNIT=in1,FILE=path(1:len_path)//'Ei_Arp',STATUS='unknown',IOSTAT=ios)     

             ! Check file status 
             CALL check_file (ios, 'Ei_Arp')

             READ(in1,*)dum_int 
             DO i = 1,levels_Arp
                READ(in1,*)gi_Arp(i),EiJ_Arp(i)
             ENDDO
             
             CLOSE(in1)

             ! Conversions [eV] -> [J]
             EiJ_Arp = EiJ_Arp*ue

             ! The levels of Arp alreadty account for the formation enthalpy
             ! (hence the ionization and energy of Arp are set to zero)
             ! Ionization energy and formation enthalpy of Arp ([J] and [J/Kg], respectively)
             Eion   = 0.d0
             hf_Arp = 0.d0
               
             ! Conversions [J] -> [J/kg]
             ei_Arp = EiJ_Arp*una/mm_Arp
            
             ! Reading core angular momentum jc for Ar electronic levels
             OPEN(UNIT=in1,FILE=path(1:len_path)//'jc',STATUS='unknown',IOSTAT=ios)     
             
             ! Check file status 
             CALL check_file (ios, 'jc')

             DO i = 1,levels_Ar
                READ(in1,*)jc(i)
             ENDDO
             
             CLOSE(in1)
             
             ! Arp level involved for heavy particle impact ionization (Ar(i)
             ! and Arp must have the same value of the core angular momentum)
             DO i = 1,levels_Ar 
                IF(jc(i).EQ.jc_Arp1) THEN
                  Arp_lev(i) = 1
                ELSE 
                  Arp_lev(i) = 2
                ENDIF
             ENDDO
              
             ! Ionization potential of Ar levels (computed)
             DO i = 1,levels_Ar
                Im(i) = EiJ_Arp(Arp_lev(i)) - EiJ_Ar(i)
             ENDDO
            
             ! Reading temperature file
             OPEN(UNIT=in1,FILE=path(1:len_path)//'T',STATUS='unknown',IOSTAT=ios)

             ! Check file status 
             CALL check_file (ios, 'T')

             READ(in1,*)nb_points
             READ(in1,*)Tmin
             READ(in1,*)Tmax
             READ(in1,*)deltaT

             ALLOCATE(tvec(nb_points))

             DO i = 1,nb_points 
                READ(in1,*)tvec(i)
             ENDDO

             CLOSE(in1)
             
             ! Reading rate coefficients files for collisional and radiative processes 
             ! Rate coefficients for electron impact excitation
             ! Ar(i) + em <=> Ar(j) + em ; j > i
             IF (flag_EExc.EQV..TRUE.) THEN

                ! Allocation
                ALLOCATE(Cij(levels_Ar,levels_Ar,nb_points))

                ! Initialization
                Cij = 0.d0

                DO i = 1,levels_Ar - 1

                   DO j = i + 1,levels_Ar 

                      ! Filename creaction
                      WRITE(i_char,30)i - 1
                      WRITE(j_char,30)j - 1

                      i_char   = ADJUSTL(i_char)
                      length_i = LEN_TRIM(i_char)
                      j_char   = ADJUSTL(j_char)
                      length_j = LEN_TRIM(j_char)

                      filename = 'C'//i_char(1:length_i)//'_'//j_char(1:length_j)
                      len_file = LEN_TRIM(filename)

                      OPEN(UNIT=in1,FILE=path(1:len_path)//filename(1:len_file),STATUS='unknown',IOSTAT=ios)

                      ! Check file status 
                      CALL check_file (ios, filename(1:len_file))

                      ! Read data
                      DO s = 1,nb_points 
                         READ(in1,*)dum_double
                         IF (dum_double.GT.1.d-125) Cij(i,j,s) = dum_double
                      ENDDO

                      CLOSE(in1)

                   ENDDO

                ENDDO 

                ! Conversion [m^3/s] -> [m^3/mol/s]
                Cij = una*Cij

                WRITE(*,20)solver(1:length),':: Argon CR library -> rate coefficients for electron impact excitation read'
                PRINT*

             ENDIF

             ! Rate coefficients for heavy-particle impact excitation
             ! Ar(i) + Ar(1) <=> Ar(j) + Ar(1) ; j > i
             IF (flag_HExc.EQV..TRUE.) THEN

                ! Allocation
                ALLOCATE(Kij(levels_Ar,levels_Ar,nb_points))   

                ! Initialization
                Kij = 0.d0

                DO i = 1,levels_Ar - 1

                   DO j = i + 1,levels_Ar 

                      ! Filename creaction
                      WRITE(i_char,30)i - 1
                      WRITE(j_char,30)j - 1

                     i_char   = ADJUSTL(i_char)
                     length_i = LEN_TRIM(i_char)
                     j_char   = ADJUSTL(j_char)
                     length_j = LEN_TRIM(j_char)

                     filename = 'K'//i_char(1:length_i)//'_'//j_char(1:length_j)
                     len_file = LEN_TRIM(filename)
 
                     OPEN(UNIT=in1,FILE=path(1:len_path)//filename(1:len_file),STATUS='unknown',IOSTAT=ios)

                     ! Check file status 
                     CALL check_file (ios, filename(1:len_file))

                     ! Read data
                     DO s = 1,nb_points 
                        READ(in1,*)Kij(i,j,s)
                     ENDDO

                     CLOSE(in1)

                   ENDDO

                ENDDO 

                ! Conversion [m^3/s] -> [m^3/mol/s]
                Kij = una*Kij

                ! Correction factor applied to the rate coefficients for the heavy-particle impact excitation from the ground state
                ! (Drawin's model overpreditcs the cross-section values - for more info see the PhD thesis of M. G. Kapper 2009)
                i = 1
                DO j = 1,levels_Ar
                   DO s = 1,nb_points 
                      Kij(i,j,s) = 0.085d0*Kij(i,j,s)
                   ENDDO
                ENDDO

                WRITE(*,20)solver(1:length),':: Argon CR library -> rate coefficients for heavy-particle impact excitation read'
                PRINT*
                
             ENDIF

             ! Rate coefficients for electron impact ionization 
             ! Ar(i) + em <=> Arp + em + em
             IF (flag_EI.EQV..TRUE.) THEN

                ! Allocation
                ALLOCATE(Si(levels_Ar,nb_points)) 

                ! Initialization
                Si = 0.d0

                DO i = 1,levels_Ar

                   ! Filename creaction
                   WRITE(i_char,30)i - 1 

                   i_char   = ADJUSTL(i_char)
                   length_i = LEN_TRIM(i_char) 

                   filename = 'S'//i_char(1:length_i)
                   len_file = LEN_TRIM(filename)

                   OPEN(UNIT=in1,FILE=path(1:len_path)//filename(1:len_file),STATUS='unknown',IOSTAT=ios)

                   ! Check file status 
                   CALL check_file (ios, filename(1:len_file))

                   ! Read data
                   DO s = 1,nb_points 
                      READ(in1,*)Si(i,s)
                   ENDDO

                   CLOSE(in1) 

                ENDDO             

                ! Conversion [m^3/s] -> [m^3/mol/s]
                Si = una*Si

                WRITE(*,20)solver(1:length),':: Argon CR library -> rate coefficients for electron impact ionization read'
                PRINT*

             ENDIF

             ! Rate coefficients for heavy-particle impact ionization
             ! Ar(i) + Ar(1) <=> Arp + Ar(1) + em  ; Arp and Ar(i) must have 
             ! the same value for the core angular momentum jc
             IF (flag_HI.EQV..TRUE.) THEN

                ! Allocation
                ALLOCATE(Vi(levels_Ar,nb_points)) 

                ! Initialization
                Vi = 0.d0

                DO i = 1,levels_Ar

                   ! Filename creaction
                   WRITE(i_char,30)i - 1 

                   i_char   = ADJUSTL(i_char)
                   length_i = LEN_TRIM(i_char) 

                   filename = 'V'//i_char(1:length_i)
                   len_file = LEN_TRIM(filename)

                   OPEN(UNIT=in1,FILE=path(1:len_path)//filename(1:len_file),STATUS='unknown',IOSTAT=ios)

                   ! Check file status 
                   CALL check_file (ios, filename(1:len_file))

                   ! Read data
                   DO s = 1,nb_points 
                      READ(in1,*)Vi(i,s)
                   ENDDO

                   CLOSE(in1) 

                ENDDO

                ! Conversion [m^3/s] -> [m^3/mol/s]
                Vi = una*Vi 

                WRITE(*,20)solver(1:length),':: Argon CR library -> rate coefficients for heavy-particle impact ionization read'
                PRINT*

             ENDIF

             ! Electron-neutral collision rate [m^3/s]
             ALLOCATE(ken(nb_points),kin(nb_points))
             OPEN(UNIT=in1,FILE=path(1:len_path)//'ken',STATUS='unknown',IOSTAT=ios)

             ! Check file status 
             CALL check_file (ios, 'ken')

             ! Read data
             DO i = 1,nb_points 
                READ(in1,*)ken(i)
             ENDDO

             CLOSE(in1) 

             ! Ion-neutral collision rate [m^3/s]
             OPEN(UNIT=in1,FILE=path(1:len_path)//'kin',STATUS='unknown',IOSTAT=ios)

             ! Check file status 
             CALL check_file (ios, 'kin')

             ! Read data
             DO i = 1,nb_points 
                READ(in1,*)kin(i)
             ENDDO

             CLOSE(in1) 

           !-----------------------------------------------!
           ! Use of multi-temperature model (not obtained from the original CR model)
           CASE('MT')
             WRITE(*,20)solver(1:length),':: Argon CR library -> MT model is being used'
             PRINT*

             ! Number of electronic levels Ar and Arp atoms (only for partition function evaluation)
             levels_Ar  = 31
             levels_Arp = 2 

             ! Number of electronic levels used during the computation
             lev_Ar  = 1
             lev_Arp = 1  

             ! Posistions of em, Ar and Arp in the species list 
             pos_em  = 1
             pos_Ar  = 2
             pos_Arp = 3

             ! Position of free electron temperature in the temperature vector
             pos_Te = 2

             mass(pos_em)  = mm_em
             mass(pos_Ar)  = mm_Ar
             mass(pos_Arp) = mm_Arp

             ! Model name
             model = 'MT'

             ! Path to data directory
             path     = in_path
             len_path = LEN_TRIM(path)  

             ! Position of free electron temperature in the temperature vector
             pos_Te = 2

             ! Read the transfer file (processes to be accounted for)
             OPEN(UNIT=in1,FILE=path(1:len_path)//'kin_scheme',STATUS='unknown',IOSTAT=ios)

             ! Check file status 
             CALL check_file (ios, 'kin_scheme')

             ! Heavy-particle impact ionization
             READ(in1,*)dum_char
             READ(in1,*)flag_HI

             ! Electron impact ionization
             READ(in1,*)dum_char
             READ(in1,*)flag_EI

             CLOSE(in1)

             ! Read fit coefficients for reaction rate coefficients (data given in [cm^3/mol/s])
             ! Heavy-particle impact ionization
             ! Ar + Ar <=> Arp + em
             OPEN(UNIT=in1,FILE=path(1:len_path)//'kf_HI',STATUS='unknown',IOSTAT=ios)

             ! Check file status 
             CALL check_file (ios, 'kf_HI')

             READ(in1,*)ahi_kf, bhi_kf,chi_kf

             CLOSE(in1)

             ! Electron impact ionization
             ! Ar + em <=> Arp + em + em
             OPEN(UNIT=in1,FILE=path(1:len_path)//'kf_EI',STATUS='unknown',IOSTAT=ios)

             ! Check file status 
             CALL check_file (ios, 'kf_EI')

             READ(in1,*)aei_kf, bei_kf,cei_kf

             CLOSE(in1)

             ! Unit conversion [cm^3/mol/s] -> [m^3/mol/s]
             ahi_kf = 1.d-6*ahi_kf 
             aei_kf = 1.d-6*aei_kf 

             ! Memory allocation
             ALLOCATE(mm_species(nb_ns))
             ALLOCATE(gi_Ar(levels_Ar),gi_Arp(levels_Arp))
             ALLOCATE(EiJ_Ar(levels_Ar),EiJ_Arp(levels_Arp))

             ! Vector storing the species molecular masses 
             mm_species = mass 
 
             ! Reading temperature file
             OPEN(UNIT=in1,FILE=path(1:len_path)//'T',STATUS='unknown',IOSTAT=ios)

             ! Check file status 
             CALL check_file (ios, 'T')

             READ(in1,*)nb_points
             READ(in1,*)Tmin
             READ(in1,*)Tmax
             READ(in1,*)deltaT

             ALLOCATE(tvec(nb_points),ken(nb_points))

             DO i = 1,nb_points 
                READ(in1,*)tvec(i)
             ENDDO

             CLOSE(in1)

             ! Electron-neutral collision rate [m^3/s]
             OPEN(UNIT=in1,FILE=path(1:len_path)//'ken',STATUS='unknown',IOSTAT=ios)

             ! Check file status 
             CALL check_file (ios, 'ken')

             ! Read data
             DO i = 1,nb_points 
                READ(in1,*)ken(i)
             ENDDO

             CLOSE(in1) 

             ! Open and read the energy level file for Ar
             OPEN(UNIT=in1,FILE=path(1:len_path)//'Ei_Ar',STATUS='unknown',IOSTAT=ios)     

             ! Check file status 
             CALL check_file (ios, 'Ei_Ar')

             READ(in1,*)dum_int 
             DO i = 1,levels_Ar
                READ(in1,*)gi_Ar(i),EiJ_Ar(i)
             ENDDO
             
             CLOSE(in1)

             ! Conversions [eV] -> [J], [J] -> [J/kg]
             EiJ_Ar = EiJ_Ar*ue

             OPEN(UNIT=in1,FILE=path(1:len_path)//'Ei_Arp',STATUS='unknown',IOSTAT=ios)     

             ! Check file status 
             CALL check_file (ios, 'Ei_Arp')

             READ(in1,*)dum_int 
             DO i = 1,levels_Arp
                READ(in1,*)gi_Arp(i),EiJ_Arp(i)
             ENDDO
          
             CLOSE(in1)

             ! Conversions [eV] -> [J], [J] -> [J/kg]
             EiJ_Arp = EiJ_Arp*ue

             ! Ionization energy and formation enthalpy of Arp ([J] and [J/Kg], respectively)
             Eion   = EiJ_Arp(1)
             hf_Arp = Eion*una/mm_Arp 

             ! Energy levels measured from the ground state of the Arp atom
             EiJ_Arp = EiJ_Arp - Eion

             !! TO ERASE !!
             gi_Ar  = 0.d0
             gi_Arp = 0.d0
             gi_Ar(1) = 1.d0
             gi_Arp(1)= 1.d0
             EiJ_Ar  = 0.d0
             EiJ_Arp = 0.d0

           CASE DEFAULT
             WRITE(*,10)solver(1:length)//' in "mod_argon_CR_initialize_CFD"...'
             PRINT*
             WRITE(*,10)'not implemented yet'
             PRINT*
             STOP     

         END SELECT

10  FORMAT(A)
20  FORMAT(A,A)
30  FORMAT(I3)

    END SUBROUTINE initialize

    !----------------------------------------------------!
    ! This subroutine shuts down the Argon CR thermodynamic library
    SUBROUTINE finalize ()

       ! Deallocate arrays       
       IF (ALLOCATED(mm_species))         DEALLOCATE(mm_species)
       IF (ALLOCATED(Ri))                 DEALLOCATE(Ri)
       IF (ALLOCATED(gi_Ar))              DEALLOCATE(gi_Ar)
       IF (ALLOCATED(gi_Arp))             DEALLOCATE(gi_Arp)
       IF (ALLOCATED(EiJ_Ar))             DEALLOCATE(EiJ_Ar)
       IF (ALLOCATED(EiJ_Arp))            DEALLOCATE(EiJ_Arp)
       IF (ALLOCATED(ei_Ar))              DEALLOCATE(ei_Ar)
       IF (ALLOCATED(ei_Arp))             DEALLOCATE(ei_Arp)
       IF (ALLOCATED(Im))                 DEALLOCATE(Im)
       IF (ALLOCATED(jc))                 DEALLOCATE(jc)
       IF (ALLOCATED(Si))                 DEALLOCATE(Si)
       IF (ALLOCATED(Vi))                 DEALLOCATE(Vi)
       IF (ALLOCATED(Cij))                DEALLOCATE(Cij)
       IF (ALLOCATED(Kij))                DEALLOCATE(Kij)
       IF (ALLOCATED(ken))                DEALLOCATE(ken)
       IF (ALLOCATED(kin))                DEALLOCATE(kin)

    END SUBROUTINE finalize

    !----------------------------------------------------!
    ! This subroutine check if the number of species is correct
    SUBROUTINE check_nb_species (ns, nspecies) 

      INTEGER, INTENT(IN) :: ns, nspecies

      IF (nspecies.NE.ns) THEN 
         WRITE(*,10)solver(1:LEN_TRIM(solver))//' :: number of species set to uncorrect value'
         PRINT*
         WRITE(*,10)solver(1:LEN_TRIM(solver))//' :: in mod_argon_CR_initialize_CFD ..'
         PRINT*
         STOP
      ENDIF 
        
10  FORMAT(A)

    END SUBROUTINE check_nb_species

    !----------------------------------------------------!
    ! This subroutine checks if a given file has been opened correctly 
    SUBROUTINE check_file (ios, filename)

      INTEGER, INTENT(IN) :: ios
      CHARACTER*(*), INTENT(IN) :: filename

      IF (ios.NE.0) THEN 
         WRITE(*,10)solver(1:LEN_TRIM(solver))//':: '//'"'//filename(1:LEN_TRIM(filename))//'"'//& 
                  & ' file not found, in mod_argon_CR_initialize_CFD ..'
         PRINT*
         STOP
      ENDIF

10  FORMAT(A)

    END SUBROUTINE check_file

  END MODULE mod_argon_CR_initialize_CFD
!------------------------------------------------------------------------------!
