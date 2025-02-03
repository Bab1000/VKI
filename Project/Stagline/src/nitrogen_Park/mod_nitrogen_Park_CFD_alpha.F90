!------------------------------------------------------------------------------!
! This module provides functions and subroutines dealing with thermodynamic, transport properties 
! and source terms for the N-N2 system when using the Park multi-temperature model.
! The initialization provided in this module should be use only when the library is interfaced with CFD codes
  MODULE mod_nitrogen_Park_initialize_CFD

    IMPLICIT NONE

    INTEGER, SAVE :: posTr, posTv
    INTEGER, SAVE :: pos_N, pos_N2
    INTEGER, SAVE :: table_points
    INTEGER, SAVE :: nb_ns, nb_tvib, nb_trot, nb_te, nb_dim, nb_eq, nb_temp, nb_inter
    REAL(KIND=8), PARAMETER :: Tmin_kin = 250.d0
    REAL(KIND=8), PARAMETER :: Tmin = 100.d0
    REAL(KIND=8), PARAMETER :: Tmax = 90000.d0
    REAL(KIND=8), PARAMETER :: dT = 5.d0
    REAL(KIND=8), SAVE :: ukb, una, urg, ue, uh, upi
    REAL(KIND=8), SAVE :: gn, mm_N, mm_N2, Rn, Rn2, mu, ed_n2, hf_n, theta_rot, theta_vib, fac_Q, fac_keq, fac_exp, & 
                        & m_ratio, ov_m_ratio
    REAL(KIND=8), SAVE :: c1, c2, c3, c4
    REAL(KIND=8), SAVE :: cv_rot
    REAL(KIND=8), SAVE :: fac_mu, fac_Diff
    REAL(KIND=8), DIMENSION(2), SAVE :: mw_a, mw_b
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: sqrt_mi, sqrt_mij
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: T_table
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: cv_tr 
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: q11, q12, q22, bs, cs, Ri, mi
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, SAVE :: Omega_11_table, Omega_12_table, Omega_22_table
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, SAVE :: Astar_table, Bstar_table, Cstar_table 
    CHARACTER*80, SAVE :: solver
    LOGICAL, SAVE :: diss_imp_N2 = .FALSE., vt_imp_mol = .FALSE.

    ! Subroutine for the initialization of the Park N-N2 library 
    CONTAINS

      !------------------------------------------------------!
      ! This subroutine initializes the Park N-N2 library
      SUBROUTINE initialize (in_solver, in_mixture, in_transf, in_reaction, in_path, ns, ntrot, ntvib, nte, neq, ndim, mass)
 
         INTEGER :: i, j, t, ij
         INTEGER :: pos1, pos2, pos3, pos4
         INTEGER :: comp, length
         REAL(KIND=8) :: lnt, lnt2, lnt3
         REAL(KIND=8) :: Omega_11, Omega_12, Omega_22, Astar, Bstar, Cstar

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

         posTr = 1 + nb_trot
         posTv = 1 + nb_trot + nb_tvib

         ! Solver name
         solver = in_solver
         length = LEN_TRIM(solver)

         WRITE(*,20)solver(1:length)//':: Nitrogen Park library -> initialization'
         PRINT*

         ! Default values of species location
         pos_N  = 1
         pos_N2 = 2

         ! Flags for chemical reactions and energy transfer terms
         SELECT CASE (in_mixture)

           CASE('complete')
             diss_imp_N2 = .TRUE.
             WRITE(*,20)solver(1:length),':: Nitrogen Park -> N2 + N2 = 2N + N2 reaction also considered'
             PRINT*

           CASE('reduced','none')
             WRITE(*,20)solver(1:length),':: Nitrogen Park _> only N2 + N = 3N reaction considered'
             PRINT*

           CASE('nitrogen')
             pos_N2 = 1
             WRITE(*,20)solver(1:length),':: Nitrogen Park _> only N2 species considered'
             PRINT*

           CASE('eq')
             WRITE(*,20)solver(1:length),':: Nitrogen Park _> equilibrium conditions - no finite rate chemistry'
             PRINT*

           CASE DEFAULT 
             WRITE(*,20)solver(1:length),':: Nitrogen Park -> in mod_initialize_CFD, error in mixture selection'
             PRINT* 
             STOP

         END SELECT 

         ! Flags for energy transfer mechanism
         SELECT CASE (in_transf)

           CASE('complete')
             vt_imp_mol = .TRUE.
             WRITE(*,20)solver(1:length),':: Nitrogen Park -> N2-N2 inelastic collisions also considered'
             PRINT*

           CASE('reduced','none')
             WRITE(*,20)solver(1:length),':: Nitrogen Park -> only N-N2 inelastic collisions considered'
             PRINT*

           CASE('nitrogen')
             WRITE(*,20)solver(1:length),':: Nitrogen Park -> only N2-N2 inelastic collisions considered'
             PRINT*

           CASE('eq')
             WRITE(*,20)solver(1:length),':: Nitrogen Park _> equilibrium conditions - no finite rate energy trasnfer'
             PRINT*

           CASE DEFAULT 
             WRITE(*,20)solver(1:length),':: Nitrogen Park -> in mod_initialize_CFD, error in transfer mechanism selection'
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

         ! Data of N2 molecule
         ! Characteristic vibrational and rotational temperatures  
         theta_rot = 2.88d0
         theta_vib = 3392.7d0

         ! Formation enthalpy (per unit mass) and dissociation energy (in J/mol)     
         ed_n2 = 113272.d0*ukb
         hf_n = (113272.d0*ukb/2.d0)*una/mm_N 
         
         IF (in_mixture.EQ.'nitrogen') hf_n = 0.d0

         ! Dissociation reaction N2 + N = N + N + N
         c1 = 0.049816d0
         c2 = -1.60d0
         c3 = 113272.d0 

         ! Dissociation reaction N2 + N2 = N2 + N + N (only the pre-exponential factor needs to be set)
         c4 = 0.0116237d0

         ! VT interactions (parameter for Millikan and White formula)
         ! N2-N  interaction
         mw_a(pos_N) = 180.D0
         mw_b(pos_N) = 0.0262D0

         ! N2-N2 interaction
         mw_a(pos_N2) = 221.D0 
         mw_b(pos_N2) = 0.0290D0

         ! Common factors
         fac_Q   = (2.d0*upi*ukb/(uh**2.d0*una))**1.5
         fac_exp = 1.d-6*una
         fac_keq = 1.d0/una*gn**2*(mm_N/mm_N2)**1.5

         ! Unit conversion 
         c1 = c1*fac_exp
         c4 = c4*fac_exp

         ! Translational and rotational specific heats (constant volume)
         mass(pos_N)  = mm_N
         mass(pos_N2) = mm_N2

         ALLOCATE(cv_tr(nb_ns))

         DO i = 1,nb_ns 
            cv_tr(i) = 1.5d0/mass(i) 
         ENDDO
         cv_rot = 1.d0/mass(pos_N2)
   
         cv_tr  = cv_tr*urg
         cv_rot = cv_rot*urg

         ! Curve fits for collision integrals
         nb_inter = (nb_ns + 1)*nb_ns/2
         comp     = 4*nb_inter
         ALLOCATE(q11(comp), q12(comp), q22(comp), bs(comp), cs(comp), Ri(nb_ns), mi(nb_ns))  

         Ri(pos_N)  = Rn
         Ri(pos_N2) = Rn2

         mi(pos_N)  = mm_N
         mi(pos_N2) = mm_N2

         IF (in_mixture.NE.'nitrogen') THEN
 
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

           ! N-N2 interaction
           q11(5)  = -1.51789d-3
           q11(6) = 2.30672d-2
           q11(7) = -4.22508d-1 
           q11(8) = 5.46824d0

           q22(5)  = -3.51846d-3 
           q22(6) = 7.14265d-2 
           q22(7) = -7.70251d-1
           q22(8) = 6.34262d0

           bs(5)  = -1.61616d-3
           bs(6) = 4.59611d-2 
           bs(7) = -3.90469d-1 
           bs(8) = 1.13622d0

           ! N2-N2 interaction
           q11(9)  = -5.04358d-3
           q11(10) = 1.08079d-1
           q11(11) = -9.98634d-1 
           q11(12) = 6.80089d0

           q22(9)  = -6.99424d-3 
           q22(10) = 1.54912d-1 
           q22(11) = -1.33680d0  
           q22(12) = 7.66252d0

           bs(9)  = -1.83418d-3
           bs(10) = 4.44907d-2
           bs(11) = -3.35033d-1
           bs(12) = 8.97329d-1

           ! No thermal diffusion
           cs = 0.d0

         ELSE 

           ! N2-N2 interaction
           q11(1) = -5.04358d-3
           q11(2) = 1.08079d-1
           q11(3) = -9.98634d-1 
           q11(4) = 6.80089d0

           q22(1) = -6.99424d-3 
           q22(2) = 1.54912d-1 
           q22(3) = -1.33680d0  
           q22(4) = 7.66252d0

           bs(1) = -1.83418d-3
           bs(2) = 4.44907d-2
           bs(3) = -3.35033d-1
           bs(4) = 8.97329d-1

           ! No thermal diffusion
           cs = 0.d0

         ENDIF

         ! Generate look-up tables for collision integrals
         table_points = INT((Tmax - Tmin)/dT) + 1
         ALLOCATE(T_table(table_points))
         ALLOCATE(Omega_11_table(table_points,nb_inter))
         ALLOCATE(Omega_12_table(table_points,nb_inter))
         ALLOCATE(Omega_22_table(table_points,nb_inter))
         ALLOCATE(Astar_table(table_points,nb_inter))          
         ALLOCATE(Bstar_table(table_points,nb_inter))
         ALLOCATE(Cstar_table(table_points,nb_inter))
         ALLOCATE(sqrt_mi(nb_ns),sqrt_mij(nb_inter))

         ! Common factors for viscosity and binary diffusion coefficients
         fac_mu   = (5.d0/16.d0)*SQRT(upi*ukb/una)
         fac_Diff = (16.d0/3.d0)/SQRT(2.d0*upi*urg)

         ! Square roots of molecular masses and reduced molecular masses
         DO i = 1,nb_ns 

            sqrt_mi(i) = DSQRT(mi(i))

            DO j = i,nb_ns 

               ij = ((i - 1)*(2*nb_ns - i) + 2*j)/2 

               sqrt_mij(ij) = SQRT(mi(i)*mi(j)/(mi(i) + mi(j))) 

            ENDDO

         ENDDO

         ! Temperature table
         T_table(1) = Tmin
         DO i = 2,table_points
            T_table(i) = T_table(i - 1) + dT
         ENDDO

         WRITE(*,20)solver(1:length),':: Nitrogen Park -> generating tables for collision integrals'
         PRINT* 

         ! Collision integral look-up tables
         DO t = 1,table_points 

            lnt  = DLOG(T_table(t))
            lnt2 = lnt*lnt
            lnt3 = lnt2*lnt

            DO i = 1,nb_ns 

               DO j = i,nb_ns 

                  ! Index
                  ij = ((i - 1)*(2*nb_ns - i) + 2*j)/2 

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
                  Omega_11_table(t,ij) = Omega_11
                  Omega_12_table(t,ij) = Omega_12
                  Omega_22_table(t,ij) = Omega_22
                  Astar_table(t,ij)    = Astar
                  Bstar_table(t,ij)    = Bstar
                  Cstar_table(t,ij)    = Cstar

               ENDDO


            ENDDO

         ENDDO

20  FORMAT(A,A)

    END SUBROUTINE initialize

    !------------------------------------------------------!
    ! This subroutine de-allocates vectors and nullifies eventually used pointers to function subroutines. 
    SUBROUTINE finalize ()

      USE mod_function_pointer_Park

      ! Data de-allocation
      IF (ALLOCATED(q11))                             DEALLOCATE(q11)
      IF (ALLOCATED(q12))                             DEALLOCATE(q12)
      IF (ALLOCATED(q22))                             DEALLOCATE(q22)
      IF (ALLOCATED(bs))                              DEALLOCATE(bs)
      IF (ALLOCATED(cs))                              DEALLOCATE(cs)
      IF (ALLOCATED(Ri))                              DEALLOCATE(Ri) 
      IF (ALLOCATED(mi))                              DEALLOCATE(mi)
      IF (ALLOCATED(sqrt_mi))                         DEALLOCATE(sqrt_mi)
      IF (ALLOCATED(sqrt_mij))                        DEALLOCATE(sqrt_mij) 
      IF (ALLOCATED(T_table))                         DEALLOCATE(T_table)
      IF (ALLOCATED(Omega_11_table))                  DEALLOCATE(Omega_11_table)
      IF (ALLOCATED(Omega_12_table))                  DEALLOCATE(Omega_12_table)
      IF (ALLOCATED(Omega_22_table))                  DEALLOCATE(Omega_22_table)
      IF (ALLOCATED(Astar_table))                     DEALLOCATE(Astar_table)
      IF (ALLOCATED(Bstar_table))                     DEALLOCATE(Bstar_table)
      IF (ALLOCATED(Cstar_table))                     DEALLOCATE(Cstar_table)

      ! Nullify pointers
      IF (ASSOCIATED(get_temp_Park))                  NULLIFY(get_temp_Park)
      IF (ASSOCIATED(get_species_cv_Park))            NULLIFY(get_species_cv_Park)
      IF (ASSOCIATED(get_species_energy_cv_Park))     NULLIFY(get_species_energy_cv_Park)
      IF (ASSOCIATED(get_source_term))                NULLIFY(get_source_term)
      IF (ASSOCIATED(get_source_term_Jac))            NULLIFY(get_source_term_Jac)

    END SUBROUTINE finalize 

  END MODULE mod_nitrogen_Park_initialize_CFD
!------------------------------------------------------------------------------!
