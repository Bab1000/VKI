!------------------------------------------------------------------------------!
! This module provides functions and subroutines dealing with thermodynamic, transport properties 
! and source terms for the N-N2 system when using the DSMC multi-temperature model.
! The initialization provided in this module should be use only when the library is interfaced with CFD codes
  MODULE mod_nitrogen_DSMC_initialize_CFD

    IMPLICIT NONE

    INTEGER, SAVE :: pos_N, pos_N2, posTr, posTv
    INTEGER, SAVE :: nb_ns, nb_tvib, nb_trot, nb_te, nb_dim, nb_eq, nb_temp
    REAL(KIND=8), SAVE :: ukb, una, urg, ue, uh, upi
    REAL(KIND=8), SAVE :: par1, par2
    REAL(KIND=8), SAVE :: gn, mm_N, mm_N2, Rn, Rn2, mu, ed_n2, hf_N, theta_rot, theta_vib, cv_rot, fac_Q
    REAL(KIND=8), SAVE :: T_ref, exp_visc, eta_ref, Sc_ref
    REAL(KIND=8), SAVE :: ov_alpha, Zr_inf, Zrv, T_star
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: cv_tr, Ri, mm 
    CHARACTER*80, SAVE :: solver

    ! Subroutine for the initialization of the DSMC N-N2 library 
    CONTAINS

      !------------------------------------------------------!
      ! This subroutine initializes the DSMC N-N2 library
      SUBROUTINE initialize (in_solver, in_mixture, in_transf, in_reaction, in_path, ns, ntrot, ntvib, nte, neq, ndim, mass)
 
         INTEGER :: i, length
         REAL(KIND=8) :: tmp

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

         WRITE(*,20)solver(1:length),':: Nitrogen DSMC library -> initialization'
         PRINT*

         ! Number of species
         IF (nb_ns.EQ.1) THEN

            WRITE(*,20)solver(1:length),':: Nitrogen DSMC library -> pure N2 mixture'
            PRINT*

            pos_N  = 1
            pos_N2 = 1

         ELSEIF (nb_ns.EQ.2) THEN

            WRITE(*,20)solver(1:length),':: Nitrogen DSMC library -> N-N2 mixture'
            PRINT* 

            pos_N  = 1
            pos_N2 = 2

         ELSE

            WRITE(*,20)solver(1:length),':: Nitrogen DSMC library -> error in number of species'
            PRINT* 
            STOP

         ENDIF

         ! Physical constants
         ukb = 1.380658d-23   
         una = 6.0221367d23   
         urg = ukb*una
         ue  = 1.602191d-19
         uh  = 6.626075d-34
         upi = 3.14159265d0 
   
         ! Parameters for Parker's formula for Zrc
         par1 = 0.25d0*upi**2 + upi
         par2 = 0.5d0*upi**1.5

         ! Various parameters 
         T_ref    = 273.d0
         eta_ref  = 1.656d-5
         Sc_ref   = 1.34d0
         exp_visc = 0.74d0  
         ov_alpha = 4.d0/8.333d0
         Zr_inf   = 23.5d0
         Zrv      = 120.d0
         T_star   = 91.5d0

         ! Molecular masses of N and N2 species
         mm_N  = 14.0067d-3
         mm_N2 = 28.0134d-3 
  
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
         hf_N = (113272.d0*ukb/2.d0)*una/mm_N 

         IF (nb_ns.EQ.1) hf_N = 0.d0

         ! Molecular masses
         mass(pos_N)  = mm_N
         mass(pos_N2) = mm_N2

         ! Common factor for translational partition function 
         fac_Q = (2.d0*upi*ukb/uh**2.d0/una)**1.5

         ! Translational and rotational specific heats (constant volume)
         ALLOCATE(cv_tr(nb_ns), Ri(nb_ns), mm(nb_ns))

         mm = mass
         DO i = 1,nb_ns 
            tmp       = 1.d0/mm(i)
            cv_tr(i)  = 1.5d0*tmp 
            Ri(i)     = tmp 
         ENDDO
         cv_rot = tmp

         cv_tr  = cv_tr*urg
         cv_rot = cv_rot*urg
         Ri     = Ri*urg

20  FORMAT(A,A)

    END SUBROUTINE initialize

    !------------------------------------------------------!
    ! This subroutine de-allocates vectors and nullifies eventually used pointers to function subroutines. 
    SUBROUTINE finalize ()

      USE mod_function_pointer_DSMC

      ! Memory de-allocation
      IF (ALLOCATED(mm))                              DEALLOCATE(mm)
      IF (ALLOCATED(Ri))                              DEALLOCATE(Ri)
      IF (ALLOCATED(cv_tr))                           DEALLOCATE(cv_tr)

      ! Nullify pointers
      IF (ASSOCIATED(get_temp_DSMC))                  NULLIFY(get_temp_DSMC)
      IF (ASSOCIATED(get_species_cv_DSMC))            NULLIFY(get_species_cv_DSMC)
      IF (ASSOCIATED(get_species_energy_cv_DSMC))     NULLIFY(get_species_energy_cv_DSMC)

    END SUBROUTINE finalize 

  END MODULE mod_nitrogen_DSMC_initialize_CFD
!------------------------------------------------------------------------------!
