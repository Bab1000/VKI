!------------------------------------------------------------------------------!
! This module provides functions and subroutines dealing with thermodynamic, transport properties 
! and source terms for a polyatomic gas according to the Whan-Chang and Ulhenbeck (WCU) theory.
! The initialization provided in this module should be used only when the library is interfaced with CFD codes
  MODULE mod_polyat_gas_wcu_initialize_CFD

    IMPLICIT NONE

    INTEGER, SAVE :: nb_ns, nb_tvib, nb_trot, nb_te, nb_dim, nb_eq, nb_temp
    INTEGER, SAVE :: table_points
    REAL(KIND=8), SAVE :: ukb, una, urg, ue, uh, upi
    REAL(KIND=8), SAVE :: cv_tr, mass, Rgas, d0, sigma_hs
    REAL(KIND=8), SAVE :: Tmin, Tmax, dT
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: T_table, cv_table, e_table, mu_table, kappa_table, lambda_table 
    CHARACTER*80, SAVE :: solver

    ! Subroutine for the initialization of the polyat_gas_wcu library  
    CONTAINS

      !------------------------------------------------------!
      ! This subroutine initializes the polyat_gas_wcu library
      SUBROUTINE initialize (in_solver, in_mixture, in_transf, in_reaction, in_path, ns, ntrot, ntvib, nte, neq, ndim, mm)
 
         INTEGER, PARAMETER :: in1 = 10, in2 = 11, in3 = 12
         INTEGER :: i, ios, len_path, length
         REAL(KIND=8) :: tmp1, tmp2, tmp3
         CHARACTER*200 :: path

         INTEGER, INTENT(IN) :: ns, ntrot, ntvib, nte, ndim, neq
         CHARACTER*(*), INTENT(IN) :: in_solver, in_mixture, in_transf, in_reaction, in_path
         REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: mm

         ! Common data
         nb_ns   = ns
         nb_trot = ntrot
         nb_tvib = ntvib
         nb_te   = nte
         nb_dim  = ndim
         nb_eq   = neq 
         nb_temp = 1 + ntvib + ntrot

         ! Solver name
         solver = in_solver
         length = LEN_TRIM(solver)

         ! Path to date files
         path     = in_path
         len_path = LEN_TRIM(path)

         ! Physical constants
         ukb = 1.380658d-23   
         una = 6.0221367d23   
         urg = ukb*una
         ue  = 1.602191d-19
         uh  = 6.626075d-34
         upi = 3.14159265d0 

         WRITE(*,20)solver(1:length),':: Polyat_gas_WCU library -> initialization'
         PRINT*

         ! Read gas data (mass and diameter)
         WRITE(*,20)solver(1:length),':: Polyat_gas_WCU library -> reading gas data'
         PRINT*

         OPEN(UNIT=in1,FILE=path(1:len_path)//'species_data',STATUS='unknown',IOSTAT=ios)

         CALL check_file (ios, 'species_data')

         ! Mass and diamter
         READ(in1,*)mass
         READ(in1,*)d0

         CLOSE(in1)
        
         sigma_hs = 0.25d0*d0**2
         Rgas  = ukb/mass
         mm    = mass*una
         cv_tr = 1.5d0*Rgas
        
         ! Read and load look-up tables for thermodynamics and transport 
         WRITE(*,20)solver(1:length),':: Polyat_gas_WCU library -> loading look-up tables'
         PRINT*

         ! Temperature table
         OPEN(UNIT=in1,FILE=path(1:len_path)//'T',STATUS='unknown',IOSTAT=ios)

         CALL check_file (ios, 'T')

         READ(in1,*)table_points
         READ(in1,*)Tmin
         READ(in1,*)Tmax
         READ(in1,*)dT

         ! Allocate memory for look-up tables
         ALLOCATE(T_table(table_points))
         ALLOCATE(e_table(table_points),cv_table(table_points))
         ALLOCATE(mu_table(table_points),kappa_table(table_points),lambda_table(table_points))

         DO i = 1,table_points
            READ(in1,*)T_table(i)
         ENDDO

         CLOSE(in1)

         ! Internal energy and specific heat (while reading the translational contribution is added)
         OPEN(UNIT=in1,FILE=path(1:len_path)//'eint',STATUS='unknown',IOSTAT=ios)
         CALL check_file (ios, 'eint')
         OPEN(UNIT=in2,FILE=path(1:len_path)//'cvint',STATUS='unknown',IOSTAT=ios)
         CALL check_file (ios, 'cvint')

         DO i = 1,table_points
            READ(in1,*)tmp1
            READ(in2,*)tmp2
            tmp3 = cv_tr*T_table(i)
            e_table(i)  = tmp3 + tmp1
            cv_table(i) = cv_tr + tmp2
         ENDDO

         CLOSE(in1)
         CLOSE(in2)
         
         ! Shear and bulk viscosity, and total thermal conductivity
         OPEN(UNIT=in1,FILE=path(1:len_path)//'eta',STATUS='unknown',IOSTAT=ios)
         CALL check_file (ios, 'eta')
         OPEN(UNIT=in2,FILE=path(1:len_path)//'kappa',STATUS='unknown',IOSTAT=ios)
         CALL check_file (ios, 'kappa')
         OPEN(UNIT=in3,FILE=path(1:len_path)//'lambda',STATUS='unknown',IOSTAT=ios)
         CALL check_file (ios, 'lambda')

         DO i = 1,table_points
            READ(in1,*)mu_table(i)
            READ(in2,*)kappa_table(i)
            READ(in3,*)lambda_table(i)
         ENDDO

         CLOSE(in1)
         CLOSE(in2)
         CLOSE(in3)
 
20  FORMAT(A,A)

    END SUBROUTINE initialize

    !----------------------------------------------------!
    ! This subroutine checks if a given file has been opened correctly 
    SUBROUTINE check_file (ios, filename)

      INTEGER, INTENT(IN) :: ios
      CHARACTER*(*), INTENT(IN) :: filename

      IF (ios.NE.0) THEN 
          WRITE(*,10)solver(1:LEN_TRIM(solver))//':: '//'"'//filename(1:LEN_TRIM(filename))//'"'//& 
                 & ' file not found, in mod_Polyat_gas_WCU_initialize_CFD ..'
          PRINT*
          STOP
       ENDIF

10  FORMAT(A)

    END SUBROUTINE check_file

    !------------------------------------------------------!
    ! This subroutine de-allocates vectors and nullifies eventually used pointers to function subroutines. 
    SUBROUTINE finalize ()

      ! Memory de-allocation
      IF (ALLOCATED(T_table))                         DEALLOCATE(T_table)
      IF (ALLOCATED(e_table))                         DEALLOCATE(e_table)
      IF (ALLOCATED(cv_table))                        DEALLOCATE(cv_table)
      IF (ALLOCATED(mu_table))                        DEALLOCATE(mu_table)
      IF (ALLOCATED(kappa_table))                     DEALLOCATE(kappa_table)
      IF (ALLOCATED(lambda_table))                    DEALLOCATE(lambda_table)

    END SUBROUTINE finalize 

  END MODULE mod_polyat_gas_wcu_initialize_CFD
!------------------------------------------------------------------------------!
