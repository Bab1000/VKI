!------------------------------------------------------!
! This program solves the Euler and Navier-Stokes equations by using different physical models (number of dimensions of the domain, 
! thermodynamics etc., etc. ). Space discretization is performed by means of Finite Volume Method (cell centred formulation) (FVMCC).
! Time integration is performed by means of schemes provided in input by the user.    
!------------------------------------------------------!
!> @author
!> Alessandro Munafò
!>
!>@brief
!>Main program
!>
!>\b DESCRIPTION: \n 
!> This program solves the Euler and Navier-Stokes equations by using different physical models (number of dimensions of the domain, 
!! thermodynamics etc., etc. ). 
!! Space discretization is performed by means of the Finite Volume Method (cell centred formulation - FVMCC).
!! Time integration is performed by means of schemes specified in input by the user. 
!! When the number of govering equations to be solved is high computations can performed faster by using the 
!! OPEN BLAS library. In this last case the number of threads must be provided as input argument to the main.
!! A local installation of the OPEN BLAS library (compiled with the option USE_OPENMP=1) is needed. 
!! The environment variable OPENMP_NUM_THREADS must be also set in order to run the code with the use of the OPEN BLAS library.  
  PROGRAM solver_fvmcc

#include "config.h"

    IMPLICIT NONE


    INTEGER :: dum
    INTEGER :: argc
    INTEGER :: nb_threads
    CHARACTER*32 :: argv
    LOGICAL :: output_only

    EXTERNAL domain_mesh, initialize, write_sol
    EXTERNAL compute_flowfield 
    EXTERNAL close_solver
    EXTERNAL post_process_solution
    EXTERNAL write_transport_fluxes 

    output_only = .false.

#ifdef GOTO_BLAS 
    PRINT*
    WRITE(*,10)'solver_fvmcc_F90:: multi-thread running option selected' 
    PRINT*
    argc = iargc() 
    IF (argc.LT.1) THEN
       WRITE(*,10)'solver_fvmcc_F90:: number of threads not specified'
       PRINT*
       STOP
    ELSE 
       CALL getarg(1,argv)
       READ(argv,20)dum
       nb_threads = INT(dum)
    ENDIF

    ! Initialize OPEN_BLAS with the number of threads specified 
    CALL openblas_set_num_threads(nb_threads)

    WRITE(*,30)'solver_fvmcc_F90:: number of threads is ',nb_threads    
    PRINT*
#else
    PRINT* 
    WRITE(*,10)'solver_fvmcc_F90:: no multi-thread running' 
    PRINT*
    argc = iargc()
    if (argc .gt. 0) then
        call getarg(1,argv)
        READ(argv,20)dum
        if (dum .eq. 1) then
            write(*,*) 'only updating output files...'
            output_only = .true.
        end if
    end if
#endif

    ! Read the input file
    CALL read_input ()

    ! Write input date (for checking and debugging purposes)
    CALL write_input ()
    
    ! Read and load the mesh 
    CALL domain_mesh ()
    
    ! Solution initialization (allocate memory and set procedure pointers)
    CALL initialize ()

    ! Compute solution 
    if (.not. output_only) then
        CALL compute_flowfield ()
    end if

    ! Write solution output file(s)
    CALL write_sol ()

    ! Solution post-processing (only if specified)
    CALL post_process_solution ()

    ! Write transport fluxes (only for viscous calculations)
    !  CALL write_transport_fluxes ()

    ! Data de-allocation (arrays and procedure pointers)
    CALL close_solver ()

10  FORMAT(A)
20  FORMAT(I4)
30  FORMAT(A,I4)

    ! Subroutines for input file reading and writing 
    CONTAINS 

    !------------------------------------------------------!
    !> This subroutine reads the input file and fills physical model and other vectors containing 
    !! common data using during the simulation.

    SUBROUTINE read_input () 

      USE mod_general_data,           ONLY: read_rate, write_rate, delta, p_in, mu_ref, T_ref, exp_visc, & 
                                          & flag_extra_data, flag_pres_elec, in_dim_id, library_name,     & 
                                          & library_data_path, viscosity_law, inter_cfl, cfl_file,        & 
                                          & xi_tol, Vdi_tol, simulation_id, state_model, ios_dir,         &
                                          & log_cfl, log_cfl_file, alpha_ser, adapt_cfl, gamma_ser, ser_cfl,  &
                                             &  ablation_library_name,  in_mixture 
      USE mod_domain_boundary
      USE mod_numerics_data 
      USE mod_physical_model
      USE mod_radiation,              ONLY: flag_radiation, coupling_period, relax, relax_factor, flag_precursor, &
                                      & flag_photochemistry,  wall_Rad
      USE mod_general_ablation,       ONLY: Flag_ablation

#ifdef CARBONABLA
      USE mod_ablation_data,          ONLY:gamma_rec, O_prob, O2_prob, N_prob, C_prob, wall_emiss, phi_pyro, exp_mdot
#endif

      IMPLICIT NONE

      INTEGER, PARAMETER :: in1 = 10
      INTEGER :: ios
      INTEGER :: i, j, m, n, sub
      INTEGER :: in_dim, in_ns, in_temp, in_tvib, in_te, in_trot
      REAL(KIND=8) :: in_gamma, in_rgas, in_pr
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: values, inter
      CHARACTER*80, PARAMETER :: empty = 'none'
      CHARACTER*80, ALLOCATABLE, DIMENSION(:) :: quantity
      CHARACTER*80 :: in_phys, in_reaction, in_transf, in_library, in_viscosity_law
      CHARACTER*80 :: line, var
      CHARACTER*80 :: in_ablation_library
      LOGICAL :: check, in_diss

      ! Initialization
      simulation_id = ''
      ! Thermodynamic library path to data files 
      library_data_path = empty

      ! Thermodynamic library (default option)
      in_library = empty 

      ! Mixture, chemical reaction kinetic scheme and energy transfer processes (default option)
      in_mixture  = empty
      in_reaction = empty
      in_transf   = empty
      state_model = empty

      ! Gas constants (default option)
      in_gamma = 1.4d0
      in_rgas  = 287.06d0
      in_pr    = 0.72d0

      ! Number of species and temperatures (default option)
      in_ns   = 1
      in_temp = 1
      in_tvib = 0
      in_trot = 0
      in_te   = 0

      ! Dissipation phenomena (default option)
      in_diss = .FALSE.
      check   = .FALSE.

      ! Diffusive flux Jacobian
      diff_flux_Jac = 'analytical'

      ! Boundary condition Jacobian (only for implicit schemes)
      bc_Jac = 'analytical'  

      ! Viscosity law and related parameters (default values)
      in_viscosity_law = 'Sutherland'
      mu_ref   = 1.458205d-6
      T_ref    = 110.333d0
      exp_visc = 0.d0

      ! Jacobian approximation (default option)
      flux_Jac_approx = 'Yoon_Jameson'

      ! Output update rate (default option)
      write_rate = 1000

      ! interactive CFL and VNN numbers (disabled as default option) 
      inter_cfl = .FALSE. 
      read_rate  = 1000

      ! Interactive file for CFL and VNN numbers
      cfl_file = empty
     
            ! History log CFL and file
      log_cfl = .FALSE.
      log_cfl_file = empty

      ! Adaptive log CFL 
      adapt_cfl = .FALSE.

      ! SER CFL
      ser_cfl = .FALSE.
      gamma_ser = 0.01
      alpha_ser = 1.0

      ! Ablation library (empty as default option)
      in_ablation_library = empty 

      WRITE(*,20)'solver_fvmcc_F90:: Reading input file'
      PRINT*

      ios_dir = .TRUE.

      OPEN(UNIT=in1,FILE='input', STATUS='old',IOSTAT=ios)

      IF(ios.NE.0) THEN
          PRINT*, 'Reading input from /input/ directory'
      OPEN(UNIT=in1,FILE='../input/input',STATUS='old',IOSTAT=ios)
          ios_dir = .FALSE.
      ENDIF
     
      IF (ios.NE.0) THEN
         PRINT*
         WRITE(*,20)'solver_fvmcc_F90:: "input" file not found'
         PRINT*
         STOP
      ENDIF
     
      line = ''
      DO WHILE (line(1:4)/='Stop')
          

         ! Physical model type and number of dimensions
         READ(in1,20)line
         IF (line=='Simulation_Name')             READ(in1,*)simulation_id
         IF (line=='Physical_model')              READ(in1,*)in_phys 
         IF (line=='Dissipation_phenomena')       READ(in1,*)in_diss
         IF (line=='Dimensions')                  READ(in1,*)in_dim

         ! Specific heat ratio, gas constant and Prandtl number (calorically perfect gas model only) 
         IF (line=='Specific_heat_ratio')         READ(in1,*)in_gamma
         IF (line=='Gas_constant')                READ(in1,*)in_rgas
         IF (line=='Prandtl_number')              READ(in1,*)in_pr
         
         ! Viscosity law
         IF (line=='Viscosity_law') THEN

             READ(in1,*)in_viscosity_law 

             viscosity_law = in_viscosity_law             

             SELECT CASE(viscosity_law)

               ! Sutherland's law for viscosity
               CASE('Sutherland')
                 READ(in1,*)mu_ref
                 READ(in1,*)T_ref

               ! Inverse power (interaction potential law)
               CASE('Inv_power')
                 READ(in1,*)mu_ref
                 READ(in1,*)T_ref
                 READ(in1,*)exp_visc

               CASE DEFAULT
                 PRINT*
                 WRITE(*,20)'solver_fvmcc_F90:: error in setting the viscosity law'
                 PRINT*
                 WRITE(*,20)'Please, modify the input file'
                 PRINT*
                 STOP

             END SELECT
 
         ENDIF

         ! More info's on the physics (nonequilibrium flows only) 
         IF (in_phys/='pg') THEN
            IF (line=='Number_of_species')            READ(in1,*)in_ns
            IF (line=='Number_of_vib_temp')           READ(in1,*)in_tvib
            IF (line=='Number_of_el_temp')            READ(in1,*)in_te
            IF (line=='Number_of_rot_temp')           READ(in1,*)in_trot
            IF (line=='Mixture')                      READ(in1,20)in_mixture
            IF (line=='State_Model')                  READ(in1,20)state_model
            IF (line=='Transfer_mechanism')           READ(in1,20)in_transf
            IF (line=='Reaction')                     READ(in1,20)in_reaction
            IF (line=='Thermodynamic_library_path')   READ(in1,20)library_data_path
            IF (line=='Thermodynamic_library') THEN
               READ(in1,20)in_library
               library_name = in_library
            ENDIF
         ENDIF 

         ! Setting the physical model parameters
         IF (in_phys=='pg') THEN
            library_name = empty  
            physical_model = set_phys_model (phys = in_phys, diss = in_diss, ndim = in_dim, ns = in_ns, ntvib = in_tvib, & 
                                       & nte = in_te, ntrot = in_trot, mixture = in_mixture, transf = in_transf,     & 
                                       & reaction = in_reaction, library = in_library, pr = in_pr, rgas = in_rgas,   & 
                                       & gamma = in_gamma)            
         ELSE 
            physical_model = set_phys_model (phys = in_phys, diss = in_diss, ndim = in_dim, ns = in_ns, ntvib = in_tvib, & 
                                       & nte = in_te, ntrot = in_trot, mixture = in_mixture, transf = in_transf,     & 
                                       & reaction = in_reaction, library = in_library, pr = in_pr, rgas = in_rgas,   & 
                                       & gamma = in_gamma)
         ENDIF

         ! Radiation coupling option
         IF (line=='Radiation_coupling') READ(in1,*) flag_radiation
         IF (line=='Precursor_effect') READ(in1,*) flag_precursor
         IF (line=='Photochemistry_source') READ(in1,*) flag_photochemistry
         IF (line=='Radiation_coupling_period') READ(in1,*) coupling_period
         IF (line=='Radiation_relaxation') READ(in1,*) relax
         IF (line=='Radiation_relaxation_factor') READ(in1,*) relax_factor
       
         ! Reading source terms
         IF (line=='Source_terms') THEN
            READ(in1,*)nb_source_terms
            ALLOCATE(source_terms(nb_source_terms),compute_source_term_Jac(nb_source_terms))
            READ(in1,20) (source_terms(i), i = 1,nb_source_terms)
            DO i = 1,nb_source_terms
               compute_source_term_Jac(i) = 'analytical'
            ENDDO
         ENDIF 

         ! Read options for source term Jacobian (numerical or analytical - default is analytical)
         IF ((line=='Source_Jac').AND.(nb_source_terms.NE.0)) THEN
            READ(in1,20) (compute_source_term_Jac(i), i = 1,nb_source_terms)
         ENDIF

         ! Read options for boundary condition Jacobian (numerical or analytical - default is analytical)
         IF (line=='Bound_Jac') THEN
            READ(in1,20)(bc_Jac(i), i = 1,2)
         ENDIF

         IF (line=='Flux_splitter')                READ(in1,*)flux_splitter
         IF (line=='Flux_Jac')                     READ(in1,*)flux_Jac_approx
         IF (line=='Diff_flux_Jac')                READ(in1,*)diff_flux_Jac
         IF (line=='Time_discretization') THEN
            READ(in1,*)time_disc
            READ(in1,*)nb_stages
            IF (nb_stages.GT.1) THEN
               ALLOCATE(stage_coeff(nb_stages))
               DO i = 1,nb_stages
                 READ(in1,*)stage_coeff(i)
               ENDDO
            ENDIF

            ! Flag for source term Jacobian
            IF (time_disc.EQ.'fe') flag_Jac = .FALSE.

            ! Flag for fully implicit scheme  
            IF (time_disc.EQ.'fi') flag_full_Impl = .TRUE.

         ENDIF

         ! Flag for stagnation line
         IF (line=='Stag_line')                    READ(in1,*)flag_stag_line

         ! Stagnation line flow geometry (0 - cylinder 1 - sphere)
         IF (line=='Stag_line_geom')               READ(in1,*)stag_line_geom

         IF (line=='CFL_number')                   READ(in1,*)cfl
         IF (line=='VNN_number')                   READ(in1,*)vnn
         IF (line=='Time_step')                    READ(in1,*)time_step
         IF (line=='Polynomial_reconstruction')    READ(in1,*)poly_rec
         IF (line=='Limiter_function')             READ(in1,*)fun_limiter
         IF (line=='Presence_electrons')           READ(in1,*)flag_pres_elec

         ! Flag for interactive cfl
         IF ((line=='Inter_CFL').OR.(line=='inter_cfl')) THEN
            READ(in1,*) inter_CFL
            READ(in1,20) cfl_file
            READ(in1,*) read_rate
         ENDIF
         ! Flag for history log cfl
         IF ((line=='Log_CFL').OR.(line=='log_cfl')) THEN
             READ(in1,*) log_cfl
             IF(log_cfl) THEN
                 READ(in1,20) log_cfl_file
                 CALL HISTORY_LOG_CFL ()
             END IF
         ENDIF
         
         IF ((line=='Adaptive_CFL').OR.(line=='adaptive_cfl')) THEN
             READ(in1,*) adapt_cfl
         ENDIF
         

         ! Flag for SER CFL
         IF ((line=='SER_CFL').OR.(line=='ser_cfl')) THEN
            ser_cfl = .TRUE.
            READ(in1,*) gamma_ser
            READ(in1,*) alpha_ser
         ENDIF
         
         ! Flag for activation of metrics in the flux computations
         IF ((line=='Metrics').OR.(line=='metrics')) READ(in1,*) flag_metrics

         ! Tolerance on species mole fractions
         IF (line=='xi_tol') READ(in1,*)xi_tol

         ! Tolerance on species diffusion velocities (used for analytical viscous Jacobians)
         IF (line=='Vdi_tol') READ(in1,*)Vdi_tol

         ! Reading the stop condition options
         IF (line=='Stop_condition') THEN 
            READ(in1,*) line
            SELECT CASE(line)

              CASE('Residual')
                stop_id = 1
                READ(in1,*)lim_res

              CASE('Time')
                stop_id = 2
                READ(in1,*)time_limit

              CASE('Iterations')
                stop_id = 3
                READ(in1,*)nb_itmax

              CASE('Residual_or_Iterations')
                stop_id = 4
                READ(in1,*)lim_res
                READ(in1,*)nb_itmax

  

              CASE DEFAULT
                PRINT*
                WRITE(*,20)'solver_fvmcc_F90:: error in setting the stop condition'
                PRINT*
                WRITE(*,20)'Please, modify the input file'
                PRINT*
                STOP
  
            END SELECT 
         ENDIF

         ! Restart options
         IF (line=='Restart')   READ(in1,*)flag_restart
         IF (line=='RestartPG') READ(in1,*)flag_restart_pg
         IF (line=='RestartEQ') READ(in1,*)flag_restart_eq

         ! Solution output file update
         IF (line=='Output_update') READ(in1,*)write_rate

         ! Output option (printing extra-data in output files)
         IF (line=='Print_extra_data')  READ(in1,*)flag_extra_data

         ! Reading boundary condition kinds and boundary condition data 
         IF (line=='Boundary_conditions') THEN
            READ(in1,*)nb_bc
            ALLOCATE(bc(nb_bc),name_bc(nb_bc))
            READ(in1,20) (name_bc(i), i = 1,nb_bc)
            READ(in1,20) (line, i =1,3 )
            bc = name_bc            

            ! Allocation of vector containing boundary data infos
            ALLOCATE(boundary_data(nb_bc))
            ALLOCATE(inlet(nb_bc*get_nb_eq (physical_model)))

            ! Species mass fractions
            IF ((get_nb_species (physical_model).GE.1).AND.(in_phys.EQ.'neq')) THEN
               ALLOCATE(yi_bound(nb_bc*get_nb_species(physical_model)))
               yi_bound = 0.d0
            ENDIF

            ! Species densities
            IF ((get_nb_species (physical_model).GE.1).AND.(in_phys.EQ.'neq')) THEN 
               ALLOCATE(rhoi_bound(nb_bc*get_nb_species(physical_model)))
               rhoi_bound = 0.d0
            ENDIF

            ! Rotational temperatures 
            IF (get_nb_trot (physical_model).GT.0) THEN
               ALLOCATE(Tr_bound(nb_bc*get_nb_trot(physical_model)))
               Tr_bound = 0.d0
            ENDIF

            ! Vibrational temperatures
            IF (get_nb_tvib (physical_model).GT.0) THEN
               ALLOCATE(Tv_bound(nb_bc*get_nb_tvib(physical_model)))
               Tv_bound = 0.d0
            ENDIF

            ! Inlet temperature vector
            ALLOCATE(inlet_temp(nb_bc*(get_nb_eq (physical_model) - get_nb_dim (physical_model) - get_nb_species (physical_model))))
            inlet_temp = 0.d0

            DO i = 1,nb_bc
      
               SELECT CASE (name_bc(i))

                 ! Subsonic inlet (total conditons for nonequilibrium flows)
                 CASE('sub_in_Tot_neq')
                   n = get_nb_species (physical_model) + 2
                   m = n - (get_nb_species (physical_model) + get_nb_tvib (physical_model) + get_nb_te (physical_model)) + &
              &        presence_tvib (physical_model)  + presence_te (physical_model) + 1

                 ! Subsonic inlet 
                 CASE ('sub_in')
                   n = get_nb_eq (physical_model) - 1
                   m = n - (get_nb_species (physical_model) + get_nb_tvib (physical_model) + get_nb_te (physical_model)) + &
              &        presence_tvib (physical_model)  + presence_te (physical_model) + 1

                 ! Subsonic inlet: radial velocity and tangential velocity derivative (Klomfass formulation),
                 ! temperature and mass fractions imposed 
                 CASE ('sub_in_vel_equi_comp')
                   n = get_nb_species (physical_model) + 4
                   m = n - get_nb_species (physical_model) + 1 


                 ! Subsonic outlet
                 CASE ('sub_out')
                   n = 1
                   m = n 

                 ! Supersonic inlet
                 CASE ('sup_in')
                   n = get_nb_eq (physical_model) 
                   m = n - (get_nb_species (physical_model) + get_nb_tvib (physical_model) + get_nb_te (physical_model)) + &
              &        presence_tvib (physical_model)  + presence_te (physical_model) + 1

                 ! No slip wall (isothermal case)
                 CASE ('no_slip_iso')
                   n = 1 
                   m = n 

                 ! No slip adiabatic wall
                 CASE ('no_slip_adiabatic')
                   n = 0 
                   m = n 

                 ! No slip ablation w/o imposing temeprature (Surface Energy Balance is solved)
                 CASE ('no_slip_SEB_ablation')
                   n = 0 
                   m = n 

                 ! No slip wall imposing the temperature (Surface Energy Balance is not solved) 
                 CASE ('no_slip_isothermal_ablation')
                   n = 1 
                   m = n 
                  
                CASE ('no_slip_iso_General_Ablation')
                   Flag_ablation = .TRUE.
                   n = 1 
                   m = n 
                
               CASE ('no_slip_seb_General_Ablation')
                   Flag_ablation = .TRUE.
                   n = 0
                   m = n 
                
                 ! Boundary conditions for which no data are imposed
                 CASE DEFAULT 
                   n = 0
                   m = n

               END SELECT 

               ! Reading data
               IF (n.GT.0) THEN 
                 
                  ALLOCATE(quantity(m))
                  ALLOCATE(values(n))
                  
                  READ(in1,*)quantity(1:m)

                  DO j = 1,n 
                     READ(in1,*)values(j)
                  ENDDO

                  ! Creation of boundary condition subroutine names
                  DO j = 1,m 
                     bc(i) = bc(i)(1:LEN_TRIM(bc(i)))//'_'//quantity(j) 
                  ENDDO

                  ! Setting boundary condition data
                  CALL set_bc_value (i, quantity, values, in_phys)
                  
                  DEALLOCATE(quantity)
                  DEALLOCATE(values)
                
               ELSEIF ((n.EQ.0).AND.(bc(i)(1:LEN_TRIM(bc(i))).EQ.'no_slip_adiabatic')) THEN 

                  bc(i) = bc(i)(1:LEN_TRIM(bc(i))) 

               ELSEIF ((n.EQ.0).AND.(bc(i)(1:LEN_TRIM(bc(i))).EQ.'no_slip_SEB_ablation')) THEN 

                  bc(i) = bc(i)(1:LEN_TRIM(bc(i))) 
                
               ENDIF

            ENDDO

         ENDIF

         ablation_library_name = in_ablation_library
         IF (line=='Ablation_library') THEN
            READ(in1,20) in_ablation_library
            ablation_library_name = in_ablation_library
            IF (ablation_library_name=='Custom') THEN
               read(in1,*) O_prob
               read(in1,*) O2_prob
               read(in1,*) N_prob
               read(in1,*) C_prob
            ELSEIF (ablation_library_name=='Custom_nitri') THEN
               read(in1,*) N_prob
               read(in1,*) C_prob
            ELSEIF (ablation_library_name=='Nitri_probability_rebuild') THEN
               read(in1,*) exp_mdot  
            ENDIF
            ELSE
         ENDIF

         IF (line=='Nitrogen_recombination_coefficient')  READ(in1,*) gamma_rec 

         IF (line=='Wall_emissivity')  READ(in1,*) wall_emiss  

         IF (line=='Pyrolysis/Ablation mass injection ratio')  READ(in1,*) phi_pyro

         ! Reading the data relative to the initial field
         IF (flag_restart.EQV..FALSE.) THEN

            IF (line=='Initial_field_var') READ(in1,*)var
            IF (var=='x') THEN 
               in_dim_id = 1
            ELSEIF (var=='y') THEN 
               in_dim_id = 2  
            ENDIF
 
            IF (line=='Number_of_subintervals') THEN
               READ(in1,*)sub
               ALLOCATE(inter(sub))
               READ(in1,*)(inter(i), i = 1,sub)
               
               ALLOCATE(p_in(get_nb_eq(physical_model),sub))

               READ(in1,20)(line, i = 1,2)
              
               ! Initial field
               DO i = 1,get_nb_eq(physical_model)
                  READ(in1,*)p_in(i,1:sub)
               ENDDO
               
               ! Increments 
               ALLOCATE(delta(sub))
               DO i = 1,sub 
                  delta(i) = inter(i)
               ENDDO
               DEALLOCATE(inter)

            ENDIF
         ENDIF

      END DO
      
      CLOSE(in1)

20  FORMAT(A)
    
    END SUBROUTINE read_input 

    !------------------------------------------------------!
    ! This subroutine assigns the physical data read in the input file for the boundary conditions.
    !------------------------------------------------------!
    !> @author
    !> Alessandro Munafò
    !>
    !>@brief
    !> Assigns phisical data for BC
    !>
    !>\b DESCRIPTION: 
    !>
    !> This subroutine assigns the physical data read in the input file for the boundary conditions.
    !------------------------------------------------------!
    SUBROUTINE set_bc_value (nb, label, values, model_name)

      USE mod_domain_boundary
      USE mod_physical_model
  
      IMPLICIT NONE    
  
      INTEGER :: i, j, k, m, neq, ns, ntrot, ntvib, nte, nt, pos
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: yi_in, rhoi_in, Tr_in, Tv_in, sup_inlet
  
      INTEGER, INTENT(IN) :: nb
      CHARACTER(*), INTENT(IN) :: model_name
      CHARACTER*(*), DIMENSION(:), INTENT(IN) :: label
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: values 
  
      ! Useful data   
      neq   = get_nb_eq (physical_model)
      ns    = get_nb_species (physical_model)
      ntrot = get_nb_trot (physical_model)
      ntvib = get_nb_tvib (physical_model)
      nte   = get_nb_te (physical_model)
      nt    = 1 + ntrot + ntvib + nte
       
      ! Size of the array
      k = 0    
      m = SIZE(label)
     
      ! Setting the boundary data according to what specified in the input file 
      CALL set_boundary_id (nb, boundary_data(nb))
      
      DO i = 1,m

         SELECT CASE (label(i))
  
           ! Species mass fractions at inlet
           CASE ('yiin')
             ALLOCATE(yi_in(ns))
             DO j = 1,ns
                yi_in (j) = values(k + j)
             ENDDO
             k = k + ns
  
           ! Species densities at inlet
           CASE ('rhoiin')
             ALLOCATE(rhoi_in(ns))
             DO j = 1,ns
                rhoi_in (j) = values(k + j)
             ENDDO
             k = k + ns 
  
           ! Rotational temperature(s) at inlet
           CASE ('Trotin')
             ALLOCATE(Tr_in(ntrot))
             DO j = 1,ntrot
                Tr_in(j) = values(k + j)
             ENDDO
             k = k + ntrot
  
           ! Vibrational temperature(s) at inlet
           CASE ('Tvin')
             ALLOCATE(Tv_in(ntvib))
             DO j = 1,ntvib
                Tv_in(j) = values(k + j)
             ENDDO
             k = k + ntvib
             
           ! Static temperature at inlet
           CASE ('Tin')
             CALL set_boundary_Tin (values(k + 1), boundary_data(nb))
             k = k + 1 
  
           ! Free electron electronic temperature at inlet
           CASE ('Tein') 
             CALL set_boundary_Tein (values(SIZE(values)), boundary_data(nb))
  
           ! Static pressure at inlet
           CASE ('pin')
             CALL set_boundary_pin (values(k + 1), boundary_data(nb))       
             k = k + 1  
  
           ! Static pressure at outlet
           CASE ('pout')
             CALL set_boundary_pout (values(k + 1), boundary_data(nb))
             k = k + 1  

           ! Static temperature at outlet
           CASE ('Tout')
             CALL set_boundary_Tout (values(k + 1), boundary_data(nb))
             k = k + 1  
  
           ! Static density at inlet
           CASE ('rhoin')
             CALL set_boundary_rhoin (values(k + 1), boundary_data(nb)) 
             k = k + 1
  
           ! Total pressure at inlet
           CASE ('p0')
             CALL set_boundary_p0 (values(k + 1), boundary_data(nb))  
             k = k + 1  
  
           ! Total temperature at inlet
           CASE ('T0')
             CALL set_boundary_T0 (values(k + 1), boundary_data(nb)) 
             k = k + 1 
   
           ! Derivative of the "physical" velocity (normal to the stagnation
           ! line) in the direction normal to the stagnation line
           CASE ('dv_dyin')
             CALL set_boundary_dv_dyin (values(k + 1), boundary_data(nb)) 
             k = k + 1 

           ! Flow angle at inlet
           CASE ('alpha')
             CALL set_boundary_alpha (values(k + 1), boundary_data(nb))
             k = k + 1  
    
           ! x-component of velocity at inlet
           CASE ('uin')
             CALL set_boundary_uin (values(k + 1), boundary_data(nb))
             k = k + 1  
  
           ! y-component of velocity at inlet
           CASE ('vin')
             CALL set_boundary_vin (values(k + 1), boundary_data(nb))
             k = k + 1  
  
           ! Wall temperature 
           CASE ('Twall')
             CALL set_boundary_Twall (values(k + 1), boundary_data(nb))
             k = k + 1  
  
           CASE DEFAULT
             WRITE(*,10)'Error in boundary condition, modify the input file'
             PRINT*
             STOP

         END SELECT 
      ENDDO
  
      ! Supersonic inlet boundary condition (primitive variables stored in a unique vector)
      IF (name_bc(nb)=='sup_in') THEN
  
         bc(nb) = ''  
         bc(nb) = 'sup_in'
         pos = (nb - 1)*neq
         DO i = 1,neq
            inlet(pos + i) = values(i)
         ENDDO
 
         pos = (nb - 1)*nt       

         ! Translational temperarure
         inlet_temp(pos + 1) = get_boundary_Tin (boundary_data(nb)) 
  
         ! Rotational temperatures
         IF (ALLOCATED(Tr_in).AND.(ntrot.GT.0).AND.(model_name.EQ.'neq')) THEN
            DO i = 1,ntrot
               inlet_temp(pos + 1 + i) = Tr_in(i)
            ENDDO
         ENDIF
            
         ! Vibrational temperatures
         IF (ALLOCATED(Tv_in).AND.(ntvib.GT.0).AND.(model_name.EQ.'neq')) THEN
            DO i = 1,ntvib
               inlet_temp(pos + 1 + ntrot + i) = Tv_in(i)
            ENDDO
         ENDIF
       
         ! Free-electron electronic temperature
         IF ((nte.NE.0).AND.(model_name.EQ.'neq')) THEN 
            inlet_temp(pos + 1 + ntrot + ntvib + 1) = get_boundary_Tein (boundary_data(nb)) 
         ENDIF 
 
      ENDIF 
          
      ! Subsonic inlet boundary condition
      IF (name_bc(nb)=='sub_in') THEN
         
         pos = (nb - 1)*nt       

         ! Translational temperarure
         inlet_temp(pos + 1) = get_boundary_Tin (boundary_data(nb)) 
  
         ! Rotational temperatures
         IF (ALLOCATED(Tr_in).AND.(ntrot.GT.0).AND.(model_name.EQ.'neq')) THEN
            DO i = 1,ntrot
               inlet_temp(pos + 1 + i) = Tr_in(i)
            ENDDO
         ENDIF
            
         ! Vibrational temperatures
         IF (ALLOCATED(Tv_in).AND.(ntvib.GT.0).AND.(model_name.EQ.'neq')) THEN
            DO i = 1,ntvib
               inlet_temp(pos + 1 + ntrot + i) = Tv_in(i)
            ENDDO
         ENDIF
       
         ! Free-electron electronic temperature
         IF ((nte.NE.0).AND.(model_name.EQ.'neq')) THEN 
            inlet_temp(pos + 1 + ntrot + ntvib + 1) = get_boundary_Tein (boundary_data(nb)) 
         ENDIF 
  
      ELSEIF (name_bc(nb)=='sub_in_vel_equi_comp') THEN
         bc(nb) = ''  
         bc(nb) = 'sub_in_vel_equi_comp'




      ENDIF

      ! Boundary species mass fractions
      IF ((ALLOCATED(yi_in)).AND.(ns.GE.1).AND.(model_name.EQ.'neq')) THEN 
         pos = ns*(nb - 1)
         DO i = 1,ns
            yi_bound(pos + i) = yi_in(i)
         ENDDO
      ENDIF
  
      ! Boundary species densities
      IF ((ALLOCATED(rhoi_in)).AND.(ns.GE.1).AND.(model_name.EQ.'neq')) THEN 
         pos = ns*(nb - 1)
         DO i = 1,ns
            rhoi_bound(pos + i) = rhoi_in(i)
         ENDDO
      ENDIF
         
      ! Boundary rotational temperatures
      IF ((ALLOCATED(Tr_in)).AND.(ntrot.GT.0).AND.(model_name.EQ.'neq')) THEN 
         pos = ntrot*(nb - 1)
         DO i = 1,ntrot
            Tr_bound(pos + i) = Tr_in(i)
         ENDDO
      ENDIF
  
      ! Boundary vibrational temperatures
      IF (ALLOCATED(Tv_in).AND.(ntvib.GT.0).AND.(model_name.EQ.'neq')) THEN 
         pos = ntvib*(nb - 1)
         DO i = 1,ntvib
            Tv_bound(pos + i) = Tv_in(i)
         ENDDO
      ENDIF
  
      ! Vector de-allocation
      IF (ALLOCATED(yi_in))          DEALLOCATE(yi_in)
      IF (ALLOCATED(rhoi_in))        DEALLOCATE(rhoi_in)
      IF (ALLOCATED(Tr_in))          DEALLOCATE(Tr_in)
      IF (ALLOCATED(Tv_in))          DEALLOCATE(Tv_in)

10  FORMAT(A)

    END SUBROUTINE set_bc_value 

    !------------------------------------------------------!
    ! This subroutine writes the input file data for checking purposes (input_data in FVMCC_F90/output).
    !------------------------------------------------------!
    !> @author
    !> Alessandro Munafò
    !>
    !>@brief
    !> Re-writes input data for check
    !>
    !>\b DESCRIPTION: 
    !>
    !> This subroutine writes the input file data for checking purposes (\e "input_data" in \e "FVMCC_F90/output").
    !------------------------------------------------------!
    SUBROUTINE write_input ()    

      USE mod_general_data,       ONLY: flag_extra_data, library_name, library_data_path, & 
                                      & xi_tol, Vdi_tol, simulation_id, ios_dir
      USE mod_numerics_data
      USE mod_domain_boundary
      USE mod_physical_model
  
      IMPLICIT NONE
  
      INTEGER :: i, j, neq, ns, ntrot, ntvib, nte, ntemp, ndim 
      INTEGER, PARAMETER :: out1 = 10
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: inlet_data, species_dens, species_mass_frac, temp_rot, temp_vib
   
      ! Useful data
      neq   = get_nb_eq (physical_model)
      ns    = get_nb_species (physical_model)
      ntrot = get_nb_trot (physical_model)
      ntvib = get_nb_tvib (physical_model)
      nte   = get_nb_te (physical_model)
      ntemp = get_nb_temp (physical_model)
      ndim  = get_nb_dim (physical_model)
     
      IF (ntrot.GT.0) ALLOCATE(temp_rot(ntrot))
      IF (ntvib.GT.0) ALLOCATE(temp_vib(ntvib))
 
      IF (ios_dir) THEN
          OPEN(UNIT=out1,FILE='./'//TRIM(simulation_id)//'_input_data',STATUS='unknown')
      ELSE
          OPEN(UNIT=out1,FILE='../output/'//TRIM(simulation_id)//'_input_data',STATUS='unknown')
      ENDIF
      WRITE(out1,20)'Input data specified by the user: '                                               
      WRITE(out1,*)
      WRITE(out1,30)'simulation name', simulation_id
      WRITE(out1,*)
      WRITE(out1,30)'physical model: ',get_name (physical_model)                                         
      WRITE(out1,*)
      WRITE(out1,30)'thermodynamic library: ',library_name                                           
      WRITE(out1,*)
      WRITE(out1,30)'thermodynamic library data path: ',library_data_path(1:LEN_TRIM(library_data_path)) 
      WRITE(out1,*)
      WRITE(out1,60)'gamma: ',get_gamma (physical_model)                                             
      WRITE(out1,*)
      WRITE(out1,60)'rgas: ',get_rgas (physical_model)                                               
      WRITE(out1,*)
      WRITE(out1,60)'Prandtl number: ',get_prandtl (physical_model)                                  
      WRITE(out1,*)
      WRITE(out1,30)'Mixture name: ',get_mixture (physical_model)                                         
      WRITE(out1,*)  
      WRITE(out1,30)'Reaction name: ',get_reaction (physical_model)                                       
      WRITE(out1,*)
      WRITE(out1,30)'Transfer mechanism name: ',get_transfer (physical_model)                             
      WRITE(out1,*)
      WRITE(out1,50)'Number of species: ',ns                                                       
      WRITE(out1,*)
      WRITE(out1,50)'Number of rotational temperatures: ',ntrot                                    
      WRITE(out1,*)
      WRITE(out1,50)'Number of vibrational temperatures: ',ntvib                                   
      WRITE(out1,*)
      WRITE(out1,50)'Number of free-electron electronic temperatures: ',nte                        
      WRITE(out1,*)
      WRITE(out1,50)'Number of temperatures: ',ntemp                                               
      WRITE(out1,*)
      WRITE(out1,50)'Number of dimensions: ',ndim                                               
      WRITE(out1,*) 
      WRITE(out1,50)'Number of conservation equations: ',neq                                       
      WRITE(out1,*)  
      WRITE(out1,40)'Inclusion of dissipation phenomena: ',get_diss_flag (physical_model)                 
      WRITE(out1,*) 
      WRITE(out1,40)'Flag stagnation line: ',flag_stag_line                 
      WRITE(out1,*)  
      IF (flag_stag_line.EQV..TRUE.) THEN
         IF (stag_line_geom.EQ.0) THEN
            WRITE(out1,90)'Stagnation line geometry: ',stag_line_geom,' - cylinder'                 
            WRITE(out1,*)  
         ELSEIF (stag_line_geom.EQ.1) THEN
            WRITE(out1,90)'Stagnation line geometry: ',stag_line_geom,' - sphere'                 
            WRITE(out1,*)  
         ENDIF
      ENDIF
      WRITE(out1,30)'Flux splitter: ',flux_splitter                                                   
      WRITE(out1,*)
      IF (flag_full_Impl.EQV..TRUE.) THEN
         WRITE(out1,30)'Flux Jacobian approximation: ',flux_Jac_approx                                                   
         WRITE(out1,*) 
      ENDIF
      WRITE(out1,30)'Active metrics: ',flag_metrics
      WRITE(out1,30)'Polynomial reconstruction: ',poly_rec                                            
      WRITE(out1,*) 
      WRITE(out1,30)'Limiter function: ',fun_limiter                                                  
      WRITE(out1,*) 
      WRITE(out1,30)'Time discretization scheme: ',time_disc                                          
      WRITE(out1,*) 
      WRITE(out1,50)'Number of stages: ',nb_stages                                          
      WRITE(out1,*)
      IF (nb_stages.GT.1) THEN
         DO i = 1,nb_stages 
            WRITE(out1,80)'Coefficient of stage ',i,' :',stage_coeff(i)
         ENDDO
         WRITE(out1,*)
      ENDIF
      WRITE(out1,40)'Restart from a previous solution: ',flag_restart                                 
      WRITE(out1,*) 
      WRITE(out1,60)'CFL number: ',cfl                                                            
      WRITE(out1,*)  
      WRITE(out1,60)'VNN number: ',vnn                                                            
      WRITE(out1,*)  
      WRITE(out1,60)'Time step: ',time_step                                                       
      WRITE(out1,*) 
      WRITE(out1,60)'Tolerance on molar fraction ',xi_tol
      WRITE(out1,*)
      WRITE(out1,60)'Tolerance on diffusion velocity (only for analytical viscous Jacobians) ',Vdi_tol
      WRITE(out1,*)
      WRITE(out1,45)'Number of source terms: ',nb_source_terms                                       
      WRITE(out1,*)

      IF (nb_source_terms.GT.0) THEN 
         DO i = 1,nb_source_terms 
            WRITE(out1,55)'Source term ',i,' is: ', source_terms(i)     
         WRITE(out1,*) 
         ENDDO
      ENDIF
      WRITE(out1,20)'Stop condition'                                                                    
      WRITE(out1,*) 
      SELECT CASE (stop_id)
  
        CASE (1) 
          WRITE(out1,60)'Residual value: ',lim_res                                              
          WRITE(out1,*) 
  
        CASE (2)
          WRITE(out1,60)'Time level: ',time_limit                                               
          WRITE(out1,*) 
  
        CASE (3) 
          WRITE(out1,'(A,I8)')'Number of iterations: ',nb_itmax                                           
          WRITE(out1,*) 
  
      END SELECT 
  
      WRITE(out1,40)'Print extra data: ',flag_extra_data                                              
      WRITE(out1,*)
  
      WRITE(out1,45)'Number of boundary conditions: ',nb_bc                                          
      WRITE(out1,*)
      DO i = 1,nb_bc  
         WRITE(out1,55)'Boundary condition on side ', i,' is: ',name_bc(i)        
         WRITE(out1,*) 
      ENDDO

      WRITE(out1,20)'Boundary condition data: '                                                         
      WRITE(out1,*) 
  
      DO i = 1,nb_bc 

         WRITE(out1,55)'bc ',i,' - ',name_bc(i)                                                          
         WRITE(out1,*) 
  
         SELECT CASE (get_name (physical_model))
  
           ! Calorically perfect gas case
           CASE ('pg')
  
             SELECT CASE (name_bc(i)) 
        
               ! Subsonic inlet  
               CASE ('sub_in')
  
                 WRITE(out1,70)'Inlet x-velocity component: ',get_boundary_uin (boundary_data(i))                   
                 WRITE(out1,*)
                 IF (ndim.EQ.2) THEN 
                    WRITE(out1,70)'Inlet y-velocity component: ',get_boundary_vin (boundary_data(i))                
                    WRITE(out1,*)
                 ENDIF
                 WRITE(out1,70)'Inlet static density:  ',get_boundary_rhoin (boundary_data(i))                       
                 WRITE(out1,*)
                 WRITE(out1,70)'Inlet static pressure: ',get_boundary_pin (boundary_data(i))                        
                 WRITE(out1,*)
                 WRITE(out1,70)'Inlet static temperature: ',get_boundary_Tin (boundary_data(i))                     
                 WRITE(out1,*)
                 WRITE(out1,70)'Inlet total pressure: ',get_boundary_p0 (boundary_data(i))                          
                 WRITE(out1,*)
                 WRITE(out1,70)'Inlet total temperature: ',get_boundary_T0 (boundary_data(i))                       
                 WRITE(out1,*)
  
               ! Subsonic outlet
               CASE ('sub_out')
                 WRITE(out1,70)'Outlet static pressure: ',get_boundary_pout (boundary_data(i))                      
                 WRITE(out1,*)
                 WRITE(out1,70)'Outlet static temperature: ',get_boundary_Tout (boundary_data(i))                   
                 WRITE(out1,*)
                  
               ! Supersonic inlet
               CASE ('sup_in')
                  ALLOCATE(inlet_data(neq))
  
                  inlet_data = get_boundary_inlet (neq, boundary_data(i))
                  DO j = 1,NEQ
                     WRITE(out1,80)'Inlet variable ',j,': ',inlet_data(j)                                    
                     WRITE(out1,*)
                  ENDDO
  
                  DEALLOCATE(inlet_data) 
  
               ! No slip wall (isothermal)
               CASE ('no_slip_iso')
                 WRITE(out1,70)'Wall temperature: ',get_boundary_Twall(boundary_data(i))                           
                 WRITE(out1,*)
  
               CASE DEFAULT 
                 WRITE(out1,20)'No data specified'                                                                         
                 WRITE(out1,*) 
  
             END SELECT 
  
           ! Nonequilibrium case 
           CASE DEFAULT
              
             SELECT CASE (name_bc(i))
  
               ! Subsonic inlet
               CASE ('sub_in')
  
                 IF (ALLOCATED(yi_bound).EQV..TRUE.) THEN
                    ALLOCATE(species_mass_frac(ns))
                    species_mass_frac = get_boundary_yi(ns, boundary_data(i))
                    WRITE(out1,80)('Mass fraction at inlet of species ',j,': ',species_mass_frac(j), j = 1,ns)
                 ENDIF 
                 IF (ALLOCATED(rhoi_bound).EQV..TRUE.) THEN
                    ALLOCATE(species_dens(ns))
                    species_dens = get_boundary_rhoi(ns, boundary_data(i))
                    WRITE(out1,80)('Density at inlet of species ',j,': ',species_dens(j),  j = 1,ns)    
                 ENDIF
                 WRITE(out1,*)
  
                 WRITE(out1,70)'Inlet x-velocity component: ',get_boundary_uin (boundary_data(i))                   
                 WRITE(out1,*)
                 IF (ndim.EQ.2) THEN
                    WRITE(out1,70)'Inlet y-velocity component: ',get_boundary_vin (boundary_data(i))                
                    WRITE(out1,*)
                 ENDIF
  
                 WRITE(out1,70)'Static temperature at inlet: ',get_boundary_Tin (boundary_data(i))                  
                 WRITE(out1,*)
               
                 IF (presence_trot (physical_model).EQ.1) THEN
                    temp_rot = get_boundary_Tr(ntrot, boundary_data(i))
                    WRITE(out1,80) ('Rotational temperature ',j,' at inlet: ',temp_rot(j), j = 1,ntrot)              
                    WRITE(out1,*) 
                 ENDIF
  
                 IF (presence_tvib (physical_model).EQ.1) THEN
                    temp_vib = get_boundary_Tv(ntvib, boundary_data(i))
                    WRITE(out1,80) ('Vibrational temperature ',j,' at inlet: ',temp_vib(j), j = 1,ntvib)              
                    WRITE(out1,*) 
                 ENDIF
  
                 IF (presence_te (physical_model).EQ.1) THEN 
                    WRITE(out1,70) 'Free-electron electronic temperature at inlet: ',& 
                                            & get_boundary_Tein (boundary_data(i)) 
                    WRITE(out1,*)
                 ENDIF
  
               ! Subsonic outlet 
               CASE ('sub_out')
                 WRITE(out1,70)'Outlet static pressure: ',get_boundary_pout (boundary_data(i))                       
                 WRITE(out1,*)
                 WRITE(out1,70)'Outlet static temperature: ',get_boundary_Tout (boundary_data(i))                    
                 WRITE(out1,*)
  
               ! Supersonic inlet
               CASE ('sup_in')
                 ALLOCATE(inlet_data(neq))
  
                 inlet_data = get_boundary_inlet (neq, boundary_data(i))
                 DO j = 1,neq
                    WRITE(out1,80)'Inlet variable ',j,': ',inlet_data(j)                                      
                    WRITE(out1,*)
                 ENDDO
  
                 DEALLOCATE(inlet_data)
  
               ! No slip wall (isothermal)
               CASE ('no_slip_iso')
                 WRITE(out1,70)'Wall temperature: ',get_boundary_Twall (boundary_data(i))                            
                 WRITE(out1,*)

               CASE ('no_slip_iso_General_Ablation')
                 WRITE(out1,70)'Wall temperature: ',get_boundary_Twall (boundary_data(i))                            
                 WRITE(out1,*)
  
               CASE DEFAULT 
                 WRITE(out1,20)'No data specified'                                                                          
                 WRITE(out1,*) 
  
             END SELECT
  
         END SELECT 
      ENDDO
  
      WRITE(out1,20)'Boundary condition subroutine names: '                                                                
      WRITE(out1,*)

      DO i = 1,nb_bc
         WRITE(out1,'(A,I1,A,A)')'Boundary condition ',i,' subroutine name: ',bc(i)                                           
         WRITE(out1,*)
      ENDDO
  
      CLOSE(out1)
  
      ! Vector de-allocation 
      IF (ALLOCATED(inlet_data))          DEALLOCATE(inlet_data)
      IF (ALLOCATED(species_dens))        DEALLOCATE(species_dens)
      IF (ALLOCATED(species_mass_frac))   DEALLOCATE(species_mass_frac)
      IF (ALLOCATED(temp_rot))            DEALLOCATE(temp_rot)
      IF (ALLOCATED(temp_vib))            DEALLOCATE(temp_vib)

20  FORMAT(A)
30  FORMAT(A,A)
40  FORMAT(A,L)
45  FORMAT(A,I2)
50  FORMAT(A,I4)
55  FORMAT(A,I2,A,A)
60  FORMAT(A,F15.8)
70  FORMAT(A,F25.10)
80  FORMAT(A,I4,A,E25.10)
90  FORMAT(A,I2,A)

    END SUBROUTINE write_input    
    !------------------------------------------------------!
    SUBROUTINE HISTORY_LOG_CFL ()

   
    USE mod_general_data,           ONLY: log_cfl_file, log_cfl_number, log_cfl_ite, log_cfl_file_lines

    INTEGER, PARAMETER :: in1 = 1
    INTEGER :: lines, i

    OPEN(UNIT=in1,FILE=log_cfl_file(1:LEN_TRIM(log_cfl_file)),STATUS='unknown')

    READ(in1,*) log_cfl_file_lines
    
    ALLOCATE(log_cfl_number(log_cfl_file_lines))
    ALLOCATE(log_cfl_ite(log_cfl_file_lines))

    DO i = 1,log_cfl_file_lines
    	READ(in1,*) log_cfl_number(i), log_cfl_ite(i)
    ENDDO
    CLOSE(in1)
    END SUBROUTINE HISTORY_LOG_CFL
    !------------------------------------------------------!
  END PROGRAM solver_fvmcc
!------------------------------------------------------------------------------!
