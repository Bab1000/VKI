!------------------------------------------------------!
! This module contains general data.
!------------------------------------------------------!
!> @author
!> Alessandro Munafò
!>
!>@brief
!>Contains general data:w
!>
!>\b DESCRIPTION: 
!>
!> This module contains general data.
!------------------------------------------------------!
  MODULE mod_general_data

     IMPLICIT NONE

     INTEGER, SAVE :: in_dim_id         !< id of the coordinate where sub-division is applied for initial field
     INTEGER, SAVE :: nb_cells          !< number of physical cells 
     INTEGER, SAVE :: nb_nodes          !< number of nodes
     INTEGER, SAVE :: nb_ghost_cells    !< number of ghost cells
     INTEGER, SAVE :: nb_tot_cells      !< number of cells (physical + ghost)
     INTEGER, SAVE :: nb_rows           !< number of rows in the mesh 
     INTEGER, SAVE :: nb_cols           !< number of columns in the mesh 
     INTEGER, SAVE :: nb_cell_rows      !< number of cells in each row of the mesh 
     INTEGER, SAVE :: nb_eq             !< number of equations 
     INTEGER, SAVE :: nb_prop           !< number of physical properties 
     INTEGER, SAVE :: nb_ns             !< number of species
     INTEGER, SAVE :: nb_ele            !< number of elements
     INTEGER, SAVE :: nb_dim            !< number of dimensions 
     INTEGER, SAVE :: nb_trot = 0       !< number of rotational temperatures
     INTEGER, SAVE :: nb_tvib = 0       !< number of vibrational temperatures
     INTEGER, SAVE :: nb_te   = 0       !< number of free-electron (electronic) temperatures
     INTEGER, SAVE :: nb_int_temp = 0   !< number of internal temperatures
     INTEGER, SAVE :: nb_temp           !< number of temperatures
     INTEGER, SAVE :: pos_em = 0        !< position of the free-electron species within the species list
     INTEGER, SAVE :: pos_rho           !< position of the density in the primitive variable vector
     INTEGER, SAVE :: pos_u             !< position of the x-component of velocity in the primitive variable vector
     INTEGER, SAVE :: pos_v             !< position of the y-component of velocity in the primitive variable vector
     INTEGER, SAVE :: pos_p             !< position of the pressure in the primitive variable vector
     INTEGER, SAVE :: pos_T             !< position of the temperature in the primitive variable vector
     INTEGER, SAVE :: pos_Tk            !< position of the first internal temperature in the primitive variable vector
     INTEGER, SAVE :: pos_Te            !< position of the free-electron (electronic) temperature in the primitive variable vector
     INTEGER, SAVE :: pos_rhou          !< position of the x-component of momentum density in the conservative variable vector
     INTEGER, SAVE :: pos_rhov          !< position of the y-component of momentum density in the conservative variable vector  
     INTEGER, SAVE :: pos_rhoE          !< position of the total energy density in the conservative variable vector
     INTEGER, SAVE :: pos_rhoek         !< position of the nonequilibrium internal energy density in the conservative variable vector 
     INTEGER, SAVE :: pos_rhose         !< position of the free-electron pseudo-entropy in the conservative variable vector 
     INTEGER, SAVE :: pos_u_cell        !< position of the x-component of velocity in the physical property vector 
     INTEGER, SAVE :: pos_v_cell        !< position of the y-component of velocity in the physical property vector 
     INTEGER, SAVE :: pos_h0_cell       !< position of the specific total enthalpy in the physical property vector
     INTEGER, SAVE :: pos_c_cell        !< position of the frozen speed of sound in the physical property vector
     INTEGER, SAVE :: pos_pres_cell     !< position of the pressure in the physical property vector
     INTEGER, SAVE :: pos_ek_cell       !< position of the kinetic energy per unit mass in the physical property vector  
     INTEGER, SAVE :: pos_rho_cell      !< position of the density in the physical property vector
     INTEGER, SAVE :: pos_T_cell        !< position of the translational temperature in the physical property vector
     INTEGER, SAVE :: pos_alpha_cell    !< position of the \f$ \alpha \f$ factor in the physical property vector
     INTEGER, SAVE :: pos_beta_cell     !< position of the \f$ \beta \f$ factor in the physical property vector 
     INTEGER, SAVE :: pos_betak_cell    !< position of the first \f$ \beta_k \f$ factor in the physical property vector
     INTEGER, SAVE :: pos_chi_cell      !< position of the first species thermal diffusion ratio in the physical property vector
     INTEGER, SAVE :: pos_Di_cell       !< position of the first species diffusion coefficient \f$ D_i \f$ in the physical property vector
     INTEGER, SAVE :: pos_gamma_cell    !< position of the frozen specific heat ratio in the physical property vector 
     INTEGER, SAVE :: pos_ei_cell       !< position of the first species specific energy in the physical property vector
     INTEGER, SAVE :: pos_eki_cell      !< position of the first entry of the nonequilibrium energy matrix in the physical property vector 
     INTEGER, SAVE :: pos_lambda_cell   !< position of the first component of the thermal conductivity in the physical property vector
     INTEGER, SAVE :: pos_mu_cell       !< position of the dynamic viscosity in the physical property vector
     INTEGER, SAVE :: pos_kappa_cell    !< position of the bull viscosity in the physical property vector
     INTEGER, SAVE :: start_u_phys      !< start location of physical cell conservative variables
     INTEGER, SAVE :: start_prop_phys   !< start location of physical cell physical properties 
     INTEGER, SAVE :: finish_u_phys     !< finish location of physical cell conservative variables
     INTEGER, SAVE :: finish_prop_phys  !< finish location of physical cell physical properties 
     INTEGER, SAVE :: write_rate        !< output file writing rate
     INTEGER, SAVE :: read_rate         !< CFL and VNN file reading rate
     INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: nb_cell_cols      !< number of cells in each column of the mesh 
     REAL(KIND=8), PARAMETER :: eta = 1.d-5                        !< factor used for numerical Jacobians
     REAL(KIND=8), PARAMETER :: gamma_e = 5.d0/3.d0                !< specific heat ratio of free-electron gas
     REAL(KIND=8), PARAMETER :: gamma_e_m1 = gamma_e - 1.d0        !< useful common factor for free-electron gas
     REAL(KIND=8), PARAMETER :: ov_gamma_e_m1 = 1.d0/gamma_e_m1    !< useful common factor for free-electron gas
     REAL(KIND=8), SAVE :: gamma                                   !< specific heat ratio
     REAL(KIND=8), SAVE :: Pr                                      !< Prandtl number
     REAL(KIND=8), SAVE :: R_gas                                   !< specific gas constant 
     REAL(KIND=8), SAVE :: cv_gas                                  !< constant volume specific heat           
     REAL(KIND=8), SAVE :: cp_gas                                  !< constant pressure specific heat
     REAL(KIND=8), SAVE :: d_ref                                   !< reference value for diameter 
     REAL(KIND=8), SAVE :: mu_ref                                  !< reference value for dynamic viscosity 
     REAL(KIND=8), SAVE :: T_ref                                   !< reference value for temperature 
     REAL(KIND=8), SAVE :: exp_visc                                !< exponent of the viscosity law 
     REAL(KIND=8), SAVE :: Re                                      !< free-electron specifici gas constrant
     REAL(KIND=8), SAVE :: p_inf                                   !< free-stream pressure 
     REAL(KIND=8), SAVE :: V_inf                                   !< free-stream velocity
     REAL(KIND=8), SAVE :: xi_tol = 1.d-15                         !< tolerance on species mole fractions
     REAL(KIND=8), SAVE :: Vdi_tol = 1.d-12                        !< tolerance on species diffusion velocity
     REAL(KIND=8), SAVE :: dx_vol
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: cell_prop    !< physical properties of all cells 
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: u            !< conservative variables of all cells   
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: u_old        !< conservative variables of all cells (only for multi-stage schemes)
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: du           !< variation of conservative variables between successive time-steps
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: dres         !< variation of conservative variables for the residual
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: ru           !< rhs residual in terms of conservative variables for all cells
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: u_vec        !< conservative variables (working vector)
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: phys_prop    !< physical properties (working vector)
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: p            !< primitive variables
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: delta        !< vector of increments for initial field  
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: geom_source  !< geometrical factor for all cells (for geometrical-like source terms)
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: s_vec        !< source term (working vector) 
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), SAVE :: js_mat     !< source term Jacobian (working matrix)
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: volumes      !< cell volumes
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: x            !< x coordinate of cell nodes
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: y            !< y coordinate of cell nodes
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: xc           !< x coordinate of cell centroids
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: yc           !< y coordinate of cell centroids
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: Ri           !< species specific gas constant
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: mi           !< species molar masses
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), SAVE :: p_in       !< initial distribution of primitive variables
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), SAVE :: iden       !< identity matrix
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: vel          !< velocity components (working vector) 
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: yi           !< species mass fractions \f$ y_i \f$ (working vector)
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: xi           !< species mole fractions \f$ x_i \f$ (working vector)
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: rhoi         !< species densities \f$ \rho_i \f$ (working vector)
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: chi          !< species thermal diffusion ratios \f$ \chi_i \f$ (working vector)
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: Di           !< species diffusion coefficients \f$ D_i \f$ (working vector)
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: Dij_mat    !< species diffusion coefficients \f$ D_ij \f$ (working matrix)
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: diff_driv    !< species diffusion driving forces \f$ d_i \f$ (working vector)
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: Ji           !< species diffusion fluxes (working vector)
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: ei           !< species specific energy \f$ e_i \f$ (working vector)
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: hi           !< species specific enthalpy \f$ h_i \f$ (working vector)
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: eiT          !< component of species specific energy in equilibtrium with translation of heavy-particles \f$ {\tilde{e}}_i \f$ (working vector)
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: eistar       !< component of species specific energy in non-equilibtrium with translation of heavy-particles \f$ e^{*}_i = e_i - {\tilde{e}}_i \f$  (working vector) 
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: epsi         !< \f$ \varepsilon_i = R_i T - \alpha e^{*}_i, \quad i \neq e \f$ 
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: sigmai       !< \f$ \sigma_i = e_i - R_i T/\alpha, \quad i \neq e\f$
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: eintk        !< specific non-equilibrium internal energy (working vector)
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: beta         !< beta factor for translational and nonequilibrium internal energy (working vector)
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: betak        !< beta factor for non-equilibrium internal energy \f$ \beta_k \f$ (working vector)
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: ov_betak     !< reciprocal beta factor for non-equilibrium internal \f$ 1/\beta_k \f$ energy (working vector) 
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: temp         !< temperatures (working vector)
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: rho_eint     !< energy densities (working vector)
     REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: lambda_vec   !< thermal conductivity components (working vector)

     REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: log_cfl_number  !< History Log CFL number
     INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: log_cfl_ite          !< History Log CFL iterations increment
     INTEGER, SAVE :: log_cfl_file_lines                              !< History Log CFL lines file
     INTEGER, SAVE :: log_cfl_int
     
     REAL(KIND=8) :: gamma_ser, alpha_ser  !< Used in the SER CFL strategy

     CHARACTER*80,  SAVE :: solver_name           !< name of the Finite volume solver
     CHARACTER*80,  SAVE :: viscosity_law         !< viscosity law
     CHARACTER*80,  SAVE :: library_name          !< name of the library providing thermodynamics, transport properties and source terms 
     CHARACTER*80,  SAVE :: model_name            !< name of the physical model 
     CHARACTER*80,  SAVE :: state_model          !< name of the state model
     CHARACTER*80,  SAVE :: simulation_id        !< name of the simulation
     CHARACTER*80, SAVE :: in_mixture             !< Mixture name  
     CHARACTER*200, SAVE :: library_data_path     !< path to the dir storiong fit coefficients to be used by the library 
     CHARACTER*200, SAVE :: ablation_library_name !< name of the library providing ablation reaction rates 
     CHARACTER*200, SAVE :: cfl_file              !< path to the CFL and VNN number file
     CHARACTER*200, SAVE :: log_cfl_file         !< path to the CFL and VNN number history log file
     LOGICAL, SAVE :: inter_cfl       = .FALSE.   !< flag for indicating if the CFL and VNN numbers are changed interactively
     LOGICAL, SAVE :: log_cfl         = .FALSE.  !< flag for indicating if the CFL and VNN numbers are changed in function of a history log
     LOGICAL, SAVE :: adapt_cfl        = .FALSE.

     LOGICAL, SAVE :: ser_cfl         = .FALse.   !< flag indicating if CFL is controlled by the switched evolution relaxation strategy
     LOGICAL, SAVE :: flag_pres_elec  = .FALSE.   !< flag for indicating if free-electrons are present within the mixture 
     LOGICAL, SAVE :: flag_extra_data = .FALSE.   !< flag for printing extra-data in the solution output file
     LOGICAL, SAVE :: flag_diss       = .FALSE.   !< flag indicating is dissipation phenomena are accounted for
     LOGICAL, SAVE :: flag_post_CR    = .FALSE.   !< flag for post-processsing (only for nonequilibrioum flows)
     LOGICAL, SAVE :: ios_dir         = .TRUE.   !< flag for the input/output directory
     

    REAL(KIND=8), SAVE :: MachIn
    INTEGER, SAVE :: shockIndex
    REAL(KIND=8), SAVE :: shockPosition

  END MODULE mod_general_data
!------------------------------------------------------------------------------!

