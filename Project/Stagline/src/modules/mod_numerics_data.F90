!------------------------------------------------------------------------------!
!> This module stores variables representing options for the discretization of 
!! the governing equations.
!! Working vectors used during computations are also allocated and stored in order to save 
!! computational time. 
  MODULE mod_numerics_data

    IMPLICIT NONE

    INTEGER, SAVE :: nb_source_terms = 0                          !< number of source terms
    INTEGER, SAVE :: nb_inv_source_terms = 0                      !< number of inviscid source terms (used only for stagnation-line flows)
    INTEGER, SAVE :: nb_itmax  = 0                                !< maximum number of iterations
    INTEGER, SAVE :: stop_id   = 0                                !< option for stopping computation
    INTEGER, SAVE :: nb_stages = 1                                !< number of stages of the time-integration method (default is 1)
    INTEGER, SAVE :: stag_line_geom = 0                           !< stagnation-line flow geometry (0 = cylinder, 1 = sphere)
    REAL(KIND=8), SAVE :: cfl = 0.d0                              !< CFL (Courant-Friedrich-Lewi) number
    REAL(KIND=8), SAVE :: vnn = 0.d0                              !< Von-Neumann number 
    REAL(KIND=8), SAVE :: time_step = 0.d0                        !< time-step value
    REAL(KIND=8), SAVE :: time_limit = 0.d0                       !< time-level final value (for time-accurate calculations)
    REAL(KIND=8), SAVE :: lim_res    = 0.d0                       !< lower bound for the residual for convergence testing
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: stage_coeff  !< stage coefficients of the time-integration scheme
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: fc                 !< convective flux
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: fd                 !< diffusive flux
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: fl                 !< left component of numerical flux in fvs schemes
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: fr                 !< right component of numerical flux in fvs schemes 
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: u_left             !< conservative variables of reconstructed left state
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: u_right            !< conservative variables of reconstructed right state 
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ull                !< conservative variables of left-left state
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ul                 !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ur                 !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: urr                !< conservative variables of right-right state
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: prim_left          !< primitive variables of left state 
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: prim_right         !< primitive variables of right state
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: prop_left          !< physical properties of reconstructed left state
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: prop_right         !< physical properties of reconstructed right state
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: prop_ll            !< physical properties of left-left state 
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: prop_l             !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: prop_r             !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: prop_rr            !< physical properties of right-right state
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: s                  !< source term 
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: aux                !< auxiliaty vector 
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: deltau             !< conservative variable jump
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: lambda             !< eigenvalue vector 
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: diss               !< dissipation term in Roe's approximate Riemann solver
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: rhoil              !< left state species densities
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: rhoir              !< right state species densities
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: xil                !< left state species mole fractions
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: xir                !< right state species mole fractions
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: cons_Roe           !< conservarive properties of Roe's averaged state
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: phys_data_Roe      !< physical properties of Roe's averaged state
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: js_line            !< array for row storage of source term Jacobiab
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: left_eig         !< left eigenvector matrix
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: right_eig        !< right eigenvector matrix  
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: absA_eig         !< absolute eigenvalue diagonal matrix
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: js               !< source term Jacobian
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: inv_js           !< inverse of the source term Jacobian
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: jsl              !< left state source term Jacobian 
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: jsr              !< right state source term Jacobian 
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: jfl              !< left state inviscid flux Jacobian             
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: jfr              !< right state inviscid flux Jacobian 
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: jfdl             !< left state diffusive flux Jacobian
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: jfdr             !< right state diffusive flux Jacobian
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: j_ghost          !< Jacobian for implicit boundary condition
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: central          !< auxiliary matrix 
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: left_bound       !< left ghost cell Jacobian
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: right_bound      !< right ghost cell Jacobian
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Ad               !< working matrix for diffusive flux Jacobian (stagnation line flows)
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Bd               !< working matrix for diffusive flux Jacobian (stagnation line flows)
    CHARACTER*80, SAVE :: flux_splitter    !< numerical flux
    CHARACTER*80, SAVE :: time_disc        !< time-discretization scheme
    CHARACTER*80, SAVE :: poly_rec         !< polynomial re-construction
    CHARACTER*80, SAVE :: fun_limiter      !< slope limiter function
    CHARACTER*80, SAVE :: flux_Jac_approx  !< approximation for numerical flux Jacobian 
    CHARACTER*80, SAVE :: diff_flux_Jac    !< variable specifying how to compute the diffusive flux Jacobian ("analytical" or "numerical")
    CHARACTER*80, ALLOCATABLE, DIMENSION(:), SAVE :: source_terms             !< source term names
    CHARACTER*80, ALLOCATABLE, DIMENSION(:), SAVE :: compute_source_term_Jac  !< option for computing the source term  Jacobian ("analytical" of "numerical")  
    CHARACTER*80, DIMENSION(2), SAVE :: bc_Jac                                 !< variable specifying how to compute the boundary conditions Jacobian ("analytical" or "numerical")
    LOGICAL, SAVE :: flag_stag_line = .FALSE.   !< flag indicating if the computation is for stagnation-line flows
    LOGICAL, SAVE :: flag_restart = .FALSE.     !< flag for solution restart from a previous solution 
    LOGICAL, SAVE :: flag_restart_pg = .FALSE.  !< flag for solution restart from a calorically perfect gas solution
    LOGICAL, SAVE :: flag_restart_eq = .FALSE.  !< flag for solution restart from an equilibrium solution 
    LOGICAL, SAVE :: flag_Jac     = .TRUE.      !< flag indicating if Jacobians need to be computed 
    LOGICAL, SAVE :: flag_full_Impl  = .FALSE.  !< flag indicating if the time-integration method is fully implicit
    LOGICAL, SAVE :: flag_metrics    = .FALSE.  !< flag indicating if metrics is used for the flux linear (2nd order) reconstruction (A. Turchi)

  END MODULE mod_numerics_data
!------------------------------------------------------------------------------!

