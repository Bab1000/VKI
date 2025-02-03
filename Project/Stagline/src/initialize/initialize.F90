!------------------------------------------------------------------------------!
!> This subroutine sets the initial field. 
  SUBROUTINE initialize () 

#include"../config.h"

    USE mod_general_data           
    USE mod_numerics_data
    USE mod_constants
#ifdef GOTO_BLAS
    USE mod_algebra,                  ONLY: initialize_GOTO_BLAS, initialize_Thomas
#else 
    USE mod_algebra,                  ONLY: initialize_Thomas
#endif
    USE mod_neq_function_pointer,     ONLY: library_initialize, library_get_Ri, library_get_el_data, &
                                          & library_number_of_elements       
    USE mod_physical_model
    USE mod_radiation,                ONLY: initialize_radiation,  wall_Rad

    USE mod_general_ablation,         ONLY: mdot_i

    USE mod_surface_properties,       ONLY: rho_i_surface, T_surface

    IMPLICIT NONE
     
    INTEGER :: i, j, k
    INTEGER :: nb_sub
    REAL(KIND=8) :: pos
    CHARACTER*80 :: mixture, reaction, transf
    EXTERNAL :: set_procedure_pointer, set_library_pointer, physical_data, & 
              & initialize_solution, read_solution 

    ! Solver name
    solver_name = 'solver_fvmcc_F90'

    ! Physical model name
    model_name = get_name (physical_model)
  
    ! Flag for dissipation phenomena
    flag_diss = get_diss_flag (physical_model)

    ! Number of conservation equations
    nb_eq = get_nb_eq (physical_model)
  
    ! Number of dimensions
    nb_dim = get_nb_dim (physical_model)
  
    ! Number of species
    nb_ns = get_nb_species (physical_model)
  
    ! Number of vibrational temperatures
    nb_tvib = get_nb_tvib (physical_model)
  
    ! Number of free-electron electronic temperatures
    nb_te = get_nb_te (physical_model)
  
    ! Number of rotational temperatures
    nb_trot = get_nb_trot (physical_model)
  
    ! Number of tempertures
    nb_temp = get_nb_temp (physical_model) 
  
    ! Number of internal temperatures (rotational and vibrational only)
    nb_int_temp = nb_temp - (nb_te + 1)
 
    ! Number of subintervals
    nb_sub = SIZE(p_in)

    MachIn = 0.0D0

    wall_Rad = 0.0D0
    
    ! Vector containing quantities used when computing fluxes
    IF (model_name.EQ.'pg') THEN
 
       ! Specific heat ratio (calorically perfect gas flows)
       gamma = get_gamma(physical_model)

       ! Prandtl number  
       Pr = get_prandtl(physical_model)

       ! Specific gas constant
       r_gas = get_rgas(physical_model)

       ! Constant volume and pressure specific heats
       cv_gas = R_gas/(gamma - 1.d0)
       cp_gas = cv_gas + R_gas

       ! Number of physical quantities to be stored per each cell 
       nb_prop = 0
  
       ! Pointers to physical data
       ! Velocity components
       pos_u_cell = 1
       nb_prop    = nb_prop + 1

       IF (nb_dim.EQ.2) pos_v_cell = pos_u_cell + 1 

       ! Total enthalpy
       pos_h0_cell = pos_u_cell + nb_dim
       nb_prop     = nb_prop + nb_dim

       ! Speed of sound
       pos_c_cell = pos_h0_cell + 1
       nb_prop    = nb_prop + 1
      
       ! Kinetic energy
       pos_ek_cell = pos_c_cell + 1
       nb_prop     = nb_prop + 1

       ! Pressure
       pos_pres_cell = pos_ek_cell + 1
       nb_prop       = nb_prop + 1

       ! Density 
       pos_rho_cell = pos_pres_cell + 1
       nb_prop      = nb_prop + 1

       ! Temperature 
       pos_T_cell = pos_rho_cell + 1
       nb_prop    = nb_prop + 1

       ! Dynamic viscosity and thermal conductivity (only for viscous computations)
       IF (flag_diss.EQV..TRUE.) THEN

          ! Dynamic viscosity 
          pos_mu_cell = pos_T_cell + 1 
          nb_prop     = nb_prop + 1

          ! Bulk viscosity 
          pos_kappa_cell = pos_mu_cell + 1 
          nb_prop        = nb_prop + 1

          ! Thermal conductivity 
          pos_lambda_cell = pos_kappa_cell + 1 
          nb_prop         = nb_prop + 1
 
       ELSE 

          pos_mu_cell = nb_prop 

       ENDIF
       
    ELSE 
 
       ! Specific gas constant (for restarting from calorically perfect gas solutions)
       r_gas = get_rgas(physical_model)
 
       ! Number of physical quantities to be stored per each cell 
       nb_prop = 0

       ! Pointers to physical data
       ! Temperature(s)
       pos_T_cell = 1
       nb_prop    = nb_prop + nb_temp

       ! Velocity component(s) 
       pos_u_cell = pos_T_cell + nb_temp
       nb_prop    = nb_prop + nb_dim

       IF (nb_dim.EQ.2) pos_v_cell = pos_u_cell + 1 

       ! Total enthalpy 
       pos_h0_cell = pos_u_cell + nb_dim
       nb_prop     = nb_prop + 1

       ! Frozen speed of sound 
       pos_c_cell = pos_h0_cell + 1
       nb_prop    = nb_prop + 1

       ! Frozen specific heat ratio
       pos_gamma_cell = pos_c_cell + 1
       nb_prop        = nb_prop + 1

       ! Pressure 
       pos_pres_cell = pos_gamma_cell + 1
       nb_prop       = nb_prop + 1

       ! Density
       pos_rho_cell = pos_pres_cell + 1
       nb_prop      = nb_prop + 1

       ! Kinetic energy 
       pos_ek_cell = pos_rho_cell + 1
       nb_prop     = nb_prop + 1

       ! Alpha factor  
       pos_alpha_cell = pos_ek_cell + 1
       nb_prop        = nb_prop + 1
       IF (nb_te.EQ.1) nb_prop = nb_prop+1

       ! Beta factor 
       pos_beta_cell = pos_alpha_cell + 1
       nb_prop       = nb_prop + 1

       ! Betak factor(s)
       IF ((nb_int_temp+nb_te).GE.1) THEN

          pos_betak_cell = pos_beta_cell + 1
          nb_prop        = nb_prop + nb_int_temp + nb_te

       ELSE 

          pos_betak_cell = pos_beta_cell

       ENDIF

       ! Species internal energy (total value)
       pos_ei_cell  = pos_betak_cell + MAX(1,nb_int_temp+nb_te)
       nb_prop      = nb_prop + nb_ns

       ! Species internal energy (nonequilibrium components)
       IF ((nb_int_temp+nb_te).GE.1) THEN

          pos_eki_cell = pos_ei_cell + nb_ns
          nb_prop      = nb_prop + nb_ns*(nb_int_temp+nb_te)

       ELSE 

          pos_eki_cell = pos_ei_cell 

       ENDIF
  
       ! Species diffusion coefficients, thermal diffusion ratios, mixture dynamic viscosity and thermal conductivity components 
       ! (translational and internal components)
       IF (flag_diss.EQV..TRUE.) THEN

          ! Species diffusion coefficients
          IF ((nb_int_temp+nb_te).GE.1) THEN
             pos_Di_cell = pos_eki_cell + nb_ns*(nb_int_temp+nb_te)
          ELSE
             pos_Di_cell = pos_ei_cell + nb_ns
          ENDIF
          nb_prop = nb_prop + nb_ns
          
          ! Species thermal diffusion ratios
          pos_chi_cell = pos_Di_cell + nb_ns
          nb_prop      = nb_prop + nb_ns

          ! Mixture dynamic viscosity 
          pos_mu_cell = pos_chi_cell + nb_ns 
          nb_prop     = nb_prop + 1

          ! Mixture bulk viscosity 
          pos_kappa_cell = pos_mu_cell + 1 
          nb_prop        = nb_prop + 1

          ! Thermal conductivity (translational and internal components)
          pos_lambda_cell = pos_kappa_cell + 1 
          nb_prop         = nb_prop + nb_temp
         
       ELSE 

          pos_mu_cell = nb_prop

       ENDIF

    ENDIF
    
    ! Useful indices
    ! Starting and finishing location of conservative variables and physical properties for physical cells
    start_u_phys    = 2*nb_eq
    start_prop_phys = 2*nb_prop

    finish_u_phys    = (nb_cells + 2)*nb_eq
    finish_prop_phys = (nb_cells + 2)*nb_prop

    ! Primitive variable indices
    pos_rho = 1
    pos_u   = nb_ns + 1
    pos_v   = nb_ns + nb_dim 
    pos_p   = nb_ns + nb_dim + 1
    pos_T   = pos_v + 1
    pos_Tk  = pos_T
    pos_Te  = pos_T
    IF (nb_int_temp.GT.0) pos_Tk = pos_Tk + 1
    IF (nb_te.EQ.1)       pos_Te = pos_Tk + MAX(1,nb_int_temp)

    ! Conservative variable indices
    pos_rhou  = pos_u 
    pos_rhov  = pos_v
    pos_rhoE  = pos_T
    pos_rhoek = pos_rhoE 
    pos_rhose = pos_rhoE
    IF ((nb_int_temp+nb_te).GT.0) pos_rhoek = pos_rhoek + 1
    IF (nb_te.EQ.1)       pos_rhose = pos_rhoek + MAX(0,nb_int_temp)
 
    ! Allocation and initialization of solution common vectors
    ALLOCATE(u(nb_eq*nb_tot_cells), p(nb_eq*nb_tot_cells), du(nb_eq*nb_tot_cells), &
            & dres((nb_eq*nb_tot_cells)))
    ALLOCATE(cell_prop(nb_prop*nb_tot_cells))
   
    IF (nb_stages.GT.1) ALLOCATE(u_old(nb_eq*nb_tot_cells))

    ! Identity matrix 
    ALLOCATE(iden(nb_eq,nb_eq))
    iden = 0.d0
    DO i = 1,nb_eq 
       iden(i,i) = 1.d0
    ENDDO

    ! Eigensystem and boundary condition data
    ALLOCATE(deltau(nb_eq), lambda(nb_eq), diss(nb_eq), u_vec(nb_eq), phys_prop(nb_prop))
    ALLOCATE(right_eig(nb_eq,nb_eq), left_eig(nb_eq,nb_eq), absA_eig(nb_eq,nb_eq))
    ALLOCATE(epsi(nb_ns), sigmai(nb_ns), rhoi(nb_ns), yi(nb_ns), xi(nb_ns))
    ALLOCATE(ei(nb_ns*nb_temp))
    ALLOCATE(hi(nb_ns*nb_temp))
    ALLOCATE(eiT(nb_ns),eistar(nb_ns))
    ALLOCATE(beta(nb_temp))
    ALLOCATE(rho_eint(nb_temp), temp(nb_temp))
    ALLOCATE(eintk(MAX(1,nb_int_temp+nb_te)),betak(MAX(1,nb_int_temp+nb_te)),ov_betak(MAX(1,nb_int_temp+nb_te)))
    ALLOCATE(vel(nb_dim))
   
    ! Ablation properties
    ALLOCATE(mdot_i(nb_ns))
    mdot_i = 0.0D0

    ! Surface properties
    ALLOCATE(rho_i_surface(nb_ns))
    ALLOCATE(T_surface(nb_temp))

    ! Physical properties and conservative variables 
    ALLOCATE(u_left(nb_eq), u_right(nb_eq))
    ALLOCATE(prim_left(nb_eq), prim_right(nb_eq))
    ALLOCATE(prop_left(nb_prop), prop_right(nb_prop))
 
    ! Additional working vectors for nonequilibrium inviscid and viscous flows
    IF (model_name.EQ.'neq') THEN
       ALLOCATE(rhoil(nb_ns), rhoir(nb_ns))
       ALLOCATE(xil(nb_ns), xir(nb_ns))
       IF (flag_diss.EQV..TRUE.) THEN
          ALLOCATE(Di(nb_ns),chi(nb_ns))
          ALLOCATE(Dij_mat(nb_ns,nb_ns))
          ALLOCATE(Ji(nb_dim*nb_ns))
          ALLOCATE(diff_driv(nb_dim*nb_ns))
          ALLOCATE(lambda_vec(nb_temp))
       ENDIF
    ENDIF 

    ! Data for time-intergarion
    ALLOCATE(fc(nb_eq), fd(nb_eq), fl(nb_eq), fr(nb_eq))  
    ALLOCATE(s(nb_eq), s_vec(nb_eq))
    ALLOCATE(js_mat(nb_eq,nb_eq))
    ALLOCATE(js(nb_eq,nb_eq), inv_js(nb_eq,nb_eq))
    ALLOCATE(aux(nb_eq))
    ALLOCATE(ull(nb_eq), ul(nb_eq), ur(nb_eq), urr(nb_eq))
    ALLOCATE(prop_ll(nb_prop), prop_l(nb_prop), prop_r(nb_prop), prop_rr(nb_prop))
    ALLOCATE(jsl(nb_eq,nb_eq),jsr(nb_eq,nb_eq))

    ! Source term Jacobian (stored in a one dimensional vector)
    ALLOCATE(js_line((nb_ns + nb_temp)*(nb_ns + nb_temp - 1)))

    ! Pre-allocation of lhs matrices for the solution of block tridiagonal system (fully implicit scheme)
    IF (flag_full_Impl.EQV..TRUE.) THEN

        CALL initialize_Thomas(nb_eq, nb_cells)

        ALLOCATE(ru(nb_eq*nb_tot_cells)) 
        ALLOCATE(left_bound(nb_eq,nb_eq), right_bound(nb_eq,nb_eq))
        ALLOCATE(central(nb_eq,nb_eq))
        ALLOCATE(jfl(nb_eq,nb_eq), jfr(nb_eq,nb_eq))
        ALLOCATE(jfdl(nb_eq,nb_eq), jfdr(nb_eq,nb_eq))
        ALLOCATE(j_ghost(nb_eq,nb_eq))
        IF ((flag_stag_line.EQV..TRUE.).AND.(flag_diss.EQV..TRUE.)) THEN
           ALLOCATE(Ad(nb_eq,nb_eq),Bd(nb_eq,nb_eq))
        ENDIF

    ENDIF

    ! Data for the GOTO_BLAS library
#ifdef GOTO_BLAS
    CALL initialize_GOTO_BLAS (nb_eq)
#endif

    ! Data for Roe's average state
    IF ((flux_splitter.EQ.'roe').OR.(flux_splitter.EQ.'hlle').OR.& 
      & (flux_splitter.EQ.'Roe').OR.(flux_splitter.EQ.'HLLE').OR.&
      & (flux_splitter.EQ.'modified_SW').OR.(flux_splitter.EQ.'MODIFIED_SW')) THEN
       ALLOCATE(cons_Roe(nb_eq))
       ALLOCATE(phys_data_Roe(nb_prop))
       cons_Roe      = 0.d0
       phys_data_Roe = 0.d0
    ENDIF

    ! Initialization
    u = 0.d0
    cell_prop = 0.d0

    ! Thermodynamic library inizialization 
    IF (model_name.EQ.'neq') THEN 

       ! Set library pointers
       CALL set_library_pointer()

       mixture  = get_mixture  (physical_model)
       reaction = get_reaction (physical_model)
       transf   = get_transfer (physical_model)

       ! Mixture, reaction and transfer input files specified for the library
       CALL library_initialize (nb_ns, nb_trot, nb_tvib, nb_te, nb_temp, nb_eq, nb_dim,     & 
                              & solver_name, mixture, reaction, state_model, transf,    &
                              & library_data_path, xi_tol)

       ALLOCATE(Ri(nb_ns),mi(nb_ns))
       CALL library_get_Ri (Ri)
       
       nb_ele = 0
       CALL library_number_of_elements(nb_ele)
   
       ! Species molar masses
       DO i = 1,nb_ns 
          mi(i) = urg/Ri(i)
       ENDDO

       ! Default position of free electron species (if present)
       pos_em = 1
       IF (flag_pres_elec.EQV..TRUE.) CALL library_get_el_data (pos_em, Re)

    ENDIF

    ! Pre-processor. In this stage subroutines and/or functions pointers are initialized. This step is 
    ! performed in order not to use SELECT CASE statements that can slow down program execution.  
    CALL set_procedure_pointer ()

    ! Solution initialization
    IF ((flag_restart.EQV..FALSE.).AND.(flag_restart_pg.EQV..FALSE.).AND.(flag_restart_eq.EQV..FALSE.)) THEN 
      
       SELECT CASE(in_dim_id)
         
         ! Variation over x (1D or 1D stagnation-line flows)
         CASE(1)

           IF (delta(1).GT.xc(1)) delta(1) = xc(1)

           DO i = 1,nb_tot_cells
  
              DO j = 1,SIZE(delta,1)
                 IF (xc(i).GE.delta(j)) THEN
                    k = j
                 ENDIF
              ENDDO
  
              ! Initialization of primitive variables
              DO j = 1,nb_eq
                 p((i - 1)*nb_eq + j) = p_in(j,k)
              ENDDO

           ENDDO
 
         ! Variation over y (or r) (2D flows only)
         CASE(2)
           
            IF (delta(1).GT.yc(1)) delta(1) = yc(1)

            DO i = 1,nb_tot_cells
  
               DO j = 1,SIZE(delta,1)
                  IF (yc(i).GE.delta(j)) THEN
                     k = j
                  ENDIF
               ENDDO
  
               ! Initialization of primitive variables
               DO j = 1,nb_eq
                  p((i - 1)*nb_eq + j) = p_in(j,k)
               ENDDO

            ENDDO

       END SELECT

    ! Restart from a previous solution
    ELSEIF ((flag_restart.EQV..TRUE.).OR.(flag_restart_pg.EQV..TRUE.).OR.(flag_restart_eq.EQV..TRUE.)) THEN

       ! Read primitive variables from a previsous solution
       CALL read_solution()

    ELSE 

       WRITE(*,10)'In "initialize.F90", error in setting the flag_restart option ...'
       PRINT*

    ENDIF

    CALL initialize_radiation()

    ! Initialize conservative variables and physical properties of all cells
    CALL initialize_solution()
    CALL physical_data()
    CALL computeMachInlet()
   
    ! De-allocation of no longer needed vectors (when solution is started from scratch)
    IF (ALLOCATED(delta)) DEALLOCATE(delta)
    IF (ALLOCATED(p_in))  DEALLOCATE(p_in)
    IF (ALLOCATED(p))     DEALLOCATE(p)
    
10 FORMAT(A)

  END SUBROUTINE initialize 
!------------------------------------------------------------------------------!
  SUBROUTINE computeMachInlet ()
  
  USE mod_general_data,          ONLY: nb_ns, nb_eq, pos_u, pos_T, MachIn, nb_temp
  USE mod_domain_boundary,       ONLY: boundary, boundary_data, get_boundary_inlet
  USE mod_neq_function_pointer,  ONLY: library_get_frozen_sound_speed
  
 
  REAL(KIND=8), DIMENSION(nb_eq) :: inlet_data, prim
  TYPE(boundary) :: bound
  REAL(KIND=8) :: c

  INTEGER :: id
    
  id = 2 
 
  bound      = boundary_data(id)
  inlet_data = get_boundary_inlet(nb_eq, bound)

  CALL library_get_frozen_sound_speed (inlet_data(1:nb_ns), inlet_data(pos_T:pos_T + nb_temp -1 ), c)

  MachIn = abs(inlet_data(pos_u))/c


  END SUBROUTINE computeMachInlet
