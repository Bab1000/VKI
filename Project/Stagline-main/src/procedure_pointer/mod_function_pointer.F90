!------------------------------------------------------------------------------!
!> This module provides the support for function pointers to be used in computations of 
!! calorically perfect gas or nonequilibrium flows.
  MODULE mod_function_pointer

    IMPLICIT NONE

    ! Abstract interface definitions 
    ABSTRACT INTERFACE

      ! Interface for subroutine computing conservative variables from primitive variables 
      SUBROUTINE prim_to_cons (prim, cons)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prim
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cons
      END SUBROUTINE prim_to_cons 
 
      ! Interface for subroutine computing primitive variables from conservative variables  
      SUBROUTINE cons_to_prim (cons, prim)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: prim
      END SUBROUTINE cons_to_prim

      ! Interface for subroutine computing physical properties from conservative variables
      SUBROUTINE cons_to_phys (cons, phys)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: phys 
      END SUBROUTINE 

      ! Interface for subroutine computing physical properties and conservative variables from 
      ! primitive variables
      SUBROUTINE prim_to_cons_phys(prim, cons, phys)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prim
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cons, phys
      END SUBROUTINE prim_to_cons_phys

      ! Interface for time-step for 1D flows
      SUBROUTINE time_step_1D (vol, cell_data, dt)
        REAL(KIND=8), INTENT(IN) :: vol
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cell_data
        REAL(KIND=8), INTENT(OUT) :: dt
      END SUBROUTINE time_step_1D

      ! Interface for time-step for 1D stagnation line flows
      SUBROUTINE time_step_1D_SL (vol, cell_data, dt)
        REAL(KIND=8), INTENT(IN) :: vol
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cell_data
        REAL(KIND=8), INTENT(OUT) :: dt
      END SUBROUTINE time_step_1D_SL

      ! Interface for 1D numerical flux 
      SUBROUTINE flux_splitter_1D (nx, left_data, right_data, cons_l, cons_r, fc)
        REAL(KIND=8), INTENT(IN) :: nx
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data, right_data, cons_l, cons_r
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: fc
      END SUBROUTINE flux_splitter_1D 

      ! Interface for 1D numerical flux Jacobian 
      SUBROUTINE flux_splitter_Jac_1D (nx, left_data, right_data, cons_l, cons_r, jfcl, jfcr)
        REAL(KIND=8), INTENT(IN) :: nx
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data, right_data, cons_l, cons_r
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfcl, jfcr
      END SUBROUTINE flux_splitter_Jac_1D

      ! Interface for 1D stagnation line numerical flux 
      SUBROUTINE flux_splitter_1D_SL (nx, vol_l, vol_r, left_data, right_data, cons_l, cons_r, fc)
        REAL(KIND=8), INTENT(IN) :: nx
        REAL(KIND=8), INTENT(IN) :: vol_l, vol_r
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data, right_data, cons_l, cons_r
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: fc
      END SUBROUTINE flux_splitter_1D_SL 

      ! Interface for 1D stagnation line numerical flux Jacobian 
      SUBROUTINE flux_splitter_Jac_1D_SL (nx, vol_l, vol_r, left_data, right_data, cons_l, cons_r, jfcl, jfcr)
        REAL(KIND=8), INTENT(IN) :: nx
        REAL(KIND=8), INTENT(IN) :: vol_l, vol_r
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data, right_data, cons_l, cons_r
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfcl, jfcr
      END SUBROUTINE flux_splitter_Jac_1D_SL  

      ! Interface for 1D inviscid flux Jacobian (nonequilibrium flow only)
      SUBROUTINE flux_Jacobian_neq_1D (nx, cons, phys_data, a)
        REAL(KIND=8), INTENT(IN) :: nx
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons, phys_data
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: a
      END SUBROUTINE flux_Jacobian_neq_1D

      ! Interface for 1D stagnataion line inviscid flux Jacobian (nonequilibrium flow only)
      SUBROUTINE flux_Jacobian_neq_1D_SL (nx, cons, phys_data, a)
        REAL(KIND=8), INTENT(IN) :: nx
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons, phys_data
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: a
      END SUBROUTINE flux_Jacobian_neq_1D_SL

      ! Interface for 1D eigensystem (nonequilibrium flow only) 
      SUBROUTINE eigensystem_neq_1D (nx, cons, phys_data, lambda, right, left)
        REAL(KIND=8), INTENT(IN) :: nx
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons, phys_data
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: lambda
        REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: right, left
      END SUBROUTINE eigensystem_neq_1D  

      ! Interface for 1D stagnation line eigensystem (nonequilibrium flow only) 
      SUBROUTINE eigensystem_neq_1D_SL (nx, cons, phys_data, lambda, right, left)
        REAL(KIND=8), INTENT(IN) :: nx
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons, phys_data
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: lambda
        REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: right, left
      END SUBROUTINE eigensystem_neq_1D_SL  

      ! Interface for 1D Roe average state (nonequilibrium flow only)
      SUBROUTINE roe_avgState_neq_1D (u_left, u_right, left_data, right_data, u_Roe, data_Roe)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left, u_right, left_data, right_data
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_Roe, data_Roe
      END SUBROUTINE roe_avgState_neq_1D 

      ! Interface for 1D diffusive flux 
      SUBROUTINE diffusive_flux_1D (voll, volr, left_data, right_data, u_left, u_right, fd)
        REAL(KIND=8), INTENT(IN) :: voll, volr
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left, u_right
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data, right_data
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: fd
      END SUBROUTINE diffusive_flux_1D 

      ! Interface for numerical 1D diffusive flux and diffusive flux Jacobians
      SUBROUTINE diffusive_flux_1D_Jac (voll, volr, left_data, right_data, u_left, u_right, fd, jfdl, jfdr)
        REAL(KIND=8), INTENT(IN) :: voll, volr
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left, u_right
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data, right_data
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: fd
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfdl, jfdr
      END SUBROUTINE diffusive_flux_1D_Jac 

      ! Interface for 1D stagnation line diffusive flux 
      SUBROUTINE diffusive_flux_1D_SL (rl, rr, voll, volr, left_data, right_data, u_left, u_right, fd)
        REAL(KIND=8), INTENT(IN) :: rl, rr
        REAL(KIND=8), INTENT(IN) :: voll, volr
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left, u_right
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data, right_data
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: fd
      END SUBROUTINE diffusive_flux_1D_SL 

      ! Interface for 1D stagnation line diffusive flux 
      SUBROUTINE diffusive_flux_1D_SL_Jac (rl, rr, voll, volr, left_data, right_data, u_left, u_right, fd, jfdl, jfdr)
        REAL(KIND=8), INTENT(IN) :: rl, rr
        REAL(KIND=8), INTENT(IN) :: voll, volr
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left, u_right
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data, right_data
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: fd
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfdl, jfdr
      END SUBROUTINE diffusive_flux_1D_SL_Jac       

      ! Interface for 1D spectral radius 
      SUBROUTINE spectral_radius_1D (phys_data, sp)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data
        REAL(KIND=8), INTENT(OUT) :: sp 
      END SUBROUTINE spectral_radius_1D

      ! Interface for 1D stagnation line spectral radius 
      SUBROUTINE spectral_radius_1D_SL (phys_data, sp)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data
        REAL(KIND=8), INTENT(OUT) :: sp 
      END SUBROUTINE spectral_radius_1D_SL

      ! Interface for 1D pre-conditioning velocity
      SUBROUTINE prec_vel_1D (vol, rho, mu, V, delta_p, c, Vp)
        REAL(KIND=8), INTENT(IN) :: vol, rho, mu, V, delta_p, c
        REAL(KIND=8), INTENT(OUT) :: Vp
      END SUBROUTINE prec_vel_1D 

      ! Interface for 1D stagnation line pre-conditioning velocity
      SUBROUTINE prec_vel_1D_SL (vol, rho, mu, V, delta_p, c, Vp)
        REAL(KIND=8), INTENT(IN) :: vol, rho, mu, V, delta_p, c
        REAL(KIND=8), INTENT(OUT) :: Vp
      END SUBROUTINE prec_vel_1D_SL 
      
      ! Interface for 1D stagnation line stress tensor 
      SUBROUTINE stress_tensor_1D_SL (mu, r, u, v, du_dr, dv_dr, tau_rr, tau_rt, tau_tt)
        REAL(KIND=8), INTENT(IN) :: mu, r, u, v, du_dr, dv_dr
        REAL(KIND=8), INTENT(OUT) :: tau_rr, tau_rt, tau_tt
      END SUBROUTINE stress_tensor_1D_SL

      ! Interface for computation of tranport properties (for calorically perfect gas flow only) 
      SUBROUTINE transp_Coeff (T, mu, lambda) 
        REAL(KIND=8), INTENT(IN) :: T
        REAL(KIND=8), INTENT(OUT) :: mu, lambda
      END SUBROUTINE transp_Coeff 

      ! Interface for 1D polynomial reconstruction
      SUBROUTINE muscl_1D (ull, ul, ur, urr, prop_ll, prop_l, prop_r, prop_rr, u_left, u_right, prop_left, prop_right)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: ull, ul, ur, urr,  prop_ll, prop_l, prop_r, prop_rr
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) ::  u_left, u_right, prop_left, prop_right
      END SUBROUTINE muscl_1D 

      ! Interface for 1D stagnation line polynomial reconstruction modified to consider also the case w/ metrics (A. Turchi)
      SUBROUTINE muscl_1D_SL (vol_ll, vol_l, vol_r, vol_rr, ull, ul, ur, urr, prop_ll, prop_l, prop_r, prop_rr, u_left,&
                                & u_right, prop_left, prop_right)
        REAL(KIND=8), INTENT(IN) :: vol_ll, vol_l, vol_r, vol_rr
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: ull, ul, ur, urr,  prop_ll, prop_l, prop_r, prop_rr
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) ::  u_left, u_right, prop_left, prop_right
      END SUBROUTINE muscl_1D_SL 

      ! Interface for 1D boundary condition (explicit or point implicit time integration scheme - 1st order accuracy in space)
      SUBROUTINE bc_expl_1D_1st (bound_id, phys_data, u_phys, ghost_data, u_ghost)
        INTEGER, INTENT(IN) :: bound_id
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys, phys_data
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost, ghost_data
      END SUBROUTINE bc_expl_1D_1st 

      ! Interface for 1D boundary condition (explicit or point implicit time integration scheme - 2nd order accuracy in space)
      SUBROUTINE bc_expl_1D_2nd (bound_id, phys_data1, u_phys1, phys_data2, u_phys2, ghost_data1, u_ghost1, ghost_data2, u_ghost2)
        INTEGER, INTENT(IN) :: bound_id
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys1, phys_data1, u_phys2, phys_data2 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost1, ghost_data1, u_ghost2, ghost_data2
      END SUBROUTINE bc_expl_1D_2nd 

      ! Interface for 1D boundary condition (fully implicit time integration scheme - 1st order accuracy in space)
      SUBROUTINE bc_impl_1D_1st (bound_id, phys_data, u_phys, ghost_data, u_ghost, ju_ghost)
        INTEGER, INTENT(IN) :: bound_id
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys, phys_data
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost, ghost_data
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: ju_ghost
      END SUBROUTINE bc_impl_1D_1st

      ! Interface for 1D boundary condition (fully implicit time integration scheme - 2nd order accuracy in space)
      SUBROUTINE bc_impl_1D_2nd (bound_id, phys_data1, u_phys1, phys_data2, u_phys2, & 
               &                 ghost_data1, u_ghost1, ghost_data2, u_ghost2, ju_ghost)
        INTEGER, INTENT(IN) :: bound_id
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys1, phys_data1, u_phys2, phys_data2
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost1, ghost_data1, u_ghost2, ghost_data2
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: ju_ghost
      END SUBROUTINE bc_impl_1D_2nd

      ! Interface for 1D stagnation line boundary condition (explicit or point implicit time integration scheme - 
      ! 1st order accuracy in space)
      SUBROUTINE bc_expl_1D_SL_1st (bound_id, phys_data, u_phys, ghost_data, u_ghost)
        INTEGER, INTENT(IN) :: bound_id
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys, phys_data
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost, ghost_data
      END SUBROUTINE bc_expl_1D_SL_1st 

      ! Interface for 1D stagtaion line boundary condition (explicit or point implicit time integration scheme - 
      ! 2nd order accuracy in space)
      SUBROUTINE bc_expl_1D_SL_2nd (bound_id, phys_data1, u_phys1, phys_data2, u_phys2, ghost_data1, u_ghost1, & 
                                  & ghost_data2, u_ghost2)
        INTEGER, INTENT(IN) :: bound_id
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys1, phys_data1, u_phys2, phys_data2 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost1, ghost_data1, u_ghost2, ghost_data2
      END SUBROUTINE bc_expl_1D_SL_2nd 

      ! Interface for 1D stagnation line boundary condition (fully implicit time integration scheme - 1st order accuracy in space)
      SUBROUTINE bc_impl_1D_SL_1st (bound_id, phys_data, u_phys, ghost_data, u_ghost, ju_ghost)
        INTEGER, INTENT(IN) :: bound_id
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys, phys_data
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost, ghost_data
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: ju_ghost
      END SUBROUTINE bc_impl_1D_SL_1st

      ! Interface for 1D stagnation line boundary condition (fully implicit time integration scheme - 2nd order accuracy in space)
      SUBROUTINE bc_impl_1D_SL_2nd (bound_id, phys_data1, u_phys1, phys_data2, u_phys2, & 
                                 &  ghost_data1, u_ghost1, ghost_data2, u_ghost2, ju_ghost)
        INTEGER, INTENT(IN) :: bound_id
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys1, phys_data1, u_phys2, phys_data2
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost1, ghost_data1, u_ghost2, ghost_data2
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: ju_ghost
      END SUBROUTINE bc_impl_1D_SL_2nd 

      ! Interface for source term (no Jacobian provided in output).
      SUBROUTINE source_term (cell_id, prop, cons, s)
        INTEGER, INTENT(IN) :: cell_id
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prop, cons 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s
      END SUBROUTINE source_term

      ! Interface for source term (Jacobian provided in output).
      SUBROUTINE source_term_Jac (source_id, cell_id, prop, cons, s, js)
        INTEGER, INTENT(IN) :: source_id, cell_id
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prop, cons 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: js 
      END SUBROUTINE source_term_Jac

      ! Interface for invisicd source term for stagnation line flows (no Jacobian provided in output).
      SUBROUTINE inv_source_term_SL (cell_id, r_cell, prop_cell, cons_cell, s)
        INTEGER, INTENT(IN) :: cell_id
        REAL(KIND=8), INTENT(IN) :: r_cell
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prop_cell, cons_cell 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s
      END SUBROUTINE inv_source_term_SL

      ! Interface for invisicd source term for stagnation line flows (Jacobian provided in output).
      SUBROUTINE inv_source_term_Jac_SL (cell_id, source_id, r_cell, prop_cell, cons_cell, s, js)
        INTEGER, INTENT(IN) :: cell_id, source_id
        REAL(KIND=8), INTENT(IN) :: r_cell
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prop_cell, cons_cell 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: js 
      END SUBROUTINE inv_source_term_Jac_SL

      ! Interface for diffusive source term for stagnation line flows (no Jacobian provided in output).
      SUBROUTINE diff_source_term_SL (r_int, vol_left, vol_cell, vol_right, prop_left, prop_cell, prop_right, & 
                                    & u_left, u_cell, u_right, s)
        REAL(KIND=8), INTENT(IN) :: r_int
        REAL(KIND=8), INTENT(IN) :: vol_left, vol_cell, vol_right
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left, u_cell, u_right
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prop_left, prop_cell, prop_right
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s
      END SUBROUTINE diff_source_term_SL

      ! Interface for diffusive source term and related Jacobians for 1D stagnation line flows (Jacobians provided in output).
      SUBROUTINE diff_source_term_SL_Jac (r_int, vol_left, vol_cell, vol_right, prop_left, prop_cell, prop_right, & 
                                        & u_left, u_cell, u_right, s, js_left, js_cell, js_right)
        REAL(KIND=8), INTENT(IN) :: r_int
        REAL(KIND=8), INTENT(IN) :: vol_left, vol_cell, vol_right
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left, u_cell, u_right
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prop_left, prop_cell, prop_right
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: js_left, js_cell, js_right 
      END SUBROUTINE diff_source_term_SL_Jac

      ! Interface for physical data 
      SUBROUTINE phys_data () 
      END SUBROUTINE phys_data

      ! Interface for single-step time integration  
      SUBROUTINE time_integration()
      END SUBROUTINE time_integration 

      ! Interface for limiter function
      FUNCTION muscl_limiter (r)
        REAL(KIND=8), INTENT(IN) :: r
        REAL(KIND=8) :: muscl_limiter
      END FUNCTION muscl_limiter
 
    END INTERFACE

    ! Type definition (pointers to subroutines for source terms - no Jacobian in output)
    TYPE source
         PROCEDURE(source_term), POINTER, NOPASS :: source 
    END TYPE source 

    ! Type definition (pointers to subroutines for source terms - Jacobian in output)
    TYPE source_Jac
         PROCEDURE(source_term_Jac), POINTER, NOPASS :: source_Jac 
    END TYPE source_Jac

    ! Type definition (pointers to subroutines for source terms - no Jacobian in output)
    TYPE inv_source_SL
         PROCEDURE(inv_source_term_SL), POINTER, NOPASS :: source 
    END TYPE inv_source_SL 

    ! Type definition (pointers to subroutines for source terms - no Jacobian in output)
    TYPE inv_source_Jac_SL
         PROCEDURE(inv_source_term_Jac_SL), POINTER, NOPASS :: source_Jac 
    END TYPE inv_source_Jac_SL 

    ! Type definition (pointers to subroutines for boundary condition - 
    ! (explicit or point implicit scheme - 1st order accuracy in space))
    TYPE ghost_state_Expl_1D_1st
         PROCEDURE(bc_expl_1D_1st), POINTER, NOPASS :: bc_procedure 
    END TYPE ghost_state_Expl_1D_1st  

    ! Type definition (pointers to subroutines for boundary condition - 
    ! (explicit or point implicit scheme - 2nd order accuracy in space))
    TYPE ghost_state_Expl_1D_2nd
         PROCEDURE(bc_expl_1D_2nd), POINTER, NOPASS :: bc_procedure 
    END TYPE ghost_state_Expl_1D_2nd  

    ! Type definition (pointers to subroutines for boundary condition - 
    ! (fully implicit scheme - 1st order accuracy in space))
    TYPE ghost_state_Impl_1D_1st
         PROCEDURE(bc_Impl_1D_1st), POINTER, NOPASS :: bc_procedure 
    END TYPE ghost_state_Impl_1D_1st  

    ! Type definition (pointers to subroutines for boundary condition - 
    ! (fully implicit scheme - 2nd order accuracy in space))
    TYPE ghost_state_Impl_1D_2nd
         PROCEDURE(bc_impl_1D_2nd), POINTER, NOPASS :: bc_procedure 
    END TYPE ghost_state_Impl_1D_2nd  

    ! Type definition (pointers to subroutines for boundary condition - 
    ! (explicit or point implicit scheme - 1st order accuracy in space))
    TYPE ghost_state_Expl_1D_SL_1st
         PROCEDURE(bc_expl_1D_SL_1st), POINTER, NOPASS :: bc_procedure 
    END TYPE ghost_state_Expl_1D_SL_1st

    ! Type definition (pointers to subroutines for boundary condition - 
    ! (explicit or point implicit scheme - 2nd order accuracy in space))
    TYPE ghost_state_Expl_1D_SL_2nd
         PROCEDURE(bc_expl_1D_SL_2nd), POINTER, NOPASS :: bc_procedure 
    END TYPE ghost_state_Expl_1D_SL_2nd

    ! Type definition (pointers to subroutines for boundary condition - 
    ! (fully implicit scheme - 2nd order accuracy in space))
    TYPE ghost_state_Impl_1D_SL_1st
         PROCEDURE(bc_Impl_1D_SL_1st), POINTER, NOPASS :: bc_procedure 
    END TYPE ghost_state_Impl_1D_SL_1st

    ! Type definition (pointers to subroutines for boundary condition - 
    ! (fully implicit scheme - 2nd order accuracy in space))
    TYPE ghost_state_Impl_1D_SL_2nd
         PROCEDURE(bc_Impl_1D_SL_2nd), POINTER, NOPASS :: bc_procedure 
    END TYPE ghost_state_Impl_1D_SL_2nd

    ! Function/subroutine pointer definitions

    ! Subroutine for computing conservative to primitive variables 
    PROCEDURE(prim_to_cons), POINTER, SAVE :: get_cons_from_prim
    PROCEDURE(prim_to_cons) :: prim_to_cons_1D, prim_to_cons_neq_1D
    PROCEDURE(prim_to_cons) :: prim_to_cons_1D_SL, prim_to_cons_neq_1D_SL

    ! Subroutine for computing primitive variables from conservative variables 
    PROCEDURE(cons_to_prim), POINTER, SAVE :: get_prim_from_cons
    PROCEDURE(cons_to_prim) :: cons_to_prim_1D, cons_to_prim_neq_1D
    PROCEDURE(cons_to_prim) :: cons_to_prim_1D_SL, cons_to_prim_neq_1D_SL

    ! Subroutine for computing conservative variables from physical variables
    PROCEDURE(cons_to_phys), POINTER, SAVE :: get_phys_from_cons
    PROCEDURE(cons_to_phys) :: cons_to_phys_1D_Eu, cons_to_phys_1D_Ns
    PROCEDURE(cons_to_phys) :: cons_to_phys_1D_SL_Eu, cons_to_phys_1D_SL_Ns
    PROCEDURE(cons_to_phys) :: cons_to_phys_neq_1D_Eu, cons_to_phys_neq_1D_Ns
    PROCEDURE(cons_to_phys) :: cons_to_phys_neq_1D_SL_Eu, cons_to_phys_neq_1D_SL_Ns

    ! Subroutine for computing conservative variables and physical properties from primitive variables 
    PROCEDURE(prim_to_cons_phys), POINTER, SAVE :: get_cons_phys_from_prim
    PROCEDURE(prim_to_cons_phys) :: prim_to_cons_phys_1D_Eu, prim_to_cons_phys_1D_Ns
    PROCEDURE(prim_to_cons_phys) :: prim_to_cons_phys_1D_SL_Eu, prim_to_cons_phys_1D_SL_Ns
    PROCEDURE(prim_to_cons_phys) :: prim_to_cons_phys_neq_1D_Eu, prim_to_cons_phys_neq_1D_Ns
    PROCEDURE(prim_to_cons_phys) :: prim_to_cons_phys_neq_1D_SL_Eu, prim_to_cons_phys_neq_1D_SL_Ns

    ! 1D flow time-step
    PROCEDURE(time_step_1D), POINTER, SAVE :: compute_time_step_1D
    PROCEDURE(time_step_1D) :: time_step_time_acc_1D
    PROCEDURE(time_step_1D) :: inv_time_step_1D, visc_time_step_1D, inv_visc_time_step_1D

    ! 1D stagnation line flow time-step
    PROCEDURE(time_step_1D_SL), POINTER, SAVE :: compute_time_step_1D_SL
    PROCEDURE(time_step_1D_SL) :: time_step_time_acc_1D_SL
    PROCEDURE(time_step_1D_SL) :: inv_time_step_1D_SL, visc_time_step_1D_SL, inv_visc_time_step_1D_SL

    ! 1D inviscid flux
    PROCEDURE(flux_splitter_1D), POINTER, SAVE :: conv_flux_1D
    PROCEDURE(flux_splitter_1D) :: van_Leer_1D, roe_1D, hlle_1D, Steger_Warming_1D
    PROCEDURE(flux_splitter_1D) :: van_Leer_neq_1D, roe_neq_1D, hlle_neq_1D, Steger_Warming_neq1T_1D, & 
                                 & Steger_Warming_neqNT_1D, Steger_Warming_neqNT_Te_1D

    ! 1D inviscid flux Jacobians
    PROCEDURE(flux_splitter_Jac_1D), POINTER, SAVE :: conv_flux_Jac_1D
    PROCEDURE(flux_splitter_Jac_1D) :: conv_flux_num_Jac_1D
    PROCEDURE(flux_splitter_Jac_1D) :: Yoon_Jameson_Jac_1D, Pos_Neg_split_Jac_1D   
    PROCEDURE(flux_splitter_Jac_1D) :: Yoon_Jameson_Jac_neq_1D, Pos_Neg_split_Jac_neq1T_1D, & 
                                     & Pos_Neg_split_Jac_neqNT_1D, Pos_Neg_split_Jac_neqNT_Te_1D 

    ! 1D stagnation line inviscid flux
    PROCEDURE(flux_splitter_1D_SL), POINTER, SAVE :: conv_flux_1D_SL
    PROCEDURE(flux_splitter_1D_SL) :: van_Leer_1D_SL, roe_1D_SL, Steger_Warming_1D_SL
    PROCEDURE(flux_splitter_1D_SL) :: van_Leer_neq_1D_SL, roe_neq_1D_SL, Steger_Warming_neq1T_1D_SL, & 
                                    & Steger_Warming_neqNT_1D_SL, Steger_Warming_neqNT_Te_1D_SL,     & 
                                    & ausmPC_1D_SL, ausmPC_neq_1D_SL, ausmP_up_as_neq_1D_SL,         &
                                    & ausmP_up_as_CD_neq_1D_SL,StegerWarming_neq_1D_SL,              &
                                    & LDFSS2_neq_1D_SL, modified_SW_neq_1D_SL, ausmP_up2_neq_1D_SL

     ! 1D stagnation line inviscid flux Jacobians
    PROCEDURE(flux_splitter_Jac_1D_SL), POINTER, SAVE :: conv_flux_Jac_1D_SL
    PROCEDURE(flux_splitter_Jac_1D_SL) :: conv_flux_num_Jac_1D_SL    
    PROCEDURE(flux_splitter_Jac_1D_SL) :: Yoon_Jameson_Jac_1D_SL, Pos_Neg_split_Jac_1D_SL
    PROCEDURE(flux_splitter_Jac_1D_SL) :: Yoon_Jameson_Jac_neq_1D_SL, Pos_Neg_split_Jac_neq1T_1D_SL, & 
                                        & Pos_Neg_split_Jac_neqNT_1D_SL, Pos_Neg_split_Jac_neqNT_Te_1D_SL, &
                                        & Pos_Neg_split_Jac_general_1D_SL 

    ! 1D eigensystem (nonequilibrium flow only)
    PROCEDURE(eigensystem_neq_1D), POINTER, SAVE :: eigsys_neq_1D
    PROCEDURE(eigensystem_neq_1D) :: eigensystem_neq1T_1D, eigensystem_neqNT_1D, & 
                                   & eigensystem_neqNT_Te_1D
 
    ! 1D stagnation line eigensystem (nonequilibrium flow only)
    PROCEDURE(eigensystem_neq_1D_SL), POINTER, SAVE :: eigsys_neq_1D_SL
    PROCEDURE(eigensystem_neq_1D_SL) :: eigensystem_neq1T_1D_SL, eigensystem_neqNT_1D_SL, & 
                                      & eigensystem_neqNT_Te_1D_SL

    ! 1D Roe's average state (nonequilibrium flow only)
    PROCEDURE(roe_avgState_neq_1D), POINTER, SAVE :: roe_avg_neq_1D
    PROCEDURE(roe_avgState_neq_1D) :: roe_avgState_neq1T_1D, roe_avgState_neqNT_1D, & 
                                    & roe_avgState_neqNT_Te_1D

    ! 1D stagnation line Roe's average state (nonequilibrium flow only)
    PROCEDURE(roe_avgState_neq_1D), POINTER, SAVE :: roe_avg_neq_1D_SL
    PROCEDURE(roe_avgState_neq_1D) :: roe_avgState_neq1T_1D_SL, roe_avgState_neqNT_1D_SL, & 
                                    & roe_avgState_neqNT_Te_1D_SL

    ! 1D inviscid flux Jacobian (nonequilibrium flow only)
    PROCEDURE(flux_Jacobian_neq_1D), POINTER, SAVE :: inv_flux_Jac_neq_1D 
    PROCEDURE(flux_Jacobian_neq_1D) :: inviscid_flux_Jac_neq1T_1D, inviscid_flux_Jac_neqNT_1D, & 
                                     & inviscid_flux_Jac_neqNT_Te_1D 
   
    ! 1D stagnation line inviscid flux Jacobian (nonequilibrium flow only)
    PROCEDURE(flux_Jacobian_neq_1D_SL), POINTER, SAVE :: inv_flux_Jac_neq_1D_SL 
    PROCEDURE(flux_Jacobian_neq_1D_SL) :: inviscid_flux_Jac_neq1T_1D_SL, inviscid_flux_Jac_neqNT_1D_SL, & 
                                        & inviscid_flux_Jac_neqNT_Te_1D_SL

    ! 1D Diffusive flux
    PROCEDURE(diffusive_flux_1D), POINTER, SAVE :: diff_flux_1D
    PROCEDURE(diffusive_flux_1D) :: ns_flux_1D, null_ns_flux_1D
    PROCEDURE(diffusive_flux_1D) :: ns_flux_neq1T_1D, ns_flux_neqNT_1D, ns_flux_neqNT_Te_1D, & 
                                  & null_ns_flux_neq_1D

    ! 1D diffusive flux and related Jacobians
    PROCEDURE(diffusive_flux_1D_Jac), POINTER, SAVE :: diff_flux_1D_Jac
    PROCEDURE(diffusive_flux_1D_Jac) :: diff_flux_num_Jac_1D
    PROCEDURE(diffusive_flux_1D_Jac) :: ns_flux_1D_Jac, null_ns_flux_1D_Jac
    PROCEDURE(diffusive_flux_1D_Jac) :: ns_flux_neq1T_1D_Jac, ns_flux_neqNT_1D_Jac, ns_flux_neqNT_Te_1D_Jac, & 
                                      & null_ns_flux_neq_1D_Jac

    ! 1D stagnation line diffusive flux
    PROCEDURE(diffusive_flux_1D_SL), POINTER, SAVE :: diff_flux_1D_SL
    PROCEDURE(diffusive_flux_1D_SL) :: null_ns_flux_1D_SL
    PROCEDURE(diffusive_flux_1D_SL) :: ns_flux_1D_SL_cyl, ns_flux_1D_SL_sph 
    PROCEDURE(diffusive_flux_1D_SL) :: ns_flux_1D_SL_sph_metr
    PROCEDURE(diffusive_flux_1D_SL) :: ns_flux_neq1T_1D_SL_cyl, ns_flux_neq1T_1D_SL_sph 
    PROCEDURE(diffusive_flux_1D_SL) :: ns_flux_neq1T_1D_SL_sph_metr
    PROCEDURE(diffusive_flux_1D_SL) :: ns_flux_neqNT_1D_SL_cyl, ns_flux_neqNT_1D_SL_sph, &
                                     & ns_flux_neqNT_1D_SL_sph_metr
    PROCEDURE(diffusive_flux_1D_SL) :: ns_flux_neqNT_Te_1D_SL_cyl, ns_flux_neqNT_Te_1D_SL_sph, &
                                     & ns_flux_neqNT_Te_1D_SL_sph_metr 

    ! 1D stagnation line diffusive flux and related Jacobian
    PROCEDURE(diffusive_flux_1D_SL_Jac), POINTER, SAVE :: diff_flux_Jac_1D_SL
    PROCEDURE(diffusive_flux_1D_SL_Jac) :: diff_flux_num_Jac_1D_SL
    PROCEDURE(diffusive_flux_1D_SL_Jac) :: null_ns_flux_1D_SL_Jac
    PROCEDURE(diffusive_flux_1D_SL_Jac) :: ns_flux_1D_SL_cyl_Jac, ns_flux_1D_SL_sph_Jac
    PROCEDURE(diffusive_flux_1D_SL_Jac) :: ns_flux_1D_SL_sph_Jac_metr
    PROCEDURE(diffusive_flux_1D_SL_Jac) :: ns_flux_neq1T_1D_SL_cyl_Jac, ns_flux_neq1T_1D_SL_sph_Jac 
    PROCEDURE(diffusive_flux_1D_SL_Jac) :: ns_flux_neq1T_1D_SL_sph_Jac_metr
    PROCEDURE(diffusive_flux_1D_SL_Jac) :: ns_flux_neqNT_1D_SL_cyl_Jac, ns_flux_neqNT_1D_SL_sph_Jac
    PROCEDURE(diffusive_flux_1D_SL_Jac) :: ns_flux_neqNT_Te_1D_SL_cyl_Jac, ns_flux_neqNT_Te_1D_SL_sph_Jac

    ! 1D spectral radius 
    PROCEDURE(spectral_radius_1D), POINTER, SAVE :: get_inv_spectral_radius_1D
    PROCEDURE(spectral_radius_1D), POINTER, SAVE :: get_visc_spectral_radius_1D 
    PROCEDURE(spectral_radius_1D) :: null_visc_spectral_radius_1D
    PROCEDURE(spectral_radius_1D) :: inv_spectral_radius_1D, visc_spectral_radius_1D, visc_spectral_radius_neq_1D

    ! 1D preconditioning velocity
    PROCEDURE(prec_vel_1D), POINTER, SAVE :: get_prec_vel_1D
    PROCEDURE(prec_vel_1D) :: inv_prec_vel_1D, visc_prec_vel_1D

    ! 1D stagnation line pre-conditioning velocity
    PROCEDURE(prec_vel_1D_SL), POINTER, SAVE :: get_prec_vel_1D_SL
    PROCEDURE(prec_vel_1D_SL) :: inv_prec_vel_1D_SL, visc_prec_vel_1D_SL

    ! 1D stagnation line spectral radius
    PROCEDURE(spectral_radius_1D_SL), POINTER, SAVE :: get_inv_spectral_radius_1D_SL
    PROCEDURE(spectral_radius_1D_SL), POINTER, SAVE :: get_visc_spectral_radius_1D_SL
    PROCEDURE(spectral_radius_1D_SL) :: null_visc_spectral_radius_1D_SL
    PROCEDURE(spectral_radius_1D_SL) :: inv_spectral_radius_1D_SL, visc_spectral_radius_1D_SL, visc_spectral_radius_neq_1D_SL

    ! 1D stagnation line stress tensor
    PROCEDURE(stress_tensor_1D_SL), POINTER, SAVE :: get_stress_tensor_1D_SL
    PROCEDURE(stress_tensor_1D_SL) :: stress_tensor_1D_SL_cyl, stress_tensor_1D_SL_sph

    ! Transport coefficients (calorically perfect gas flow only)
    PROCEDURE(transp_Coeff), POINTER, SAVE :: get_transpCoeff
    PROCEDURE(transp_Coeff) :: transpCoeff_sutherland, transpCoeff_inv_power

    ! 1D Polynomial reconstruction
    PROCEDURE(muscl_1D), POINTER, SAVE :: rec_1D
    PROCEDURE(muscl_1D) :: poly_rec_1D
    PROCEDURE(muscl_1D) :: poly_rec_neq_1D
    PROCEDURE(muscl_1D) :: null_poly_rec_1D

    ! 1D stagnation line polynomial reconstruction
    PROCEDURE(muscl_1D_SL), POINTER, SAVE :: rec_1D_SL
    PROCEDURE(muscl_1D_SL) :: poly_rec_1D_SL
    PROCEDURE(muscl_1D_SL) :: poly_rec_1D_SL_metr
    PROCEDURE(muscl_1D_SL) :: poly_rec_neq_1D_SL
    PROCEDURE(muscl_1D_SL) :: poly_rec_neq_1D_SL_metr
    PROCEDURE(muscl_1D_SL) :: null_poly_rec_1D_SL

    ! 1D Boundary conditions (explicit or point implicit scheme - 1st order accuracy in space)
    PROCEDURE(bc_expl_1D_1st) :: sup_out_Expl, sub_in_rhoin_pin_1DExpl, sup_in_1DExpl, & 
                               & sub_out_pout_1DExpl 
    PROCEDURE(bc_expl_1D_1st) :: sup_in_neq_1DExpl, sub_in_neq_rhoiin_T_1DExpl,        & 
                               & sub_out_neq_Tout_1DExpl

    ! 1D Boundary conditions (explicit or point implicit scheme - 2nd order accuracy in space)
    PROCEDURE(bc_expl_1D_2nd) :: sup_out_1DExpl_2nd, sub_in_rhoin_pin_1DExpl_2nd,       & 
                               & sup_in_1DExpl_2nd, sub_out_pout_1DExpl_2nd
    PROCEDURE(bc_expl_1D_2nd) :: sup_in_neq_1DExpl_2nd, sub_in_neq_rhoiin_T_1DExpl_2nd, &  
                               & sup_out_neq_1DExpl_2nd 

    ! 1D Boundary conditions (fully implicit scheme - 1st order accuracy in space)
    PROCEDURE(bc_impl_1D_1st) :: sup_out_1DImpl, sub_in_rhoin_pin_1DImpl, sub_out_pout_1DImpl,  & 
                               & sup_in_1DImpl
    PROCEDURE(bc_impl_1D_1st) :: sub_in_neq_rhoiin_T_1DImpl, sup_out_neq_1DImpl, sup_in_neq_1DImpl, & 
                               & sub_out_neq_Tout_1DImpl

    ! 1D Boundary conditions (fully implicit scheme - 2nd order accuracy in space)
    PROCEDURE(bc_impl_1D_2nd) :: sup_out_1DImpl_2nd, sub_in_rhoin_pin_1DImpl_2nd,             & 
                               & sup_in_1DImpl_2nd, sub_out_pout_1DImpl_2nd
    PROCEDURE(bc_impl_1D_2nd) :: sup_out_neq_1DImpl_2nd, sub_in_neq_rhoiin_T_1DImpl_2nd, sup_in_neq_1DImpl_2nd

    ! 1D stagnation line boundary conditions (explicit or point implicit scheme - 1st order accuracy in space)
    PROCEDURE(bc_expl_1D_SL_1st) :: sup_in_1D_SL_Expl, slip_wall_1D_SL_Expl, no_slip_iso_Twall_1D_SL_Expl
    PROCEDURE(bc_expl_1D_SL_1st) :: slip_wall_neq_1D_SL_Expl, sup_in_neq_1D_SL_Expl, no_slip_iso_Twall_neq_1D_SL_Expl
    PROCEDURE(bc_expl_1D_SL_1st) :: sub_in_1D_SL_Expl 
    PROCEDURE(bc_expl_1D_SL_1st) :: sub_in_vel_neq_1D_SL_Expl
    PROCEDURE(bc_expl_1D_SL_1st) :: no_slip_SEB_ablation_neq_1D_SL_Expl
    PROCEDURE(bc_expl_1D_SL_1st) :: no_slip_isothermal_ablation_neq_1D_SL_Expl
    PROCEDURE(bc_expl_1D_SL_1st) :: no_slip_adiabatic_neq_1D_SL_Expl
    PROCEDURE(bc_expl_1D_SL_1st) :: no_slip_isothermal_general_ablation_neq_1D_SL_Expl
    PROCEDURE(bc_expl_1D_SL_1st) :: no_slip_SEB_general_ablation_neq_1D_SL_Expl

    ! 1D stagnation line boundary conditions (explicit or point implicit scheme - 2nd order accuracy in space)
    PROCEDURE(bc_expl_1D_SL_2nd) :: sup_in_1D_SL_Expl_2nd, slip_wall_1D_SL_Expl_2nd, no_slip_iso_Twall_1D_SL_Expl_2nd
    PROCEDURE(bc_expl_1D_SL_2nd) ::  sub_in_1D_SL_Expl_2nd   
    PROCEDURE(bc_expl_1D_SL_2nd) :: sup_in_neq_1D_SL_Expl_2nd, slip_wall_neq_1D_SL_Expl_2nd, & 
                                  & no_slip_iso_Twall_neq_1D_SL_Expl_2nd
    PROCEDURE(bc_expl_1D_SL_2nd) :: no_slip_adiabatic_neq_1D_SL_Expl_2nd
    PROCEDURE(bc_expl_1D_SL_2nd) :: sub_in_vel_neq_1D_SL_Expl_2nd
    PROCEDURE(bc_expl_1D_SL_2nd) :: no_slip_SEB_ablation_neq_1D_SL_Expl_2nd
    PROCEDURE(bc_expl_1D_SL_2nd) :: no_slip_isothermal_ablation_neq_1D_SL_Expl_2nd
    PROCEDURE(bc_expl_1D_SL_2nd) :: no_slip_isothermal_general_ablation_neq_1D_SL_Expl_2nd
 
    ! 1D stagnation line boundary conditions (fully implicit scheme - 1st order accuracy in space)
    PROCEDURE(bc_impl_1D_SL_1st) :: sup_in_1D_SL_Impl, slip_wall_1D_SL_Impl, no_slip_iso_Twall_1D_SL_Impl 
    PROCEDURE(bc_impl_1D_SL_1st) :: sup_in_neq_1D_SL_Impl

    ! 1D stagnation line boundary conditions (fully implicit scheme - 2nd order accuracy in space)
    PROCEDURE(bc_impl_1D_SL_2nd) :: sup_in_1D_SL_Impl_2nd, slip_wall_1D_SL_Impl_2nd, no_slip_iso_Twall_1D_SL_Impl_2nd
    PROCEDURE(bc_impl_1D_SL_2nd) :: sup_in_neq_1D_SL_Impl_2nd

    ! Evolve solution 
    PROCEDURE(time_integration), POINTER, SAVE :: evolve_solution
    PROCEDURE(time_integration) :: fe_1D, si_1D, pi_1D, fi_1D
    PROCEDURE(time_integration) :: fe_1D_SL, si_1D_SL, fi_1D_SL

    ! Limiter
    PROCEDURE(muscl_limiter), POINTER, SAVE :: limiter
    PROCEDURE(muscl_limiter) :: van_leer, minimod, superbee, koren, mc, charm, hcus, & 
                              & hquick, van_albada1, van_albada2, ospre, no_limiter

    ! Source term(s) and source term Jacobian(s)
    PROCEDURE(source_term) :: source_term_quasi1D
    PROCEDURE(source_term) :: source_term_neq_kinetics, source_term_neq_el_kinetics, source_term_neq_quasi1D
    PROCEDURE(source_term_Jac) :: source_term_num_Jac
    PROCEDURE(source_term_Jac) :: source_term_quasi1D_Jac
    PROCEDURE(source_term_Jac) :: source_term_neq_kinetics_Jac, source_term_neq_el_kinetics_Jac, & 
                                & source_term_neq_quasi1D_Jac, source_term_neq_el_quasi1D_Jac

    PROCEDURE(inv_source_term_SL) :: source_term_inv_1D_SL_cyl, source_term_inv_1D_SL_sph
    PROCEDURE(inv_source_term_SL) :: source_term_neq_kinetics_SL, source_term_neq_el_kinetics_SL
    PROCEDURE(inv_source_term_SL) :: source_term_inv_neq_1D_SL_cyl, source_term_inv_neq_1D_SL_sph
    PROCEDURE(inv_source_term_Jac_SL) :: source_term_inv_num_Jac_SL 
    PROCEDURE(inv_source_term_Jac_SL) :: source_term_inv_1D_SL_cyl_Jac, source_term_inv_1D_SL_sph_Jac  
    PROCEDURE(inv_source_term_Jac_SL) :: source_term_inv_neq_1D_SL_cyl_Jac, source_term_inv_neq_1D_SL_sph_Jac
    PROCEDURE(inv_source_term_Jac_SL) :: source_term_neq_kinetics_Jac_SL, source_term_neq_el_kinetics_Jac_SL

    ! Diffusive source term for stagnation line flows 
    PROCEDURE(diff_source_term_SL), POINTER, SAVE :: get_diff_source_term_SL
    PROCEDURE(diff_source_term_SL) :: null_diff_source_term_SL
    PROCEDURE(diff_source_term_SL) :: source_term_diff_1D_SL_cyl, source_term_diff_1D_SL_sph
    PROCEDURE(diff_source_term_SL) :: source_term_diff_neq1T_1D_SL_cyl, source_term_diff_neq1T_1D_SL_sph
    PROCEDURE(diff_source_term_SL) :: source_term_diff_neqNT_1D_SL_cyl, source_term_diff_neqNT_1D_SL_sph
    PROCEDURE(diff_source_term_SL) :: source_term_diff_neqNT_Te_1D_SL_cyl, source_term_diff_neqNT_Te_1D_SL_sph

    PROCEDURE(diff_source_term_SL_Jac), POINTER, SAVE :: get_diff_source_term_SL_Jac
    PROCEDURE(diff_source_term_SL_Jac) :: null_diff_source_term_SL_Jac
    PROCEDURE(diff_source_term_SL_Jac) :: source_term_diff_num_Jac_SL
    PROCEDURE(diff_source_term_SL_Jac) :: source_term_diff_1D_SL_cyl_Jac, source_term_diff_1D_SL_sph_Jac
    PROCEDURE(diff_source_term_SL_Jac) :: source_term_diff_neq1T_1D_SL_cyl_Jac, source_term_diff_neq1T_1D_SL_sph_Jac
    PROCEDURE(diff_source_term_SL_Jac) :: source_term_diff_neqNT_1D_SL_cyl_Jac, source_term_diff_neqNT_1D_SL_sph_Jac
    PROCEDURE(diff_source_term_SL_Jac) :: source_term_diff_neqNT_Te_1D_SL_cyl_Jac, source_term_diff_neqNT_Te_1D_SL_sph_Jac

    TYPE(source), ALLOCATABLE, DIMENSION(:), SAVE :: get_source_term
    TYPE(source_Jac), ALLOCATABLE, DIMENSION(:), SAVE :: get_source_term_Jac

    TYPE(inv_source_SL), ALLOCATABLE, DIMENSION(:), SAVE :: get_inv_source_term_SL
    TYPE(inv_source_Jac_SL), ALLOCATABLE, DIMENSION(:), SAVE :: get_inv_source_term_Jac_SL

    TYPE(ghost_state_Expl_1D_1st), ALLOCATABLE, DIMENSION(:), SAVE :: get_ghost_state_Expl_1D_1st
    TYPE(ghost_state_Expl_1D_2nd), ALLOCATABLE, DIMENSION(:), SAVE :: get_ghost_state_Expl_1D_2nd
    TYPE(ghost_state_Impl_1D_1st), ALLOCATABLE, DIMENSION(:), SAVE :: get_ghost_state_Impl_1D_1st
    TYPE(ghost_state_Impl_1D_2nd), ALLOCATABLE, DIMENSION(:), SAVE :: get_ghost_state_Impl_1D_2nd

    TYPE(ghost_state_Expl_1D_SL_1st), ALLOCATABLE, DIMENSION(:), SAVE :: get_ghost_state_Expl_1D_SL_1st
    TYPE(ghost_state_Expl_1D_SL_2nd), ALLOCATABLE, DIMENSION(:), SAVE :: get_ghost_state_Expl_1D_SL_2nd
    TYPE(ghost_state_Impl_1D_SL_1st), ALLOCATABLE, DIMENSION(:), SAVE :: get_ghost_state_Impl_1D_SL_1st
    TYPE(ghost_state_Impl_1D_SL_2nd), ALLOCATABLE, DIMENSION(:), SAVE :: get_ghost_state_Impl_1D_SL_2nd

  END MODULE mod_function_pointer
!------------------------------------------------------------------------------!
