!------------------------------------------------------------------------------!
!> This subroutine computes the solution variation in case of use of a fully implicit time-integration scheme scheme for 1D flows.
  SUBROUTINE fi_1D()

#include "../config.h"

    USE mod_general_data,            ONLY: nb_cells, nb_eq, nb_prop, start_u_phys, finish_u_phys, start_prop_phys,       & 
                                         & volumes, iden, s_vec, js_mat, u, ru, du, dres, cell_prop
    USE mod_numerics_data,           ONLY: nb_source_terms, fc, fd, s, u_right, u_left, ull, ul, ur, urr, prop_left,     &
                                         & prop_right, prop_ll, prop_l, prop_r, prop_rr, js, jfl, jfr, jfdl, jfdr,       & 
                                         & j_ghost, central, left_bound, right_bound
    USE mod_function_pointer,        ONLY: conv_flux_1D, conv_flux_Jac_1D, diff_flux_1D_Jac, compute_time_step_1D,       & 
                                         & rec_1D, get_source_term_Jac
#ifdef GOTO_BLAS 
    USE mod_algebra,                 ONLY: tridiag_block_solver_GOTO_BLAS, Ai_mat, Bi_mat, Ci_mat
#else
    USE mod_algebra,                 ONLY: tridiag_block_solver, Ai_mat, Bi_mat, Ci_mat
#endif

    IMPLICIT NONE

    INTEGER :: i, j, k
    INTEGER :: pos1, pos2
    INTEGER :: l, ll, r, rr, start_l, start_r
    INTEGER :: ll_prop, l_prop, r_prop, rr_prop
    INTEGER :: ll_eq, l_eq, r_eq, rr_eq
    REAL(KIND=8) :: tmp1, tmp2
    REAL(KIND=8) :: dt, vol_ov_dt, ov_dt
    REAL(KIND=8) :: vol_l, vol_r

    INTERFACE 
      SUBROUTINE apply_bc_1D_Impl(bound_id, phys_data1, u_phys1, phys_data2, u_phys2, & 
                                & ghost_data1, u_ghost1, ghost_data2, u_ghost2, jb)
        INTEGER, INTENT(IN) :: bound_id
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys1, phys_data1, u_phys2, phys_data2 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost1, ghost_data1, u_ghost2, ghost_data2
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jb
      END SUBROUTINE apply_bc_1D_Impl
    END INTERFACE  

    ! Boundary conditions
    ! Left boundary - boundary condition 1

    ! Physical cells
    ! Conservative variables
    DO i = 1,nb_eq
       ur(i)  = u(start_u_phys + i)
       urr(i) = u(start_u_phys + nb_eq + i)
    ENDDO

    ! Physical properties
    DO i = 1,nb_prop
       prop_r(i)  = cell_prop(start_prop_phys + i)
       prop_rr(i) = cell_prop(start_prop_phys + nb_prop + i)
    ENDDO
    
    CALL apply_bc_1D_Impl(1, prop_r, ur, prop_rr, urr, prop_l, ul, prop_ll, ull, j_ghost)
       
    ! Fill conservative variable and physical variable solution vectors
    ! Conservative variables
    DO i = 1,nb_eq
       u(i) = ull(i)
       u(nb_eq + i) = ul(i)
    ENDDO

    ! Physical properties
    DO i = 1,nb_prop
       cell_prop(i) = prop_ll(i)
       cell_prop(nb_prop + i) = prop_l(i)
    ENDDO 

    ! Implicit left boundary condition
    ! Convective flux Jacobians
    vol_l = volumes(2)
    vol_r = volumes(3)
    CALL conv_flux_Jac_1D (1.d0, prop_l, prop_r, ul, ur, jfl, jfr)
    
    ! Diffusive flux and Jacobians
    CALL diff_flux_1D_Jac (vol_l, vol_r, prop_l, prop_r, ul, ur, fd, jfdl, jfdr)
    
    ! Implicit left boundary condition 
#ifdef GOTO_BLAS 
    left_bound = jfr - jfdr 
    js_mat     = jfl - jfdl
    CALL DGEMM('n', 'n', nb_eq, nb_eq, nb_eq, -1.d0, js_mat, nb_eq, j_ghost, nb_eq, -1.d0, left_bound, nb_eq)
#else
    left_bound = 0.d0

    DO j = 1,nb_eq
       DO i = 1,nb_eq
          tmp1 = 0.d0
          DO k = 1,nb_eq
             tmp1 = tmp1 + (jfl(i,k) - jfdl(i,k))*j_ghost(k,j)
          ENDDO
          left_bound(i,j) = - tmp1 - (jfr(i,j) - jfdr(i,j))
       ENDDO
    ENDDO
#endif
   
    ! Right boundary - boundary condition 2
    ! Physical cells
    ! Conservative variables
    DO i = 1,nb_eq
       ul(i)  = u((nb_cells + 1)*nb_eq + i)
       ull(i) = u(nb_cells*nb_eq + i)
    ENDDO

    ! Physical properties
    DO i = 1,nb_prop
       prop_l(i)  = cell_prop((nb_cells + 1)*nb_prop + i)
       prop_ll(i) = cell_prop(nb_cells*nb_prop + i)
    ENDDO

    CALL apply_bc_1D_Impl(2, prop_l, ul, prop_ll, ull, prop_r, ur, prop_rr, urr, j_ghost) 
    
    ! Fill conservative variable and physical variable solution vectors
    ! Conservative variables
    DO i = 1,nb_eq
       u((nb_cells + 2)*nb_eq + i) = ur(i)
       u((nb_cells + 3)*nb_eq + i) = urr(i)
    ENDDO

    ! Physical properties
    DO i = 1,nb_prop
       cell_prop((nb_cells + 2)*nb_prop + i) = prop_r(i)
       cell_prop((nb_cells + 3)*nb_prop + i) = prop_rr(i)
    ENDDO

    ! Implicit right boundary condition 
    ! Convective flux Jacobians 
    vol_l = volumes(nb_cells + 2)
    vol_r = volumes(nb_cells + 3)
    CALL conv_flux_Jac_1D (1.d0, prop_l, prop_r, ul, ur, jfl, jfr)
    
    ! Diffusive flux and Jacobians
    CALL diff_flux_1D_Jac (vol_l, vol_r, prop_l, prop_r, ul, ur, fd, jfdl, jfdr)

#ifdef GOTO_BLAS
    right_bound = 0.d0
    js_mat      = jfr - jfdr
    CALL DGEMM('n', 'n', nb_eq, nb_eq, nb_eq, 1.d0, js_mat, nb_eq, j_ghost, nb_eq, 1.d0, right_bound, nb_eq)
#else
    right_bound = 0.d0

    DO j = 1,nb_eq
       DO i = 1,nb_eq
          tmp1 = 0.d0
          DO k = 1,nb_eq
             tmp1 = tmp1 + (jfr(i,k) - jfdr(i,k))*j_ghost(k,j)
          ENDDO
          right_bound(i,j) = tmp1 
       ENDDO
    ENDDO
#endif

    ! Initialization for loop over cells
    ll = 1
    l  = ll + 1
    r  = l + 1
    rr = r + 1
       
    ll_eq = (ll - 1)*nb_eq
    l_eq  = (l - 1)*nb_eq
    r_eq  = (r - 1)*nb_eq
    rr_eq = (rr - 1)*nb_eq 

    ll_prop = (ll - 1)*nb_prop
    l_prop  = (l - 1)*nb_prop
    r_prop  = (r - 1)*nb_prop
    rr_prop = (rr - 1)*nb_prop

    DO i = 1,nb_eq 
       ull(i) = u(ll_eq + i) 
       ul(i)  = u(l_eq + i)
       ur(i)  = u(r_eq + i)        
       urr(i) = u(rr_eq + i) 
    ENDDO       

    DO i = 1,nb_prop
       prop_ll(i) = cell_prop(ll_prop + i) 
       prop_l(i)  = cell_prop(l_prop + i)
       prop_r(i)  = cell_prop(r_prop + i)
       prop_rr(i) = cell_prop(rr_prop + i) 
    ENDDO      
 
    ! Volumes of left and right states
    vol_l = volumes(l) 
    vol_r = volumes(r) 
    
    ! Polynomial re-construction
    CALL rec_1D (ull, ul, ur, urr, prop_ll, prop_l, prop_r, prop_rr, u_left, u_right, prop_left, prop_right)
     
    ! Convective flux
    CALL conv_flux_1D (1.d0, prop_left, prop_right, u_left, u_right, fc) 
 
    ! Diffusive flux 
    CALL diff_flux_1D_Jac (vol_l, vol_r, prop_l, prop_r, ul, ur, fd, jfdl, jfdr)
   
    ! Rhs residual initialization
    start_l = l*nb_eq
    DO i = 1,nb_eq 
       ru(start_l + i) = fc(i) - fd(i)
    ENDDO

    ! Data exchange (left -> left-left, right -> left, right-right -> right)
    ! Indices 
    ll = l
    l  = r
    r  = rr
    rr = r + 1      
 
    ! Volumes 
    vol_l = vol_r
    vol_r = volumes(r)

    ! Conservative variables
    ull = ul
    ul  = ur
    ur  = urr

    ! Physical properties
    prop_ll = prop_l 
    prop_l  = prop_r 
    prop_r  = prop_rr
       
    start_l = (l - 1)*nb_eq
    start_r = start_l + nb_eq
    rr_eq   = (rr - 1)*nb_eq 
    rr_prop = (rr - 1)*nb_prop

    ! Initialization of Ai and Bi matrices
    DO j = 1,nb_eq 
       DO i = 1,nb_eq
          Ai_mat(i,j) = 0.d0  
          Bi_mat(i,j) = 0.d0
       ENDDO
    ENDDO

    ! Loop over physical cells
    DO i = 1,nb_cells

       ! Right-right state
       ! Conservative variables
       DO j = 1,nb_eq
          urr(j) = u(rr_eq + j)
       ENDDO

       ! Physical properties
       DO j = 1,nb_prop 
          prop_rr(j) = cell_prop(rr_prop + j)
       ENDDO
       
       ! Polynomial re-construction
       CALL rec_1D (ull, ul, ur, urr, prop_ll, prop_l, prop_r, prop_rr, u_left, u_right, prop_left, prop_right)

       ! Convective flux
       CALL conv_flux_1D (1.d0, prop_left, prop_right, u_left, u_right, fc) 

       ! Convective flux Jacobians 
       CALL conv_flux_Jac_1D (1.d0, prop_l, prop_r, ul, ur, jfl, jfr)

       ! Diffusive flux and Jacobians
       CALL diff_flux_1D_Jac (vol_l, vol_r, prop_l, prop_r, ul, ur, fd, jfdl, jfdr)

       ! Source term(s) and Jacobian(s)
       s  = 0.d0
       js = 0.d0
       DO j = 1,nb_source_terms 
          CALL get_source_term_Jac(j)%source_Jac(j, i, prop_l, ul, s_vec, js_mat)
#ifdef GOTO_BLAS
          CALL DAXPY (nb_eq, 1.d0, s_vec, 1, s, 1)
#else 
          s  = s + s_vec
#endif
          js = js + js_mat
       ENDDO  

       ! Time-step  
       CALL compute_time_step_1D (vol_l, prop_l, dt) 
       vol_ov_dt = vol_l/dt       
       ov_dt = 1.d0 / dt

       ! Fill the block tridiagonal system
       central = iden*vol_ov_dt - js*vol_l

       DO j = 1,nb_eq 
          tmp1 = fc(j) - fd(j)
          ru(start_l + j) = ru(start_l + j) - tmp1 + s(j)*vol_l
          ru(start_r + j) = tmp1
          pos1 = (i - 1)*nb_eq + j
          pos2 = pos1 + nb_eq
          DO k = 1,nb_eq
             tmp1 = jfl(k,j) - jfdl(k,j)
             tmp2 = jfr(k,j) - jfdr(k,j)
             Ai_mat(k,pos1) = Ai_mat(k,pos1) + central(k,j) + tmp1
             Ai_mat(k,pos2) = - tmp2
             Bi_mat(k,pos2) = - tmp1
             Ci_mat(k,pos1) =   tmp2 
          ENDDO
       ENDDO     
 
       ! Data exchange (left -> left-left, right -> left, right-right -> right)
       ! Indices
       ll = l
       l  = r
       r  = rr
       rr = rr + 1

       start_l = start_r 
       start_r = start_r + nb_eq
       rr_eq   = rr_eq + nb_eq
       rr_prop = rr_prop + nb_prop

       ! Volumes 
       vol_l = vol_r
       vol_r = volumes(r)

       ! Conservative variables 
       ull = ul
       ul  = ur
       ur  = urr
 
       ! Physical properties      
       prop_ll = prop_l 
       prop_l  = prop_r 
       prop_r  = prop_rr 

    ENDDO

    ! Implicit boundary conditions
    ! Left boundary - boundary condition 1
    DO j = 1,nb_eq 
       DO i = 1,nb_eq
          Ai_mat(i,j) = Ai_mat(i,j) + left_bound(i,j)  
       ENDDO
    ENDDO

    ! Right boundary - boundary condition 2
    pos1 = (nb_cells - 1)*nb_eq
    DO j = 1,nb_eq 
       DO i = 1,nb_eq
          Ai_mat(i,pos1 + j) = Ai_mat(i,pos1 + j) + right_bound(i,j)
       ENDDO
    ENDDO 

    ! Solve the block tridiagonal system
    pos1 = start_u_phys + 1
    pos2 = finish_u_phys
#ifdef GOTO_BLAS
    CALL tridiag_block_solver_GOTO_BLAS (nb_eq, nb_cells, ru(pos1:pos2), du(pos1:pos2))
#else
    CALL tridiag_block_solver (nb_eq, nb_cells, ru(pos1:pos2), du(pos1:pos2))
#endif

    ! Sign change
    du = -du
    
    dres = du * ov_dt
    
  END SUBROUTINE fi_1D
!------------------------------------------------------------------------------!
