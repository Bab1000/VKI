!------------------------------------------------------------------------------!
!> This subroutine computes the solution variation in case of use of a source implicit time-integration scheme for 1D flows.
  SUBROUTINE si_1D ()

#include "../config.h"

    USE mod_general_data,            ONLY: nb_cells, nb_eq, nb_prop, start_u_phys, start_prop_phys, volumes, s_vec,  & 
                                         & js_mat, iden, u, du, dres, cell_prop
    USE mod_numerics_data,           ONLY: nb_source_terms, fc, fd, u_left, u_right, ull, ul, ur, urr, prop_left,    &   
                                         & prop_right,  prop_ll, prop_l, prop_r, prop_rr, aux, s, js, inv_js
    USE mod_function_pointer,        ONLY: conv_flux_1D, diff_flux_1D, rec_1D, compute_time_step_1D, get_source_term_Jac
#ifdef GOTO_BLAS 
    USE mod_algebra,                 ONLY: inv_matrix_GOTO_BLAS
#else 
    USE mod_algebra,                 ONLY: inv_matrix
#endif 

    IMPLICIT NONE

    INTEGER :: i, j, k
    INTEGER :: l, ll, r, rr, start_l, start_r
    INTEGER :: ll_prop, l_prop, r_prop, rr_prop
    INTEGER :: ll_eq, l_eq, r_eq, rr_eq
    REAL(KIND=8) :: tmp, tmp2
    REAL(KIND=8) :: dt, vol_ov_dt, ov_dt
    REAL(KIND=8) :: vol_l, vol_r

    INTERFACE 
      SUBROUTINE apply_bc_1D_Expl(bound_id, phys_data1, u_phys1, phys_data2, u_phys2, & 
                                & ghost_data1, u_ghost1, ghost_data2, u_ghost2)
        INTEGER, INTENT(IN) :: bound_id
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys1, phys_data1, u_phys2, phys_data2 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost1, ghost_data1, u_ghost2, ghost_data2
      END SUBROUTINE apply_bc_1D_Expl
    END INTERFACE 

    ! Initialization
    du = 0.d0
    dres = 0.d0
    
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
    
    CALL apply_bc_1D_Expl(1, prop_r, ur, prop_rr, urr, prop_l, ul, prop_ll, ull)
         
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

    CALL apply_bc_1D_Expl(2, prop_l, ul, prop_ll, ull, prop_r, ur, prop_rr, urr) 

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
    CALL diff_flux_1D (vol_l, vol_r, prop_l, prop_r, ul, ur, fd)
 
    ! Rhs residual initialization
    start_l = l*nb_eq
    DO i = 1,nb_eq 
       du(start_l + i) = - (fc(i) - fd(i))
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
 
       ! Diffusive flux 
       CALL diff_flux_1D (vol_l, vol_r, prop_l, prop_r, ul, ur, fd) 

       ! Time-step 
       CALL compute_time_step_1D (vol_l, prop_l, dt)
       vol_ov_dt = vol_l/dt

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
       
       ! Lhs matrix (inclusion of the source term Jacobian(s))
       js = iden*vol_ov_dt - js*vol_l

       ! Lhs matrix inversion
#ifdef GOTO_BLAS
       CALL inv_matrix_GOTO_BLAS (nb_eq, js)   

       DO j = 1,nb_eq 
          tmp  = fc(j) - fd(j)
          aux(j)          = du(start_l + j) + tmp - s(j)*vol_l
          du(start_r + j) = - tmp
       ENDDO      

       CALL DGEMV ('n', nb_eq, nb_eq, 1.d0, js, nb_eq, aux, 1, 0.d0, du(start_l + 1:start_r), 1)
#else 
       CALL inv_matrix (nb_eq, js, inv_js)

       ! Rhs residual update (matrix multiplication is performed in the loop)
       aux = 0.d0
       DO j = 1,nb_eq 
          tmp = fc(j) - fd(j)
          tmp2 = du(start_l + j) + tmp - s(j)*vol_l
          du(start_r + j) = du(start_r + j) - tmp
          DO k = 1,nb_eq 
             aux(k) = aux(k) + inv_js(k,j)*tmp2  
          ENDDO
       ENDDO
      
       du(start_l + 1:start_r) = aux
#endif
       
       dres(start_l + 1:start_r) = du(start_l + 1:start_r)*ov_dt

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

  END SUBROUTINE si_1D
!------------------------------------------------------------------------------!
