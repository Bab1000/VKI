!------------------------------------------------------------------------------!
!> This subroutine computes the convective flux by means of the HLLE
!! approximate Riemann solver for nonequilibrium gas flows.
  SUBROUTINE hlle_neq_1D (nx, left_data, right_data, u_left, u_right, f)

#include "../config.h"
 
    USE mod_general_data,             ONLY: nb_ns, nb_eq, nb_temp, pos_u_cell, pos_pres_cell, pos_c_cell, & 
                                          & pos_rhou, pos_rhoE
    USE mod_numerics_data,            ONLY: deltau, cons_Roe, phys_data_Roe, fl, fr
    USE mod_function_pointer,         ONLY: roe_avg_neq_1D

    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8) :: tmp1, tmp2, tmp3
    REAL(KIND=8) :: bl, br, b_plus, b_minus, diff
    REAL(KIND=8) :: c, u, un
    REAL(KIND=8) :: cl, cr, pl, pr, ul, ur, unl, unr 
    
    REAL(KIND=8), INTENT(IN) :: nx                        !< normal to the cell interface
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: f          !< numerical flux

    ! Speed of sound, velocity and pressure of left and right states
    cl = left_data(pos_c_cell)
    ul = left_data(pos_u_cell)
    pl = left_data(pos_pres_cell)

    cr = right_data(pos_c_cell)
    ur = right_data(pos_u_cell)
    pr = right_data(pos_pres_cell)

    ! Normal velocity of left and right states
    unl = ul*nx 
    unr = ur*nx

    ! Conservative variable jump (ur - ul) and inviscid fluxes
    ! of left and right states
    DO i = 1,nb_eq 
       tmp1 = u_left(i)
       tmp2 = u_right(i)
       deltau(i) = tmp2 - tmp1
       fl(i)     = tmp1*unl 
       fr(i)     = tmp2*unr 
    ENDDO

    ! Global momentum and energy equations 
    fl(pos_rhou) = fl(pos_rhou) + pl*nx
    fl(pos_rhoE) = fl(pos_rhoE) + pl*unl

    fr(pos_rhou) = fr(pos_rhou) + pr*nx
    fr(pos_rhoE) = fr(pos_rhoE) + pr*unr

    ! Compute Roe's averaged state
    CALL roe_avg_neq_1D (u_left, u_right, left_data, right_data, cons_Roe, phys_data_Roe)

    c  = phys_data_Roe(pos_c_cell)
    u  = phys_data_Roe(pos_u_cell)
    un = u*nx

    ! Wave speeds for the HLLE approximate Riemann solver 
    ! (see the original paper of Einfeldt for more details - JCP vol. 92 1991 pp 273)
    bl = MIN(un - c,unl - cl)
    br = MAX(un + c,unr + cr)
    b_plus  = MAX(br,0.d0)
    b_minus = MIN(bl,0.d0)
    diff    = 1.d0/(b_plus - b_minus)

    ! Numerical flux 
#ifdef GOTO_BLAS
    tmp1 = b_plus*diff
    tmp2 = - b_minus*diff
    tmp3 = tmp1*b_minus
 
    CALL DSCAL (nb_eq, tmp3, deltau, 1)
    CALL DCOPY (nb_eq, deltau, 1, f, 1)
    CALL DAXPY (nb_eq, tmp1, fl, 1, f, 1)
    CALL DAXPY (nb_eq, tmp2, fr, 1, f, 1)
#else 
    f = (b_plus*fl - b_minus*fr + b_plus*b_minus*deltau)*diff
#endif

  END SUBROUTINE hlle_neq_1D 
!------------------------------------------------------------------------------!
