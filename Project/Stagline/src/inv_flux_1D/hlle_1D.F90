!------------------------------------------------------------------------------!
!> This subroutine computes the convective flux by means of the HLLE approximate Riemann solver 
!! for 1D calorically perfect gas flows.
  SUBROUTINE hlle_1D (nx, left_data, right_data, u_left, u_right, f) 

    USE mod_general_data,     ONLY: nb_eq, pos_u_cell, pos_pres_cell, pos_c_cell
    USE mod_numerics_data,    ONLY: fl, fr, deltau, phys_data_Roe

    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8) :: tmp1, tmp2
    REAL(KIND=8) :: c, u, un, diff
    REAL(KIND=8) :: cl, cr, pl, pr, ul, ur, unl, unr
    REAL(KIND=8) :: b_plus, b_minus, bl, br

    REAL(KIND=8), INTENT(IN) :: nx                        !< normal to the cell interface
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: f          !< numerical flux

    ! Interface for external subroutine roe_avgState_1D 
    INTERFACE 
      SUBROUTINE roe_avgState_1D (left_data, right_data, roe_data) 
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data, right_data
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: roe_data
      END SUBROUTINE roe_avgState_1D
    END INTERFACE

    ! Physical data of left and right states
    ul = left_data(pos_u_cell)
    pl = left_data(pos_pres_cell)
    cl = left_data(pos_c_cell)

    ur = right_data(pos_u_cell)
    pr = right_data(pos_pres_cell)
    cr = right_data(pos_c_cell)

    ! Normal velocities of left and right states
    unl = ul*nx
    unr = ur*nx 

    ! Conservative variable jump (ur - ul) and inviscid physical fluxes of left and right states
    DO i = 1,nb_eq 
       tmp1 = u_left(i) 
       tmp2 = u_right(i)
       deltau(i) = tmp2 - tmp1
       fl(i)     = tmp1*unl
       fr(i)     = tmp2*unr
    ENDDO

    ! Adding the pressure contribution to the global momentum and energy equations
    fl(2) = fl(2) + pl*nx
    fl(3) = fl(3) + pl*unl

    fr(2) = fr(2) + pr*nx
    fr(3) = fr(3) + pr*unr

    ! Roe's averaged state
    CALL roe_avgState_1D (left_data, right_data, phys_data_Roe)  

    u  = phys_data_Roe(pos_u_cell)
    c  = phys_data_Roe(pos_c_cell)
    un = u*nx

    ! Wave speeds for the HLLE approximate Riemann solver 
    ! (see the original paper of Einfeldt for more details - JCP vol. 92 1991 pp 273)
    bl = MIN(un - c,unl - cl)
    br = MAX(un + c,unr + cr)
    b_plus  = MAX(br,0.d0)
    b_minus = MIN(bl,0.d0)
    diff    = 1.d0/(b_plus - b_minus)

    ! Numerical flux 
    f = (b_plus*fl - b_minus*fr + b_plus*b_minus*deltau)*diff

  END SUBROUTINE hlle_1D
!------------------------------------------------------------------------------!
