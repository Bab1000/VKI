!------------------------------------------------------------------------------!
!> This subroutine computes the convective flux by means of the approximate Riemann solver 
!! proposed by Roe for a 1D stagnation line calorically perfect gas flows. The Harten-Hyman entropy fix is used
!! in order to remove unphysical expansion shocks.
  SUBROUTINE roe_1D_SL (nx, vol_l, vol_r, left_data, right_data, u_left, u_right, f)

    USE mod_general_data,         ONLY: pos_u_cell, pos_v_cell, pos_pres_cell, pos_c_cell, pos_rhou, pos_rhoe, nb_eq
    USE mod_numerics_data,        ONLY: deltau, lambda, right_eig, left_eig, absA_eig, diss, phys_data_Roe

    IMPLICIT NONE

    INTEGER :: i, j, k
    REAL(KIND=8) :: tmp1, tmp2, tmp3
    REAL(KIND=8) :: eig, eps  
    REAL(KIND=8) :: cl, cr, pl, pr, ul, ur, unl, unr

    REAL(KIND=8), INTENT(IN) :: nx                        !< normal to the cell interface
    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: f          !< numerical flux

    ! Explicit interfaces 
    INTERFACE
      ! Subroutine eigensystem_1D_SL 
      SUBROUTINE eigensystem_1D_SL(nx, phys_data, lambda, right, left)
        IMPLICIT NONE
        REAL(KIND=8), INTENT(IN) :: nx 
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: lambda
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: right, left
      END SUBROUTINE eigensystem_1D_SL

      ! Subroutine roe_avgState_1D_SL 
      SUBROUTINE roe_avgState_1D_SL(left_data, right_data, roe_data) 
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data, right_data
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: roe_data
      END SUBROUTINE roe_avgState_1D_SL
    END INTERFACE

    ! Data of left and right states
    ul = left_data(pos_u_cell)
    pl = left_data(pos_pres_cell)
    cl = left_data(pos_c_cell)

    ur = right_data(pos_u_cell)
    pr = right_data(pos_pres_cell)
    cr = right_data(pos_c_cell)
  
    unl = ul*nx
    unr = ur*nx 

    ! Conservative variable jump (ur - ul)
    DO i = 1,nb_eq 
       deltau(i) = u_right(i) - u_left(i)
    ENDDO

    ! Average flux 
    DO i = 1,nb_eq
       f(i) = 0.5d0*(u_left(i)*unl + u_right(i)*unr)
    ENDDO

    ! Adding the pressure contribution to the radial momentum and energy equations
    f(pos_rhou) = f(pos_rhou) + 0.5d0*(pl + pr)*nx
    f(pos_rhoE) = f(pos_rhoE) + 0.5d0*(pl*unl + pr*unr)

    ! Roe's averaged state
    CALL roe_avgState_1D_SL (left_data, right_data, phys_data_Roe) 
 
    ! Eigensystem 
    CALL eigensystem_1D_SL (nx, phys_data_Roe, lambda, right_eig, left_eig)

    ! Entropy fix (Harten-Hyman - 1983)
    ! First eigenvalue
    eig  = lambda(1)
    tmp1 = eig - unl
    tmp2 = unr - eig
    tmp3 = ABS(eig)
    eps  = MAX(0.d0,tmp1,tmp2)

    IF (tmp3.LT.eps) THEN
       lambda(1) = eps
    ELSE 
       lambda(1) = tmp3
    ENDIF 

    ! Second eigenvalue
    eig  = lambda(2)
    tmp1 = eig   - (unl - cl)
    tmp2 = (unr - cr) - eig
    tmp3 = ABS(eig)
    eps  = MAX(0.d0,tmp1,tmp2)

    IF (tmp3.LT.eps) THEN
       lambda(2) = eps
    ELSE 
       lambda(2) = tmp3
    ENDIF 

    ! Third eigenvalue
    lambda(3) = lambda(1)

    ! Fourth eigenvalue
    eig  = lambda(4)
    tmp1 = eig   - (unl + cl)
    tmp2 = (unr + cr) - eig
    tmp3 = ABS(eig)
    eps  = MAX(0.d0,tmp1,tmp2)

    IF (tmp3.LT.eps) THEN
       lambda(4) = eps
    ELSE 
       lambda(4) = tmp3
    ENDIF 
    
    ! Dissipation term
    ! Matrix product R*|Lambda|
    DO j = 1,nb_eq 
       tmp1 = lambda(j)
       DO i = 1,nb_eq 
          right_eig(i,j) = right_eig(i,j)*tmp1
       ENDDO
    ENDDO

    ! Matrix product R*|Lambda|*L
    DO j = 1,nb_eq 
       DO i = 1,nb_eq 
          absA_eig(i,j) = 0.d0
       ENDDO
    ENDDO
    
    DO k = 1,nb_eq 
       DO j = 1,nb_eq
          tmp1 = left_eig(k,j)
          DO i = 1,nb_eq 
             absA_eig(i,j) = absA_eig(i,j) + right_eig(i,k)*tmp1 
          ENDDO
       ENDDO
    ENDDO

    ! Dissipation term D = 0.5*R*|Lambda|*L*DU
    diss = 0.d0
    DO j = 1,nb_eq 
       tmp1 = deltau(j)
       DO i = 1,nb_eq 
          diss(i) = diss(i) + absA_eig(i,j)*tmp1
       ENDDO
    ENDDO  

    ! Numerical flux 
    f = f - 0.5d0*diss
     
  END SUBROUTINE roe_1D_SL
!------------------------------------------------------------------------------!
