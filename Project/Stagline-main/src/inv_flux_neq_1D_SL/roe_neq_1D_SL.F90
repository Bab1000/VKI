!------------------------------------------------------------------------------!
!> This subroutine computes the convective flux by means of the approximate Riemann solver 
!! proposed by Roe for 1D stagnation line nonequilibrium gas flows.    
  SUBROUTINE roe_neq_1D_SL (nx, vol_l, vol_r, left_data, right_data, u_left, u_right, f)

#include"../config.h"
 
    USE mod_general_data,             ONLY: nb_ns, nb_eq, nb_temp, pos_u_cell, pos_pres_cell, pos_c_cell,   & 
                                          & pos_rhou, pos_rhov, pos_rhoE
    USE mod_numerics_data,            ONLY: deltau, lambda, left_eig, right_eig, absA_eig, diss, cons_Roe,  & 
                                          & phys_data_Roe
    USE mod_function_pointer,         ONLY: eigsys_neq_1D_SL, roe_avg_neq_1D_SL
    USE mod_neq_1D_SL 

    IMPLICIT NONE

    INTEGER :: i, j, k
    REAL(KIND=8) :: eps, eig
    REAL(KIND=8) :: tmp1, tmp2, tmp3
    REAL(KIND=8) :: cl, cr, pl, pr, ul, ur, unl, unr 
    
    REAL(KIND=8), INTENT(IN) :: nx                        !< normal to the cell interface
    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
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

    ! Normal velocities of left and right states
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

    ! Radial momentum and energy equations (pressure contribution is added)
    f(pos_rhou) = f(pos_rhou) + 0.5d0*(pl + pr)*nx
    f(pos_rhoE) = f(pos_rhoE) + 0.5d0*(pl*unl + pr*unr)

    ! Compute Roe's averaged state
    CALL roe_avg_neq_1D_SL (u_left, u_right, left_data, right_data, cons_Roe, phys_data_Roe)

    ! Eigensystem 
    CALL eigsys_neq_1D_SL (nx, cons_Roe, phys_data_Roe, lambda, right_eig, left_eig)

    ! Entropy fix (Harten-Hyman - 1983)
    ! Species continuity equations
    eig  = lambda(1)
    tmp1 = eig - unl
    tmp2 = unr - eig
    tmp3 = ABS(eig)
    eps  = MAX(0.d0,tmp1,tmp2)
 
    IF (tmp3.LE.eps) THEN
       DO i = 1,nb_ns 
          lambda(i) = eps
       ENDDO
    ELSE 
       DO i = 1,nb_ns 
          lambda(i) = tmp3
       ENDDO
    ENDIF

    ! Radial momentum equation
    eig  = lambda(pos_rhou)
    tmp1 = eig   - (unl - cl)
    tmp2 = (unr - cr) - eig
    tmp3 = ABS(eig)
    eps  = MAX(0.d0,tmp1,tmp2)

    IF (tmp3.LE.eps) THEN
       lambda(pos_rhou) = eps
    ELSE 
       lambda(pos_rhou) = tmp3
    ENDIF

    ! Circumferential momentum equation
    lambda(pos_rhov) = lambda(nb_ns)

    ! Global energy equation
    eig  = lambda(pos_rhoE)
    tmp1 = eig   - (unl + cl)
    tmp2 = (unr + cr) - eig
    tmp3 = ABS(eig)
    eps  = MAX(0.d0,tmp1,tmp2)

    IF (tmp3.LE.eps) THEN
       lambda(pos_rhoE) = eps
    ELSE 
       lambda(pos_rhoE) = tmp3
    ENDIF

    ! Internal energy equations and free electron entropy equation (if present)
    IF (nb_temp.GT.1) THEN

       eig  = lambda(pos_rhoE + 1)
       tmp1 = eig - unl
       tmp2 = unr - eig
       tmp3 = ABS(eig)
       eps  = MAX(0.d0,tmp1,tmp2)

       IF (tmp3.LE.eps) THEN
         DO i = 2,nb_temp
            lambda(pos_rhoE + i-1) = eps
         ENDDO
       ELSE 
         DO i = 2,nb_temp
            lambda(pos_rhoE + i-1) = tmp3
         ENDDO
       ENDIF

    ENDIF

#ifdef GOTO_BLAS
    ! Matrix product R*|Lambda|
    DO j = 1,nb_eq 
       tmp1 = lambda(j)
       DO i = 1,nb_eq 
          right_eig(i,j) = right_eig(i,j)*tmp1
       ENDDO
    ENDDO

    ! Matrix product R*|Lambda|*L
    CALL DGEMM('n', 'n', nb_eq, nb_eq, nb_eq, 1.d0, right_eig, nb_eq, left_eig, nb_eq, 0.d0, absA_eig, nb_eq) 

    ! Numerical flux F = 0.5*(FL + FR) - 0.5*R*|Lambda|*L*DU
    CALL DGEMV ('n', nb_eq, nb_eq, -0.5d0, absA_eig, nb_eq, deltau, 1, 1.d0, f, 1) 

#else

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
  
    ! Subtraction of dissipatiion term
    f = f - 0.5d0*diss
#endif

  END SUBROUTINE roe_neq_1D_SL 
!------------------------------------------------------------------------------!
