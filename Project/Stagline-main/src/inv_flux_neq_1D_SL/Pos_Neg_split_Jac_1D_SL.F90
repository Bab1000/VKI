!------------------------------------------------------------------------------!
!> This subrotine provides the numerical flux Jacobian for 1D stagnation line nonequilibrium gas flows
!! according to the positive and negative split of the inviscid flux Jacobian.
!! In this case the flow is characterized by N temperatures plus a separated temperature Te for free electrons.
  SUBROUTINE Pos_Neg_split_Jac_general_1D_SL (nx, vol_l, vol_r, left_data, right_data, u_left, u_right, jfl, jfr)

    USE mod_general_data, ONLY: nb_ns, pos_ei_cell, pos_h0_cell, pos_rho_cell, pos_T_cell, pos_u_cell,  &
                              & pos_v_cell, yi, Ri, nb_temp, pos_pres_cell, pos_beta_cell, nb_eq, pos_c_cell
    USE mod_function_pointer,         ONLY: eigsys_neq_1D_SL
    USE mod_numerics_data,            ONLY: deltau, lambda, left_eig, right_eig, absA_eig

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: nx                        !< normal to the cell interface
    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8) :: cl, cr
    REAL(KIND=8) :: tmp1
    REAL(KIND=8), PARAMETER :: eps_zero = 0.0D0
    REAL(KIND=8) :: eps
    REAL(KIND=8), DIMENSION(nb_eq, nb_eq) :: absA_eig_positive, absA_eig_negative
    REAL(KIND=8), DIMENSION(nb_eq) :: lambda_minus, lambda_plus, lambda_positive, lambda_negative 
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfl      !< numerical flux Jacobian with respect to the left state 
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfr      !< numerical flux Jacobian with respect to the right state
    
    INTEGER :: i, j, k
    REAL(KIND=8) :: a, a2, ov_a2, u, v, rho, p, H, l1, l2, l3
    REAL(KIND=8), DIMENSION(nb_temp) :: beta, em
    REAL(KIND=8), DIMENSION(nb_ns) :: gam
    REAL(KIND=8), DIMENSION(nb_ns+2+nb_temp, nb_ns+2+nb_temp) :: L, R


        ! Speed of sound, velocity and pressure of left and right states
    cl = left_data(pos_c_cell)
    cr = right_data(pos_c_cell)
    
   ! To compute A+ = R Lambda+ L
    CALL eigsys_neq_1D_SL (nx,u_left,left_data, lambda, right_eig, left_eig)

   eps = eps_zero* cl

   DO i = 1, nb_eq
   lambda_positive(i) = 0.5D0* ( lambda(i) +  sqrt(lambda(i)**2 + eps**2))
   ENDDO      
    ! Matrix product R*|Lambda|
    DO j = 1,nb_eq 
       tmp1 = lambda_positive(j)
       DO i = 1,nb_eq 
          right_eig(i,j) = right_eig(i,j)*tmp1
       ENDDO
    ENDDO

    ! Matrix product R*|Lambda|*L
    DO j = 1,nb_eq 
       DO i = 1,nb_eq 
          absA_eig_positive(i,j) = 0.d0
       ENDDO
    ENDDO
    
    DO k = 1,nb_eq 
       DO j = 1,nb_eq
          tmp1 = left_eig(k,j)
        DO i = 1,nb_eq 
             absA_eig_positive(i,j) = absA_eig_positive(i,j) + right_eig(i,k)*tmp1
         ENDDO
       ENDDO
    ENDDO
  
    jfl = absA_eig_positive
    

! To compute A- = R Lambda- L
    lambda = 0.D0
    right_eig =0.D0
    left_eig = 0.D0  
    
    
    CALL eigsys_neq_1D_SL (nx,u_right,right_data, lambda, right_eig, left_eig)

    eps = eps_zero* cr 
   
   DO i = 1, nb_eq
   lambda_negative(i) = 0.5D0* ( lambda(i) -  sqrt(lambda(i)**2+eps**2))
   ENDDO
    
    ! Matrix product R*|Lambda|
    DO j = 1,nb_eq 
       tmp1 = lambda_negative(j)
       DO i = 1,nb_eq 
          right_eig(i,j) = right_eig(i,j)*tmp1
       ENDDO
    ENDDO

    ! Matrix product R*|Lambda|*L
    DO j = 1,nb_eq 
       DO i = 1,nb_eq 
          absA_eig_negative(i,j) = 0.d0
       ENDDO
    ENDDO
    
    DO k = 1,nb_eq 
       DO j = 1,nb_eq
          tmp1 = left_eig(k,j)
        DO i = 1,nb_eq 
             absA_eig_negative(i,j) = absA_eig_negative(i,j) + right_eig(i,k)*tmp1
         ENDDO
       ENDDO
    ENDDO

   jfr = absA_eig_negative

 
  END SUBROUTINE Pos_Neg_split_Jac_general_1D_SL
!------------------------------------------------------------------------------!
