!------------------------------------------------------------------------------!
!> This subroutine computes the convective flux by means of the approximate Riemann solver 
!! proposed by Candler by using a Modified Steger-Warming with pressure sensor
!! for 1D stagnation line nonequilibrium gas flows.    
!#define DEBUG
  SUBROUTINE modified_SW_neq_1D_SL (nx, vol_l, vol_r, left_data, right_data, u_left, u_right, f)

#include"../config.h"

 
    USE mod_general_data,             ONLY: nb_ns, nb_eq, nb_temp, pos_u_cell, pos_pres_cell, pos_c_cell,   & 
                                          & pos_rhou, pos_rhov, pos_rhoE, nb_prop
    USE mod_numerics_data,            ONLY: deltau, lambda, left_eig, right_eig, absA_eig, diss, cons_Roe,  & 
                                          & phys_data_Roe
    USE mod_function_pointer,         ONLY: eigsys_neq_1D_SL
    USE mod_neq_1D_SL 

    IMPLICIT NONE
    
    INTEGER :: j, i, k
    REAL(KIND=8), PARAMETER :: eps_zero = 0.0D0
    REAL(KIND=8) :: tmp1
    REAL(KIND=8) :: cl, cr, pl, pr, ul, ur, unl, unr 
    REAL(KIND=8) :: average_speed_sound, eps, omega_interface
    REAL(KIND=8), DIMENSION(nb_eq) :: U_minus, U_plus
    REAL(KIND=8), DIMENSION(nb_eq) :: lambda_minus, lambda_plus, lambda_positive, lambda_negative, f_P, f_M
    REAL(KIND=8), DIMENSION(nb_eq, nb_eq) :: absA_eig_positive, absA_eig_negative 
    REAL(KIND=8), DIMENSION(nb_eq) :: average_u
    REAL(KIND=8), DIMENSION(nb_prop) :: data_plus, data_minus
    
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
    
    CALL pressure_weight_function(pr, pl, omega_interface) 
  
    U_plus = (1.D0- omega_interface)* u_right +  omega_interface*u_left
    data_plus = (1.D0- omega_interface)* right_data +  omega_interface*left_data
   
   ! To compute A+ = R Lambda+ L
    CALL eigsys_neq_1D_SL (nx,U_plus,data_plus, lambda, right_eig, left_eig)

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


   DO i = 1,nb_eq
    f_P(i) = 0.0D0
   ENDDO    

    DO j = 1,nb_eq 
      tmp1 = u_left(j)
       DO i = 1,nb_eq 
       f_P(i) =f_P(i)+ absA_eig_positive(i, j)*tmp1
        ENDDO
    ENDDO

! To compute A- = R Lambda- L
    lambda = 0.D0
    right_eig =0.D0
    left_eig = 0.D0  
    
    U_minus = (1.D0 - omega_interface)* u_left + omega_interface* u_right
    data_minus = (1.D0 - omega_interface)* left_data + omega_interface* right_data
    
    CALL eigsys_neq_1D_SL (nx,U_minus,data_minus, lambda, right_eig, left_eig)

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

    DO i = 1,nb_eq
    f_M(i)= 0.0D0
    ENDDO

    DO j = 1,nb_eq 
       tmp1 = u_right(j)
       DO i = 1,nb_eq 
       f_M(i) = f_M(i)+ absA_eig_negative(i, j)*tmp1  
        ENDDO
    ENDDO
    
  
     DO i = 1,nb_eq
       f(i) = f_P(i)+f_M(i)
     ENDDO

!      DO i = 1,nb_eq
!      WRITE(*,*) f(i), f_P(i), f_M(i)
!     ENDDO


  END SUBROUTINE modified_SW_neq_1D_SL
!------------------------------------------------------------------------------!
