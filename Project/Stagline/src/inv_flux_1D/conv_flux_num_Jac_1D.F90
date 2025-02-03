!------------------------------------------------------------------------------!
!> This subroutine computes the left and right state Jacobian of the selected numerical 
!! flux by means of numerical differentiation for 1D flows.
  SUBROUTINE conv_flux_num_Jac_1D (nx, left_data, right_data, u_left, u_right, jfl, jfr)

    USE mod_general_data,              ONLY: nb_eq, nb_prop, eta
    USE mod_function_pointer,          ONLY: conv_flux_1D, get_phys_from_cons 

    IMPLICIT NONE

    INTEGER :: i, j 
    REAL(KIND=8) :: dul, dur, ov_dul, ov_dur, ul, ur
    REAL(KIND=8), DIMENSION(nb_eq) :: u_left_pert, u_right_pert, f, fl_pert, fr_pert
    REAL(KIND=8), DIMENSION(nb_prop) :: left_data_pert, right_data_pert

    REAL(KIND=8), INTENT(IN) :: nx                        !< normal to the cell interface
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfl      !< numerical flux Jacobian with respect to the left state 
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfr      !< numerical flux Jacobian with respect to the right state

    ! Inviscid flux
    CALL conv_flux_1D(nx, left_data, right_data, u_left, u_right, f)

    ! Loop for left and righ state flux Jacobians by means of numerical differentiation 
    DO j = 1,nb_eq 

       ! Perturb left and right state conservative variables 
       DO i = 1,nb_eq 
          u_left_pert(i)  = u_left(i)
          u_right_pert(i) = u_right(i)
       ENDDO

       ul  = u_left(j)
       dul = eta*MAX(ABS(ul),1.d-20)*SIGN(1.d0,ul) 
       u_left_pert(j) = ul + dul

       ur  = u_right(j)
       dur = eta*MAX(ABS(ur),1.d-20)*SIGN(1.d0,ur) 
       u_right_pert(j) = ur + dur

       ! Get the physical properties corresponding to the perturbed left and righ state conservative variables
       CALL get_phys_from_cons(u_left_pert, left_data_pert) 
       CALL get_phys_from_cons(u_right_pert, right_data_pert) 
 
       ! Get the inviscid flux corresponding to the perturbed left and right state conservative variables
       CALL conv_flux_1D(nx, left_data_pert, right_data, u_left_pert, u_right, fl_pert)
       CALL conv_flux_1D(nx, left_data, right_data_pert, u_left, u_right_pert, fr_pert)

       ov_dul = 1.d0/dul
       ov_dur = 1.d0/dur
       DO i = 1,nb_eq 
          jfl(i,j) = (fl_pert(i) - f(i))*ov_dul
          jfr(i,j) = (fr_pert(i) - f(i))*ov_dur
       ENDDO
      
    ENDDO

  END SUBROUTINE conv_flux_num_Jac_1D
!------------------------------------------------------------------------------!
