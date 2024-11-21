!------------------------------------------------------------------------------!
! This subroutine computes the source term and the related Jacobian by means of numerical differentiation 
! for calorically perfect gas flows
  SUBROUTINE source_term_num_Jac (source_id, cell_id, cell_data, cons, s, js)

    USE mod_general_data,                    ONLY: nb_eq, nb_prop, eta
    USE mod_function_pointer,                ONLY: get_source_term, get_phys_from_cons

    INTEGER :: i, j
    REAL(KIND=8) :: du, ov_du, u
    REAL(KIND=8), DIMENSION(nb_eq) :: s_pert, cons_pert
    REAL(KIND=8), DIMENSION(nb_prop) :: cell_data_pert

    INTEGER, INTENT(IN) :: source_id, cell_id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cell_data, cons
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: js

    ! Source term 
    CALL get_source_term(source_id)%source(cell_id, cell_data, cons, s)

    ! Loop for source term Jacobian by means of numerical differentiation
    DO j = 1,nb_eq

       ! Perturb conservative variables 
       DO i = 1,nb_eq 
          cons_pert(i) = cons(i)
       ENDDO

       u  = cons(j)
       du = eta*MAX(ABS(u),1.d-20)*SIGN(1.d0,u)
    
       cons_pert(j) = u + du

       ! Get the physical properties corresponding to the perturbed conservative variables
       CALL get_phys_from_cons(cons_pert, cell_data_pert) 

       ! Get the source term corresponding to the perturbed conservative variables
       CALL get_source_term(source_id)%source(cell_id, cell_data_pert, cons_pert, s_pert)

       ! Numerical differentiation
       ov_du = 1.d0/du
       DO i = 1,nb_eq
          js(i,j) = (s_pert(i) - s(i))*ov_du
       ENDDO

    ENDDO

  END SUBROUTINE source_term_num_Jac
!------------------------------------------------------------------------------!
