!------------------------------------------------------------------------------!
! This subroutine computes the diffusive source term and the related Jacobians by means of numerical differentiation 
! for 1D stagnation line flows (calorically perfect gas or nonequilibriumf flows).
  SUBROUTINE source_term_diff_num_Jac_SL (r_c, vol_l, vol_c, vol_r, prop_left, prop_cell, prop_right, & 
                                        & u_left, u_cell, u_right, s, jsl, jsc, jsr)

    USE mod_general_data,                 ONLY: nb_eq, nb_prop, eta
    USE mod_function_pointer,             ONLY: get_diff_source_term_SL, get_phys_from_cons

    IMPLICIT NONE

    INTEGER :: i, j
    REAL(KIND=8) :: tmp
    REAL(KIND=8) :: ul, uc, ur, dul, duc, dur, ov_dul, ov_duc, ov_dur
    REAL(KIND=8), DIMENSION(nb_eq) :: sl, sc, sr
    REAL(KIND=8), DIMENSION(nb_eq) :: u_left_pert, u_cell_pert, u_right_pert
    REAL(KIND=8), DIMENSION(nb_prop) :: prop_left_pert, prop_cell_pert, prop_right_pert

    REAL(KIND=8), INTENT(IN) :: r_c
    REAL(KIND=8), INTENT(IN) :: vol_l, vol_c, vol_r
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left, u_cell, u_right
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prop_left, prop_cell, prop_right
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jsl, jsc, jsr
   
    ! Source term
    CALL get_diff_source_term_SL(r_c, vol_l, vol_c, vol_r, prop_left, prop_cell, prop_right, & 
                               & u_left, u_cell, u_right, s)   

    ! Loop for left, center and right state flux Jacobians by means of numerical differentiation
    DO j = 1,nb_eq

       ! Perturb left, center and right state conservative variables 
       DO i = 1,nb_eq 
          u_left_pert(i)  = u_left(i)
          u_cell_pert(i)  = u_cell(i)
          u_right_pert(i) = u_right(i)
       ENDDO

       ul  = u_left(j)
       dul = eta*MAX(ABS(ul),1.d-20)*SIGN(1.d0,ul) 
       u_left_pert(j) = ul + dul

       uc  = u_cell(j)
       duc = eta*MAX(ABS(uc),1.d-20)*SIGN(1.d0,uc) 
       u_cell_pert(j) = uc + duc 

       ur  = u_right(j)
       dur = eta*MAX(ABS(ur),1.d-20)*SIGN(1.d0,ur) 
       u_right_pert(j) = ur + dur

       ! Get the physical properties corresponding to the perturbed left, center and righ state conservative variables
       CALL get_phys_from_cons(u_left_pert, prop_left_pert) 
       CALL get_phys_from_cons(u_cell_pert, prop_cell_pert)
       CALL get_phys_from_cons(u_right_pert, prop_right_pert) 

       ! Get the diffusive source term corresponding to the perturbed left, center and right state conservative variables
       CALL get_diff_source_term_SL(r_c, vol_l, vol_c, vol_r, prop_left_pert, prop_cell, prop_right, & 
                                  & u_left_pert, u_cell, u_right, sl)
       CALL get_diff_source_term_SL(r_c, vol_l, vol_c, vol_r, prop_left, prop_cell_pert, prop_right, & 
                                  & u_left, u_cell_pert, u_right, sc)
       CALL get_diff_source_term_SL(r_c, vol_l, vol_c, vol_r, prop_left, prop_cell, prop_right_pert, & 
                                  & u_left, u_cell, u_right_pert, sr) 

       ov_dul = 1.d0/dul 
       ov_duc = 1.d0/duc
       ov_dur = 1.d0/dur
       DO i = 1,nb_eq
          tmp = s(i) 
          jsl(i,j) = (sl(i) - tmp)*ov_dul
          jsc(i,j) = (sc(i) - tmp)*ov_duc
          jsr(i,j) = (sr(i) - tmp)*ov_dur
       ENDDO

    ENDDO   

  END SUBROUTINE source_term_diff_num_Jac_SL
!------------------------------------------------------------------------------!
