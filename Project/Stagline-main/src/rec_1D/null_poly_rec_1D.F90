!------------------------------------------------------------------------------!
!> This subroutine applies a null 1D polynomial reconstruction for 1D calorically perfect gas flows.
!! It is used for first order space accurate solutions.
  SUBROUTINE null_poly_rec_1D (ull, ul, ur, urr, prop_ll, prop_l, prop_r, prop_rr, & 
                            &  u_left, u_right, prop_left, prop_right)

    IMPLICIT NONE
  
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: ull         !< conservative variables of left-left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: ul          !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: ur          !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: urr         !< conservative variables of right-right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: prop_ll     !< physical properties of left-left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: prop_l      !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: prop_r      !< physical properties of right state     
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: prop_rr     !< physical properties of right-right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_left      !< conservative variables of reconstructed left state 
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_right     !< conservative variables of reconstructed right state 
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: prop_left   !< physical properties of reconstructed left state 
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: prop_right  !< physical properties of reconstructed right state  
  
    ! Left state 
    u_left = ul
    prop_left = prop_l

    ! Right state 
    u_right = ur
    prop_right = prop_r

  END SUBROUTINE null_poly_rec_1D
!------------------------------------------------------------------------------!
