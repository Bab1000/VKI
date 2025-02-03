!------------------------------------------------------------------------------!
!> This subroutine computes the diffusive source term and the related Jacobians for 1D stagnation line flows 
!! for nonequilibrium flows (sphere case - N temperatures).
  SUBROUTINE source_term_diff_neqNT_1D_SL_sph_Jac(r_c, vol_l, vol_c, vol_r, prop_left, prop_cell, prop_right, & 
                                              & u_left, u_cell, u_right, s, js_left, js_cell, js_right)

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: r_c                        !< radial location of cell centroid
    REAL(KIND=8), INTENT(IN) :: vol_l                      !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_c                      !< volume of central state
    REAL(KIND=8), INTENT(IN) :: vol_r                      !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left       !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_cell       !< conservative variables of central state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right      !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prop_left    !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prop_cell    !< physical properties of central state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prop_right   !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s           !< source term
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: js_left   !< source term Jacobian with respect to the left state
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: js_cell   !< source term Jacobian with respect to the central state
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: js_right  !< source term Jacobian with respect to the right state 

    s = 0.d0
    js_left  = 0.d0
    js_cell  = 0.d0
    js_right = 0.d0
    
    PRINT*,'In "source_term_diff_neqNT_1D_SL_sph_Jac.F90",implementation to be finshed!'
    STOP

  END SUBROUTINE source_term_diff_neqNT_1D_SL_sph_Jac
!------------------------------------------------------------------------------!
