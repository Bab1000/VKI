!------------------------------------------------------------------------------!
!> This subroutine sets the diffusive flux to zero when solving the 1D Euler equations.
  SUBROUTINE null_ns_flux_1D (vol_l, vol_r, left_data, right_data, u_left, u_right, fd)

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: fd         !< diffusive flux

    ! Diffusive flux components set to zero 
    fd = 0.d0
 
  END SUBROUTINE null_ns_flux_1D 
!------------------------------------------------------------------------------!

