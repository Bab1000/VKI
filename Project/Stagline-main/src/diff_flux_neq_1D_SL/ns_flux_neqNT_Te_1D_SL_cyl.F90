!------------------------------------------------------------------------------!
!> This subroutine computes the diffusive flux when solving the 1D stagnation line equations 
!! for nonequilibrium flows (cylinder case - N temperatures and one temperature for free electrons Te).
  SUBROUTINE ns_flux_neqNT_Te_1D_SL_cyl (r_l, r_r, vol_l, vol_r, left_data, right_data, u_left, u_right, fd)

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: r_l                       !< radial position of left state
    REAL(KIND=8), INTENT(IN) :: r_r                       !< radial position of rigth state
    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state 
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: fd         !< diffusive flux

    ! Diffusive flux set to zero 
    fd = 0.d0
 
    PRINT*,'In "ns_flux_neqNT_Te_1D_SL_cyl.F90", implementation to be finshed!'
    STOP

  END SUBROUTINE ns_flux_neqNT_Te_1D_SL_cyl
!------------------------------------------------------------------------------!

