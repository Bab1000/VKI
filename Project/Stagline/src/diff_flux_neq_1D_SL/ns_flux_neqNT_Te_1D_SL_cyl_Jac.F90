!------------------------------------------------------------------------------!
!> This subroutine computes the diffusive flux and the related Jacobians when solving the 1D stagnation line equations 
!! for nonequilibrium flows (cylinder case - N temperatures and one temperature for free electrons Te).
  SUBROUTINE ns_flux_neqNT_Te_1D_SL_cyl_Jac (r_l, r_r, vol_l, vol_r, left_data, right_data, u_left, u_right, fd, & 
                                           & jfdl, jfdr)

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
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfdl     !< diffusive flux Jacobian with respect to the left state
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfdr     !< diffusive flux Jacobian with respect to the right state

    ! Diffusive flux and Jacobians set to zero 
    fd = 0.d0
    jfdl = 0.d0
    jfdr = 0.d0
 
    PRINT*,'In "ns_flux_neqNT_Te_1D_SL_cyl_Jac.F90", implementation to be finshed'
    STOP

  END SUBROUTINE ns_flux_neqNT_Te_1D_SL_cyl_Jac
!------------------------------------------------------------------------------!

