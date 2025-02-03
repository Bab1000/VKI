!------------------------------------------------------------------------------!
!> This subroutine computes the diffusive flux and diffusive flux Jacobians for the 1D Navier-Stokes equations 
!! for nonequilibrium flows when using an NT_Te model (diffusive flux Jacobians are provided in output). 
  SUBROUTINE ns_flux_neqNT_Te_1D_Jac (vol_l, vol_r, left_data, right_data, u_left, u_right, fd, jfdl, jfdr)

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: fd         !< diffusive flux
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfdl     !< diffusive flux Jacobian with respect to the left state
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfdr     !< diffusive flux Jacobian with respect to the right state

    fd = 0.d0
    jfdl = 0.d0
    jfdr = 0.d0

    PRINT*,'in ns_flux_neqNT_Te_1D_Jac.F90'
    PRINT*,'not implemented yet'
    STOP

  END SUBROUTINE ns_flux_neqNT_Te_1D_Jac 
!------------------------------------------------------------------------------!
