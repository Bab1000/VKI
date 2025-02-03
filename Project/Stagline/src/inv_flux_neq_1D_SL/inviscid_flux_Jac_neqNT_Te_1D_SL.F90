!------------------------------------------------------------------------------!
!> This subroutine computes the inviscid flux Jacobian for the 1D Euler stagnation line equations for nonequilibrium flows.
!! In this case the flow is characterized by N temperatures plus a separated temperature Te for free electrons.
  SUBROUTINE inviscid_flux_Jac_neqNT_Te_1D_SL (nx, cons, phys_data, a)

    IMPLICIT NONE
 
    REAL(KIND=8), INTENT(IN) :: nx                       !> normal to the cell interface
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons       !> conservative variables
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data  !> physical properties
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: a       !> inviscid flux Jacobian

    a = 0.d0

    PRINT*
    WRITE(*,'(A)')'In "inviscid_flux_Jac_neqNT_Te_1D_SL.F90", not implemented yet...'
    PRINT*
    STOP

  END SUBROUTINE inviscid_flux_Jac_neqNT_Te_1D_SL
!------------------------------------------------------------------------------!
