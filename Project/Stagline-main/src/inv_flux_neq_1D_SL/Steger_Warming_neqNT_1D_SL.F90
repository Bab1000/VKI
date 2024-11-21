!------------------------------------------------------------------------------!
!> This subroutine computes the numerical flux according to the Steger-Warming flux vector splitting 
!! for 1D stagnation line nonequilibrium gas flows. In this case the flow is characterized by N temperatures. 
!! However, no separated temperature exist for the free electrons.
  SUBROUTINE Steger_Warming_neqNT_1D_SL (nx, vol_l, vol_r, left_data, right_data, u_left, u_right, f) 

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: nx                        !< normal to the cell interface
    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: f          !< numerical flux

    f = 0.d0

    PRINT*
    WRITE(*,'(A)')'In "Steger_Warming_neqNT_1D_SL.F90", not implemented yet...'
    PRINT*
    STOP

  END SUBROUTINE Steger_Warming_neqNT_1D_SL
!------------------------------------------------------------------------------!
