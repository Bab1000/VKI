!------------------------------------------------------------------------------!
!> This subroutine provides the stress tensor components for 1D stagnation line flow (sphere case).
  SUBROUTINE stress_tensor_1D_SL_sph (mu, r, u, v, du_dr, dv_dr, tau_rr, tau_rt, tau_tt)

    IMPLICIT NONE

    REAL(KIND=8), PARAMETER :: coeff = 4.d0/3.d0
    REAL(KIND=8) :: vel_sum_ov_r 

    REAL(KIND=8), INTENT(IN) :: mu       !< dynamic viscosity
    REAL(KIND=8), INTENT(IN) :: r        !< radial location
    REAL(KIND=8), INTENT(IN) :: u        !< radial velocity
    REAL(KIND=8), INTENT(IN) :: v        !< circumferential velocity
    REAL(KIND=8), INTENT(IN) :: du_dr    !< radial velocity gradient
    REAL(KIND=8), INTENT(IN) :: dv_dr    !< circumferential velocity gradient 
    REAL(KIND=8), INTENT(OUT) :: tau_rr  !< rr viscous stress tensor component 
    REAL(KIND=8), INTENT(OUT) :: tau_rt  !< r-theta viscous stress tensor component
    REAL(KIND=8), INTENT(OUT) :: tau_tt  !< theta-theta viscous stress tensor component

    ! Useful quantity
    vel_sum_ov_r = (u + v)/r

    ! Stres tensor components
    tau_rr = mu*coeff*(du_dr - vel_sum_ov_r)
    tau_rt = mu*(dv_dr - vel_sum_ov_r)
    tau_tt = - 0.5d0*tau_rr
  
  END SUBROUTINE stress_tensor_1D_SL_sph
!------------------------------------------------------------------------------!
