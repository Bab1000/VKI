!------------------------------------------------------------------------------!
!> This subroutine sets the time-step to the value specified in input by the user 
!! for time-accurate solutions
   SUBROUTINE time_step_time_acc_1D_SL (vol, phys_prop, dt)

    USE mod_numerics_data,         ONLY: time_step

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: vol                        !< cell volume
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_prop    !< physical properties
    REAL(KIND=8), INTENT(OUT) :: dt                        !< time-step

    ! Time-step
    dt = time_step

   END SUBROUTINE time_step_time_acc_1D_SL
!------------------------------------------------------------------------------!
!> This subroutine computes the time-step for 1D stagnation line inviscid flows.
  SUBROUTINE inv_time_step_1D_SL (vol, phys_prop, dt)

    USE mod_general_data,          ONLY: pos_u_cell, pos_c_cell
    USE mod_numerics_data,         ONLY: cfl

    IMPLICIT NONE

    REAL(KIND=8) :: u, c

    REAL(KIND=8), INTENT(IN) :: vol                        !< cell volume
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_prop    !< physical properties
    REAL(KIND=8), INTENT(OUT) :: dt                        !< time-step

    u = phys_prop(pos_u_cell)
    c = phys_prop(pos_c_cell)

    ! Time-step
    dt = cfl*vol/(ABS(u) + c)

  END SUBROUTINE inv_time_step_1D_SL 
!------------------------------------------------------------------------------!
!> This subroutine computes the time-step for 1D stagnation line viscous flows.
  SUBROUTINE visc_time_step_1D_SL (vol, phys_prop, dt)

    USE mod_numerics_data,         ONLY: vnn
    USE mod_function_pointer,      ONLY: get_visc_spectral_radius_1D_SL

    IMPLICIT NONE

    REAL(KIND=8) :: sp

    REAL(KIND=8), INTENT(IN) :: vol                        !< cell volume
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_prop    !< physical properties
    REAL(KIND=8), INTENT(OUT) :: dt                        !< time-step

    ! Compute the diffusive flux Jacobian spectral radius
    CALL get_visc_spectral_radius_1D_SL (phys_prop, sp)

    ! Time-step
    dt = 0.5d0*vnn*vol**2/sp

  END SUBROUTINE visc_time_step_1D_SL
!------------------------------------------------------------------------------!
!> This subroutine computes the time-step for 1D stagnation line viscous flows by combining the 
!! contributions of the inviscid and viscous spectral radii. 
  SUBROUTINE inv_visc_time_step_1D_SL (vol, phys_prop, dt)

    USE mod_general_data,          ONLY: pos_u_cell, pos_c_cell
    USE mod_numerics_data,         ONLY: cfl
    USE mod_function_pointer,      ONLY: get_visc_spectral_radius_1D_SL

    IMPLICIT NONE

    REAL(KIND=8) :: u, c
    REAL(KIND=8) :: sp_inv, sp_visc

    REAL(KIND=8), INTENT(IN) :: vol                        !< cell volume
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_prop    !< physical properties
    REAL(KIND=8), INTENT(OUT) :: dt                        !< time-step

    ! Velocity and speed of sound
    u = phys_prop(pos_u_cell)
    c = phys_prop(pos_c_cell)

    ! Inviscid spectral radius
    sp_inv = ABS(u) + c

    ! Viscous spectral radius
    CALL get_visc_spectral_radius_1D_SL (phys_prop, sp_visc)

    ! Time-step    
    dt = cfl*vol/(sp_inv + sp_visc/vol)

  END SUBROUTINE inv_visc_time_step_1D_SL
!------------------------------------------------------------------------------!
