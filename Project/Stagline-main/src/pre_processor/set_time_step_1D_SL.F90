!------------------------------------------------------------------------------!
!> This subroutine associates the pointers for the subroutine computing the time-step for 1D stagnation line flows.
  SUBROUTINE set_time_step_1D_SL

    USE mod_general_data,        ONLY: flag_diss
    USE mod_numerics_data,       ONLY: time_step, cfl
    USE mod_function_pointer

    IMPLICIT NONE

    ! Time-accurate solution (time-step specified in input)
    IF (time_step.NE.0.d0) THEN

       compute_time_step_1D_SL => time_step_time_acc_1D_SL 

    ! Time-step computed based on the CFL number
    ELSE

       IF (cfl.EQ.0.d0) THEN
          WRITE(*,10)'In "set_time_step_1D_SL.F90", CFL number value not set!'
          WRITE(*,10)'Please change the input file'
          STOP
       ENDIF

       ! Inviscid flow
       IF (flag_diss.EQV..FALSE.) THEN

          compute_time_step_1D_SL => inv_time_step_1D_SL

       ! Viscous flow
       ELSEIF (flag_diss.EQV..TRUE.) THEN

          compute_time_step_1D_SL => inv_visc_time_step_1D_SL

       ENDIF

    ENDIF

10 FORMAT(A)

  END SUBROUTINE set_time_step_1D_SL
!------------------------------------------------------------------------------!
