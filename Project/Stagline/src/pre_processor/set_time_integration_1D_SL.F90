!-----------------------------------------------------------------------------!
!> This subroutine associates a subroutine pointer according to the time-integration scheme specified in the input file 
!! for 1D stagnation line flows. 
  SUBROUTINE set_time_integration_1D_SL ()

#include "../config.h"

    USE mod_function_pointer
    USE mod_numerics_data,        ONLY: time_disc

    IMPLICIT NONE

    SELECT CASE(time_disc)

      ! Fully explicit time-integration scheme
      CASE('fe')
        evolve_solution => fe_1D_SL

  ! Source-implicit time-integration scheme 
      CASE('si')
        evolve_solution => si_1D_SL
 
      ! Fully implicit time-integration scheme 
      CASE('fi') 
        evolve_solution => fi_1D_SL     

      CASE DEFAULT 
        WRITE(*,10)'In "set_time_integration_1D_SL.F90", error in the selection of time-integration method...'
        PRINT* 
        STOP

    END SELECT 

10 FORMAT(A)

  END SUBROUTINE set_time_integration_1D_SL  
!------------------------------------------------------------------------------!
