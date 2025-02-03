!------------------------------------------------------------------------------!
!> This subroutine sets the limiter for polynomial re-construction (2nd order accuracy in space)
  SUBROUTINE set_limiter ()

    USE mod_numerics_data,          ONLY: fun_limiter
    USE mod_function_pointer

    IMPLICIT NONE

    ! Limiter selection
    SELECT CASE (fun_limiter)

      CASE ('van_leer')
        limiter => van_leer

      CASE ('minimod')
         limiter => minimod

      CASE ('superbee')
        limiter => superbee

      CASE ('koren')
        limiter => koren

      CASE ('mc')
        limiter => mc

      CASE ('charm')
        limiter => charm

      CASE ('hcus')
        limiter => hcus

      CASE ('hquick')
        limiter => hquick

      CASE ('van_albada1')
        limiter => van_albada1

      CASE ('van_albada2')
        limiter => van_albada2

      CASE ('ospre')
        limiter => ospre

      CASE ('none')
        limiter => no_limiter

      CASE DEFAULT 
        WRITE(*,10)'In set_limiter.F90, error in limiter selection'
        PRINT*
        STOP 

    END SELECT 

10 FORMAT(A)

   END SUBROUTINE set_limiter
!------------------------------------------------------------------------------!
