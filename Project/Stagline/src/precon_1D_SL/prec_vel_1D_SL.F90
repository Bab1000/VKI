!-------------------------------------------------------------------------------!
!> This subroutine computes the preconditioning velocity for 1D stagnation line inviscid flows.
   SUBROUTINE inv_prec_vel_1D_SL (vol, rho, mu, V, delta_p, c, Vp)

     IMPLICIT NONE

     REAL(KIND=8), INTENT(IN) :: vol           !< cell volume
     REAL(KIND=8), INTENT(IN) :: rho           !< density 
     REAL(KIND=8), INTENT(IN) :: mu            !< viscosity (not used)
     REAL(KIND=8), INTENT(IN) :: V             !< velocity magnitude
     REAL(KIND=8), INTENT(IN) :: delta_p       !< pressure difference 
     REAL(KIND=8), INTENT(IN) :: c             !< speed of sound
     REAL(KIND=8), INTENT(OUT) :: Vp           !< inviscid pre-conditioning velocity

     ! Inviscid pre-conditioning velocity
     Vp = MIN(MAX(V,DSQRT(delta_p/rho)),c)

   END SUBROUTINE inv_prec_vel_1D_SL
!------------------------------------------------------------------------------!
!> This subroutine computes the preconditioning velocity for 1D stagnation line viscous flows.
   SUBROUTINE visc_prec_vel_1D_SL (vol, rho, mu, V, delta_p, c, Vp)

     IMPLICIT NONE

     REAL(KIND=8), INTENT(IN) :: vol           !< cell volume
     REAL(KIND=8), INTENT(IN) :: rho           !< density 
     REAL(KIND=8), INTENT(IN) :: mu            !< viscosity
     REAL(KIND=8), INTENT(IN) :: V             !< velocity magnitude
     REAL(KIND=8), INTENT(IN) :: delta_p       !< pressure difference 
     REAL(KIND=8), INTENT(IN) :: c             !< speed of sound
     REAL(KIND=8), INTENT(OUT) :: Vp           !< viscous pre-conditioning velocity

     ! Viscous pre-conditioning velocity
     Vp = MIN(MAX(V,DSQRT(delta_p/rho),mu/(rho*vol)),c)

   END SUBROUTINE visc_prec_vel_1D_SL
!------------------------------------------------------------------------------!
