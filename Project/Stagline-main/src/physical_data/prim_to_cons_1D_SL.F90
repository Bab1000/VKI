!------------------------------------------------------------------------------!
!> This subroutine computes the set of conservative variable vector from the set of primitive variable vector
!! for 1D stagnation line calorically perfect gas flows.
  SUBROUTINE prim_to_cons_1D_SL (prim, cons)
  
    USE mod_general_data,      ONLY: gamma

    IMPLICIT NONE

    REAL(KIND=8) :: gamma_minus1
    REAL(KIND=8) :: u, v, rho, p, v2 
  
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prim   !< vector of primitive variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cons  !< vector of conservative variables
 
    ! Useful quantity
    gamma_minus1 = gamma - 1.d0
    
    rho = prim(1)
    u   = prim(2)
    v   = prim(3)
    p   = prim(4)
    v2  = u*u
    
    cons(1) = rho
    cons(2) = rho*u
    cons(3) = rho*v
    cons(4) = p/gamma_minus1 + 0.5d0*rho*v2
      
  END SUBROUTINE prim_to_cons_1D_SL
!------------------------------------------------------------------------------!
