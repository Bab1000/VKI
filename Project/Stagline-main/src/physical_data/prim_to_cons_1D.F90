!------------------------------------------------------------------------------!
!> This subroutine computes the set of conservative variable vector from the set of primitive variable vector
!! for 1D calorically perfect gas flows.
  SUBROUTINE prim_to_cons_1D (prim, cons)
  
    USE mod_general_data,      ONLY: gamma

    IMPLICIT NONE

    REAL(KIND=8) :: gamma_minus1
    REAL(KIND=8) :: u, rho, p, v2 
  
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prim   !< vector of primitive variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cons  !< vector of conservative variables
 
    ! Useful quantity
    gamma_minus1 = gamma - 1.d0
    
    rho = prim(1)
    u   = prim(2)
    p   = prim(3)
    v2  = u*u
    
    cons(1) = rho
    cons(2) = rho*u
    cons(3) = p/gamma_minus1 + 0.5d0*rho*v2
      
  END SUBROUTINE prim_to_cons_1D
!------------------------------------------------------------------------------!
