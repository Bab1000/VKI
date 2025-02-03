!------------------------------------------------------------------------------!
!> This subroutione computes the set of conservative variable vector from the set of primitive variable vector
!! for 1D calorically perfect gas flows.
  SUBROUTINE cons_to_prim_1D (cons, prim)
  
    USE mod_general_data,     ONLY: gamma

    IMPLICIT NONE

    REAL(KIND=8) :: gamma_minus1
    REAL(KIND=8) :: rho, rhou, rhoE, rhov2
  
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons   !< vector of conservative variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: prim  !< vector of primitive variables

    ! Useful quantity
    gamma_minus1 = gamma - 1.d0
        
    rho   = cons(1)
    rhou  = cons(2)
    rhoE  = cons(3)
    rhov2 = rhou*rhou
      
    prim(1) = rho
    prim(2) = rhou/rho
    prim(3) = gamma_minus1*(rhoE - 0.5d0*rhov2/rho)
    
  END SUBROUTINE cons_to_prim_1D
!------------------------------------------------------------------------------!
