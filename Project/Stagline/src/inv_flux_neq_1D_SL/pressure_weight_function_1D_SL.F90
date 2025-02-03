!------------------------------------------------------------------------------!
!> This subroutine computes the pressure-dependend sensor function for the modifiend Steger-Warming
!> described by Candler on the Hypersonic Nonequilibrium Flows (Chapter 5 by Candler)

  SUBROUTINE pressure_weight_function(p_right, p_left, omega_interface) 

    IMPLICIT NONE
    
    REAL(KIND=8) :: Delta_P, g_square_dP
    REAL(KIND=8), PARAMETER :: g =0.5D0
   
    REAL(KIND=8), INTENT(IN)  :: p_left, p_right 
    REAL(KIND=8), INTENT(OUT) :: omega_interface
  
    g_square_dP = (g* Delta_P(p_left,p_right))**2+ 1.0D0

    omega_interface = 1.D0-0.5D0*(1.D0/g_square_dP)
     

  END SUBROUTINE pressure_weight_function

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
  FUNCTION Delta_P (p_left, p_right)

    REAL(KIND=8), INTENT(IN)  :: p_left, p_right 
    REAL(KIND=8) :: Delta_P
   
  
  Delta_P  = abs(p_right-p_left)/min(p_left,p_right)

  END FUNCTION Delta_P

