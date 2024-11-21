!------------------------------------------------------------------------------!
!> This module provides subroutines for variable transformation and eigensystem for 2D calorically perfect gas flows.
  MODULE mod_pg_2D

    IMPLICIT NONE
  
    REAL(KIND=8), SAVE :: gamma, gamma_minus1, gamma_plus1, gamma_minus3, R_gas   
 
    ! Subroutine and functions for handling data.
    CONTAINS 

      !------------------------------------------------------! 
      !> This subroutine initializes the specific heat ratios and the specific gas constants.
      !! It has to be called once before the use of the other subroutines provided in the present module
      SUBROUTINE initialize_pg_2D(in_gamma, in_R_gas) 

        REAL(KIND=8), INTENT(IN) :: in_gamma, in_R_gas

        ! Specific heat ratio and gas constant
        gamma = in_gamma
        gamma_plus1  = gamma + 1.d0
        gamma_minus1 = gamma - 1.d0
        gamma_minus3 = gamma - 3.d0
         
        R_gas = in_R_gas 

      END SUBROUTINE initialize_pg_2D

      !------------------------------------------------------!
      !> This subroutine computes the set of conservative variable vector from the set of primitive variable vector
      !! for 2D calorically perfect gas flows.
      SUBROUTINE prim_to_cons_2D (prim, cons)
  
        REAL(KIND=8) :: u, v, rho, p, v2 
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prim
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cons
 
        rho = prim(1)
        u   = prim(2)
        v   = prim(3)
        p   = prim(4)
        v2  = u*u + v*v
    
        cons(1) = rho
        cons(2) = rho*u
        cons(3) = rho*v
        cons(4) = p/gamma_minus1 + 0.5d0*rho*v2
    
      END SUBROUTINE prim_to_cons_2D

      !------------------------------------------------------!
      !> This subroutione computes the set of conservative variable vector from the set of primitive variable vector
      !! for 2D calorically perfect gas flows.
      SUBROUTINE cons_to_prim_2D (cons, prim)
  
        REAL(KIND=8) :: rho, rhou, rhov, rhoE, rhov2
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: prim
    
        rho   = cons(1)
        rhou  = cons(2)
        rhov  = cons(3)
        rhoE  = cons(4)
        rhov2 = rhou*rhou + rhov*rhov
      
        prim(1) = rho
        prim(2) = rhou/rho
        prim(3) = rhov/rho
        prim(4) = gamma_minus1*(rhoE - 0.5d0*rhov2/rho)
    
      END SUBROUTINE cons_to_prim_2D

      !----------------------------------------------------!
      !> This subroutine computes the inviscid flux for the 2D Euler equations
      !! along the direction specified by the vector n = (nx, ny).
      SUBROUTINE inv_flux_2D (nx, ny, prim, f)
  
        REAL(KIND=8) :: rho, u, v, vn, p, v2, h0, tmp 
  
        REAL(KIND=8), INTENT(IN) :: nx, ny
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prim
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: f 
  
        rho = prim(1)
        u   = prim(2)
        v   = prim(3)
        p   = prim(4)
  
        vn = u*nx + v*ny
        v2 = u*u + v*v
        h0 = p*gamma/(rho*gamma_minus1) + 0.5d0*v2
  
        tmp  = rho*vn
        f(1) = tmp
        f(2) = tmp*u + p*nx
        f(2) = tmp*v + p*ny
        f(4) = tmp*h0
            
      END SUBROUTINE inv_flux_2D

      !----------------------------------------------------!
      !> This subroutine compitues the inviscid flux Jacobians for the 2D Euler equations.
      SUBROUTINE flux_Jacobian_2D (nx, ny, prim, j)
 
        REAL(KIND=8), INTENT(IN) :: nx, ny
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prim
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: j

        j = 0.d0

      END SUBROUTINE flux_Jacobian_2D

  END MODULE mod_pg_2D
!------------------------------------------------------------------------------!
