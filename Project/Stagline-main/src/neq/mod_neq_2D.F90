!------------------------------------------------------------------------------!
!> This module provides implementation for subroutine used when computing 2D nonequilibrium flows.  
  MODULE mod_neq_2D

    USE mod_general_data,     ONLY: nb_ns, nb_dim, nb_temp, nb_int_temp, nb_te, nb_eq, pos_em,  & 
                                  & pos_u, pos_v, pos_T, pos_Tk, pos_Te, pos_rhou, pos_rhov,    & 
                                  & pos_rhoE, pos_rhoek, pos_rhose, gamma_e, gamma_e_m1,        & 
                                  & ov_gamma_e_m1, Ri
  
    IMPLICIT NONE
  
    ! Subroutines and functions for handling data.
    CONTAINS 
  
      !------------------------------------------------------!
      !> This subroutine computes the conservative variable vector from the primitive variable vector 
      !! for 2D nonequilibrium flows.
      SUBROUTINE prim_to_cons_neq_2D (prim, cons)
   
        USE mod_neq_function_pointer,    ONLY: library_get_energy_densities
  
        INTEGER :: i 
        REAL(KIND=8) :: rho, u, v, ek, tmp
   
        REAL(KIND=8), DIMENSION(nb_temp) :: temp, rho_e  
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prim
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cons
  
        ! Species densities
        rho = 0.d0   
        DO i = 1,nb_ns 
           tmp = prim(i)
           cons(i) = tmp
           rho     = rho + tmp
        ENDDO
  
        ! x and y velocity components and kinetic energy per unit mass
        u  = prim(pos_u)
        v  = prim(pos_v)
        ek = 0.5d0*(u*u + v*v)
  
        ! Momentum densities
        cons(pos_rhou) = rho*u
        cons(pos_rhov) = rho*v
  
        ! Temperature vector
        DO i = 1,nb_temp
           temp(i) = prim(pos_T + i - 1)
        ENDDO
   
        ! Total energy per unit volume (total and needed components)
        CALL library_get_energy_densities (prim(1:nb_ns), temp, rho_e)
   
        ! Remaining components of conservative variable vector
        DO i = 1,nb_temp
           cons(pos_rhoE + i - 1) = rho_e(i)
        ENDDO
  
        ! Adding the kinetic energy contribution (rhoE)
        cons(pos_rhoE) = cons(pos_rhoE) + rho*ek 
        
      END SUBROUTINE prim_to_cons_neq_2D
  
      !------------------------------------------------------!
      !> This subroutine computes the primitive variable vector from the conservative variable vector 
      !! for 2D nonequilibrium flows.
      SUBROUTINE cons_to_prim_neq_2D (cons, prim)
  
        USE mod_neq_function_pointer,    ONLY: library_get_temperatures
  
        INTEGER :: i
        REAL(KIND=8) :: rho, rho_ek, u, v, tmp
        REAL(KIND=8), DIMENSION(nb_temp) :: temp, rho_eint
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: prim
  
        ! Species densities 
        rho = 0.d0
        DO i = 1,nb_ns 
           tmp = cons(i)
           prim(i) = tmp
           rho     = rho + tmp
        ENDDO
  
        ! Velocity components and kinetic energy density
        u      = cons(pos_rhou)/rho
        v      = cons(pos_rhov)/rho
        rho_ek = 0.5d0*rho*(u*u + v*v)
  
        prim(pos_u) = u
        prim(pos_v) = v

        ! Energy densities
        rho_eint(1) = cons(pos_rhoE) - rho_ek 
        DO i = 1,nb_int_temp + nb_te
           rho_eint(pos_rhoE + i) = cons(pos_rhoE + i)
        ENDDO
             
        ! Computing the temperatures
        CALL library_get_temperatures (prim(1:nb_ns), rho_eint, temp) 
  
        ! Temperatures 
        DO i = 1,nb_temp
           prim(pos_T + i - 1) = temp(i)
        ENDDO
  
      END SUBROUTINE cons_to_prim_neq_2D
      
  END MODULE mod_neq_2D
!------------------------------------------------------------------------------!
