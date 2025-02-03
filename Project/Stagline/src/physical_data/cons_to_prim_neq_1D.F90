!------------------------------------------------------------------------------!
!> This subroutine computes the primitive variable vector from the conservative variable vector 
!! for 1D nonequilibrium flows.
  SUBROUTINE cons_to_prim_neq_1D (cons, prim)
   
    USE mod_general_data,            ONLY: nb_ns, nb_temp, nb_int_temp, nb_te, pos_u, pos_T, pos_rhou, & 
                                         & pos_rhoE, temp, rho_eint 
    USE mod_neq_function_pointer,    ONLY: library_get_temperatures
  
    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8) :: rho, rho_ek, u, tmp
  
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons   !< vector of conservative variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: prim  !< vector of primitive variables
  
    ! Species densities 
    rho = 0.d0
    DO i = 1,nb_ns 
       tmp = cons(i)
       prim(i) = tmp
       rho     = rho + tmp
    ENDDO
  
    ! Velocity and kinetic energy density
    u      = cons(pos_rhou)/rho
    rho_ek = 0.5d0*rho*u*u
  
    ! Velocity
    prim(pos_u) = u

    ! Energy densities
    rho_eint(1) = cons(pos_rhoE) - rho_ek 
    DO i = 1,nb_int_temp + nb_te
       rho_eint(i + 1) = cons(pos_rhoE + i)
    ENDDO
             
    ! Computing the temperatures
    CALL library_get_temperatures (prim(1:nb_ns), rho_eint, temp)   

    ! Temperatures 
    DO i = 1,nb_temp
       prim(pos_T + i - 1) = temp(i)
    ENDDO
  
  END SUBROUTINE cons_to_prim_neq_1D
!------------------------------------------------------------------------------!
