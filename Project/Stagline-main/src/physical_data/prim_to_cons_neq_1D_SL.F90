!------------------------------------------------------------------------------!
!> This subroutine computes the conservative variable vector from the primitive variable vector 
!! for 1D stagnation line nonequilibrium flows.
  SUBROUTINE prim_to_cons_neq_1D_SL (prim, cons)
   
    USE mod_general_data,            ONLY: nb_ns, nb_temp, pos_u, pos_v, pos_T, pos_rhou, pos_rhov, & 
                                         & pos_rhoE, temp, rho_eint
    USE mod_neq_function_pointer,    ONLY: library_get_energy_densities

    IMPLICIT NONE
  
    INTEGER :: i
    REAL(KIND=8) :: rho, ek, u, v, tmp
  
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prim   !< vector of primitive variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cons  !< vector of conservative variables
  
    ! Species densities
    rho = 0.d0   
    DO i = 1,nb_ns 
       tmp = prim(i)
       cons(i) = tmp
       rho     = rho + tmp
    ENDDO
  
    ! Radial velocity and kinetic energy per unit mass
    u  = prim(pos_u)
    ek = 0.5d0*u**2
  
    ! Circumferential velocity 
    v = prim(pos_v)

    ! Radial momentum density
    cons(pos_rhou) = rho*u
  
    ! Circumferential momentum density 
    cons(pos_rhov) = rho*v

    ! Temperature vector
    DO i = 1,nb_temp
       temp(i) = prim(pos_T + i - 1)
    ENDDO
   
    ! Total energy per unit volume (total and needed components)
    CALL library_get_energy_densities (prim(1:nb_ns), temp, rho_eint)
    
    ! Remaining components of conservative variable vector
    DO i = 1,nb_temp
       cons(pos_rhoE + i - 1) = rho_eint(i)
    ENDDO
  
    ! Adding the kinetic energy contribution (rhoE)
    cons(pos_rhoE) = cons(pos_rhoE) + rho*ek 
  
  END SUBROUTINE prim_to_cons_neq_1D_SL
!------------------------------------------------------------------------------!
