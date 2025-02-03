!------------------------------------------------------------------------------!
!> This subroutine computes the conservative and physical variables from the primitive variables for 
!! 1D nonequilibrium flows gas flows when solving the Euler equations
  SUBROUTINE prim_to_cons_phys_neq_1D_Eu (prim, cons, phys_prop)

    USE mod_general_data,                 ONLY: nb_ns, nb_temp, nb_int_temp, pos_u, pos_T, pos_rhou, pos_rhoE, & 
                                              & pos_u_cell, pos_T_cell, pos_alpha_cell, pos_beta_cell,         & 
                                              & pos_pres_cell, pos_c_cell, pos_gamma_cell, pos_rho_cell,       &
                                              & pos_h0_cell, pos_ei_cell, pos_ek_cell, rhoi, ei, beta, temp,   & 
                                              & rho_eint
    USE mod_neq_function_pointer,         ONLY: library_get_thermodynamic_data, library_get_energy_densities

    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8) :: tmp
    REAL(KIND=8) :: rho
    REAL(KIND=8) :: alpha, c, gamma, ek, p, u

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prim        !< vector of primitive variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cons       !< vector of conservative variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: phys_prop  !< vector of physical properties 

    ! Species and mixture densities
    rho = 0.d0
    DO i = 1,nb_ns 
       tmp = prim(i)
       rhoi(i) = tmp
       rho     = rho + tmp
    ENDDO

    ! Velocity and kinetic energy per unit mass
    u  = prim(pos_u)
    ek = 0.5d0*u**2

    ! Temperature(s) 
    DO i = 1,nb_temp 
       temp(i) = prim(pos_T + i - 1)
    ENDDO 

    ! Compute energy densities
    CALL library_get_energy_densities (rhoi, temp, rho_eint)

    ! Compute thermodynamic data
    CALL library_get_thermodynamic_data (rho, rhoi, temp, c, gamma, p, alpha, beta, ei)

    ! Fill the vector of conservative variables
    DO i = 1,nb_ns 
       cons(i) = rhoi(i)
    ENDDO

    cons(pos_rhou) = rho*u

    DO i = 1,nb_temp
       cons(pos_rhoE + i - 1) = rho_eint(i)
    ENDDO
    cons(pos_rhoE) = cons(pos_rhoE) + rho*ek

    ! Fill the vector of physical properties
    DO i = 1,nb_temp
       phys_prop(pos_T_cell + i - 1) = temp(i)
    ENDDO

    phys_prop(pos_rho_cell)   = rho
    phys_prop(pos_u_cell)     = u
    phys_prop(pos_h0_cell)    = (rho_eint(1) + p)/rho + ek
    phys_prop(pos_c_cell)     = c
    phys_prop(pos_gamma_cell) = gamma
    phys_prop(pos_pres_cell)  = p
    phys_prop(pos_ek_cell)    = ek
    phys_prop(pos_alpha_cell) = alpha

    DO i = 1,nb_int_temp + 1
       phys_prop(pos_beta_cell + i - 1) = beta(i)
    ENDDO
    
    DO i = 1,nb_ns + nb_ns*nb_int_temp  
       phys_prop(pos_ei_cell + i - 1) = ei(i) 
    ENDDO 

  END SUBROUTINE prim_to_cons_phys_neq_1D_Eu 
!------------------------------------------------------------------------------!
!> This subroutine computes the conservative and physical variables from the primitive variables for 
!! 1D nonequilibrium flows gas flows when solving the Navie-Stokes equations
  SUBROUTINE prim_to_cons_phys_neq_1D_Ns (prim, cons, phys_prop)

    USE mod_general_data,                 ONLY: nb_ns, nb_temp, nb_int_temp, pos_u, pos_T, pos_rhou, pos_rhoE,      & 
                                              & pos_u_cell, pos_T_cell, pos_alpha_cell, pos_beta_cell,              & 
                                              & pos_pres_cell, pos_c_cell, pos_gamma_cell, pos_rho_cell,            &
                                              & pos_h0_cell, pos_ei_cell, pos_ek_cell, pos_mu_cell, pos_kappa_cell, &
                                              & pos_Di_cell, pos_lambda_cell, pos_chi_cell, rhoi, ei, beta, temp,   & 
                                              & rho_eint, Di, chi, xi, lambda_vec
    USE mod_neq_function_pointer,         ONLY: library_get_thermodynamic_data, library_get_energy_densities,    &
                                              & library_get_transpCoeff, library_get_molar_fractions

    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8) :: tmp
    REAL(KIND=8) :: rho
    REAL(KIND=8) :: alpha, c, gamma, ek, kappa, mu, p, u

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prim        !< vector of primitive variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cons       !< vector of conservative variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: phys_prop  !< vector of physical properties

    ! Species and mixture densities
    rho = 0.d0
    DO i = 1,nb_ns 
       tmp = prim(i)
       rhoi(i) = tmp
       rho     = rho + tmp
    ENDDO

    ! Velocity and kinetic energy per unit mass
    u  = prim(pos_u)
    ek = 0.5d0*u**2

    ! Temperature(s) 
    DO i = 1,nb_temp 
       temp(i) = prim(pos_T + i - 1)
    ENDDO 

    ! Compute energy densities and species molar fractions
    CALL library_get_energy_densities (rhoi, temp, rho_eint)
    CALL library_get_molar_fractions(rhoi, xi)

    ! Compute thermodynamic data
    CALL library_get_thermodynamic_data (rho, rhoi, temp, c, gamma, p, alpha, beta, ei)

    ! Compute transport coefficients
    CALL library_get_transpCoeff(p, xi, temp, mu, kappa, lambda_vec, Di, chi)  

    ! Fill the vector of conservative variables
    DO i = 1,nb_ns 
       cons(i) = rhoi(i)
    ENDDO

    cons(pos_rhou) = rho*u

    DO i = 1,nb_temp
       cons(pos_rhoE + i - 1) = rho_eint(i)
    ENDDO
    cons(pos_rhoE) = cons(pos_rhoE) + rho*ek

    ! Fill the vector of physical properties
    DO i = 1,nb_temp
       phys_prop(pos_T_cell + i - 1) = temp(i)
    ENDDO

    phys_prop(pos_rho_cell)   = rho
    phys_prop(pos_u_cell)     = u
    phys_prop(pos_h0_cell)    = (rho_eint(1) + p)/rho + ek
    phys_prop(pos_c_cell)     = c
    phys_prop(pos_gamma_cell) = gamma
    phys_prop(pos_pres_cell)  = p
    phys_prop(pos_ek_cell)    = ek
    phys_prop(pos_alpha_cell) = alpha

    DO i = 1,nb_int_temp + 1
       phys_prop(pos_beta_cell + i - 1) = beta(i)
    ENDDO
    
    DO i = 1,nb_ns + nb_ns*nb_int_temp  
       phys_prop(pos_ei_cell + i - 1) = ei(i) 
    ENDDO 

    DO i = 1,nb_ns 
       phys_prop(pos_Di_cell + i - 1) = Di(i) 
    ENDDO

    DO i = 1,nb_ns 
       phys_prop(pos_chi_cell + i - 1) = chi(i) 
    ENDDO

    phys_prop(pos_mu_cell) = mu
    phys_prop(pos_kappa_cell) = kappa

    DO i = 1,nb_temp
       phys_prop(pos_lambda_cell + i - 1) = lambda_vec(i) 
    ENDDO

  END SUBROUTINE prim_to_cons_phys_neq_1D_Ns 
!------------------------------------------------------------------------------!

