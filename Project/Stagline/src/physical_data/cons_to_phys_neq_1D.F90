!------------------------------------------------------------------------------!
!> This subroutine computes the physical properties from conservative variables for 1D nonequilibrium gas flows
!! when solving the Euler equations. 
  SUBROUTINE cons_to_phys_neq_1D_Eu (cons, phys_prop)

    USE mod_general_data,            ONLY: nb_ns, nb_temp, nb_int_temp, pos_rhou, pos_rhoE, pos_u_cell, & 
                                         & pos_h0_cell, pos_c_cell, pos_gamma_cell, pos_pres_cell,      & 
                                         & pos_rho_cell, pos_ek_cell, pos_alpha_cell, pos_beta_cell,    & 
                                         & pos_ei_cell, pos_T_cell, rho_eint, rhoi, temp, ei, beta
    USE mod_neq_function_pointer,    ONLY: library_get_data 

    IMPLICIT NONE

    INTEGER :: i, k
    REAL(KIND=8) :: tmp
    REAL(KIND=8) :: rho, rhoE
    REAL(KIND=8) :: alpha, c, gamma, ek, p, u

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons        !< vector of conservative variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: phys_prop  !< vector of physical properties

    ! Mixture and species densities 
    rho = 0.d0
    DO i = 1,nb_ns 
       tmp = cons(i)
       rhoi(i) = tmp
       rho     = rho + tmp 
    ENDDO

    ! Velocity and kinetic energy per unit mass
    u  = cons(pos_rhou)/rho
    ek = 0.5d0*u**2 

    ! Total Internal energy density
    rhoE = cons(pos_rhoE)
    rho_eint(1) = rhoE - rho*ek

    ! Internal energy densities
    DO k = 1,nb_temp - 1
       rho_eint(1 + k) = cons(pos_rhoE + k)
    ENDDO 

    ! Compute temperatures and other thermodynamic data
    CALL library_get_data(rho, rhoi, rho_eint, temp, c, gamma, p, alpha, beta, ei)

    ! Fill vector of physical properties 
    DO i = 1,nb_temp
       phys_prop(pos_T_cell + i - 1) = temp(i)
    ENDDO
      
    phys_prop(pos_u_cell)     = u
    phys_prop(pos_h0_cell)    = (rhoE + p)/rho
    phys_prop(pos_c_cell)     = c
    phys_prop(pos_gamma_cell) = gamma
    phys_prop(pos_pres_cell)  = p
    phys_prop(pos_rho_cell)   = rho
    phys_prop(pos_ek_cell)    = ek
    phys_prop(pos_alpha_cell) = alpha
 
    DO i = 1,nb_int_temp + 1
       phys_prop(pos_beta_cell + i - 1)  = beta(i) 
    ENDDO

    DO i = 1,nb_ns + nb_ns*nb_int_temp
       phys_prop(pos_ei_cell + i - 1) = ei(i)  
    ENDDO

  END SUBROUTINE cons_to_phys_neq_1D_Eu
!------------------------------------------------------------------------------!
!> This subroutine computes the physical properties from conservative variables for 1D nonequilibrium gas flows
!! when solving the Navier-Stokes equations. 
  SUBROUTINE cons_to_phys_neq_1D_Ns (cons, phys_prop)

    USE mod_general_data,            ONLY: nb_ns, nb_temp, nb_int_temp, pos_rhou, pos_rhoE, pos_u_cell, & 
                                         & pos_h0_cell, pos_c_cell, pos_gamma_cell, pos_pres_cell,      & 
                                         & pos_rho_cell, pos_ek_cell, pos_alpha_cell, pos_beta_cell,    & 
                                         & pos_ei_cell, pos_mu_cell, pos_kappa_cell, pos_lambda_cell,   & 
                                         & pos_T_cell, pos_Di_cell, pos_chi_cell, rho_eint, rhoi, temp, & 
                                         & ei, Di, chi, beta, xi, lambda_vec
    USE mod_neq_function_pointer,    ONLY: library_get_data, library_get_transpCoeff, library_get_molar_fractions

    IMPLICIT NONE

    INTEGER :: i, k
    REAL(KIND=8) :: tmp
    REAL(KIND=8) :: rho, rhoE
    REAL(KIND=8) :: alpha, c, gamma, ek, kappa, mu, p, u

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons        !< vector of conservative variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: phys_prop  !< vector of physical properties

    ! Mixture and species densities 
    rho = 0.d0
    DO i = 1,nb_ns 
       tmp = cons(i)
       rhoi(i) = tmp
       rho     = rho + tmp 
    ENDDO

    ! Velocity and kinetic energy per unit mass
    u  = cons(pos_rhou)/rho
    ek = 0.5d0*u**2 

    ! Total Internal energy density
    rhoE = cons(pos_rhoE)
    rho_eint(1) = rhoE - rho*ek

    ! Internal energy densities
    DO k = 1,nb_temp - 1
       rho_eint(1 + k) = cons(pos_rhoE + k)
    ENDDO 

    ! Compute temperatures and other thermodynamic data
    CALL library_get_data(rho, rhoi, rho_eint, temp, c, gamma, p, alpha, beta, ei)
 
    ! Compute species molar fractions
    CALL library_get_molar_fractions(rhoi, xi)

    ! Compute species diffusion coefficients, mixture viscosity and thermal conductivity components 
    CALL library_get_transpCoeff(p, xi, temp, mu, kappa, lambda_vec, Di, chi)
 
    ! Fill vector of physical properties 
    DO i = 1,nb_temp
       phys_prop(pos_T_cell + i - 1) = temp(i)
    ENDDO
      
    phys_prop(pos_u_cell)     = u
    phys_prop(pos_h0_cell)    = (rhoE + p)/rho
    phys_prop(pos_c_cell)     = c
    phys_prop(pos_gamma_cell) = gamma
    phys_prop(pos_pres_cell)  = p
    phys_prop(pos_rho_cell)   = rho
    phys_prop(pos_ek_cell)    = ek
    phys_prop(pos_alpha_cell) = alpha
 
    DO i = 1,nb_int_temp + 1
       phys_prop(pos_beta_cell + i - 1)  = beta(i) 
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
    
  END SUBROUTINE cons_to_phys_neq_1D_Ns
!------------------------------------------------------------------------------!

