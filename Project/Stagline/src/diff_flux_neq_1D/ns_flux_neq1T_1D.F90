!------------------------------------------------------------------------------!
!> This subroutine computes the diffusive flux for the 1D Navier-Stokes equations for nonequilibrium flows
!! when using a 1T model (diffusive flux Jacobians are not provided in output).
  SUBROUTINE ns_flux_neq1T_1D (vol_l, vol_r, left_data, right_data, u_left, u_right, fd)

    USE mod_general_data,            ONLY: nb_ns, nb_dim, pos_u_cell, pos_ei_cell, pos_T_cell, pos_mu_cell,   &
                                         & pos_kappa_cell, pos_lambda_cell, pos_pres_cell, pos_rhou,          & 
                                         & pos_rhoE, rhoi, xi, Ri, ei, diff_driv, Ji 
    USE mod_neq_function_pointer,    ONLY: library_get_species_DiffFlux, library_get_molar_fractions,         & 
                                         & library_comp_tol
    USE mod_numerics_data,           ONLY: xil, xir, rhoil, rhoir

    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8), PARAMETER :: coeff   = 4.d0/3.d0
    REAL(KIND=8) :: ov_dx
    REAL(KIND=8) :: kappal, kappar, mul, mur, lambdal, lambdar, pl, pr, Tl, Tr, ul, ur
    REAL(KIND=8) :: kappa, lambda, mu, p, rhoi_Vdi, u, T, tau
    REAL(KIND=8) :: du_dx, dT_dx
    REAL(KIND=8) :: q, q_Diff, q_Fourier

    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: fd         !< diffusive flux
  
    ! Temperature, pressure, velocity, bulk viscosity, dynamic viscosity, thermal conductivity, 
    ! species densities and molar fractions of left and right states
    ! Left state
    Tl      = left_data(pos_T_cell)
    pl      = left_data(pos_pres_cell)
    ul      = left_data(pos_u_cell)
    mul     = left_data(pos_mu_cell)
    kappal  = left_data(pos_kappa_cell)
    lambdal = left_data(pos_lambda_cell) 

    DO i = 1,nb_ns 
       rhoil(i) = u_left(i)
    ENDDO
    CALL library_get_molar_fractions(rhoil, xil)   
    CALL library_comp_tol(xil)

    ! Right state
    Tr      = right_data(pos_T_cell)
    pr      = right_data(pos_pres_cell)
    ur      = right_data(pos_u_cell)
    mur     = right_data(pos_mu_cell)
    kappar  = right_data(pos_kappa_cell)
    lambdar = right_data(pos_lambda_cell)

    DO i = 1,nb_ns 
       rhoir(i) = u_right(i)
    ENDDO
    CALL library_get_molar_fractions(rhoir, xir)   
    CALL library_comp_tol(xir)

    ! Temperature, pressure, velocity, dynamic viscosity bulk viscosity and thermal conductivity at cell interface
    T   = 0.5d0*(Tl + Tr)
    p   = 0.5d0*(pl + pr)
    u   = 0.5d0*(ul + ur)
    mu  = 0.5d0*(mul + mur)
    kappa  = 0.5d0*(kappal + kappar)
    lambda = 0.5d0*(lambdal + lambdar) 
 
    ! Species densities and specific energies at cell interface 
    DO i = 1,nb_ns 
       rhoi(i) = 0.5d0*(rhoil(i) + rhoir(i))
       ei(i)   = 0.5d0*(left_data(pos_ei_cell + i - 1) + right_data(pos_ei_cell + i - 1)) 
    ENDDO

    ! Species molar fractions at cell interface
    CALL library_get_molar_fractions (rhoi, xi)
    CALL library_comp_tol(xi)

    ! Velocity, temperature, pressure and species molar fraction gradients
    ov_dx = 2.d0/(vol_l + vol_r)
    du_dx = (ur - ul)*ov_dx
    dT_dx = (Tr - Tl)*ov_dx

    ! Compute linearly independent diffusion driving forces (thermal diffusion is added)
    DO i = 1,nb_ns 
       diff_driv(i) = (xir(i) - xil(i))*ov_dx
    ENDDO
    
    ! Compute the species mass diffusion flux 
    CALL library_get_species_DiffFlux(p, T, T, xi, diff_driv, Ji)

    ! Diffusive flux components 
    ! Species continuity equations
    q_Fourier = - lambda*dT_dx
    q_Diff    = 0.d0
    DO i = 1,nb_ns
       rhoi_Vdi = Ji(i)
       fd(i)    = - rhoi_Vdi  
       q_Diff   = q_diff + rhoi_Vdi*(ei(i) + Ri(i)*T)
    ENDDO

    ! Global momentum and energy equations
    q = q_Fourier + q_Diff

    tau = (coeff*mu + kappa)*du_dx
    fd(pos_rhou) = tau
    fd(pos_rhoE) = tau*u - q

  END SUBROUTINE ns_flux_neq1T_1D 
!------------------------------------------------------------------------------!
