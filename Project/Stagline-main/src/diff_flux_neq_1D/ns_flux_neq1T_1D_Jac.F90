!------------------------------------------------------------------------------!
!> This subroutine computes the diffusive flux and diffusive flux Jacobians for the 1D Navier-Stokes equations 
!! for nonequilibrium flows when using a 1T model (diffusive flux Jacobians are provided in output). 
  SUBROUTINE ns_flux_neq1T_1D_Jac (vol_l, vol_r, left_data, right_data, u_left, u_right, fd, jfdl, jfdr)
    
    USE mod_general_data,            ONLY: nb_ns, nb_dim, pos_u_cell, pos_ei_cell, pos_beta_cell, pos_T_cell,  & 
                                         & pos_mu_cell, pos_kappa_cell, pos_lambda_cell, pos_pres_cell,        & 
                                         & pos_Di_cell, pos_rhou, pos_rhoE, Di, xi, yi, Ri, rhoi, ei,          &
                                         & diff_driv, Ji
    USE mod_neq_function_pointer,    ONLY: library_get_species_DiffFlux, library_get_molar_fractions,          & 
                                         & library_get_mass_fractions_from_molar_fractions, library_comp_tol
    USE mod_numerics_data,           ONLY: xil, xir, rhoil, rhoir

    IMPLICIT NONE

    INTEGER :: i, j
    REAL(KIND=8), PARAMETER :: coeff   = 4.d0/3.d0
    REAL(KIND=8) :: ov_dx
    REAL(KIND=8) :: coeff_mu_ov_rho, coeff_mu_ov_rho_u, coeff_mu_ov_rho_v2, lambda_ov_beta
    REAL(KIND=8) :: betal, betar, kappal, kappar, mul, mur, lambdal, lambdar, pl, pr, Tl, Tr, ul, ur
    REAL(KIND=8) :: beta, dens, ek, kappa, mu, lambda, p, rho, rhoi_Vdi, u, T, tau, mm
    REAL(KIND=8) :: q, q_Diff, q_Fourier
    REAL(KIND=8) :: du_dx, dT_dx

    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: fd         !< diffusive flux
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfdl     !< diffusive flux Jacobian with respect to the left state
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfdr     !< diffusive flux Jacobian with respect to the right state

    ! Temperature, pressure, velocity, beta factor, dynamic viscosity, bulk viscosity thermal conductivity 
    ! and species densities and  molar fractions of left and right states
    ! Left state
    Tl      = left_data(pos_T_cell)
    pl      = left_data(pos_pres_cell)
    ul      = left_data(pos_u_cell)
    mul     = left_data(pos_mu_cell)
    betal   = left_data(pos_beta_cell)
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
    betar   = right_data(pos_beta_cell)
    kappar  = right_data(pos_kappa_cell)
    lambdar = right_data(pos_lambda_cell)

    DO i = 1,nb_ns 
       rhoir(i) = u_right(i)
    ENDDO
    CALL library_get_molar_fractions(rhoir, xir)  
    CALL library_comp_tol(xir)

    ! Temperature, pressure, density, velocity, kinetic energy, beta factor, dynamic viscosity 
    ! and thermal conductivity at cell interface
    T   = 0.5d0*(Tl + Tr)
    p   = 0.5d0*(pl + pr)
    u   = 0.5d0*(ul + ur)
    mu  = 0.5d0*(mul + mur)
    beta   = 0.5d0*(betal + betar)
    kappa  = 0.5d0*(kappal + kappar)  
    lambda = 0.5d0*(lambdal + lambdar)  

    ! Species densities, specific energies and average diffusion coefficients at cell interface
    rho = 0.d0
    DO i = 1,nb_ns 
       dens    = 0.5d0*(rhoil(i) + rhoir(i))
       rhoi(i) = dens 
       rho     = rho + dens
       ei(i)   = 0.5d0*(left_data(pos_ei_cell + i - 1) + right_data(pos_ei_cell + i - 1))
       Di(i)   = 0.5d0*(left_data(pos_Di_cell + i - 1) + right_data(pos_Di_cell + i - 1))
    ENDDO

    ! Species molar fractions at cell interface
    CALL library_get_molar_fractions (rhoi, xi)
    CALL library_comp_tol(xi)

    ! Compute species mass fractions and mixture molar mass at cell interface 
    ! from species molar fractions
    CALL library_get_mass_fractions_from_molar_fractions(xi, mm, yi)

    ! Velocity and temperature gradients
    ov_dx = 2.d0/(vol_l + vol_r)
    du_dx = (ur - ul)*ov_dx
    dT_dx = (Tr - Tl)*ov_dx

    ! Compute linearly independent diffusion driving forces
    DO i = 1,nb_ns 
       diff_driv(i) = (xir(i) - xil(i))*ov_dx
    ENDDO

    ! Compute the species mass diffusion flux 
    CALL library_get_species_DiffFlux(p, T, T, xi, diff_driv, Ji)

    ! Diffusive flux components 
    ! Species continuity conservation equations
    q_Fourier = - lambda*dT_dx
    q_Diff    = 0.d0
    DO i = 1,nb_ns
       rhoi_Vdi = Ji(i)
       fd(i)    = - rhoi_Vdi  
       q_Diff   = q_diff + rhoi_Vdi*(ei(i) + Ri(i)*T)
    ENDDO

    ! Global momentum and energy equations
    q  = q_Fourier + q_Diff 

    tau = (coeff*mu + kappa)*du_dx
    fd(pos_rhou) = tau
    fd(pos_rhoE) = tau*u - q 

    ! Right state diffusive flux Jacobian
    ! Common factors
    ek = 0.5d0*u**2
    lambda_ov_beta = lambda/beta*ov_dx
    coeff_mu_ov_rho    = (coeff*mu + kappa)/rho*ov_dx
    coeff_mu_ov_rho_u  = coeff_mu_ov_rho*u
    coeff_mu_ov_rho_v2 = coeff_mu_ov_rho_u*u

    ! Column j,  j = 1,..,nb_ns 
    DO j = 1,nb_ns 
       
       DO i = 1,nb_ns 
          jfdr(i,j) = 0.d0
       ENDDO

       jfdr(pos_rhou,j) = - coeff_mu_ov_rho_u
       jfdr(pos_rhoE,j) = - coeff_mu_ov_rho_v2  + lambda_ov_beta*(ek - ei(j))

    ENDDO

    ! Column nb_ns + 1
    DO i = 1,nb_ns 
       jfdr(i,pos_rhou) = 0.d0
    ENDDO
    
    jfdr(pos_rhou,pos_rhou) = coeff_mu_ov_rho 
    jfdr(pos_rhoE,pos_rhou) = coeff_mu_ov_rho_u - lambda_ov_beta*u 

    ! Column nb_ns + 2
    DO i = 1,nb_ns 
       jfdr(i,pos_rhoE) = 0.d0
    ENDDO
    
    jfdr(pos_rhou,pos_rhoE) = 0.d0 
    jfdr(pos_rhoE,pos_rhoE) = lambda_ov_beta

    ! Left state diffusive flux Jacobian
    jfdl = - jfdr
    
    IF (nb_ns.GT.1) THEN
      PRINT*
      WRITE(*,'(A)')'In "ns_flux_neq1T_1D_Jac.F90", species difflusive flux Jacobian sub-matrix to be completed!'
      PRINT*
      STOP
    ENDIF

  END SUBROUTINE ns_flux_neq1T_1D_Jac 
!------------------------------------------------------------------------------!
