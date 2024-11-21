!------------------------------------------------------------------------------!
!> This subroutine computes the diffusive flux and diffusive flux Jacobians for the 1D Navier-Stokes equations 
!! for nonequilibrium flows when using an NT model (diffusive flux Jacobians are provided in output). 
  SUBROUTINE ns_flux_neqNT_1D_Jac (vol_l, vol_r, left_data, right_data, u_left, u_right, fd, jfdl, jfdr)

    USE mod_general_data,            ONLY: nb_ns, nb_dim, nb_temp, nb_int_temp, pos_u_cell, pos_ei_cell,   &
                                         & pos_T_cell, pos_pres_cell, pos_mu_cell, pos_kappa_cell,         & 
                                         & pos_lambda_cell, pos_beta_cell, pos_betak_cell, pos_ek_cell,    & 
                                         & pos_Di_cell, pos_rhou, pos_rhoE, pos_rhoek, rhoi, Ri, ei, eiT,  & 
                                         & xi, yi, Di, betak, diff_driv, Ji, lambda_vec
    USE mod_neq_function_pointer,    ONLY: library_get_species_DiffFlux, library_get_molar_fractions,      & 
                                         & library_get_mass_fractions_from_molar_fractions, library_comp_tol
    USE mod_numerics_data,           ONLY: xil, xir, rhoil, rhoir


    IMPLICIT NONE

    INTEGER :: i, j, k, kp
    INTEGER :: pos_i, pos_j, pos_k
    REAL(KIND=8), PARAMETER :: coeff   = 4.d0/3.d0
    REAL(KIND=8) :: ov_dx, ov_mm
    REAL(KIND=8) :: coeff_mu_ov_rho, coeff_mu_ov_rho_u, coeff_mu_ov_rho_v2, lambda_ov_beta
    REAL(KIND=8) :: betal, betar, eil, eir, kappal, kappar, mul, mur, pl, pr, Tl, Tr, ul, ur
    REAL(KIND=8) :: beta, dens, ek, fac, rhoi_Vdi, kappa, mu, p, rho, q_intk, tau, u, T, mm
    REAL(KIND=8) :: q, q_Diff, q_Diff_int, q_Fourier_tr, q_Fourier_int
    REAL(KIND=8) :: du_dx
    REAL(KIND=8), DIMENSION(nb_temp) :: dT_dx
    REAL(KIND=8), DIMENSION(nb_int_temp) :: lambdak_ov_betak
    REAL(KIND=8), DIMENSION(nb_ns) :: lambdak_ov_betak_ei

    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: fd         !< diffusive flux
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfdl     !< diffusive flux Jacobian with respect to the left state
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfdr     !< diffusive flux Jacobian with respect to the right state

    ! Temperature, pressure, velocity, dynamic viscosity, bulk viscosity, beta factor 
    ! and species densities and molar fractions of of left and right states
    ! Left state
    Tl    = left_data(pos_T_cell)
    pl    = left_data(pos_pres_cell)
    ul    = left_data(pos_u_cell)
    mul   = left_data(pos_mu_cell)
    kappal = left_data(pos_kappa_cell)
    betal  = left_data(pos_beta_cell)   

    DO i = 1,nb_ns 
       rhoil(i) = u_left(i)
    ENDDO
    CALL library_get_molar_fractions(rhoil, xil)
    CALL library_comp_tol(xil)

    ! Right state
    Tr    = right_data(pos_T_cell)
    pr    = right_data(pos_pres_cell)
    ur    = right_data(pos_u_cell)
    mur   = right_data(pos_mu_cell)
    kappar = right_data(pos_kappa_cell)
    betar  = right_data(pos_beta_cell)   

    DO i = 1,nb_ns 
       rhoir(i) = u_right(i)
    ENDDO
    CALL library_get_molar_fractions(rhoir, xir) 
    CALL library_comp_tol(xir)

    ! Temperature, pressure, velocity, dynamic viscosity, bulk viscosity, thermal conductivity components and beta factor 
    ! at cell interface
    T    = 0.5d0*(Tl + Tr)
    p    = 0.5d0*(pl + pr)
    u    = 0.5d0*(ul + ur)
    mu   = 0.5d0*(mul + mur)
    kappa = 0.5d0*(kappal + kappar)
    beta  = 0.5d0*(betal + betar)
    DO i = 1,nb_temp
       lambda_vec(i) = 0.5d0*(left_data(pos_lambda_cell + i - 1) + right_data(pos_lambda_cell + i - 1)) 
    ENDDO
 
    ! Beta factors associated to internal temperatures at cell interface
    DO k = 1,nb_int_temp
       betak(k) = 0.5d0*(left_data(pos_betak_cell + k - 1) + right_data(pos_betak_cell + k - 1))
    ENDDO
    
    ! Species specific energy data at cell interface
    DO i = 1,nb_ns + nb_ns*nb_int_temp
       ei(i) = 0.5d0*(left_data(pos_ei_cell + i - 1) + right_data(pos_ei_cell + i - 1))
    ENDDO
   
    DO i = 1,nb_ns
       pos_k = pos_ei_cell + nb_ns + nb_int_temp*(i - 1)
       eil  = 0.d0
       eir  = 0.d0
       DO k = 1,nb_int_temp
          eil = eil + left_data(pos_k + k - 1) 
          eir = eir + right_data(pos_k + k - 1) 
       ENDDO
       eiT(i) = ei(i) - 0.5d0*(eil + eir)
    ENDDO

    ! Mixture density, species densities, specific thermal equilibrium energy, and average   
    ! diffusion coefficients at cell interface 
    rho = 0.d0
    DO i = 1,nb_ns 
       dens    = 0.5d0*(rhoil(i) + rhoir(i))
       rhoi(i) = dens 
       rho     = rho + dens
       Di(i)   = 0.5d0*(left_data(pos_Di_cell + i - 1) + right_data(pos_Di_cell + i - 1))
    ENDDO

    ! Species molar fractions at cell interface
    CALL library_get_molar_fractions (rhoi, xi)
    CALL library_comp_tol(xi)

    ! Compute species mass fractions and mixture molar mass at cell interface 
    ! from species molar fractions
    CALL library_get_mass_fractions_from_molar_fractions(xi, mm, yi)
    ov_mm = 1.d0/mm

    ! Velocity and temperature gradients
    ov_dx = 2.d0/(vol_l + vol_r)
    du_dx = (ur - ul)*ov_dx
    DO i = 1,nb_temp
       dT_dx(i) = (right_data(pos_T_cell + i - 1) - left_data(pos_T_cell + i - 1))*ov_dx
    ENDDO

    ! Compute diffusion driving forces (thermal diffusion is included)
    DO i = 1,nb_ns 
       diff_driv(i) = (xir(i) - xil(i))*ov_dx
    ENDDO

    ! Compute the species mass diffusion fluxes 
    CALL library_get_species_DiffFlux(p, T, T, xi, diff_driv, Ji) 

    ! Diffusive flux components 
    ! Species continuity equations
    q_Fourier_tr  = - lambda_vec(1)*dT_dx(1)
    q_Diff        = 0.d0
    DO i = 1,nb_ns
       rhoi_Vdi = Ji(i)
       fd(i)    = - rhoi_Vdi  
       q_Diff   = q_Diff + rhoi_Vdi*(ei(i) + Ri(i)*T)
    ENDDO

    ! Fourier heat flux (components associated to internal energy)
    q_Fourier_int = 0.d0
    DO k = 1,nb_int_temp
       q_intk = - lambda_vec(k + 1)*dT_dx(k + 1)
       q_Fourier_int = q_Fourier_int + q_intk
       q_Diff_int    = 0.d0
       DO i = 1,nb_ns
          q_Diff_int = q_Diff_int + Ji(i)*ei(nb_ns + nb_int_temp*(i - 1) + k)
       ENDDO
       fd(pos_rhoek + k - 1) = - (q_intk + q_Diff_int)
    ENDDO

    ! Global momentum and energy equations
    q   = q_Fourier_tr + q_Fourier_int + q_Diff
    tau = (coeff*mu + kappa)*du_dx
    fd(pos_rhou) = tau
    fd(pos_rhoE) = tau*u - q

    ! RIght state diffusive flux Jacobians
    ! Common factors
    ek = 0.5d0*u**2
    lambda_ov_beta = lambda_vec(1)/beta*ov_dx
    coeff_mu_ov_rho    = (coeff*mu + kappa)/rho*ov_dx
    coeff_mu_ov_rho_u  = coeff_mu_ov_rho*u
    coeff_mu_ov_rho_v2 = coeff_mu_ov_rho_u*u

    DO k = 1,nb_int_temp
       lambdak_ov_betak(k) = lambda_vec(k + 1)/betak(k)*ov_dx
    ENDDO

    DO i = 1,nb_ns
       pos_i = nb_ns + nb_int_temp*(i - 1) 
       fac   = 0.d0
       DO k = 1,nb_int_temp
          fac = fac + ei(pos_i + k)*lambdak_ov_betak(k) 
       ENDDO
       lambdak_ov_betak_ei(i) = fac
    ENDDO

    ! Column j,  j = 1,..,nb_ns 
    DO j = 1,nb_ns 
       
       pos_j = nb_ns + nb_int_temp*(j - 1)

       DO i = 1,nb_ns 
          jfdr(i,j) = 0.d0
       ENDDO

       jfdr(pos_rhou,j) = - coeff_mu_ov_rho_u
       jfdr(pos_rhoE,j) = - coeff_mu_ov_rho_v2  + lambda_ov_beta*(ek - eiT(j)) - lambdak_ov_betak_ei(j)

       DO k = 1,nb_int_temp
          jfdr(pos_rhoek + k - 1,j) = - ei(pos_j + k)*lambdak_ov_betak(k)
       ENDDO

    ENDDO

    ! Column nb_ns + 1
    DO i = 1,nb_ns 
       jfdr(i,pos_rhou) = 0.d0
    ENDDO
    
    jfdr(pos_rhou,pos_rhou) = coeff_mu_ov_rho 
    jfdr(pos_rhoE,pos_rhou) = coeff_mu_ov_rho_u - lambda_ov_beta*u 

    DO k = 1,nb_int_temp
       jfdr(pos_rhoek + k - 1,pos_rhou) = 0.d0
    ENDDO

    ! Column nb_ns + 2
    DO i = 1,nb_ns 
       jfdr(i,pos_rhoE) = 0.d0
    ENDDO
    
    jfdr(pos_rhou,pos_rhoE) = 0.d0 
    jfdr(pos_rhoE,pos_rhoE) = lambda_ov_beta

    DO k = 1,nb_int_temp
       jfdr(pos_rhoek + k - 1,pos_rhoE) = 0.d0
    ENDDO

    ! Column nb_ns + 2 + k,  k = 1,..,nb_int_temp
    DO k = 1,nb_int_temp 

       pos_k = pos_rhoek + k - 1

       fac = lambdak_ov_betak(k)

       DO i = 1,nb_ns 
          jfdr(i,pos_k) = 0.d0
       ENDDO

       jfdr(pos_rhou,pos_k) = 0.d0
       jfdr(pos_rhoE,pos_k) = fac - lambda_ov_beta

       DO kp = 1,k - 1 
          jfdr(pos_rhoek + kp - 1,pos_k) = 0.d0 
       ENDDO

       jfdr(pos_k,pos_k) = fac

       DO kp = k + 1,nb_int_temp 
          jfdr(pos_rhoek + kp - 1,pos_k) = 0.d0 
       ENDDO

    ENDDO

    ! Left state diffusive flux Jacobian
    jfdl = - jfdr

    IF (nb_ns.GT.1) THEN
      PRINT*
      WRITE(*,'(A)')'In "ns_flux_neqNT_1D_Jac.F90", species difflusive flux Jacobian sub-matrix to be completed!'
      PRINT*
      STOP
    ENDIF

  END SUBROUTINE ns_flux_neqNT_1D_Jac 
!------------------------------------------------------------------------------!
