!------------------------------------------------------------------------------!
!> This subroutine computes the diffusive flux for the 1D Navier-Stokes equations for nonequilibrium flows 
!! when using a NT temperature model (diffusive flux Jacobians are not provided in output).
  SUBROUTINE ns_flux_neqNT_1D (vol_l, vol_r, left_data, right_data, u_left, u_right, fd)

    USE mod_general_data,            ONLY: nb_ns, nb_dim, nb_temp, nb_int_temp, pos_u_cell, pos_ei_cell,  &
                                         & pos_T_cell, pos_pres_cell, pos_mu_cell, pos_kappa_cell,        & 
                                         & pos_lambda_cell, pos_rhou, pos_rhoE, pos_rhoek, rhoi, xi, Ri,  & 
                                         & ei, diff_driv, Ji, lambda_vec
    USE mod_neq_function_pointer,    ONLY: library_get_species_DiffFlux, library_get_molar_fractions,     & 
                                         & library_comp_tol         
    USE mod_numerics_data,           ONLY: xil, xir, rhoil, rhoir

    IMPLICIT NONE

    INTEGER :: i, k
    REAL(KIND=8), PARAMETER :: coeff   = 4.d0/3.d0
    REAL(KIND=8) :: tmp
    REAL(KIND=8) :: ov_dx
    REAL(KIND=8) :: kappal, kappar, mul, mur, pl, pr, Tl, Tr, ul, ur
    REAL(KIND=8) :: kappa, mu, p, q_intk, tau, u, T, du_dx
    REAL(KIND=8) :: q, q_Diff, q_Diff_int, q_Fourier_tr, q_Fourier_int
    REAL(KIND=8), DIMENSION(nb_temp) :: dT_dx

    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: fd         !< diffusive flux

    ! Temperature, pressure, velocity, dynamic viscosity, bulkl viscosity, 
    ! species densities and molar fractions of left and right states
    ! Left state
    Tl   = left_data(pos_T_cell) 
    pl   = left_data(pos_pres_cell)
    ul   = left_data(pos_u_cell)
    kappal = left_data(pos_kappa_cell)
    mul    = left_data(pos_mu_cell)
   
    DO i = 1,nb_ns 
       rhoil(i) = u_left(i)
    ENDDO
    CALL library_get_molar_fractions(rhoil, xil) 
    CALL library_comp_tol(xil)

    ! Right state
    Tr   = right_data(pos_T_cell)
    pr   = right_data(pos_pres_cell)
    ur   = right_data(pos_u_cell)
    kappar = right_data(pos_kappa_cell)
    mur    = right_data(pos_mu_cell)
   
    DO i = 1,nb_ns 
       rhoir(i) = u_right(i)
    ENDDO
    CALL library_get_molar_fractions(rhoir, xir)  
    CALL library_comp_tol(xir)

    ! Temperature, pressure, velocity, dynamic viscosity, bulk viscosity and thermal conductivity components at cell interface 
    T   = 0.5d0*(Tl + Tr)
    p   = 0.5d0*(pl + pr)
    u   = 0.5d0*(ul + ur)
    mu  = 0.5d0*(mul + mur)
    kappa = 0.5d0*(kappal + kappar)
    DO i = 1,nb_temp
       lambda_vec(i) = 0.5d0*(left_data(pos_lambda_cell + i -1) + right_data(pos_lambda_cell + i -1)) 
    ENDDO

    ! Species densities at cell interface  
    DO i = 1,nb_ns 
       rhoi(i) = 0.5d0*(rhoil(i) + rhoir(i))
    ENDDO

    ! Species specific energies at cell interface 
    DO i = 1,nb_ns + nb_int_temp*nb_ns
       ei(i) = 0.5d0*(left_data(pos_ei_cell + i - 1) + right_data(pos_ei_cell + i - 1))
    ENDDO 

    ! Species molar fractions at cell interface
    CALL library_get_molar_fractions (rhoi, xi)
    CALL library_comp_tol(xi)

    ! Velocity and temperature gradients
    ov_dx = 2.d0/(vol_l + vol_r)
    du_dx = (ur - ul)*ov_dx
    DO i = 1,nb_temp
       dT_dx(i) = (right_data(pos_T_cell + i - 1) - left_data(pos_T_cell + i - 1))*ov_dx
    ENDDO
    
    ! Compute diffusion driving forces 
    DO i = 1,nb_ns 
       diff_driv(i) = (xir(i) - xil(i))*ov_dx
    ENDDO
  
    ! Compute the species mass diffusion flux 
    CALL library_get_species_DiffFlux(p, T, T, xi, diff_driv, Ji)
   
    ! Diffusive flux components 
    ! Species continuity equations
    q_Fourier_tr  = - lambda_vec(1)*dT_dx(1)
    q_Diff        = 0.d0
    DO i = 1,nb_ns
       tmp    = Ji(i) 
       fd(i)  = - tmp  
       q_Diff = q_Diff + tmp*(ei(i) + Ri(i)*T)
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
    q = q_Fourier_tr + q_Fourier_int + q_Diff

    tau = (coeff*mu + kappa)*du_dx
    fd(pos_rhou) = tau
    fd(pos_rhoE) = tau*u - q
   
  END SUBROUTINE ns_flux_neqNT_1D 
!------------------------------------------------------------------------------!
