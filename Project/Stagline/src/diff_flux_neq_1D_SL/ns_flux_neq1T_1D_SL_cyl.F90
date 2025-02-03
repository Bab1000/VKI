!------------------------------------------------------------------------------!
!> This subroutine computes the diffusive flux when solving the 1D stagnation line equations 
!! for nonequilibrium flows (cylinder case - 1 temperature).
  SUBROUTINE ns_flux_neq1T_1D_SL_cyl (r_l, r_r, vol_l, vol_r, left_data, right_data, u_left, u_right, fd)

    USE mod_general_data,            ONLY: nb_ns, nb_dim, pos_u_cell, pos_v_cell, pos_ei_cell, pos_T_cell,  & 
                                         & pos_pres_cell, pos_mu_cell, pos_lambda_cell, pos_rhou, pos_rhov, & 
                                         & pos_rhoE, rhoi, xi, Ri, ei, chi, diff_driv, Ji 
    USE mod_neq_function_pointer,    ONLY: library_get_species_DiffFlux, library_get_molar_fractions,       & 
                                         & library_comp_tol
    USE mod_numerics_data,           ONLY: xil, xir, rhoil, rhoir

    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8) :: ov_dr
    REAL(KIND=8) :: dens, r, u, v, lambda, p, mu, rhoi_Vdi, T
    REAL(KIND=8) :: tau_rr, tau_rt, tau_tt
    REAL(KIND=8) :: mul, mur, lambdal, lambdar, pl, pr,  Tl, Tr, ul, ur, vl, vr
    REAL(KIND=8) :: du_dr, dv_dr, dT_dr
    REAL(KIND=8) :: q, q_Diff, q_Fourier, q_Therm_Diff

    REAL(KIND=8), INTENT(IN) :: r_l                       !< radial position of left state
    REAL(KIND=8), INTENT(IN) :: r_r                       !< radial position of rigth state
    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state 
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: fd         !< diffusive flux

    ! Explicit interface for subroutine stress_tensor_1D_SL_cyl
    INTERFACE
      SUBROUTINE stress_tensor_1D_SL_cyl (mu, r, u, v, du_dr, dv_dr, tau_rr, tau_rt, tau_tt)
        IMPLICIT NONE
        REAL(KIND=8), INTENT(IN) :: mu, r, u, v, du_dr, dv_dr
        REAL(KIND=8), INTENT(OUT) :: tau_rr, tau_rt, tau_tt
      END SUBROUTINE stress_tensor_1D_SL_cyl
    END INTERFACE

    ! Cell interface position
    r = 0.5d0*(r_l + r_r)

    ! Temperature, pressure, velocity components, dynamic viscosity, thermal conductivity, species densities 
    ! and molar fractions of left and right states
    ! Left state
    Tl      = left_data(pos_T_cell)
    pl      = left_data(pos_pres_cell)
    ul      = left_data(pos_u_cell)
    vl      = left_data(pos_v_cell)
    mul     = left_data(pos_mu_cell)
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
    vr      = right_data(pos_v_cell)
    mur     = right_data(pos_mu_cell)
    lambdar = right_data(pos_lambda_cell)

    DO i = 1,nb_ns 
       rhoir(i) = u_right(i)
    ENDDO
    CALL library_get_molar_fractions(rhoir, xir)
    CALL library_comp_tol(xir)

    ! Temperature, pressure, velocity, dynamic viscosity and thermal conductivity at cell interface
    T   = 0.5d0*(Tl + Tr)
    p   = 0.5d0*(pl + pr)
    u   = 0.5d0*(ul + ur)
    v   = 0.5d0*(vl + vr)
    mu  = 0.5d0*(mul + mur)
    lambda = 0.5d0*(lambdal + lambdar) 
 
    ! Species densities and specific energies at cell interface  
    DO i = 1,nb_ns 
       rhoi(i) = 0.5d0*(rhoil(i) + rhoir(i))
       ei(i)   = 0.5d0*(left_data(pos_ei_cell + i - 1) + right_data(pos_ei_cell + i - 1))
    ENDDO

    ! Species molar fractions at cell interface
    CALL library_get_molar_fractions (rhoi, xi)
    CALL library_comp_tol(xi)
 
    ! Velocity, temperature and species molar fraction gradients
    ov_dr = 2.d0/(vol_l + vol_r)
    du_dr = ov_dr*(ur - ul)
    dv_dr = ov_dr*(vr - vl)
    dT_dr = ov_dr*(Tr - Tl)

    ! Compute linearly independent diffusion driving forces
    diff_driv = 0.d0
    DO i = 1,nb_ns 
       diff_driv(i) = (xir(i) - xil(i))*ov_dr
    ENDDO
    
    ! Compute the species mass diffusion velocities
    CALL library_get_species_DiffFlux(p, T, T, xi, diff_driv, Ji)

    ! Stress tensor
    CALL stress_tensor_1D_SL_cyl (mu, r, u, v, du_dr, dv_dr, tau_rr, tau_rt, tau_tt)

    ! Diffusive flux components 
    ! Species continuity equations
    q_Fourier = - lambda*dT_dr
    q_Diff    = 0.d0
    DO i = 1,nb_ns
       rhoi_Vdi = Ji(i)
       fd(i)    = - rhoi_Vdi  
       q_Diff   = q_Diff + rhoi_Vdi*(ei(i) + Ri(i)*T)
    ENDDO
    q = q_Fourier + q_Diff

    ! Radial and circumferential momentum equations
    fd(pos_rhou) = tau_rr
    fd(pos_rhov) = tau_rt

    ! Global energy equation
    fd(pos_rhoE) = tau_rr*u - q

  END SUBROUTINE ns_flux_neq1T_1D_SL_cyl
!------------------------------------------------------------------------------!
