!------------------------------------------------------------------------------!
!> This subroutine computes the diffusive source term for 1D stagnation line nonequilibrium flows 
!! (sphere case - 1 temperature).
  SUBROUTINE source_term_diff_neq1T_1D_SL_sph (r_c, vol_l, vol_c, vol_r, prop_left, prop_cell, prop_right, & 
                                             & u_left, u_cell, u_right, s)

    USE mod_general_data,            ONLY: nb_ns, nb_dim, pos_u_cell, pos_v_cell, pos_ei_cell, pos_T_cell, & 
                                         & pos_pres_cell, pos_mu_cell, pos_lambda_cell, pos_rhou,          & 
                                         & pos_rhov, pos_rhoE, rhoi, xi, Ri, ei, chi, diff_driv, Ji 
    USE mod_neq_function_pointer,    ONLY: library_get_species_DiffFlux, library_get_molar_fractions,      & 
                                         & library_comp_tol
    USE mod_numerics_data,           ONLY: xil, xir, rhoil, rhoir

    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8) :: ov_dr, ov_rc
    REAL(KIND=8) :: u, v, lambda, mu, p, rhoi_Vdi, T
    REAL(KIND=8) :: tau_rr, tau_rt, tau_tt
    REAL(KIND=8) :: pl, pr, Tl, Tr, ul, ur, vl, vr
    REAL(KIND=8) :: du_dr, dv_dr, dT_dr
    REAL(KIND=8) :: q, q_Diff, q_Fourier

    REAL(KIND=8), INTENT(IN) :: r_c                        !< radial location of cell centroid
    REAL(KIND=8), INTENT(IN) :: vol_l                      !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_c                      !< volume of central state
    REAL(KIND=8), INTENT(IN) :: vol_r                      !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left       !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_cell       !< conservative variables of central state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right      !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prop_left    !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prop_cell    !< physical properties of central state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prop_right   !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s           !< source term

    ! Explicit interface for subroutine stress_tensor_1D_SL_sph
    INTERFACE
      SUBROUTINE stress_tensor_1D_SL_sph (mu, r, u, v, du_dr, dv_dr, tau_rr, tau_rt, tau_tt)
        IMPLICIT NONE
        REAL(KIND=8), INTENT(IN) :: mu, r, u, v, du_dr, dv_dr
        REAL(KIND=8), INTENT(OUT) :: tau_rr, tau_rt, tau_tt
      END SUBROUTINE stress_tensor_1D_SL_sph
    END INTERFACE

    ! Temperature, velocity components, species densities and molar fractions of left and right states
    ! Left state
    Tl = prop_left(pos_T_cell)
    ul = prop_left(pos_u_cell)
    vl = prop_left(pos_v_cell)

    DO i = 1,nb_ns 
       rhoil(i) = u_left(i)
    ENDDO
    CALL library_get_molar_fractions(rhoil, xil)   
    CALL library_comp_tol(xil)

    ! Right state
    Tr = prop_right(pos_T_cell)
    ur = prop_right(pos_u_cell)
    vr = prop_right(pos_v_cell)

    DO i = 1,nb_ns 
       rhoir(i) = u_right(i)
    ENDDO
    CALL library_get_molar_fractions(rhoir, xir)
    CALL library_comp_tol(xir)

    ! Temperature, pressure, velocity components, species densities and specific energies
    T  = prop_cell(pos_T_cell)
    p  = prop_cell(pos_pres_cell)
    u  = prop_cell(pos_u_cell)
    v  = prop_cell(pos_v_cell)
    mu = prop_cell(pos_mu_cell)
    lambda = prop_cell(pos_lambda_cell)

    ! Species densities and specific energies of central state 
    DO i = 1,nb_ns 
       rhoi(i) = u_cell(i)
       ei(i)   = prop_cell(pos_ei_cell + i - 1)
    ENDDO

    ! Species molar fractions at cell interface
    CALL library_get_molar_fractions (rhoi, xi)
    CALL library_comp_tol(xi) 

    ! Velocity, temperature and species molar fraction gradients
    ov_rc = 1.d0/r_c
    ov_dr = 2.d0/(vol_l + 2.d0*vol_c + vol_r)
    du_dr = ov_dr*(ur - ul)
    dv_dr = ov_dr*(vr - vl)
    dT_dr = ov_dr*(Tr - Tl)

    ! Compute linearly independent diffusion driving forces
    diff_driv = 0.d0
    DO i = 1,nb_ns 
       diff_driv(i) = (xir(i) - xil(i))*ov_dr
    ENDDO
    
    ! Compute the species mass diffusion flux 
    CALL library_get_species_DiffFlux(p, T, T, xi, diff_driv, Ji)

    ! Stress tensor
    CALL stress_tensor_1D_SL_sph (mu, r_c, u, v, du_dr, dv_dr, tau_rr, tau_rt, tau_tt)

    ! Diffusive source term components 
    ! Species continuity equations
    q_Fourier = - lambda*dT_dr
    q_Diff    = 0.d0
    DO i = 1,nb_ns
       rhoi_Vdi = Ji(i)
       s(i)     = - 2.d0*rhoi_Vdi*ov_rc 
       q_Diff   = q_Diff + rhoi_Vdi*(ei(i) + Ri(i)*T)
    ENDDO
    q = q_Fourier + q_Diff

    ! Radial and circumferential momentum equations
    s(pos_rhou) = 2.d0*(tau_rr - tau_tt + tau_rt)*ov_rc
    s(pos_rhov) = (3.d0*tau_rt - tau_tt)*ov_rc

    ! Global energy equation
    s(pos_rhoE) = 2.d0*(u*(tau_rr + tau_rt) + v*tau_tt - q)*ov_rc

  END SUBROUTINE source_term_diff_neq1T_1D_SL_sph
!------------------------------------------------------------------------------!
