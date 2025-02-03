!------------------------------------------------------------------------------!
!> This subroutine computes the diffusive source term for 1D stagnation line nonequilibrium flows 
!! (cylinder case - N temperatures).
  SUBROUTINE source_term_diff_neqNT_1D_SL_cyl (r_c, vol_l, vol_c, vol_r, prop_left, prop_cell, prop_right, & 
                                             & u_left, u_cell, u_right, s)

    USE mod_general_data,            ONLY: nb_ns, nb_temp, nb_int_temp, pos_u_cell, pos_v_cell, pos_ei_cell,  & 
                                         & pos_T_cell, pos_pres_cell, pos_mu_cell, pos_lambda_cell, pos_rhou, & 
                                         & pos_rhov, pos_rhoE, pos_rhoek, rhoi, xi, Ri, ei, diff_driv, Ji,    & 
                                         & lambda_vec 
    USE mod_neq_function_pointer,    ONLY: library_get_species_DiffFlux, library_get_molar_fractions,         & 
                                         & library_comp_tol
    USE mod_numerics_data,           ONLY: xil, xir, rhoil, rhoir

    IMPLICIT NONE

    INTEGER :: i, k
    REAL(KIND=8) :: tmp
    REAL(KIND=8) :: ov_dr, ov_rc
    REAL(KIND=8) :: u, v, mu, p, q_intk, rhoi_Vdi, T
    REAL(KIND=8) :: tau_rr, tau_rt, tau_tt
    REAL(KIND=8) :: Tl, Tr, ul, ur, vl, vr
    REAL(KIND=8) :: du_dr, dv_dr, dT_dr
    REAL(KIND=8) :: q, q_Diff, q_Diff_int, q_Fourier_tr, q_Fourier_int
    REAL(KIND=8), DIMENSION(nb_temp) :: grad_T 

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

    ! Explicit interface for subroutine stress_tensor_1D_SL_cyl
    INTERFACE
      SUBROUTINE stress_tensor_1D_SL_cyl (mu, r, u, v, du_dr, dv_dr, tau_rr, tau_rt, tau_tt)
        IMPLICIT NONE
        REAL(KIND=8), INTENT(IN) :: mu, r, u, v, du_dr, dv_dr
        REAL(KIND=8), INTENT(OUT) :: tau_rr, tau_rt, tau_tt
      END SUBROUTINE stress_tensor_1D_SL_cyl
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

    ! Temperature, pressure, velocity components and thermal conductivity components of central state
    T  = prop_cell(pos_T_cell)
    p  = prop_cell(pos_pres_cell)
    u  = prop_cell(pos_u_cell)
    v  = prop_cell(pos_v_cell)
    mu = prop_cell(pos_mu_cell)
    DO i = 1,nb_temp
       lambda_vec(i) = prop_cell(pos_lambda_cell + i -1) 
    ENDDO 

    ! Species densities of central state 
    DO i = 1,nb_ns 
       rhoi(i) = u_cell(i)
    ENDDO

    ! Species specific energies of central state
    DO i = 1,nb_ns + nb_int_temp*nb_ns
       ei(i) = prop_cell(pos_ei_cell + i - 1)
    ENDDO

    ! Species mass and molar fractions at cell interface
    CALL library_get_molar_fractions (rhoi, xi)
    CALL library_comp_tol(xi)
 
    ! Velocity, temperature, pressure and species molar fraction gradients
    ov_rc = 1.d0/r_c
    ov_dr = 2.d0/(vol_l + 2.d0*vol_c + vol_r)
    du_dr = ov_dr*(ur - ul)
    dv_dr = ov_dr*(vr - vl)
    DO i = 1,nb_temp
       grad_T(i) = (prop_right(pos_T_cell + i - 1) - prop_left(pos_T_cell + i - 1))*ov_dr
    ENDDO 

    ! Compute linearly independent diffusion driving forces (thermal diffusion is added)
    diff_driv = 0.d0
    DO i = 1,nb_ns 
       diff_driv(i) = (xir(i) - xil(i))*ov_dr
    ENDDO
    
    ! Compute the species diffusion velocities 
    CALL library_get_species_DiffFlux(p, T, T, xi, diff_driv, Ji)

    ! Stress tensor
    CALL stress_tensor_1D_SL_cyl (mu, r_c, u, v, du_dr, dv_dr, tau_rr, tau_rt, tau_tt)

    ! Diffusive source term components
    ! Species continuity equations
    q_Fourier_tr = - lambda_vec(1)*grad_T(1)
    q_Diff       = 0.d0
    DO i = 1,nb_ns
       tmp    = Ji(i)
       s(i)   = - tmp*ov_rc  
       q_Diff = q_Diff + tmp*(ei(i) + Ri(i)*T)
    ENDDO

    ! Internal energy equations
    ! Fourier heat flux (components associated to internal energy)
    q_Fourier_int = 0.d0
    DO k = 1,nb_int_temp
       q_intk = - lambda_vec(k + 1)*grad_T(k + 1)
       q_Fourier_int = q_Fourier_int + q_intk
       q_Diff_int    = 0.d0
       DO i = 1,nb_ns
          q_Diff_int = q_Diff_int + Ji(i)*ei(nb_ns + nb_int_temp*(i - 1) + k)
       ENDDO
       s(pos_rhoek + k - 1) = - (q_intk + q_Diff_int)*ov_rc
    ENDDO
    q = q_Fourier_tr + q_Fourier_int + q_Diff    

    ! Radial and circumferential momentum equations
    s(pos_rhou) = (tau_rr + tau_rt - tau_tt)*ov_rc
    s(pos_rhov) = (2.d0*tau_rt - tau_tt)*ov_rc

    ! Global energy equation
    s(pos_rhoE) = (u*(tau_rr + tau_rt) + v*tau_tt - q)*ov_rc

  END SUBROUTINE source_term_diff_neqNT_1D_SL_cyl
!------------------------------------------------------------------------------!
