!------------------------------------------------------------------------------!
!> This subroutine computes the diffusive flux when solving the 1D stagnation line equations 
!! for nonequilibrium flows (sphere case - N temperatures).
  SUBROUTINE ns_flux_neqNT_1D_SL_sph (r_l, r_r, vol_l, vol_r, left_data, right_data, u_left, u_right, fd)

    USE mod_general_data,            ONLY: nb_ns, nb_dim, nb_temp, nb_int_temp, pos_u_cell, pos_v_cell,  & 
                                         & pos_ei_cell, pos_T_cell, pos_pres_cell, pos_mu_cell,          & 
                                         & pos_lambda_cell, pos_rhou, pos_rhov, pos_rhoE, pos_rhoek,     & 
                                         & rhoi, xi, Ri, ei, hi, diff_driv, Ji, lambda_vec 
    USE mod_neq_function_pointer,    ONLY: library_get_species_DiffFlux, library_get_molar_fractions,    & 
                                         & library_comp_tol
    USE mod_numerics_data,           ONLY: xil, xir, rhoil, rhoir

    IMPLICIT NONE

    INTEGER :: i, k
    REAL(KIND=8) :: tmp
    REAL(KIND=8) :: ov_dr
    REAL(KIND=8) :: r, u, v, lambda, mu, p, q_intk, rhoi_Vdi, Vdiff, T
    REAL(KIND=8) :: tau_rr, tau_rt, tau_tt
    REAL(KIND=8) :: mul, mur, pl, pr, Tl, Tr, ul, ur, vl, vr
    REAL(KIND=8) :: du_dr, dv_dr
    REAL(KIND=8) :: q, q_Diff, q_Diff_int, q_Fourier_tr, q_Fourier_int
    REAL(KIND=8), DIMENSION(nb_temp) :: grad_T

    REAL(KIND=8), INTENT(IN) :: r_l                       !< radial position of left state
    REAL(KIND=8), INTENT(IN) :: r_r                       !< radial position of rigth state
    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state 
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: fd         !< diffusive flux

    ! Explicit interface for subroutine stress_tensor_1D_SL_sph
    INTERFACE
      SUBROUTINE stress_tensor_1D_SL_sph (mu, r, u, v, du_dr, dv_dr, tau_rr, tau_rt, tau_tt)
        IMPLICIT NONE
        REAL(KIND=8), INTENT(IN) :: mu, r, u, v, du_dr, dv_dr
        REAL(KIND=8), INTENT(OUT) :: tau_rr, tau_rt, tau_tt
      END SUBROUTINE stress_tensor_1D_SL_sph
    END INTERFACE

    ! Cell interface position
    r = 0.5d0*(r_l + r_r)

    ! Temperature, pressure, velocity components, dynamic viscosity, thermal conductivity, species densities 
    ! and molar fractions of left and right states
    ! Left state
    Tl  = left_data(pos_T_cell)
    pl  = left_data(pos_pres_cell)
    ul  = left_data(pos_u_cell)
    vl  = left_data(pos_v_cell)
    mul = left_data(pos_mu_cell)

    DO i = 1,nb_ns 
       rhoil(i) = u_left(i)
    ENDDO
    CALL library_get_molar_fractions(rhoil, xil)   
    CALL library_comp_tol(xil)

    ! Right state
    Tr  = right_data(pos_T_cell)
    pr  = right_data(pos_pres_cell)
    ur  = right_data(pos_u_cell)
    vr  = right_data(pos_v_cell)
    mur = right_data(pos_mu_cell)

    DO i = 1,nb_ns 
       rhoir(i) = u_right(i)
    ENDDO
    CALL library_get_molar_fractions(rhoir, xir)   
    CALL library_comp_tol(xir)

    ! Temperature, pressure, velocity, dynamic viscosity and thermal conductivity components at cell interface
    T   = 0.5d0*(Tl + Tr)
    p   = 0.5d0*(pl + pr)
    u   = 0.5d0*(ul + ur)
    v   = 0.5d0*(vl + vr)
    mu  = 0.5d0*(mul + mur)
    DO i = 1,nb_temp
       lambda_vec(i) = 0.5d0*(left_data(pos_lambda_cell + i -1)+ right_data(pos_lambda_cell + i -1)) 
    ENDDO 

    ! Species densities at cell interface  
    DO i = 1,nb_ns 
       rhoi(i) = 0.5d0*(rhoil(i) + rhoir(i)) 
    ENDDO

    ! Specific energies at cell interface 
    DO i = 1,nb_ns + nb_int_temp*nb_ns
       ei(i) = 0.5d0*(left_data(pos_ei_cell + i - 1) + right_data(pos_ei_cell + i - 1))
    ENDDO

    ! Species molar fractions at cell interface
    CALL library_get_molar_fractions (rhoi, xi)
    CALL library_comp_tol(xi) 

    ! Velocity, temperature, pressure and species molar fraction gradients
    ov_dr = 2.d0/(vol_l + vol_r)
    du_dr = ov_dr*(ur - ul)
    dv_dr = ov_dr*(vr - vl)
    DO i = 1,nb_temp
       grad_T(i) = (right_data(pos_T_cell + i - 1) - left_data(pos_T_cell + i - 1))*ov_dr
    ENDDO

    ! Compute linearly independent diffusion driving forces (thermal diffusion is added)
    diff_driv = 0.d0
    DO i = 1,nb_ns 
       diff_driv(i) = (xir(i) - xil(i))*ov_dr
    ENDDO
    
    ! Compute the species mass diffusion flux
    CALL library_get_species_DiffFlux(p, T, T, xi, diff_driv, Ji)

    ! Stress tensor
    CALL stress_tensor_1D_SL_sph (mu, r, u, v, du_dr, dv_dr, tau_rr, tau_rt, tau_tt)

    ! Diffusive flux components 
    ! Species continuity equations
    q_Fourier_tr = - lambda_vec(1)*grad_T(1)
    q_Diff       = 0.d0
    DO i = 1,nb_ns
       tmp    = Ji(i)
       fd(i)  = - tmp  
       q_Diff = q_Diff + tmp*(ei(i) + Ri(i)*T)
    ENDDO

    ! Fourier heat flux (components associated to internal energy)
    q_Fourier_int = 0.d0
    DO k = 1,nb_int_temp
       q_intk = - lambda_vec(k + 1)*grad_T(k + 1)
       q_Fourier_int = q_Fourier_int + q_intk
       q_Diff_int    = 0.d0
       DO i = 1,nb_ns
          q_Diff_int = q_Diff_int + Ji(i)*ei(nb_ns*k + i)
       ENDDO
       fd(pos_rhoek + k - 1) = - (q_intk + q_Diff_int)
    ENDDO
    q = q_Fourier_tr + q_Fourier_int + q_Diff

    ! Radial and circumferential momentum equations
    fd(pos_rhou) = tau_rr
    fd(pos_rhov) = tau_rt

    ! Global energy equation
    fd(pos_rhoE) = tau_rr*u - q

  END SUBROUTINE ns_flux_neqNT_1D_SL_sph
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!> This subroutine computes the diffusive flux when solving the 1D stagnation line equations 
!! for nonequilibrium flows  accounting for the mesh metrics (sphere case - N temperatures).
  SUBROUTINE ns_flux_neqNT_1D_SL_sph_metr (r_l, r_r, vol_l, vol_r, left_data, right_data, u_left, u_right, fd)

    USE mod_general_data,            ONLY: nb_ns, nb_dim, nb_temp, nb_int_temp, pos_u_cell, pos_v_cell,  & 
                                         & pos_ei_cell, pos_T_cell, pos_pres_cell, pos_mu_cell,          & 
                                         & pos_lambda_cell, pos_rhou, pos_rhov, pos_rhoE, pos_rhoek,     & 
                                         & rhoi, xi, Ri, ei, diff_driv, Ji, lambda_vec 
    USE mod_neq_function_pointer,    ONLY: library_get_species_DiffFlux, library_get_molar_fractions,    & 
                                         & library_comp_tol
    USE mod_numerics_data,           ONLY: xil, xir, rhoil, rhoir

    IMPLICIT NONE

    INTEGER :: i, k
    REAL(KIND=8) :: tmp
    REAL(KIND=8) :: ov_dr
    REAL(KIND=8) :: r, u, v, lambda, mu, p, q_intk, rhoi_Vdi, Vdiff, T
    REAL(KIND=8) :: tau_rr, tau_rt, tau_tt
    REAL(KIND=8) :: mul, mur, pl, pr, Tl, Tr, ul, ur, vl, vr
    REAL(KIND=8) :: du_dr, dv_dr
    REAL(KIND=8) :: q, q_Diff, q_Diff_int, q_Fourier_tr, q_Fourier_int
    REAL(KIND=8), DIMENSION(nb_temp) :: grad_T

    REAL(KIND=8), INTENT(IN) :: r_l                       !< radial position of left state
    REAL(KIND=8), INTENT(IN) :: r_r                       !< radial position of rigth state
    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state 
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: fd         !< diffusive flux

    ! Explicit interface for subroutine stress_tensor_1D_SL_sph
    INTERFACE
      SUBROUTINE stress_tensor_1D_SL_sph (mu, r, u, v, du_dr, dv_dr, tau_rr, tau_rt, tau_tt)
        IMPLICIT NONE
        REAL(KIND=8), INTENT(IN) :: mu, r, u, v, du_dr, dv_dr
        REAL(KIND=8), INTENT(OUT) :: tau_rr, tau_rt, tau_tt
      END SUBROUTINE stress_tensor_1D_SL_sph
    END INTERFACE

       ! Cell interface position
    r = (r_l + vol_l*0.5d0)

    ! Temperature, pressure, velocity components, dynamic viscosity, thermal conductivity, species densities 
    ! and molar fractions of left and right states
    ! Left state
    Tl  = left_data(pos_T_cell)
    pl  = left_data(pos_pres_cell)
    ul  = left_data(pos_u_cell)
    vl  = left_data(pos_v_cell)
    mul = left_data(pos_mu_cell)

    DO i = 1,nb_ns 
       rhoil(i) = u_left(i)
    ENDDO
    CALL library_get_molar_fractions(rhoil, xil)   

    ! Right state
    Tr  = right_data(pos_T_cell)
    pr  = right_data(pos_pres_cell)
    ur  = right_data(pos_u_cell)
    vr  = right_data(pos_v_cell)
    mur = right_data(pos_mu_cell)

    DO i = 1,nb_ns 
       rhoir(i) = u_right(i)
    ENDDO
    CALL library_get_molar_fractions(rhoir, xir)

    ! Temperature, pressure, velocity, dynamic viscosity and thermal conductivity components at cell interface
    ! If considering non-uniform stencil (A. Turchi)
    T   = (Tl*vol_l + Tr*vol_r)/(vol_l+vol_r)
    p   = (pl*vol_l + pr*vol_r)/(vol_l+vol_r)
    u   = (ul *vol_l+ ur*vol_r)/(vol_l+vol_r)
    v   = (vl *vol_l+ vr*vol_r)/(vol_l+vol_r)
    mu  = (mul *vol_l+ mur*vol_r)/(vol_l+vol_r)
    DO i = 1,nb_temp
       lambda_vec(i) = (left_data(pos_lambda_cell + i -1) *vol_l + right_data(pos_lambda_cell + i -1) * vol_r)/(vol_l+vol_r) 
    ENDDO 

    ! Species densities at cell interface  
    DO i = 1,nb_ns 
       rhoi(i) = (rhoil(i) *vol_l + rhoir(i) *vol_r )/(vol_l+vol_r) 
    ENDDO

    ! Specific energies at cell interface 
    DO i = 1,nb_ns + nb_int_temp*nb_ns
       ei(i) = (left_data(pos_ei_cell + i - 1) *vol_l + right_data(pos_ei_cell + i - 1) *vol_r )/(vol_l+vol_r)
    ENDDO

    ! Species molar fractions at cell interface
    CALL library_get_molar_fractions (rhoi, xi)
    CALL library_comp_tol(xi) 

    ! Velocity, temperature, pressure and species molar fraction gradients
    ov_dr = 2.d0/(vol_l + vol_r)
    du_dr = ov_dr*(ur - ul)
    dv_dr = ov_dr*(vr - vl)
    DO i = 1,nb_temp
       grad_T(i) = (right_data(pos_T_cell + i - 1) - left_data(pos_T_cell + i - 1))*ov_dr
    ENDDO

    ! Compute linearly independent diffusion driving forces (thermal diffusion is added)
    diff_driv = 0.d0
    DO i = 1,nb_ns 
       diff_driv(i) = (xir(i) - xil(i))*ov_dr
    ENDDO
    
    ! Compute the species mass diffusion flux
    CALL library_get_species_DiffFlux(p, T, T, xi, diff_driv, Ji)

    ! Stress tensor
    CALL stress_tensor_1D_SL_sph (mu, r, u, v, du_dr, dv_dr, tau_rr, tau_rt, tau_tt)

    ! Diffusive flux components 
    ! Species continuity equations
    q_Fourier_tr = - lambda_vec(1)*grad_T(1)
    q_Diff       = 0.d0
    DO i = 1,nb_ns
       tmp    = Ji(i)
       fd(i)  = - tmp  
       q_Diff = q_Diff + tmp*(ei(i) + Ri(i)*T)
    ENDDO

    ! Fourier heat flux (components associated to internal energy)
    q_Fourier_int = 0.d0
    DO k = 1,nb_int_temp
       q_intk = - lambda_vec(k + 1)*grad_T(k + 1)
       q_Fourier_int = q_Fourier_int + q_intk
       q_Diff_int    = 0.d0
       DO i = 1,nb_ns
          q_Diff_int = q_Diff_int + Ji(i)*ei(nb_ns*k + i)
       ENDDO
       fd(pos_rhoek + k - 1) = - (q_intk + q_Diff_int)
    ENDDO
    q = q_Fourier_tr + q_Fourier_int + q_Diff

    ! Radial and circumferential momentum equations
    fd(pos_rhou) = tau_rr
    fd(pos_rhov) = tau_rt

    ! Global energy equation
    fd(pos_rhoE) = tau_rr*u - q

  END SUBROUTINE ns_flux_neqNT_1D_SL_sph_metr
!------------------------------------------------------------------------------!
