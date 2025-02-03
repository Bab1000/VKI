!------------------------------------------------------------------------------!
!> This subroutine computes the diffusive source term and the related Jacobians for 1D stagnation line flows 
!! for nonequilibrium flows (sphere case - N temperature and one for free electrons Te).
  SUBROUTINE source_term_diff_neqNT_Te_1D_SL_sph_Jac(r_c, vol_l, vol_c, vol_r, prop_left, prop_cell, prop_right, & 
                                                   & u_left, u_cell, u_right, s, js_left, js_cell, js_right)

     USE mod_general_data,            ONLY: nb_ns, nb_eq, nb_temp, nb_int_temp, pos_u_cell, pos_v_cell, pos_ei_cell, mi, & 
                                         & pos_rho_cell, pos_T_cell, pos_pres_cell, pos_mu_cell, pos_lambda_cell, pos_rhou, & 
                                         & pos_rhov, pos_rhoE, pos_rhoek, rhoi, xi, yi, Ri, ei, hi, diff_driv, pos_beta_cell, & 
                                         & Ji, lambda_vec, nb_te, pos_em, nb_eq, Dij_mat, Di, pos_di_cell, vdi_tol
    USE mod_neq_function_pointer,    ONLY: library_get_species_DiffFlux, library_get_molar_fractions,         & 
                                         & library_comp_tol, library_get_mass_fractions_from_molar_fractions
    USE mod_numerics_data,           ONLY: xil, xir, rhoil, rhoir, Ad, Bd

    IMPLICIT NONE

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
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: js_left   !< source term Jacobian with respect to the left state
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: js_cell   !< source term Jacobian with respect to the central state
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: js_right  !< source term Jacobian with respect to the right state 

    INTEGER :: i, j, k
    REAL(KIND=8) :: tmp, tmp1, tmp2, tmp3, tmp4, tmp5
    REAL(KIND=8) :: ov_dr, ov_rc
    REAL(KIND=8) :: dens, u, v, mu, mm, p, q_intk, rhoi_Vdi, T, Te, rho
    REAL(KIND=8) :: tau_rr, tau_rt, tau_tt
    REAL(KIND=8) :: pl, pr, Tl, Tr, ul, ur, vl, vr, ek
    REAL(KIND=8) :: du_dr, dv_dr, dT_dr
    REAL(KIND=8) :: q, q_Diff, q_Diff_int, q_Fourier_tr, q_Fourier_int
    REAL(KIND=8), DIMENSION(nb_temp) :: grad_T
    REAL(KIND=8), DIMENSION(nb_temp) :: bk
    REAL(KIND=8), DIMENSION(nb_temp) :: beta
    REAL(KIND=8), DIMENSION(nb_ns) :: sum_Dixi
    LOGICAL :: mass_diff_Jac 
    
    ! Explicit interface for subroutine stress_tensor_1D_SL_sph
    INTERFACE
      SUBROUTINE stress_tensor_1D_SL_sph (mu, r, u, v, du_dr, dv_dr, tau_rr, tau_rt, tau_tt)
        IMPLICIT NONE
        REAL(KIND=8), INTENT(IN) :: mu, r, u, v, du_dr, dv_dr
        REAL(KIND=8), INTENT(OUT) :: tau_rr, tau_rt, tau_tt
      END SUBROUTINE stress_tensor_1D_SL_sph
    END INTERFACE
    
    !PRINT*,'In "source_term_diff_neqNT_Te_1D_SL_sph_Jac.F90",implementation to be finshed!'
    !STOP
    
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
    rho = prop_cell(pos_rho_cell)
    T  = prop_cell(pos_T_cell)
    Te = prop_cell(pos_T_cell+nb_temp-1)
    p  = prop_cell(pos_pres_cell)
    u  = prop_cell(pos_u_cell)
    v  = prop_cell(pos_v_cell)
    mu = prop_cell(pos_mu_cell)
    DO i = 1,nb_temp
       lambda_vec(i) = prop_cell(pos_lambda_cell + i -1) 
       beta(i) = prop_cell(pos_beta_cell+i-1)
    ENDDO 

    ! Species data of central state 
    DO i = 1,nb_ns 
       rhoi(i) = u_cell(i)
       ei(i) = prop_cell(pos_ei_cell + i - 1)
       Di(i)   = prop_cell(pos_Di_cell + i - 1)
    ENDDO

    ! Compute the species enthalpies
    hi = ei
    DO i = 1,nb_ns
       IF (i==pos_em) THEN
         hi(i) = hi(i) + Ri(i)*Te
         hi((nb_temp-1)*nb_ns+i) = hi((nb_temp-1)*nb_ns+i) + Ri(i)*Te
       ELSE
         hi(i) = hi(i) + Ri(i)*T
       ENDIF
    ENDDO

    ! Species mass and molar fractions at cell interface
    CALL library_get_molar_fractions (rhoi, xi)
    CALL library_comp_tol(xi)
    
    ! Compute species mass fractions and mixture molar mass at cell interface 
    ! from species molar fractions
    CALL library_get_mass_fractions_from_molar_fractions(xi, mm, yi)
 
    ! Velocity, temperature, pressure and species molar fraction gradients
    ov_rc = 1.d0/r_c
    ov_dr = 2.d0/(vol_l + 2.d0*vol_c + vol_r)
    du_dr = ov_dr*(ur - ul)
    dv_dr = ov_dr*(vr - vl)
    DO i = 1,nb_temp
       grad_T(i) = (prop_right(pos_T_cell + i - 1) - prop_left(pos_T_cell + i - 1))*ov_dr
    ENDDO 

    ! Compute linearly independent diffusion driving forces
    diff_driv = 0.d0
    DO i = 1,nb_ns 
       diff_driv(i) = (xir(i) - xil(i))*ov_dr
    ENDDO
    
    ! Compute the species mass diffusion flux 
    CALL library_get_species_DiffFlux(p, T, Te, xi, diff_driv, Ji)
    
    ! Check the values of species diffusion velocities 
    mass_diff_Jac = .TRUE.
    tmp1 = T/p
    DO i = 1,nb_ns 
       IF (ABS(Ji(i)*tmp1*Ri(i)/xi(i)).LT.Vdi_tol) mass_diff_Jac = .FALSE.
    ENDDO 

    ! Stress tensor
    CALL stress_tensor_1D_SL_sph (mu, r_c, u, v, du_dr, dv_dr, tau_rr, tau_rt, tau_tt)

    ! Diffusive source term components
    ! Species continuity equations
    q_Fourier_tr = - lambda_vec(1)*grad_T(1)
    q_Diff       = 0.d0
    DO i = 1,nb_ns
       tmp    = Ji(i)
       s(i)   = -2.d0*tmp*ov_rc  
       q_Diff = q_Diff + tmp*hi(i)
    ENDDO

    ! Internal energy equations
    ! Fourier heat flux (components associated to internal energy)
    q_Fourier_int = 0.d0
    DO k = 1,nb_int_temp + nb_te
       q_intk = - lambda_vec(k + 1)*grad_T(k + 1)
       q_Fourier_int = q_Fourier_int + q_intk
       q_Diff_int    = 0.d0
       DO i = 1,nb_ns
          !q_Diff_int = q_Diff_int + Ji(i)*ei(nb_ns + nb_int_temp*(i - 1) + k)
          q_Diff_int = q_Diff_int + Ji(i)*hi(nb_ns*k + i)
       ENDDO
       s(pos_rhoek + k - 1) = -2.d0*(q_intk + q_Diff_int)*ov_rc
    ENDDO
    q = q_Fourier_tr + q_Fourier_int + q_Diff    

    ! Radial and circumferential momentum equations
    s(pos_rhou) = 2.d0*(tau_rr - tau_tt + tau_rt)*ov_rc
    s(pos_rhov) = (3.d0*tau_rt - tau_tt)*ov_rc

    ! Global energy equation
    s(pos_rhoE) = 2.d0*(u*(tau_rr + tau_rt) + v*tau_tt - q)*ov_rc

    ! Electron pressure work
    s(nb_eq) = s(nb_eq) - xi(pos_em)*p*(du_dr+2.0*(u+v)*ov_rc)
    
    ! Compute diffusion coefficient matrix 
    ! Initialization
    Dij_mat  = 0.d0
    sum_Dixi = 0.d0
    
    ! The nb_ns*nb_ns viscous Jacobian sum-matrix is computed only when 
    ! mass diffusion is occurring in order to avoid numerical problems  
    IF (mass_diff_Jac.EQV..TRUE.) THEN

      DO j = 1,nb_ns 
       
         tmp1 = Di(j)/xi(j)
         tmp2 = yi(j)*tmp1
         DO i = 1,j - 1 
            Dij_mat(i,j) = -tmp2
         ENDDO

         Dij_mat(j,j) = tmp1 - tmp2

         DO i = j + 1,nb_ns 
            Dij_mat(i,j) = -tmp2
         ENDDO

      ENDDO
   
      DO i = 1,nb_ns 
         DO k = 1,nb_ns 
            sum_Dixi(i) = sum_Dixi(i) + Dij_mat(i,k)*xi(k)
         ENDDO
      ENDDO

    ENDIF
    
    ! Jacobian Matrices
    ! A^d_s
    ! Columns 1 - nb_ns
    tmp1 = -2.0*ov_rc*ov_dr
    tmp2 = tmp1 * mu / rho * (2.0d0*u + v)
    tmp3 = tmp1 * mu / rho / 6.0d0 * (2.0d0*u + 9.0d0*v)
    tmp4 = tmp1 * mu / rho / 3.0d0 * (4.0d0*u + v)
    DO j = 1,nb_ns
        ! Species continuity
        tmp5 = tmp1 * mm / mi(j)
        bk = 0.0d0
        DO i = 1,nb_ns
            Ad(i,j) = tmp5 * yi(i) * (Dij_mat(i,j) - sum_Dixi(i))
            DO k = 1,nb_temp
                bk(k) = bk(k) + Ad(i,j)*hi((k-1)*nb_ns+i)
            END DO
        END DO
        
        tmp5 = lambda_vec(1)/beta(1)
        bk(1) = bk(1) - tmp1*tmp5*(0.5d0*u*u-ei(j))
        
        DO k = 2,nb_temp
            bk(1) = bk(1) + tmp1*(lambda_vec(k)/beta(k)-tmp5) * ei((k-1)*nb_ns+j)
            bk(k) = bk(k) + tmp1*lambda_vec(k)/beta(k) * ei((k-1)*nb_ns+j)
        END DO
        
        ! Momentum equations and total energy equation
        Ad(pos_rhou,j) = tmp2
        Ad(pos_rhov,j) = tmp3
        Ad(pos_rhoE,j) = bk(1) + tmp4
         
        ! Internal energy equations
        DO k = 1,nb_int_temp + nb_te
            Ad(pos_rhoek+k-1,j) = bk(k+1)
        END DO
    END DO
    
    ! rho*u column
    Ad(1:nb_eq,pos_rhou)  = 0.0d0
    Ad(pos_rhou,pos_rhou) = -tmp1*2.0d0*mu/rho
    Ad(pos_rhov,pos_rhou) = -tmp1*mu/rho/3.0d0
    Ad(pos_rhoE,pos_rhou) = tmp1*(lambda_vec(1)/beta(1)*u + mu/rho*4.0d0/3.0d0*(u-0.5d0*v))
    
    ! rho*v column
    Ad(1:nb_eq,pos_rhov)  = 0.0d0
    Ad(pos_rhou,pos_rhov) = -tmp1*mu/rho
    Ad(pos_rhov,pos_rhov) = -tmp1*1.5d0*mu/rho
    Ad(pos_rhoE,pos_rhov) = -tmp1*u*mu/rho
    
    ! rho*E column
    Ad(1:nb_eq,pos_rhoE)  = 0.0d0
    tmp5 = lambda_vec(1)/beta(1)
    Ad(pos_rhoE,pos_rhoE) = tmp1*tmp5
    
    ! Internal energy columns
    DO k = 1,nb_int_temp + nb_te
        Ad(1:nb_eq,pos_rhoek+k-1)  = 0.0d0
        Ad(pos_rhoE,pos_rhoek+k-1) = tmp1*(tmp5-lambda_vec(k+1)/beta(k+1))
        Ad(pos_rhoek+k-1,pos_rhoek+k-1) = -tmp1*lambda_vec(k+1)/beta(k+1)
    END DO
    
    ! Bd matrix
    tmp1 = mu/rho*ov_rc*ov_rc
    tmp2 = tmp1*(u + v)
    tmp3 = tmp1*(7.0d0*u*u + 5.0d0*u*v - 2.0d0*v*v)
    Bd = 0.0d0
    DO j = 1,nb_ns
        Bd(pos_rhou,j) = 6.0d0*tmp2
        Bd(pos_rhov,j) = 11.0d0/3.0d0*tmp2
        Bd(pos_rhoE,j) = 4.0d0/3.0d0*tmp3
    END DO
    
    Bd(pos_rhou,pos_rhou) = -tmp1*6.0d0
    Bd(pos_rhov,pos_rhou) = -tmp1*11.0d0/3.0d0
    Bd(pos_rhoE,pos_rhou) = tmp1*2.0d0/3.0d0*(14.0d0*u-5.0d0*u)
    
    Bd(pos_rhou,pos_rhov) = -tmp1*6.0d0
    Bd(pos_rhov,pos_rhov) = -tmp1*11.0d0/3.0d0
    Bd(pos_rhoE,pos_rhov) = tmp1*2.0d0/3.0d0*(4.0d0*v-5.0d0*u)
    
    ! Left, central and right state diffusive source term Jacobians
    js_left  = -Ad
    js_cell  = Bd
    js_right = Ad

  END SUBROUTINE source_term_diff_neqNT_Te_1D_SL_sph_Jac
!------------------------------------------------------------------------------!
