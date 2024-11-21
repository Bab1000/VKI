!------------------------------------------------------------------------------!
!> This subroutine computes the diffusive source term and the related Jacobians for 1D stagnation line flows 
!! for nonequilibrium flows (sphere case - 1 temperature).
  SUBROUTINE source_term_diff_neq1T_1D_SL_sph_Jac(r_c, vol_l, vol_c, vol_r, prop_left, prop_cell, prop_right, & 
                                                & u_left, u_cell, u_right, s, js_left, js_cell, js_right)

    USE mod_general_data,            ONLY: nb_ns, nb_dim, pos_u_cell, pos_v_cell, pos_ei_cell, pos_T_cell,  & 
                                         & pos_pres_cell, pos_mu_cell, pos_lambda_cell, pos_beta_cell,      & 
                                         & pos_Di_cell, pos_rho_cell, pos_rhou, pos_rhov, pos_rhoE,         & 
                                         & Vdi_tol, mi, rhoi, xi, yi, Ri, ei, hi, mi, Dij_mat, Di,          & 
                                         & diff_driv, Ji
    USE mod_neq_function_pointer,    ONLY: library_get_species_DiffFlux, library_get_molar_fractions,       & 
                                         & library_get_mass_fractions_from_molar_fractions, library_comp_tol
    USE mod_numerics_data,           ONLY: xil, xir, rhoil, rhoir, Ad, Bd

    IMPLICIT NONE

    INTEGER :: i, j, k
    REAL(KIND=8), PARAMETER :: coeff1 = 2.d0/3.d0
    REAL(KIND=8), PARAMETER :: coeff2 = 2.d0*coeff1
    REAL(KIND=8), PARAMETER :: coeff3 = 11.d0/3.d0
    REAL(KIND=8) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
    REAL(KIND=8) :: ov_dr, ov_rc, ov_drrc
    REAL(KIND=8) :: beta, ek, lambda, mu, mm, p, rho, rhoi_Vdi, u, v, uv, us, vs, vel_sum, T
    REAL(KIND=8) :: tau_rr, tau_rt, tau_tt
    REAL(KIND=8) :: pl, pr, Tl, Tr, ul, ur, vl, vr
    REAL(KIND=8) :: du_dr, dv_dr, dT_dr
    REAL(KIND=8) :: q, q_Diff, q_Fourier
    REAL(KIND=8) :: ov_rho, mu_ov_rho, m2_mu_ov_rho, m4_mu_ov_rho, m6_mu_ov_rho, lam_ov_beta
    REAL(KIND=8) :: coeff1_mu_ov_rho, coeff3_mu_ov_rho
    REAL(KIND=8), DIMENSION(nb_ns) :: sum_Diyi_ov_mi
    LOGICAL :: mass_diff_Jac 

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

    ! Density, pressure, temperature, velocity components, species densities, beta factor and specific energies
    rho = prop_cell(pos_rho_cell)
    p   = prop_cell(pos_pres_cell)
    T   = prop_cell(pos_T_cell)
    u   = prop_cell(pos_u_cell)
    v   = prop_cell(pos_v_cell)
    mu  = prop_cell(pos_mu_cell)
    beta   = prop_cell(pos_beta_cell)
    lambda = prop_cell(pos_lambda_cell)

    ! Species densities and specific energies and enthalpies of central state 
    DO i = 1,nb_ns 
       rhoi(i) = u_cell(i)
       tmp1    = prop_cell(pos_ei_cell + i - 1)
       ei(i)   = tmp1
       hi(i)   = tmp1 + Ri(i)*T
       Di(i)   = prop_cell(pos_Di_cell + i - 1)
    ENDDO

    ! Species molar fractions at cell interface
    CALL library_get_molar_fractions (rhoi, xi)
    CALL library_comp_tol(xi) 

    ! Compute species mass fractions and mixture molar mass at cell interface 
    ! from species molar fractions
    CALL library_get_mass_fractions_from_molar_fractions(xi, mm, yi)

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
    q_Fourier = - lambda*dT_dr
    q_Diff    = 0.d0
    DO i = 1,nb_ns
       rhoi_Vdi = Ji(i)
       s(i)     = - 2.d0*rhoi_Vdi*ov_rc 
       q_Diff   = q_Diff + rhoi_Vdi*hi(i)
    ENDDO
    q = q_Fourier + q_Diff

    ! Radial and circumferential momentum equations
    s(pos_rhou) = 2.d0*(tau_rr - tau_tt + tau_rt)*ov_rc
    s(pos_rhov) = (3.d0*tau_rt - tau_tt)*ov_rc

    ! Global energy equation
    s(pos_rhoE) = 2.d0*(u*(tau_rr + tau_rt) + v*tau_tt - q)*ov_rc 

    ! Compute diffusion coefficient matrix 
    ! Initialization
    Dij_mat        = 0.d0
    sum_Diyi_ov_mi = 0.d0

    ! The nb_ns*nb_ns viscous Jacobian sum-matrix is computed only when 
    ! mass diffusion is occurring in order to avoid numerical problems  
    IF (mass_diff_Jac.EQV..TRUE.) THEN

      ! Watch out: all the terms in the species diffusion coefficient matrix 
      ! are multiplied by 2 as needed by the Jacobian formulation 
      ! (this operation is performed in order to save computational time)
      DO j = 1,nb_ns 
       
         tmp1 = 2.d0*Di(j)/xi(j)
         tmp2 = yi(j)*tmp1
         DO i = 1,j - 1 
            Dij_mat(i,j) = - tmp2
         ENDDO

         Dij_mat(j,j) = tmp1 - tmp2

         DO i = j + 1,nb_ns 
            Dij_mat(i,j) = - tmp2
         ENDDO

      ENDDO
   
      DO i = 1,nb_ns 
         DO k = 1,nb_ns 
            sum_Diyi_ov_mi(i) = sum_Diyi_ov_mi(i) + Dij_mat(i,k)*yi(k)/mi(k)
         ENDDO
      ENDDO

    ENDIF

    ! Diffusive source term Jacobians
    ! Ad matrix
    ! Common factors
    us = u**2
    vs = v**2
    ek = 0.5d0*us
    uv = u*v
    vel_sum = u + v
    ov_drrc = ov_dr*ov_rc
    ov_rho = ov_drrc/rho 
    mu_ov_rho = mu*ov_rho
    m2_mu_ov_rho = 2.d0*mu_ov_rho
    coeff1_mu_ov_rho = coeff1*mu_ov_rho
    lam_ov_beta = 2.d0*ov_drrc*lambda/beta

    tmp1 = - m2_mu_ov_rho*(vel_sum + u)
    tmp2 = - mu_ov_rho*(2.d0*u + 9.d0*v)/3.d0
    tmp3 = - coeff1_mu_ov_rho*u*(vel_sum + 3.d0*u)

    ! Column j,  j = 1,..,nb_ns 
    DO j = 1,nb_ns 
 
       tmp4 = 0.d0
       tmp5 = mm*ov_drrc/mi(j)
       DO i = 1,nb_ns 
          tmp6 = yi(i)*tmp5*(Dij_mat(i,j) - mm*sum_Diyi_ov_mi(i))
          Ad(i,j) = tmp6
          tmp4    = tmp4 + tmp6*hi(i)
       ENDDO

       Ad(pos_rhou,j) = tmp1
       Ad(pos_rhov,j) = tmp2
       Ad(pos_rhoE,j) = tmp4 + tmp3 + lam_ov_beta*(ek - ei(j)) 
    
    ENDDO

    ! Column nb_ns + 1
    DO i = 1,nb_ns 
       Ad(i,pos_rhou) = 0.d0
    ENDDO
    Ad(pos_rhou,pos_rhou) = m2_mu_ov_rho + m2_mu_ov_rho
    Ad(pos_rhov,pos_rhou) = coeff1_mu_ov_rho
    Ad(pos_rhoE,pos_rhou) = coeff2*mu_ov_rho*(2.d0*u - v) - lam_ov_beta*u

    ! Column nb_ns + 2
    DO i = 1,nb_ns 
       Ad(i,pos_rhov)  = 0.d0
    ENDDO
    Ad(pos_rhou,pos_rhov) = m2_mu_ov_rho
    Ad(pos_rhov,pos_rhov) = mu_ov_rho + m2_mu_ov_rho
    Ad(pos_rhoE,pos_rhov) = m2_mu_ov_rho*u

    ! Column nb_ns + 3
    DO i = 1,nb_ns 
       Ad(i,pos_rhoE)  = 0.d0
    ENDDO
    Ad(pos_rhou,pos_rhoE) = 0.d0
    Ad(pos_rhov,pos_rhoE) = 0.d0
    Ad(pos_rhoE,pos_rhoE) = lam_ov_beta

    ! Bd matrix
    ! Common factors
    mu_ov_rho = (mu/rho)*ov_rc**2
    m6_mu_ov_rho = 6.d0*mu_ov_rho
    coeff1_mu_ov_rho = coeff1*mu_ov_rho
    coeff3_mu_ov_rho = coeff3*mu_ov_rho

    tmp1 = m6_mu_ov_rho*vel_sum
    tmp2 = coeff3_mu_ov_rho*vel_sum
    tmp3 = coeff2*mu_ov_rho*(7.d0*us + 5.d0*uv - 2.d0*vs)

    ! Column j,  j = 1,..,nb_ns
    DO j = 1,nb_ns  

       DO i = 1,nb_ns 
          Bd(i,j) = 0.d0
       ENDDO   

       Bd(pos_rhou,j) = tmp1
       Bd(pos_rhov,j) = tmp2
       Bd(pos_rhoE,j) = tmp3
 
    ENDDO

    ! Column nb_ns + 2
    DO i = 1,nb_ns
       Bd(i,pos_rhou) = 0.d0
    ENDDO
    Bd(pos_rhou,pos_rhou) = - m6_mu_ov_rho
    Bd(pos_rhov,pos_rhou) = - coeff3_mu_ov_rho
    Bd(pos_rhoE,pos_rhou) = coeff1_mu_ov_rho*(14.d0*u - 5.d0*v)

    ! Column nb_ns + 2
    DO i = 1,nb_ns 
       Bd(i,pos_rhov) = 0.d0
    ENDDO
    Bd(pos_rhou,pos_rhov) = - m6_mu_ov_rho
    Bd(pos_rhov,pos_rhov) = - coeff3_mu_ov_rho
    Bd(pos_rhoE,pos_rhov) = coeff1_mu_ov_rho*(4.d0*v - 5.d0*u)

    ! Column nb_ns + 3 
    DO i = 1,nb_ns 
       Bd(i,pos_rhoE) = 0.d0
    ENDDO
    Bd(pos_rhou,pos_rhoE) = 0.d0
    Bd(pos_rhov,pos_rhoE) = 0.d0
    Bd(pos_rhoE,pos_rhoE) = 0.d0
    
    ! Left, central and right state diffusive source term Jacobians
    js_left  = - Ad
    js_cell  = Bd
    js_right = Ad
  
  END SUBROUTINE source_term_diff_neq1T_1D_SL_sph_Jac
!------------------------------------------------------------------------------!
