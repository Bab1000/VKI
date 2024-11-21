!------------------------------------------------------------------------------!
!> This subroutine computes the diffusive flux and the related Jacobians when solving the 1D stagnation line equations 
!! for nonequilibrium flows (cylinder case - 1 temperature).
  SUBROUTINE ns_flux_neq1T_1D_SL_cyl_Jac (r_l, r_r, vol_l, vol_r, left_data, right_data, u_left, u_right, fd, jfdl, jfdr)

    USE mod_general_data,            ONLY: nb_ns, nb_dim, pos_u_cell, pos_v_cell, pos_ei_cell, pos_T_cell,  & 
                                         & pos_pres_cell, pos_mu_cell, pos_lambda_cell, pos_beta_cell,      & 
                                         & pos_Di_cell, pos_rhou, pos_rhov, pos_rhoE, Vdi_tol, mi, rhoi,    & 
                                         & xi, yi, Ri, ei, hi, mi, Dij_mat, Di, diff_driv, Ji
    USE mod_neq_function_pointer,    ONLY: library_get_species_DiffFlux, library_get_molar_fractions,       & 
                                         & library_get_mass_fractions_from_molar_fractions, library_comp_tol
    USE mod_numerics_data,           ONLY: xil, xir, rhoil, rhoir, Ad, Bd

    IMPLICIT NONE

    INTEGER :: i, j, k
    REAL(KIND=8), PARAMETER :: coeff1 = 2.d0/3.d0
    REAL(KIND=8), PARAMETER :: coeff2 = 2.d0*coeff1
    REAL(KIND=8) :: tmp1, tmp2, tmp3
    REAL(KIND=8) :: ov_dr, ov_r
    REAL(KIND=8) :: dens, beta, ek, lambda, mu, mm, p, r, rho, rhoi_Vdi, T, u, v, us
    REAL(KIND=8) :: tau_rr, tau_rt, tau_tt
    REAL(KIND=8) :: betal, betar, mul, mur, lambdal, lambdar, pl, pr, Tl, Tr, ul, ur, vl, vr
    REAL(KIND=8) :: du_dr, dv_dr, dT_dr
    REAL(KIND=8) :: q, q_Diff, q_Fourier
    REAL(KIND=8) :: lam_ov_beta, ov_rho, mu_ov_rho, u_mu_ov_rho, coeff_mu_ov_rho, coeff_u_mu_ov_rho,   & 
                  & coeff_us_mu_ov_rho, v_mu_ov_rho, us_mu_ov_rho, vel_sum, mu_ov_rho_vel_sum,         & 
                  & coeff_mu_ov_rho_vel_sum, m2_coeff_mu_ov_rho_vel_sum_u
    REAL(KIND=8), DIMENSION(nb_ns) :: sum_Diyi_ov_mi
    LOGICAL :: mass_diff_Jac

    REAL(KIND=8), INTENT(IN) :: r_l                       !< radial position of left state
    REAL(KIND=8), INTENT(IN) :: r_r                       !< radial position of rigth state
    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state 
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: fd         !< diffusive flux
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfdl     !< diffusive flux Jacobian with respect to the left state
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfdr     !< diffusive flux Jacobian with respect to the right state

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

    ! Temperature, pressure, velocity components, dynamic viscosity, thermal conductivity, beta factor, 
    ! species densities and molar fractions of left and right states
    ! Left state
    Tl      = left_data(pos_T_cell)
    pl      = left_data(pos_pres_cell)
    ul      = left_data(pos_u_cell)
    vl      = left_data(pos_v_cell)
    mul     = left_data(pos_mu_cell)
    betal   = left_data(pos_beta_cell)
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
    betar   = right_data(pos_beta_cell)
    lambdar = right_data(pos_lambda_cell)

    DO i = 1,nb_ns 
       rhoir(i) = u_right(i)
    ENDDO
    CALL library_get_molar_fractions(rhoir, xir)
    CALL library_comp_tol(xir)

    ! Temperature, pressure, velocity, dynamic viscosity, beta factor and thermal conductivity at cell interface
    T  = 0.5d0*(Tl + Tr)
    p  = 0.5d0*(pr + pr)
    u  = 0.5d0*(ul + ur)
    v  = 0.5d0*(vl + vr)
    mu = 0.5d0*(mul + mur)
    beta = 0.5d0*(betal + betar)
    lambda = 0.5d0*(lambdal + lambdar)
 
    ! Species densities, average binary diffusion coefficients and specific energies at cell interface  
    rho = 0.d0
    DO i = 1,nb_ns 
       dens = 0.5d0*(rhoil(i) + rhoir(i))
       rho  = rho + dens
       rhoi(i) = dens 
       tmp1    = 0.5d0*(left_data(pos_ei_cell + i - 1) + right_data(pos_ei_cell + i - 1))
       ei(i)   = tmp1
       hi(i)   = tmp1 + Ri(i)*T
       Di(i)   = 0.5d0*(left_data(pos_Di_cell + i - 1) + right_data(pos_Di_cell + i - 1))
    ENDDO

    ! Species molar fractions at cell interface
    CALL library_get_molar_fractions (rhoi, xi)
    CALL library_comp_tol(xi)

    ! Compute species mass fractions and mixture molar mass at cell interface 
    ! from species molar fractions
    CALL library_get_mass_fractions_from_molar_fractions(xi, mm, yi)
 
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
    
    ! Compute the species mass diffusion flux
    CALL library_get_species_DiffFlux(p, T, T, xi, diff_driv, Ji)

    ! Check the values of species diffusion velocities 
    mass_diff_Jac = .TRUE.
    tmp1 = T/p
    DO i = 1,nb_ns 
       IF (ABS(Ji(i)*tmp1*Ri(i)/xi(i)).LT.Vdi_tol) mass_diff_Jac = .FALSE.
    ENDDO

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

    ! Diffusive flux Jacobians 
    ! Ad matrix
    ! Common factors
    ! For the evaluation of the Ad and Bd matrices the radial velocity at the cell interface 
    ! is set equal to the left state velocity for sake of robustness 
    u  = ul
    us = u**2
    ek = 0.5d0*us 
    ov_rho = ov_dr/rho 
    mu_ov_rho = mu*ov_rho
    u_mu_ov_rho  = u*mu_ov_rho
    us_mu_ov_rho = u*u_mu_ov_rho
    v_mu_ov_rho  = v*mu_ov_rho
    lam_ov_beta  = ov_dr*lambda/beta
    coeff_mu_ov_rho    = coeff2*mu_ov_rho
    coeff_u_mu_ov_rho  = u*coeff_mu_ov_rho
    coeff_us_mu_ov_rho = us*coeff_mu_ov_rho

    ! Compute diffusion coefficient matrix 
    ! Initialization
    Dij_mat        = 0.d0
    sum_Diyi_ov_mi = 0.d0

    ! The nb_ns*nb_ns viscous Jacobian sum-matrix is computed only when 
    ! mass diffusion is occurring in order to avoid numerical problems  
    IF (mass_diff_Jac.EQV..TRUE.) THEN

      DO j = 1,nb_ns 
       
         tmp1 = Di(j)/xi(j)
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

    ! Column j,  j = 1,..,nb_ns 
    DO j = 1,nb_ns 
      
       tmp1 = 0.d0
       tmp2 = mm/mi(j)*ov_dr
       DO i = 1,nb_ns 
          tmp3 = yi(i)*tmp2*(Dij_mat(i,j) - mm*sum_Diyi_ov_mi(i))
          Ad(i,j) = tmp3
          tmp1    = tmp1 + tmp3*hi(i)
       ENDDO

       Ad(pos_rhou,j) = - coeff_u_mu_ov_rho
       Ad(pos_rhov,j) = - v_mu_ov_rho
       Ad(pos_rhoE,j) = tmp1 - coeff_us_mu_ov_rho + lam_ov_beta*(ek - ei(j))

    ENDDO

    ! Column nb_ns + 1
    DO i = 1,nb_ns 
       Ad(i,pos_rhou) = 0.d0
    ENDDO
    Ad(pos_rhou,pos_rhou) = coeff_mu_ov_rho
    Ad(pos_rhov,pos_rhou) = 0.d0
    Ad(pos_rhoE,pos_rhou) = coeff_u_mu_ov_rho - lam_ov_beta*u

    ! Column nb_ns + 2
    DO i = 1,nb_ns 
       Ad(i,pos_rhov) = 0.d0
    ENDDO 
    Ad(pos_rhou,pos_rhov) = 0.d0
    Ad(pos_rhov,pos_rhov) = mu_ov_rho
    Ad(pos_rhoE,pos_rhov) = 0.d0

    ! Column nb_ns + 3
    DO i = 1,nb_ns 
       Ad(i,pos_rhoE) = 0.d0
    ENDDO 
    Ad(pos_rhou,pos_rhoE) = 0.d0
    Ad(pos_rhov,pos_rhoE) = 0.d0
    Ad(pos_rhoE,pos_rhoE) = lam_ov_beta

    ! Bd matrix
    ! Common factors
    vel_sum   = (u + v)
    mu_ov_rho = mu/(rho*r) 
    coeff_mu_ov_rho = coeff1*mu_ov_rho
    mu_ov_rho_vel_sum = mu_ov_rho*vel_sum
    coeff_mu_ov_rho_vel_sum = coeff_mu_ov_rho*vel_sum     
    m2_coeff_mu_ov_rho_vel_sum_u = 2.d0*coeff_mu_ov_rho_vel_sum*u

    ! Column j,  j = 1,..,nb_ns 
    DO j = 1,nb_ns 

       DO i = 1,nb_ns 
          Bd(i,j) = 0.d0
       ENDDO

       Bd(pos_rhou,j) = coeff_mu_ov_rho_vel_sum
       Bd(pos_rhov,j) = mu_ov_rho_vel_sum
       Bd(pos_rhoE,j) = m2_coeff_mu_ov_rho_vel_sum_u 

    ENDDO

    ! Column nb_ns + 1
    DO i = 1,nb_ns 
       Bd(i,pos_rhou) = 0.d0
    ENDDO
    Bd(pos_rhou,pos_rhou) = - coeff_mu_ov_rho
    Bd(pos_rhov,pos_rhou) = - mu_ov_rho
    Bd(pos_rhoE,pos_rhou) = - coeff_mu_ov_rho*(u + vel_sum)

    ! Column nb_ns + 2
    DO i = 1,nb_ns 
       Bd(i,pos_rhov) = 0.d0
    ENDDO
    Bd(pos_rhou,pos_rhov) = - coeff_mu_ov_rho
    Bd(pos_rhov,pos_rhov) = - mu_ov_rho
    Bd(pos_rhoE,pos_rhov) = - coeff_mu_ov_rho*u

    ! Column nb_ns + 3
    DO i = 1,nb_ns 
       Bd(i,pos_rhoE) = 0.d0
    ENDDO
    Bd(pos_rhou,pos_rhoE) = 0.d0
    Bd(pos_rhov,pos_rhoE) = 0.d0
    Bd(pos_rhoE,pos_rhoE) = 0.d0  

    ! Left and right state diffusive flux Jacobians
    jfdr = Ad + Bd
    jfdl = - Ad + Bd

  END SUBROUTINE ns_flux_neq1T_1D_SL_cyl_Jac
!------------------------------------------------------------------------------!

