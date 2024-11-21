!------------------------------------------------------------------------------!
!> This subrotine provides the numerical flux Jacobian for 1D stagnation line nonequilibrium gas flows
!! according to the positive and negative split of the inviscid flux Jacobian.
!! In this case the flow is characterized by N temperatures. However, no separated temperature exist for the free electrons.
  SUBROUTINE Pos_Neg_split_Jac_neqNT_1D_SL (nx, vol_l, vol_r, left_data, right_data, u_left, u_right, jfl, jfr)

    USE mod_general_data,             ONLY: nb_ns, nb_eq, pos_alpha_cell, pos_c_cell, pos_ei_cell, pos_ek_cell, &
                                          & pos_gamma_cell, pos_h0_cell, pos_rho_cell, pos_T_cell, pos_u_cell,  &
                                          & pos_v_cell, pos_rhou, pos_rhov, pos_rhoE, yi, epsi, sigmai, Ri,     &
                                          & pos_eki_cell, eintk, nb_int_temp
    USE mod_neq_function_pointer,     ONLY: library_get_mass_fractions

    IMPLICIT NONE

    INTEGER :: i, j
    REAL(KIND=8) :: tmp1, tmp2
    REAL(KIND=8) :: alpha, eps, gamma, nx2
    REAL(KIND=8) :: c, cn, ek, h0, rho, v, u, u_vn, v_vn, vn, vn2, T, p
    REAL(KIND=8) :: alpha_ov_c2, ov_c, ov_c2, yi_epsj_ov_c2, yi_vn_ov_c, epsj_ov_c2, ov_alpha
    REAL(KIND=8) :: l1, l2, l3, l1p, l2p, l3p, l1m, l2m, l3m, l1m_u, l1p_u, l1p_ov_gamma_m1, l1m_ov_gamma_m1
    REAL(KIND=8) :: eig_diff, eig_sum, eig_diff_vn_ov_c, eig_sum_m_l1p, eig_sum_m_l1m, u_eig_sum_m_l1p,         &
                  & v_eig_sum_m_l1p, h0_eig_sum_m_l1p, u_eig_sum_m_l1m, v_eig_sum_m_l1m, h0_eig_sum_m_l1m,      &
                  & eig_sum_vn_nx, eig_diff_ov_c, eig_diff_nx_ov_c, eig_sum_vn2, &
                  & u_eig_sum_m_l1p_alpha_ov_c2, u_eig_sum_m_l1m_alpha_ov_c2, &
                  & v_eig_sum_m_l1p_alpha_ov_c2, v_eig_sum_m_l1m_alpha_ov_c2

    REAL(KIND=8), INTENT(IN) :: nx                        !< normal to the cell interface
    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfl      !< numerical flux Jacobian with respect to the left state 
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfr      !< numerical flux Jacobian with respect to the right state

    jfl = 0.d0
    jfr = 0.d0

    ! Useful quantity 
    nx2 = nx**2

    ! Left state
    ! A^+ matrix (positive split Jacobian) 
    ! Physical data 
    T     = left_data(pos_T_cell)
    u     = left_data(pos_u_cell)
    v     = left_data(pos_v_cell)
    h0    = left_data(pos_h0_cell)
    c     = left_data(pos_c_cell)
    alpha = left_data(pos_alpha_cell)
    rho   = left_data(pos_rho_cell)
    ek    = left_data(pos_ek_cell)
    gamma = left_data(pos_gamma_cell)

    ! Species mass fractions
    CALL library_get_mass_fractions (rho, u_left(1:nb_ns), yi)

    ! Eigenvalues 
    vn  = u*nx
    l1  = vn
    l2  = vn - c
    l3  = vn + c

    ! Positive split eigenvalues 
    l1p = 0.5d0*(l1 + ABS(l1))
    l2p = 0.5d0*(l2 + ABS(l2))
    l3p = 0.5d0*(l3 + ABS(l3))

    ! Common factors 
    eig_sum  = 0.5d0*(l3p + l2p)
    eig_diff = 0.5d0*(l3p - l2p)

    u_vn  = u*vn
    v_vn  = v*vn
    vn2   = vn**2
    ov_c  = 1.d0/c
    ov_c2 = ov_c/c
    ov_alpha    = 1.d0/alpha
    alpha_ov_c2 = alpha*ov_c2

    l1p_u = l1p*u

    eig_sum_m_l1p    = eig_sum - l1p
    eig_sum_vn_nx    = eig_sum*vn*nx
    eig_sum_vn2      = eig_sum*vn2

    u_eig_sum_m_l1p  = u*eig_sum_m_l1p
    v_eig_sum_m_l1p  = v*eig_sum_m_l1p
    h0_eig_sum_m_l1p = h0*eig_sum_m_l1p

    u_eig_sum_m_l1p_alpha_ov_c2 = u_eig_sum_m_l1p*alpha_ov_c2
    v_eig_sum_m_l1p_alpha_ov_c2 = v_eig_sum_m_l1p*alpha_ov_c2

    eig_diff_ov_c    = eig_diff*ov_c
    eig_diff_nx_ov_c = eig_diff_ov_c*nx
    eig_diff_vn_ov_c = eig_diff_ov_c*vn

    l1p_ov_gamma_m1  = l1p/(gamma - 1.d0)

    ! Energy related data
    DO i = 1,nb_ns
       tmp1 = Ri(i)*T
       tmp2 = left_data(pos_ei_cell + i - 1)
       epsi(i)   = tmp1 - alpha*tmp2
       sigmai(i) = tmp2 - tmp1*ov_alpha
    ENDDO

    epsi   = epsi + alpha*ek
    sigmai = sigmai + ek

    DO i = 1, nb_int_temp
      eintk(i) = 0.d0
      DO j = 1, nb_ns
        tmp1 = left_data(pos_eki_cell + j - 1 +nb_ns*(i-1))
        epsi(j) = epsi(j) + alpha*tmp1
        sigmai(j) = sigmai(j) - tmp1
        eintk(i) = eintk(i) + tmp1*yi(j)
      ENDDO
    ENDDO

    ! Column j,  j = 1,..,nb_ns 
    DO j = 1,nb_ns

       epsj_ov_c2 = epsi(j)*ov_c2

       DO i = 1,j - 1
          jfl(i,j) = yi(i)*(eig_sum_m_l1p*epsj_ov_c2 - eig_diff_vn_ov_c)
       ENDDO

       jfl(j,j) = l1p + yi(j)*(eig_sum_m_l1p*epsj_ov_c2 - eig_diff_vn_ov_c)

       DO i = j + 1,nb_ns
          jfl(i,j) = yi(i)*(eig_sum_m_l1p*epsj_ov_c2 - eig_diff_vn_ov_c)
       ENDDO

       jfl(pos_rhou,j) = epsj_ov_c2*u_eig_sum_m_l1p  + l1p_u - eig_sum_vn_nx + eig_diff_ov_c*(nx*epsi(j) - u_vn)
       jfl(pos_rhov,j) = epsj_ov_c2*v_eig_sum_m_l1p - eig_diff_ov_c*v_vn
       jfl(pos_rhoE,j) = epsj_ov_c2*h0_eig_sum_m_l1p + l1p*sigmai(j) + epsi(j)*l1p_ov_gamma_m1 &
                             & - eig_sum_vn2 + eig_diff_vn_ov_c*(epsi(j) - h0)

       DO i = 1, nb_int_temp
          jfl(pos_rhoE+i, j) = (epsj_ov_c2*eig_sum_m_l1p - eig_diff_vn_ov_c)*eintk(i)
       ENDDO

    ENDDO

    ! Column nb_ns + 1 
    tmp1 = eig_diff_nx_ov_c  - u_eig_sum_m_l1p_alpha_ov_c2
    DO i = 1,nb_ns
       jfl(i,pos_rhou) = tmp1*yi(i)
    ENDDO

    jfl(pos_rhou,pos_rhou) = - u_eig_sum_m_l1p_alpha_ov_c2*u  + eig_sum*nx2 + eig_diff_vn_ov_c*(1.d0 - alpha)
    jfl(pos_rhov,pos_rhou) = - u_eig_sum_m_l1p_alpha_ov_c2*v  + eig_diff_nx_ov_c*v
    jfl(pos_rhoE,pos_rhou) = - u_eig_sum_m_l1p_alpha_ov_c2*h0 - l1p_u + eig_sum*vn*nx + eig_diff_ov_c*(h0*nx - alpha*u_vn)

    DO i = 1, nb_int_temp
       jfl(pos_rhoE+i, pos_rhou) = (eig_diff_ov_c - eig_sum_m_l1p*alpha_ov_c2*u)*eintk(i)
    ENDDO

    ! Column nb_ns + 2
    DO i = 1,nb_ns
       jfl(i,pos_rhov) = 0.d0
    ENDDO

    jfl(pos_rhou,pos_rhov) = 0.d0
    jfl(pos_rhov,pos_rhov) = l1p
    jfl(pos_rhoE,pos_rhov) = 0.d0

    DO i = 1, nb_int_temp
       jfl(pos_rhoE+i, pos_rhov) = 0.d0
    ENDDO

    ! Column nb_ns + 3
    tmp1 = eig_sum_m_l1p*alpha_ov_c2
    DO i = 1,nb_ns
       jfl(i,pos_rhoE) = tmp1*yi(i)
    ENDDO

    jfl(pos_rhou,pos_rhoE) = u_eig_sum_m_l1p_alpha_ov_c2  + alpha*eig_diff_nx_ov_c
    jfl(pos_rhov,pos_rhoE) = v_eig_sum_m_l1p_alpha_ov_c2
    jfl(pos_rhoE,pos_rhoE) = h0_eig_sum_m_l1p*alpha_ov_c2 + l1p + eig_diff_vn_ov_c*alpha

    DO i = 1, nb_int_temp
       jfl(pos_rhoE+i, pos_rhoE) = eintk(i)*alpha_ov_c2*eig_sum_m_l1p
    ENDDO

    ! Column j = nb_ns+4,..,nb_eq
    DO j = 1, nb_int_temp

       DO i = 1, nb_ns
          jfl(i, pos_rhoE+j) = - eig_sum_m_l1p*alpha_ov_c2*yi(i)
       ENDDO

       jfl(pos_rhou, pos_rhoE+j) = - u_eig_sum_m_l1p_alpha_ov_c2 - eig_diff_ov_c*alpha 
       jfl(pos_rhov, pos_rhoE+j) = - v_eig_sum_m_l1p_alpha_ov_c2 
       jfl(pos_rhoE, pos_rhoE+j) = - h0_eig_sum_m_l1p*alpha_ov_c2 - eig_diff_vn_ov_c*alpha

       DO i = 1, j-1
          jfl(pos_rhoE+i, pos_rhoE+j) = - eig_sum_m_l1p*alpha_ov_c2*eintk(i) 
       ENDDO
       jfl(pos_rhoE+j, pos_rhoE+j) = - eig_sum_m_l1p*alpha_ov_c2*eintk(i) + l1p
       DO i = j+1, nb_int_temp
          jfl(pos_rhoE+i, pos_rhoE+j) = - eig_sum_m_l1p*alpha_ov_c2*eintk(i) 
       ENDDO

    ENDDO

    ! Right state
    ! A^- matrix (negative split Jacobian) 
    ! Physical data 
    T     = right_data(pos_T_cell)
    u     = right_data(pos_u_cell)
    v     = right_data(pos_v_cell)
    h0    = right_data(pos_h0_cell)
    c     = right_data(pos_c_cell)
    alpha = right_data(pos_alpha_cell)
    rho   = right_data(pos_rho_cell)
    ek    = right_data(pos_ek_cell)
    gamma = right_data(pos_gamma_cell)

    ! Species mass fractions
    CALL library_get_mass_fractions (rho, u_right(1:nb_ns), yi)

    ! Eigenvalues 
    vn  = u*nx
    l1  = vn
    l2  = vn - c
    l3  = vn + c

    ! Negative split eigenvalues 
    l1m = 0.5d0*(l1 - ABS(l1))
    l2m = 0.5d0*(l2 - ABS(l2))
    l3m = 0.5d0*(l3 - ABS(l3))

    ! Common factors 
    eig_sum  = 0.5d0*(l3m + l2m)
    eig_diff = 0.5d0*(l3m - l2m)

    u_vn  = u*vn
    v_vn  = v*vn
    vn2   = vn**2
    ov_c  = 1.d0/c
    ov_c2 = ov_c/c
    ov_alpha    = 1.d0/alpha
    alpha_ov_c2 = alpha*ov_c2

    l1m_u = l1m*u

    eig_sum_m_l1m    = eig_sum - l1m
    eig_sum_vn_nx    = eig_sum*vn*nx
    eig_sum_vn2      = eig_sum*vn2

    u_eig_sum_m_l1m  = u*eig_sum_m_l1m
    v_eig_sum_m_l1m  = v*eig_sum_m_l1m
    h0_eig_sum_m_l1m = h0*eig_sum_m_l1m

    u_eig_sum_m_l1m_alpha_ov_c2 = u_eig_sum_m_l1m*alpha_ov_c2
    v_eig_sum_m_l1m_alpha_ov_c2 = v_eig_sum_m_l1m*alpha_ov_c2

    eig_diff_ov_c    = eig_diff*ov_c
    eig_diff_nx_ov_c = eig_diff_ov_c*nx
    eig_diff_vn_ov_c = eig_diff_ov_c*vn

    l1m_ov_gamma_m1  = l1m/(gamma - 1.d0)

    ! Energy related data
    DO i = 1,nb_ns
       tmp1 = Ri(i)*T
       tmp2 = right_data(pos_ei_cell + i - 1)
       epsi(i)   = tmp1 - alpha*tmp2
       sigmai(i) = tmp2 - tmp1*ov_alpha
    ENDDO

    epsi   = epsi + alpha*ek
    sigmai = sigmai + ek

    DO i = 1, nb_int_temp
      eintk(i) = 0.d0
      DO j = 1, nb_ns
        tmp1 = right_data(pos_eki_cell + j - 1 +nb_ns*(i-1))
        epsi(j) = epsi(j) - alpha*tmp1
        sigmai(j) = sigmai(j) + tmp1
        eintk(i) = eintk(i) + tmp1
      ENDDO
    ENDDO

    ! Column j,  j = 1,..,nb_ns 
    DO j = 1,nb_ns

       epsj_ov_c2 = epsi(j)*ov_c2

       DO i = 1,j - 1
          jfr(i,j) = yi(i)*(eig_sum_m_l1m*epsj_ov_c2 - eig_diff_vn_ov_c)
       ENDDO

       jfr(j,j) = l1m + yi(j)*(eig_sum_m_l1m*epsj_ov_c2 - eig_diff_vn_ov_c)

       DO i = j + 1,nb_ns
          jfr(i,j) = yi(i)*(eig_sum_m_l1m*epsj_ov_c2 - eig_diff_vn_ov_c)
       ENDDO

       jfr(pos_rhou,j) = epsj_ov_c2*u_eig_sum_m_l1m  + l1m_u - eig_sum_vn_nx + eig_diff_ov_c*(nx*epsi(j) - u_vn)
       jfr(pos_rhov,j) = epsj_ov_c2*v_eig_sum_m_l1m - eig_diff_ov_c*v_vn
       jfr(pos_rhoE,j) = epsj_ov_c2*h0_eig_sum_m_l1m + l1m*sigmai(j) + epsi(j)*l1m_ov_gamma_m1 &
                              & - eig_sum_vn2 + eig_diff_vn_ov_c*(epsi(j) - h0)

       DO i = 1, nb_int_temp
          jfr(pos_rhoE+i, j) = (epsj_ov_c2*eig_sum_m_l1m - eig_diff_vn_ov_c)*eintk(i)
       ENDDO

    ENDDO

    ! Column nb_ns + 1 
    tmp1 = eig_diff_nx_ov_c - u_eig_sum_m_l1m_alpha_ov_c2
    DO i = 1,nb_ns
       jfr(i,pos_rhou) = tmp1*yi(i)
    ENDDO

    jfr(pos_rhou,pos_rhou) = - u_eig_sum_m_l1m_alpha_ov_c2*u  + eig_sum*nx2 + eig_diff_vn_ov_c*(1.d0 - alpha)
    jfr(pos_rhov,pos_rhou) = - u_eig_sum_m_l1m_alpha_ov_c2*v  + eig_diff_nx_ov_c*v
    jfr(pos_rhoE,pos_rhou) = - u_eig_sum_m_l1m_alpha_ov_c2*h0 - l1m_u + eig_sum*vn*nx + eig_diff_ov_c*(h0*nx - alpha*u_vn)

    DO i = 1, nb_int_temp
       jfr(pos_rhoE+i, pos_rhou) = (eig_diff_ov_c - eig_sum_m_l1m*alpha_ov_c2*u)*eintk(i)
    ENDDO

    ! Column nb_ns + 2
    DO i = 1,nb_ns
       jfr(i,pos_rhov) = 0.d0
    ENDDO

    jfr(pos_rhou,pos_rhov) = 0.d0
    jfr(pos_rhov,pos_rhov) = l1m
    jfr(pos_rhoE,pos_rhov) = 0.d0

    DO i = 1, nb_int_temp
       jfr(pos_rhoE+i, pos_rhov) = 0.d0
    ENDDO

    ! Column nb_ns + 3
    tmp1 = eig_sum_m_l1m*alpha_ov_c2
    DO i = 1,nb_ns
       jfr(i,pos_rhoE) = tmp1*yi(i)
    ENDDO

    jfr(pos_rhou,pos_rhoE) = u_eig_sum_m_l1m_alpha_ov_c2  + alpha*eig_diff_nx_ov_c
    jfr(pos_rhov,pos_rhoE) = v_eig_sum_m_l1m_alpha_ov_c2
    jfr(pos_rhoE,pos_rhoE) = h0_eig_sum_m_l1m*alpha_ov_c2 + l1m + eig_diff_vn_ov_c*alpha

    DO i = 1, nb_int_temp
       jfr(pos_rhoE+i, pos_rhoE) = eintk(i)*alpha_ov_c2*eig_sum_m_l1m
    ENDDO

    ! Column j = nb_ns+4,..,nb_eq
    DO j = 1, nb_int_temp

       DO i = 1, nb_ns
          jfr(i, pos_rhoE+j) = - eig_sum_m_l1m*alpha_ov_c2*yi(i)
       ENDDO

       jfr(pos_rhou, pos_rhoE+j) = - u_eig_sum_m_l1m_alpha_ov_c2 - eig_diff_ov_c*alpha     
       jfr(pos_rhov, pos_rhoE+j) = - v_eig_sum_m_l1m_alpha_ov_c2
       jfr(pos_rhoE, pos_rhoE+j) = - h0_eig_sum_m_l1m*alpha_ov_c2 - eig_diff_vn_ov_c*alpha

       DO i = 1, j-1
          jfr(pos_rhoE+i, pos_rhoE+j) = - eig_sum_m_l1m*alpha_ov_c2*eintk(i)
       ENDDO
       jfr(pos_rhoE+j, pos_rhoE+j) = - eig_sum_m_l1m*alpha_ov_c2*eintk(i) + l1m
       DO i = j+1, nb_int_temp
          jfr(pos_rhoE+i, pos_rhoE+j) = - eig_sum_m_l1m*alpha_ov_c2*eintk(i)
       ENDDO

    ENDDO

  END SUBROUTINE Pos_Neg_split_Jac_neqNT_1D_SL
!------------------------------------------------------------------------------!
