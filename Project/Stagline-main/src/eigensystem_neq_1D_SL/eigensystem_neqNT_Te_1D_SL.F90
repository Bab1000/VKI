!------------------------------------------------------------------------------!
!> This subroutine provides the eigensystem for the 1D stagnation line Euler equations for nonequilibrium flows:
!! a) Eigenvalues 
!! b) Right eigenvector matrix
!! c) Left eigenvector matrix 
!! In this case the flow is characterized by N temperatures plus a separated temperature Te for free electrons. 
  SUBROUTINE eigensystem_neqNT_Te_1D_SL (nx, cons, phys_data, eig, right_eig, left_eig) 

   USE mod_general_data,         ONLY: nb_ns, pos_u_cell, pos_v_cell, pos_c_cell, pos_ek_cell,  &
                                     & pos_alpha_cell, pos_h0_cell, pos_ei_cell, pos_rho_cell,  &
                                     & pos_T_cell, pos_rhou, pos_rhov, pos_rhoE, Ri, epsi,      &
                                     & sigmai, ei, yi, nb_eq, nb_int_temp, pos_rhoek,           &
                                     & pos_eki_cell, eintk, pos_beta_cell, pos_em, nb_temp, nb_te

    IMPLICIT NONE

    INTEGER :: i, j
    REAL(KIND=8) :: tmp1, tmp2, tmp3
    REAL(KIND=8) :: alpha, alpha_e, eps, beta_e
    REAL(KIND=8) :: c, cn, ek, h0, ov_c2, ov_rho, v_ov_rho, rho, p, T, Te, u, v, uc, vn, vvn, vcn


    REAL(KIND=8), INTENT(IN) :: nx                                       !< normal
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons                       !< conservative variables 
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data                  !< physical properties
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: eig                       !< eigenvalue vector
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: right_eig, left_eig     !< right and left eigenvector matrices

    eig = 0.d0
    right_eig = 0.d0
    left_eig  = 0.d0

    ! Physical data
    T     = phys_data(pos_T_cell)
    Te    = phys_data(pos_T_cell+nb_temp-1)
    u     = phys_data(pos_u_cell)
    v     = phys_data(pos_v_cell)
    h0    = phys_data(pos_h0_cell)
    c     = phys_data(pos_c_cell)
    rho   = phys_data(pos_rho_cell)
    ek    = phys_data(pos_ek_cell)
    alpha = phys_data(pos_alpha_cell)
    beta_e = phys_data(pos_beta_cell+nb_temp-1)
    alpha_e = cons(pos_em)*Ri(pos_em)/beta_e

    ! Common factors
    vn  = u*nx
    vvn = v*vn
    uc  = u*c
    cn  = c*nx
    vcn = uc*nx

    ov_rho = 1.d0/rho
    v_ov_rho = v*ov_rho
    DO i = 1,nb_ns
       yi(i) = cons(i)*ov_rho
       IF(i.EQ.pos_em) THEN
         tmp1  = Ri(i)*Te
       ELSE
         tmp1  = Ri(i)*T
       ENDIF
       tmp2  = phys_data(pos_ei_cell + i - 1)
       epsi(i)   = tmp1 - alpha*tmp2
       sigmai(i) = tmp2 - tmp1/alpha
    ENDDO

    epsi   = epsi + alpha*ek
    sigmai = sigmai + ek

    DO i = 1, nb_int_temp
      eintk(i) = 0.d0
      DO j = 1, nb_ns
        IF (j.NE.pos_em) THEN
          tmp1 = phys_data(pos_eki_cell + j - 1 +nb_ns*(i-1))
          epsi(j) = epsi(j) + alpha*tmp1
          sigmai(j) = sigmai(j) - tmp1
          eintk(i) = eintk(i) + tmp1*yi(j)
        ENDIF
      ENDDO
    ENDDO

    eintk(nb_temp-1) = 0.d0
    DO j = 1, nb_ns
      tmp1 = phys_data(pos_eki_cell + j - 1 +nb_ns*(nb_temp-2))
      epsi(j) = epsi(j) + (alpha-alpha_e)*tmp1
      sigmai(j) = sigmai(j) + (alpha_e/alpha-1)* tmp1
      eintk(nb_temp-1) = eintk(nb_temp-1) + tmp1*yi(j)
    ENDDO

    ! Eigenvalues 
    DO i = 1,nb_ns
       eig(i) = vn
    ENDDO
    eig(pos_rhou) = vn - c
    eig(pos_rhov) = vn
    eig(pos_rhoE) = vn + c
    DO i=1,nb_int_temp+nb_te
      eig(pos_rhoek+i-1) = vn
    ENDDO

    ! Right eigenvector matrix 
    ! Column j,  j = 1,..,nb_ns
    DO j = 1,nb_ns

       DO i = 1,j - 1
          right_eig(i,j) = 0.d0
       ENDDO

       right_eig(j,j) = 1.d0

       DO i = j + 1,nb_ns
          right_eig(i,j) = 0.d0
       ENDDO

       right_eig(pos_rhou,j) = u
       right_eig(pos_rhov,j) = v
       right_eig(pos_rhoE,j) = sigmai(j)

       DO i = 1, nb_int_temp + nb_te
         right_eig(pos_rhoE+i, j) = 0.d0
       ENDDO
    ENDDO

    ! Column nb_ns + 1
    DO i = 1,nb_ns
       right_eig(i,pos_rhou) = yi(i)
    ENDDO

    right_eig(pos_rhou,pos_rhou) = u - cn
    right_eig(pos_rhov,pos_rhou) = v
    right_eig(pos_rhoE,pos_rhou) = h0 - vcn

    DO i = 1, nb_int_temp + nb_te
      right_eig(pos_rhoE+i,pos_rhou) = eintk(i)
    ENDDO

    ! Column nb_ns + 2
    DO i = 1,nb_ns
       right_eig(i,pos_rhov) = 0.d0
    ENDDO

    right_eig(pos_rhou,pos_rhov) = 0.d0
    right_eig(pos_rhov,pos_rhov) = rho
    right_eig(pos_rhoE,pos_rhov) = 0.d0

    DO i = 1, nb_int_temp + nb_te
      right_eig(pos_rhoE+i, pos_rhov) = 0.d0
    ENDDO

    ! Column nb_ns + 3
    DO i = 1,nb_ns
       right_eig(i,pos_rhoE) = yi(i)
    ENDDO

    right_eig(pos_rhou,pos_rhoE) = u + cn
    right_eig(pos_rhov,pos_rhoE) = v
    right_eig(pos_rhoE,pos_rhoE) = h0 + vcn

    DO i = 1, nb_int_temp + nb_te
      right_eig(pos_rhoE+i, pos_rhoE) = eintk(i)
    ENDDO

    ! Column j, j = nb_ns+4,..., nb_eq-1
    DO j = 1, nb_int_temp
       DO i = 1, nb_ns
         right_eig(i, pos_rhoE+j) = 0.d0
       ENDDO

       right_eig(pos_rhou, pos_rhoE+j) = 0.d0
       right_eig(pos_rhov, pos_rhoE+j) = 0.d0
       right_eig(pos_rhoE, pos_rhoE+j) = 1.d0

       DO i = 1,j - 1
          right_eig(pos_rhoE+i, pos_rhoE+j) = 0.d0
       ENDDO
       right_eig(pos_rhoE+j, pos_rhoE+j) = 1.d0
       DO i = j + 1,nb_int_temp + nb_te
          right_eig(pos_rhoE+i, pos_rhoE+j) = 0.d0
       ENDDO
    ENDDO

    ! Column nb_eq
    DO i = 1, nb_ns
      right_eig(i, pos_rhoE+nb_temp-1) = 0.d0
    ENDDO

    right_eig(pos_rhou, pos_rhoE+nb_temp-1) = 0.d0
    right_eig(pos_rhov, pos_rhoE+nb_temp-1) = 0.d0
    right_eig(pos_rhoE, pos_rhoE+nb_temp-1) = 1.d0 - alpha_e/alpha

    DO i = 1, nb_int_temp
      right_eig(pos_rhoE+i, pos_rhoE+nb_temp-1) = 0.d0
    ENDDO

    right_eig(pos_rhoE+nb_temp-1, pos_rhoE+nb_temp-1) = 1.d0

    ! Left eigenvector matrix
    ! Column j,  j = 1,..,nb_ns
    ov_c2 = 1.d0/c**2
    tmp1  = 0.5d0*vcn*ov_c2
    DO j = 1,nb_ns

       tmp2 = epsi(j)*ov_c2
       tmp3 = 0.5d0*tmp2

       DO i = 1,j - 1
          left_eig(i,j) = - yi(i)*tmp2
       ENDDO

       left_eig(j,j) = 1.d0 - yi(j)*tmp2

       DO i = j + 1,nb_ns
          left_eig(i,j) = - yi(i)*tmp2
       ENDDO

       left_eig(pos_rhou,j) = tmp3 + tmp1
       left_eig(pos_rhov,j) = - v_ov_rho
       left_eig(pos_rhoE,j) = tmp3 - tmp1

       DO i = 1, nb_int_temp + nb_te
         left_eig(pos_rhoE+i, j) = -eintk(i)*tmp2
       ENDDO
    ENDDO

    ! Column nb_ns + 1  
    tmp1 = alpha*u*ov_c2
    tmp2 = 0.5d0*tmp1
    tmp3 = 0.5d0*nx/c

    DO i = 1,nb_ns
       left_eig(i,pos_rhou) = yi(i)*tmp1
    ENDDO

    left_eig(pos_rhou,pos_rhou) = - (tmp3 + tmp2)
    left_eig(pos_rhov,pos_rhou) = 0.d0
    left_eig(pos_rhoE,pos_rhou) = tmp3 - tmp2

    DO i = 1, nb_int_temp + nb_te
      left_eig(pos_rhoE+i, pos_rhou) = eintk(i)*tmp1
    ENDDO

    ! Column nb_ns + 2
    DO i = 1,nb_ns
       left_eig(i,pos_rhov) = 0.d0
    ENDDO

    left_eig(pos_rhou,pos_rhov) = 0.d0
    left_eig(pos_rhov,pos_rhov) = ov_rho
    left_eig(pos_rhoE,pos_rhov) = 0.d0
    DO i = 1, nb_int_temp + nb_te
      left_eig(pos_rhoE+i,pos_rhov) = 0.d0
    ENDDO

    ! Column nb_ns + 3
    tmp1 = ov_c2*alpha
    tmp2 = 0.5d0*tmp1

    DO i = 1,nb_ns
       left_eig(i,pos_rhoE) = - yi(i)*tmp1
    ENDDO

    left_eig(pos_rhou,pos_rhoE) = tmp2
    left_eig(pos_rhov,pos_rhoE) = 0.d0
    left_eig(pos_rhoE,pos_rhoE) = tmp2

    DO i = 1, nb_int_temp + nb_te
      left_eig(pos_rhoE+i, pos_rhoE) = - eintk(i)*tmp1
    ENDDO

    ! Column j, j = nb_ns+4,..., nb_eq-1
    DO j = 1, nb_int_temp
      DO i = 1, nb_ns
        left_eig(i, pos_rhoE+j) = yi(i)*tmp1
      ENDDO

      left_eig(pos_rhou, pos_rhoE+j) = - tmp2
      left_eig(pos_rhov, pos_rhoE+j) = 0.d0
      left_eig(pos_rhoE, pos_rhoE+j) = - tmp2

      DO i = 1,j - 1
         left_eig(pos_rhoE+i, pos_rhoE+j) = eintk(i)*tmp1
      ENDDO
      left_eig(pos_rhoE+j, pos_rhoE+j) = 1.d0 + eintk(j)*tmp1
      DO i = j + 1, nb_int_temp + nb_te
         left_eig(pos_rhoE+i, pos_rhoE+j) = eintk(i)*tmp1
      ENDDO
    ENDDO

    ! Column nb_eq
    tmp1 = ov_c2*(alpha_e-alpha)
    tmp2 = 0.5d0*tmp1

    DO i = 1, nb_ns
      left_eig(i, pos_rhoE+nb_temp-1) = -yi(i)*tmp1
    ENDDO

    left_eig(pos_rhou, pos_rhoE+nb_temp-1) = tmp2
    left_eig(pos_rhov, pos_rhoE+nb_temp-1) = 0.d0
    left_eig(pos_rhoE, pos_rhoE+nb_temp-1) = tmp2

    DO i = 1, nb_int_temp 
       left_eig(pos_rhoE+i, pos_rhoE+nb_temp-1) = -eintk(i)*tmp1
    ENDDO

    left_eig(pos_rhoE+nb_temp-1, pos_rhoE+nb_temp-1) = 1.d0 - eintk(nb_temp-1)*tmp1

  END SUBROUTINE eigensystem_neqNT_Te_1D_SL
!------------------------------------------------------------------------------!
