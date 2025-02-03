!------------------------------------------------------------------------------!
!> This subroutine provides the eigensystem for the 1D Euler equations for nonequilibrium flows:
!! a) Eigenvalues 
!! b) Right eigenvector matrix
!! c) Left eigenvector matrix 
!! In this case the flow is characterized by N temperatures plus a separated temperature Te for free electrons.
  SUBROUTINE eigensystem_neqNT_Te_1D (nx, cons, phys_data, eig, right_eig, left_eig) 

    USE mod_general_data,         ONLY: nb_ns, nb_int_temp, nb_te, nb_temp, pos_u_cell, pos_c_cell, pos_ek_cell,    & 
                                      & pos_alpha_cell, pos_gamma_cell, pos_betak_cell, pos_h0_cell, pos_ei_cell,   & 
                                      & pos_rho_cell, pos_T_cell, pos_em, pos_rhou, pos_rhoE, pos_rhoek, pos_rhose, & 
                                      & gamma_e, gamma_e_m1, ov_gamma_e_m1, Ri, epsi, sigmai, ei, yi, eintk, betak, & 
                                      & ov_betak
 
    IMPLICIT NONE

    INTEGER :: i, j, k, kp, pos_k
    REAL(KIND=8) :: tmp1, tmp2, tmp3
    REAL(KIND=8) :: alpha, eps, gamma, ge_m_g, ge_m1_pe_ov_rho, ge_m_g_pe_ov_rho, rho_pow_ge,  & 
                  & rho_pow_ge_m1, pe_gamma_e_ov_rhoc2, pe_ge_m1_ov_rho_pow_ge,                & 
                  & alpha_pe_gamma_e_ov_rhoc2
    REAL(KIND=8) :: c, cn, ek, eit, eistar, fe, fe_ov_c2, h0, ov_c2, ov_rho, ov_alpha,  & 
                  & rho, p, pe, se, T, Te, u, uc, vn, vcn, RiT

    REAL(KIND=8), INTENT(IN) :: nx                                       !< normal
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons                       !< conservative variables 
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data                  !< physical properties
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: eig                       !< eigenvalue vector
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: right_eig, left_eig     !< right and left eigenvector matrices

    ! Physical data
    T     = phys_data(pos_T_cell)
    Te    = phys_data(pos_T_cell + nb_temp - 1)
    u     = phys_data(pos_u_cell)
    h0    = phys_data(pos_h0_cell)
    c     = phys_data(pos_c_cell)
    rho   = phys_data(pos_rho_cell)
    ek    = phys_data(pos_ek_cell)
    alpha = phys_data(pos_alpha_cell)
    gamma = phys_data(pos_gamma_cell)

    DO i = 1,nb_ns + nb_ns*nb_int_temp
       ei(i) = phys_data(pos_ei_cell + i - 1) 
    ENDDO

    ov_rho   = 1.d0/rho
    ov_alpha = 1.d0/alpha
    DO k = 1,nb_int_temp 
       eintk(k) = cons(pos_rhoek + k - 1)*ov_rho
       tmp1     = phys_data(pos_betak_cell + k - 1) 
       betak(k)    = tmp1 
       ov_betak(k) = 1.d0/tmp1
    ENDDO

    ! Free electron pressure and specific pseudo-entropy
    se = cons(pos_rhose)*ov_rho
    pe = cons(pos_em)*Ri(pos_em)*Te

    ! Common factors
    vn    = u*nx
    uc    = u*c
    cn    = c*nx
    vcn   = uc*nx
    ov_c2 = 1.d0/c**2

    ov_alpha   = 1.d0/alpha
    rho_pow_ge = rho**gamma_e
    rho_pow_ge_m1 = rho**gamma_e_m1

    ge_m_g = gamma_e - gamma
    ge_m1_pe_ov_rho  = gamma_e_m1*pe*ov_rho
    ge_m_g_pe_ov_rho = ge_m_g*pe*ov_rho

    pe_ge_m1_ov_rho_pow_ge    = pe*gamma_e_m1/rho_pow_ge
    pe_gamma_e_ov_rhoc2       = pe*gamma_e*ov_rho*ov_c2
    alpha_pe_gamma_e_ov_rhoc2 = alpha*pe_gamma_e_ov_rhoc2

    fe = ge_m_g*ov_gamma_e_m1*(rho**gamma_e_m1)   
    fe_ov_c2 = fe*ov_c2   

    ! Free electrons
    yi(pos_em)     = cons(pos_em)*ov_rho
    epsi(pos_em)   = 0.d0  
    sigmai(pos_em) = 0.d0 

    ! Heavy-particles 
    DO i = pos_em + 1,nb_ns 
       yi(i) = cons(i)*ov_rho
       RiT   = Ri(i)*T
       eit   = ei(i)   
       eistar = 0.d0
       DO k = 1,nb_int_temp
          eistar = eistar + ei(nb_ns + (i - 1)*nb_int_temp + k)
       ENDDO   
       epsi(i)   = RiT - alpha*(eit - eistar) 
       sigmai(i) = eit - RiT*ov_alpha
    ENDDO
        
    epsi   = epsi + alpha*ek + ge_m_g_pe_ov_rho
    sigmai = sigmai + ek

    ! Eigenvalues 
    DO i = 1,nb_ns 
       eig(i) = vn
    ENDDO

    eig(pos_rhou) = vn - c 
    eig(pos_rhoE) = vn + c

    DO k = 1,nb_int_temp
       eig(pos_rhoek + k - 1) = vn
    ENDDO

    eig(pos_rhose) = vn

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
       right_eig(pos_rhoE,j) = sigmai(j)

       DO k = 1,nb_int_temp
          right_eig(pos_rhoek + k - 1,j) = ei(nb_ns + nb_int_temp*(j - 1) + k)
       ENDDO

       right_eig(pos_rhose,j) = - pe_ge_m1_ov_rho_pow_ge

    ENDDO

    ! Column nb_ns + 1
    DO i = 1,nb_ns 
       right_eig(i,pos_rhou) = yi(i)
    ENDDO
 
    right_eig(pos_rhou,pos_rhou) = u - cn
    right_eig(pos_rhoE,pos_rhou) = h0 - vcn

    DO k = 1,nb_int_temp
       right_eig(pos_rhoek + k - 1,pos_rhou) = eintk(k)
    ENDDO

    right_eig(pos_rhose,pos_rhou) = se

    ! Column nb_ns + 2
    DO i = 1,nb_ns 
       right_eig(i,pos_rhoE) = yi(i)
    ENDDO
 
    right_eig(pos_rhou,pos_rhoE) = u + cn
    right_eig(pos_rhoE,pos_rhoE) = h0 + vcn

    DO k = 1,nb_int_temp
       right_eig(pos_rhoek + k - 1,pos_rhoE) = eintk(k)
    ENDDO

    right_eig(pos_rhose,pos_rhoE) = se

    ! Column nb_ns + 2 + k,  k = 1,..,nb_int_temp 
    DO k = 1,nb_int_temp

       pos_k = pos_rhoek + k - 1
       tmp1  = betak(k)

       DO i = 1,nb_ns 
          right_eig(i,pos_k) = 0.d0
       ENDDO

       right_eig(pos_rhou,pos_k) = 0.d0 
       right_eig(pos_rhoE,pos_k) = tmp1

       DO kp = 1,k - 1
          right_eig(pos_rhoek + kp - 1,pos_k) = 0.d0
       ENDDO

       right_eig(pos_k,pos_k) = tmp1

       DO kp = k + 1,nb_int_temp
          right_eig(pos_rhoek + kp - 1,pos_k) = 0.d0
       ENDDO

       right_eig(pos_rhose,pos_k) = 0.d0

    ENDDO

    ! Column nb_ns + 2 + nb_int_temp + nb_te
    DO i = 1,nb_ns 
       right_eig(i,pos_rhose) = 0.d0
    ENDDO

    right_eig(pos_rhou,pos_rhose) = 0.d0
    right_eig(pos_rhoE,pos_rhose) = - ge_m_g*ov_gamma_e_m1/(gamma - 1.d0)

    DO k = 1,nb_int_temp
       right_eig(pos_rhoek + k - 1,pos_rhose) = 0.d0
    ENDDO

    right_eig(pos_rhose,pos_rhose) = 1.d0/rho_pow_ge_m1

    ! Left eigenvector matrix 
    ! Column i,  i = 1,..,nb_ns
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
       left_eig(pos_rhoE,j) = tmp3 - tmp1 

       DO k = 1,nb_int_temp
          left_eig(pos_rhoek + k - 1,j) = - ei(nb_ns + nb_int_temp*(j - 1) + k)*ov_betak(k)
       ENDDO

       left_eig(pos_rhose,j) = ge_m1_pe_ov_rho - epsi(j)*pe_gamma_e_ov_rhoc2

    ENDDO

    ! Column nb_ns + 1	
    tmp1 = alpha*u*ov_c2
    tmp2 = 0.5d0*tmp1
    tmp3 = 0.5d0*nx/c

    DO i = 1,nb_ns 
       left_eig(i,pos_rhou) = yi(i)*tmp1
    ENDDO
 
    left_eig(pos_rhou,pos_rhou) = - (tmp3 + tmp2)
    left_eig(pos_rhoE,pos_rhou) = tmp3 - tmp2

    DO k = 1,nb_int_temp
       left_eig(pos_rhoek + k - 1,pos_rhou) = 0.d0
    ENDDO

    left_eig(pos_rhose,pos_rhou) = alpha_pe_gamma_e_ov_rhoc2*u

    ! Column nb_ns + 2
    tmp1 = ov_c2*alpha
    tmp2 = 0.5d0*tmp1

    DO i = 1,nb_ns 
       left_eig(i,pos_rhoE) = - yi(i)*tmp1
    ENDDO
 
    left_eig(pos_rhou,pos_rhoE) = tmp2
    left_eig(pos_rhoE,pos_rhoE) = tmp2

    DO k = 1,nb_int_temp
       left_eig(pos_rhoek + k - 1,pos_rhoE) = 0.d0
    ENDDO 

    left_eig(pos_rhose,pos_rhoE) = - alpha_pe_gamma_e_ov_rhoc2 

    ! Column nb_ns + 2 + k,  k = 1,..,nb_int_temp 
    DO k = 1,nb_int_temp

       pos_k = pos_rhoek + k - 1

       DO i = 1,nb_ns 
          left_eig(i,pos_k) = yi(i)*tmp1
       ENDDO

       left_eig(pos_rhou,pos_k) = - tmp2
       left_eig(pos_rhoE,pos_k) = - tmp2

       DO kp = 1,k - 1
          left_eig(pos_rhoek + kp - 1,pos_k) = 0.d0  
       ENDDO

       left_eig(pos_k,pos_k) = ov_betak(k)
 
       DO kp = k + 1,nb_int_temp
          left_eig(pos_rhoek + kp - 1,pos_k) = 0.d0
       ENDDO

       left_eig(pos_rhose,pos_k) = alpha_pe_gamma_e_ov_rhoc2

    ENDDO

    ! Column nb_ns + 2 + nb_int_temp + nb_te 
    DO i = 1,nb_ns
       left_eig(i,pos_rhose) = - yi(i)*fe_ov_c2
    ENDDO

    tmp1 = 0.5d0*fe_ov_c2
    left_eig(pos_rhou,pos_rhose) = tmp1
    left_eig(pos_rhoE,pos_rhose) = tmp1

    DO k = 1,nb_int_temp
      left_eig(pos_rhoek + k - 1,pos_rhose) = 0.d0
    ENDDO

    left_eig(pos_rhose,pos_rhose) = rho_pow_ge_m1*(1.d0 - pe_gamma_e_ov_rhoc2*ge_m_g*ov_gamma_e_m1)

  END SUBROUTINE eigensystem_neqNT_Te_1D
!------------------------------------------------------------------------------!
