!------------------------------------------------------------------------------!
!> This subroutine provides the eigensystem for the 1D Euler equations for nonequilibrium flows:
!! a) Eigenvalues 
!! b) Right eigenvector matrix
!! c) Left eigenvector matrix 
!! In this case the flow is characterized by N temperatures (no separated electron entropy equation used).
  SUBROUTINE eigensystem_neqNT_1D (nx, cons, phys_data, eig, right_eig, left_eig) 

    USE mod_general_data,         ONLY: nb_ns, nb_int_temp, pos_u_cell, pos_c_cell, pos_ek_cell,    & 
                                      & pos_alpha_cell, pos_betak_cell, pos_h0_cell, pos_ei_cell,   & 
                                      & pos_rho_cell, pos_T_cell, pos_rhou, pos_rhoE, pos_rhoek,    & 
                                      & Ri, epsi, sigmai, ei, yi, eintk, betak, ov_betak

    IMPLICIT NONE

    INTEGER :: i, j, k, kp
    INTEGER :: pos_k
    REAL(KIND=8) :: tmp1, tmp2, tmp3
    REAL(KIND=8) :: alpha, eps, ov_alpha
    REAL(KIND=8) :: c, cn, ek, h0, ov_c2, ov_rho, rho, T, u, uc, vn, vcn

    REAL(KIND=8), INTENT(IN) :: nx                                       !< normal
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons                       !< conservative variables 
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data                  !< physical properties
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: eig                       !< eigenvalue vector
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: right_eig, left_eig     !< right and left eigenvector matrices

    ! Physical data
    T     = phys_data(pos_T_cell)
    u     = phys_data(pos_u_cell)
    h0    = phys_data(pos_h0_cell)
    c     = phys_data(pos_c_cell)
    rho   = phys_data(pos_rho_cell)
    ek    = phys_data(pos_ek_cell)
    alpha = phys_data(pos_alpha_cell)

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

    ! Common factors
    vn  = u*nx
    uc  = u*c
    cn  = c*nx
    vcn = uc*nx

    DO i = 1,nb_ns 
       yi(i) = cons(i)*ov_rho 
       tmp1  = Ri(i)*T
       tmp2  = ei(i)
       tmp3  = 0.d0
       pos_k = nb_ns + nb_int_temp*(i - 1)
       DO k = 1,nb_int_temp
          tmp3 = tmp3 + ei(pos_k + k)
       ENDDO
       epsi(i)   = tmp1 - alpha*(tmp2 - tmp3)
       sigmai(i) = tmp2 - tmp1*ov_alpha 
    ENDDO

    epsi   = epsi + alpha*ek 
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

       pos_k = nb_ns + nb_int_temp*(j - 1)
       DO k = 1,nb_int_temp
          right_eig(pos_rhoek + k - 1,j) = ei(pos_k + k)
       ENDDO

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

    ! Column nb_ns + 2
    DO i = 1,nb_ns 
       right_eig(i,pos_rhoE) = yi(i)
    ENDDO
 
    right_eig(pos_rhou,pos_rhoE) = u + cn
    right_eig(pos_rhoE,pos_rhoE) = h0 + vcn

    DO k = 1,nb_int_temp
       right_eig(pos_rhoek + k - 1,pos_rhoE) = eintk(k)
    ENDDO

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

    ENDDO

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
       left_eig(pos_rhoE,j) = tmp3 - tmp1 

       pos_k = nb_ns + nb_int_temp*(j - 1)
       DO k = 1,nb_int_temp
          left_eig(pos_rhoek + k - 1,j) = - ei(pos_k + k)*ov_betak(k)
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
    left_eig(pos_rhoE,pos_rhou) = tmp3 - tmp2

    DO k = 1,nb_int_temp
       left_eig(pos_rhoek + k - 1,pos_rhou) = 0.d0
    ENDDO

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

    ENDDO

  END SUBROUTINE eigensystem_neqNT_1D
!------------------------------------------------------------------------------!
