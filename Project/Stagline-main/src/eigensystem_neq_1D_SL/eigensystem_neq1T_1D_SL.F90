!------------------------------------------------------------------------------!
!> This subroutine provides the eigensystem for the 1D stagnation line Euler equations for nonequilibrium flows:
!! a) Eigenvalues 
!! b) Right eigenvector matrix
!! c) Left eigenvector matrix 
!! In this case the flow is characterized by one single temperature T.
  SUBROUTINE eigensystem_neq1T_1D_SL (nx, cons, phys_data, eig, right_eig, left_eig) 

    USE mod_general_data,         ONLY: nb_ns, pos_u_cell, pos_v_cell, pos_c_cell, pos_ek_cell,  & 
                                      & pos_alpha_cell, pos_h0_cell, pos_ei_cell, pos_rho_cell,  & 
                                      & pos_T_cell, pos_rhou, pos_rhov, pos_rhoE, Ri, epsi,      &
                                      & sigmai, ei, yi

    IMPLICIT NONE

    INTEGER :: i, j
    REAL(KIND=8) :: tmp1, tmp2, tmp3
    REAL(KIND=8) :: alpha, eps
    REAL(KIND=8) :: c, cn, ek, h0, ov_c2, ov_rho, v_ov_rho, rho, p, T, u, v, uc, vn, vvn, vcn

    REAL(KIND=8), INTENT(IN) :: nx                                       !< normal
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons                       !< conservative variables 
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data                  !< physical properties
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: eig                       !< eigenvalue vector
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: right_eig, left_eig     !< right and left eigenvector matrices

    ! Physical data
    T     = phys_data(pos_T_cell)
    u     = phys_data(pos_u_cell)
    v     = phys_data(pos_v_cell)
    h0    = phys_data(pos_h0_cell)
    c     = phys_data(pos_c_cell)
    rho   = phys_data(pos_rho_cell)
    ek    = phys_data(pos_ek_cell)
    alpha = phys_data(pos_alpha_cell)

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
       tmp1  = Ri(i)*T
       tmp2  = phys_data(pos_ei_cell + i - 1)
       epsi(i)   = tmp1 - alpha*tmp2 
       sigmai(i) = tmp2 - tmp1/alpha 
    ENDDO
        
    epsi   = epsi + alpha*ek 
    sigmai = sigmai + ek

    ! Eigenvalues 
    DO i = 1,nb_ns 
       eig(i) = vn
    ENDDO
    eig(pos_rhou) = vn - c 
    eig(pos_rhov) = vn 
    eig(pos_rhoE) = vn + c
   
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

    ENDDO

    ! Column nb_ns + 1
    DO i = 1,nb_ns 
       right_eig(i,pos_rhou) = yi(i)
    ENDDO
 
    right_eig(pos_rhou,pos_rhou) = u - cn
    right_eig(pos_rhov,pos_rhou) = v
    right_eig(pos_rhoE,pos_rhou) = h0 - vcn

    ! Column nb_ns + 2
    DO i = 1,nb_ns 
       right_eig(i,pos_rhov) = 0.d0
    ENDDO

    right_eig(pos_rhou,pos_rhov) = 0.d0
    right_eig(pos_rhov,pos_rhov) = rho
    right_eig(pos_rhoE,pos_rhov) = 0.d0 

    ! Column nb_ns + 3
    DO i = 1,nb_ns 
       right_eig(i,pos_rhoE) = yi(i)
    ENDDO
 
    right_eig(pos_rhou,pos_rhoE) = u + cn
    right_eig(pos_rhov,pos_rhoE) = v
    right_eig(pos_rhoE,pos_rhoE) = h0 + vcn

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

    ! Column nb_ns + 2
    DO i = 1,nb_ns 
       left_eig(i,pos_rhov) = 0.d0
    ENDDO

    left_eig(pos_rhou,pos_rhov) = 0.d0
    left_eig(pos_rhov,pos_rhov) = ov_rho
    left_eig(pos_rhoE,pos_rhov) = 0.d0 

    ! Column nb_ns + 3
    tmp1 = ov_c2*alpha
    tmp2 = 0.5d0*tmp1

    DO i = 1,nb_ns 
       left_eig(i,pos_rhoE) = - yi(i)*tmp1
    ENDDO
 
    left_eig(pos_rhou,pos_rhoE) = tmp2
    left_eig(pos_rhov,pos_rhoE) = 0.d0
    left_eig(pos_rhoE,pos_rhoE) = tmp2 

  END SUBROUTINE eigensystem_neq1T_1D_SL
!------------------------------------------------------------------------------!
