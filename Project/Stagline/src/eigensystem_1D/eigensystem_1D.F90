!------------------------------------------------------------------------------!
!> This subroutine provides the eigensystem for the 1D Euler equation:
!! a) Eigenvalues 
!! b) Right eigenvator matrix 
!! c) Left eigenvector matrix
  SUBROUTINE eigensystem_1D (nx, phys_data, eig, right_eig, left_eig) 

    USE mod_general_data,   ONLY: pos_u_cell, pos_c_cell, pos_h0_cell, pos_ek_cell, gamma

    IMPLICIT NONE

    REAL(KIND=8) :: tmp1, tmp2, tmp3
    REAL(KIND=8) :: gamma_minus1, gamma_minus1_ov_c2
    REAL(KIND=8) :: c, cn, cvn, ek, h0, ov_c2, u, uc, vn, v2

    REAL(KIND=8), INTENT(IN) :: nx                                       !< normal 
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data                  !< physical properties
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: eig                       !< eigenvalue vector
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: right_eig, left_eig     !< right and left eigenvector matrices
   
    ! Physical data 
    u  = phys_data(pos_u_cell)
    h0 = phys_data(pos_h0_cell)
    c  = phys_data(pos_c_cell)
    ek = phys_data(pos_ek_cell)

    ! Common factors  
    vn  = u*nx
    uc  = u*c
    cn  = c*nx
    cvn = vn*c
    ov_c2 = 1.d0/c**2

    gamma_minus1 = gamma - 1.d0
    gamma_minus1_ov_c2 = gamma_minus1*ov_c2

    ! Eigenvalues 
    eig(1) = vn 
    eig(2) = vn - c
    eig(3) = vn + c 

    ! Right eigenvector matrix 
    ! First column 
    right_eig(1,1) = 1.d0
    right_eig(2,1) = u
    right_eig(3,1) = ek

    ! Second column 
    right_eig(1,2) = 1.d0
    right_eig(2,2) = u - cn
    right_eig(3,2) = h0 - cvn

    ! Third column 
    right_eig(1,3) = 1.d0
    right_eig(2,3) = u + cn
    right_eig(3,3) = h0 + cvn

    ! Left eigenvector matrix
    tmp1 = gamma_minus1_ov_c2*ek
    tmp2 = 0.5d0*tmp1
    tmp3 = 0.5d0*cvn*ov_c2

    ! First column
    left_eig(1,1) = 1.d0 - tmp1
    left_eig(2,1) = tmp2 + tmp3
    left_eig(3,1) = tmp2 - tmp3

    ! Second column
    tmp1 = gamma_minus1_ov_c2*u
    tmp2 = 0.5d0*tmp1
    tmp3 = 0.5d0*cn*ov_c2 

    left_eig(1,2) = tmp1
    left_eig(2,2) = - (tmp2 + tmp3)
    left_eig(3,2) = tmp3 - tmp2

    ! Third column 
    tmp1 = gamma_minus1_ov_c2
    tmp2 = 0.5d0*tmp1

    left_eig(1,3) = - tmp1
    left_eig(2,3) = tmp2
    left_eig(3,3) = tmp2 

  END SUBROUTINE eigensystem_1D
!------------------------------------------------------------------------------!
