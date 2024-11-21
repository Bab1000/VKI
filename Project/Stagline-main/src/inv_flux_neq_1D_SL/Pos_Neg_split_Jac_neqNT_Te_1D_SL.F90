!------------------------------------------------------------------------------!
!> This subrotine provides the numerical flux Jacobian for 1D stagnation line nonequilibrium gas flows
!! according to the positive and negative split of the inviscid flux Jacobian.
!! In this case the flow is characterized by N temperatures plus a separated temperature Te for free electrons.
  SUBROUTINE Pos_Neg_split_Jac_neqNT_Te_1D_SL (nx, vol_l, vol_r, left_data, right_data, u_left, u_right, jfl, jfr)

    USE mod_general_data, ONLY: nb_ns, pos_ei_cell, pos_h0_cell, pos_rho_cell, pos_T_cell, pos_u_cell,  &
                              & pos_v_cell, yi, Ri, nb_temp, pos_pres_cell, pos_beta_cell, nb_eq

    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: nx                        !< normal to the cell interface
    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfl      !< numerical flux Jacobian with respect to the left state 
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jfr      !< numerical flux Jacobian with respect to the right state
    
    INTEGER :: i, j
    REAL(KIND=8) :: a, a2, ov_a2, u, v, rho, p, H, l1, l2, l3
    REAL(KIND=8), DIMENSION(nb_temp) :: beta, em
    REAL(KIND=8), DIMENSION(nb_ns) :: gam
    REAL(KIND=8), DIMENSION(nb_ns+2+nb_temp, nb_ns+2+nb_temp) :: L, R
    
    ! Left state
    ! A^+ matrix (positive split Jacobian) 
    u   = left_data(pos_u_cell)
    v   = left_data(pos_v_cell)
    rho = left_data(pos_rho_cell)
    p   = left_data(pos_pres_cell)
    H   = left_data(pos_h0_cell)
    yi  = u_left(1:nb_ns)/rho
    
    do i = 1,nb_temp
        em(i) = sum(yi * left_data(pos_ei_cell+(i-1)*nb_ns:pos_ei_cell+i*nb_ns))
    end do
    
    beta(1) = sum(u_left(2:nb_ns)*Ri(2:nb_ns)) / left_data(pos_beta_cell)
    beta(2:nb_temp) = -beta(1)
    beta(nb_temp) = beta(nb_temp) + u_left(1)*Ri(1) / left_data(pos_beta_cell+nb_temp-1)
    
    do i = 1,nb_ns
        gam(i) = beta(1)*(0.5d0*u*u-left_data(pos_ei_cell+i-1))
        do j = 2,nb_temp
            gam(i) = gam(i) - beta(j)*left_data(pos_ei_cell+(j-1)*nb_ns+i)
        end do
        
        if (i == 1) then
            gam(i) = gam(i) + Ri(i)*left_data(pos_T_cell+nb_temp-1)
        else
            gam(i) = gam(i) + Ri(i)*left_data(pos_T_cell)
        end if
    end do
    
    ! Sound speed
    a = (1.0d0+beta(1))*p/rho
    a2 = a*a
    ov_a2 = 1.0d0/a2
    
    ! Form the L matrix
    L = 0.0d0
    
    ! Columns 1..ns
    do j = 1,nb_ns
        L(j,j) = ov_a2
        L(nb_ns+1,j) = u*ov_a2
        L(nb_ns+3,j) = (u*u-gam(j)/beta(1))*ov_a2        
    end do
    
    ! Column ns+1
    L(nb_ns+2,nb_ns+1) = ov_a2
    
    ! Column ns+2
    L(1:nb_ns,nb_ns+2) = yi*ov_a2
    L(nb_ns+1,nb_ns+2) = (u+a)*ov_a2
    L(nb_ns+2,nb_ns+2) = v*ov_a2
    L(nb_ns+3,nb_ns+2) = (H+u*a)*ov_a2
    
    do i = 2,nb_temp
        L(nb_ns+2+i,nb_ns+2) = em(i)*ov_a2
    end do
    
    ! Column ns+3
    L(1:nb_ns,nb_ns+3) = yi*ov_a2
    L(nb_ns+1,nb_ns+3) = (u-a)*ov_a2
    L(nb_ns+2,nb_ns+3) = v*ov_a2
    L(nb_ns+3,nb_ns+3) = (H-u*a)*ov_a2
    
    do i = 2,nb_temp
        L(nb_ns+2+i,nb_ns+3) = em(i)*ov_a2
    end do
    
    ! Columns ns+4..ns+2+nt
    do j = nb_ns+4,nb_ns+2+nb_temp
        L(nb_ns+3,j) = -beta(j-nb_ns-2)/beta(1)*ov_a2
        L(j,j)       = ov_a2
    end do
    
    ! Form the R matrix
    R = 0.0d0
    
    ! Columns 1..ns
    do j = 1,nb_ns
        R(1:nb_ns,j) = -gam(j)*yi
        R(j,j) = R(j,j) + a2
        
        R(nb_ns+1,j) = -gam(j)*v
        R(nb_ns+2,j) = 0.5d0*(gam(j)-a*u)
        R(nb_ns+3,j) = 0.5d0*(gam(j)+a*u)
        
        do i = 2,nb_temp
            R(nb_ns+2+i,j) = -gam(j)*em(i)
        end do
    end do
    
    ! Column ns+1
    R(1:nb_ns,nb_ns+1) = beta(1)*u*yi
    R(nb_ns+1,nb_ns+1) = beta(1)*u*v
    R(nb_ns+2,nb_ns+1) = 0.5d0*(-beta(1)*u+a)
    R(nb_ns+3,nb_ns+1) = 0.5d0*(-beta(1)*u-a)
    
    do i = 2,nb_temp
        R(nb_ns+2+i,nb_ns+1) = beta(1)*u*em(i)
    end do
    
    ! Column ns+2
    R(nb_ns+1,nb_ns+2) = a2
    
    ! Column ns+3
    R(1:nb_ns,nb_ns+3) = -beta(1)*yi
    R(nb_ns+1,nb_ns+3) = -beta(1)*v
    R(nb_ns+2,nb_ns+3) = 0.5d0*beta(1)
    R(nb_ns+3,nb_ns+3) = 0.5d0*beta(1)
    
    do i = 2,nb_temp
        R(nb_ns+2+i,nb_ns+3) = -beta(1)*em(i)
    end do
    
    ! Columns ns+4..ns+2+nt
    do j = nb_ns+4,nb_ns+2+nb_temp
        R(1:nb_ns,j) = -beta(j-nb_ns-2)*yi
        R(nb_ns+1,j) = -beta(j-nb_ns-2)*v
        R(nb_ns+2,j) = 0.5d0*beta(j-nb_ns-2)
        R(nb_ns+3,j) = 0.5d0*beta(j-nb_ns-2)
        
        do i = 2,nb_temp
            R(nb_ns+2+i,j) = -beta(j-nb_ns-2)*em(i)
        end do
        R(j,j) = R(j,j) + a2
    end do
    
    !do i = 1,nb_eq
    !    do j = 1,nb_eq
    !        write(*,'(E15.7)',advance='no') R(i,j)
    !    end do
    !    write(*,*)
    !end do
    !stop
    
    ! Postive-split eigenvalues
    l1 = 0.5d0*(u + sqrt(u*u+0.1*a))
    l2 = 0.5d0*((u+a) + sqrt((u+a)*(u+a)+0.1*a))
    l3 = 0.5d0*((u-a) + sqrt((u-a)*(u-a)+0.1*a))
    
    ! Compute L*Lambda
    L(:,1:nb_ns+1) = l1*L(:,1:nb_ns+1)
    L(:,nb_ns+2)   = l2*L(:,nb_ns+2)
    L(:,nb_ns+3)   = l3*L(:,nb_ns+3)
    L(:,nb_ns+4:nb_ns+2+nb_temp) = l1*L(:,nb_ns+4:nb_ns+2+nb_temp) 
       
    
    ! Finally compute A+ = L*Lambda^+*R
    jfl = matmul(L, R)
    
    ! Right state
    ! A^- matrix (negative split Jacobian) 
    u   = right_data(pos_u_cell)
    v   = right_data(pos_v_cell)
    rho = right_data(pos_rho_cell)
    p   = right_data(pos_pres_cell)
    H   = right_data(pos_h0_cell)
    yi  = u_right(1:nb_ns)/rho
    
    do i = 1,nb_temp
        em(i) = sum(yi * right_data(pos_ei_cell+(i-1)*nb_ns:pos_ei_cell+i*nb_ns))
    end do
    
    beta(1) = sum(u_right(2:nb_ns)*Ri(2:nb_ns)) / right_data(pos_beta_cell)
    beta(2:nb_temp) = -beta(1)
    beta(nb_temp) = beta(nb_temp) + u_right(1)*Ri(1) / right_data(pos_beta_cell+nb_temp-1)
    
    do i = 1,nb_ns
        gam(i) = beta(1)*(0.5d0*u*u-right_data(pos_ei_cell+i-1))
        do j = 2,nb_temp
            gam(i) = gam(i) - beta(j)*right_data(pos_ei_cell+(j-1)*nb_ns+i)
        end do
        
        if (i == 1) then
            gam(i) = gam(i) + Ri(i)*right_data(pos_T_cell+nb_temp-1)
        else
            gam(i) = gam(i) + Ri(i)*right_data(pos_T_cell)
        end if
    end do
    
    ! Sound speed
    a = (1.0d0+beta(1))*p/rho
    a2 = a*a
    ov_a2 = 1.0d0/a2
    
    ! Form the L matrix
    L = 0.0d0
    
    ! Columns 1..ns
    do j = 1,nb_ns
        L(j,j) = ov_a2
        L(nb_ns+1,j) = u*ov_a2
        L(nb_ns+3,j) = (u*u-gam(j)/beta(1))*ov_a2        
    end do
    
    ! Column ns+1
    L(nb_ns+2,nb_ns+1) = ov_a2
    
    ! Column ns+2
    L(1:nb_ns,nb_ns+2) = yi*ov_a2
    L(nb_ns+1,nb_ns+2) = (u+a)*ov_a2
    L(nb_ns+2,nb_ns+2) = v*ov_a2
    L(nb_ns+3,nb_ns+2) = (H+u*a)*ov_a2
    
    do i = 2,nb_temp
        L(nb_ns+2+i,nb_ns+2) = em(i)*ov_a2
    end do
    
    ! Column ns+3
    L(1:nb_ns,nb_ns+3) = yi*ov_a2
    L(nb_ns+1,nb_ns+3) = (u-a)*ov_a2
    L(nb_ns+2,nb_ns+3) = v*ov_a2
    L(nb_ns+3,nb_ns+3) = (H-u*a)*ov_a2
    
    do i = 2,nb_temp
        L(nb_ns+2+i,nb_ns+3) = em(i)*ov_a2
    end do
    
    ! Columns ns+4..ns+2+nt
    do j = nb_ns+4,nb_ns+2+nb_temp
        L(nb_ns+3,j) = -beta(j-nb_ns-2)/beta(1)*ov_a2
        L(j,j)       = ov_a2
    end do
    
    ! Form the R matrix
    R = 0.0d0
    
    ! Columns 1..ns
    do j = 1,nb_ns
        R(1:nb_ns,j) = -gam(j)*yi
        R(j,j) = R(j,j) + a2
        
        R(nb_ns+1,j) = -gam(j)*v
        R(nb_ns+2,j) = 0.5d0*(gam(j)-a*u)
        R(nb_ns+3,j) = 0.5d0*(gam(j)+a*u)
        
        do i = 2,nb_temp
            R(nb_ns+2+i,j) = -gam(j)*em(i)
        end do
    end do
    
    ! Column ns+1
    R(1:nb_ns,nb_ns+1) = beta(1)*u*yi
    R(nb_ns+1,nb_ns+1) = beta(1)*u*v
    R(nb_ns+2,nb_ns+1) = 0.5d0*(-beta(1)*u+a)
    R(nb_ns+3,nb_ns+1) = 0.5d0*(-beta(1)*u-a)
    
    do i = 2,nb_temp
        R(nb_ns+2+i,nb_ns+1) = beta(1)*u*em(i)
    end do
    
    ! Column ns+2
    R(nb_ns+1,nb_ns+2) = a2
    
    ! Column ns+3
    R(1:nb_ns,nb_ns+3) = -beta(1)*yi
    R(nb_ns+1,nb_ns+3) = -beta(1)*v
    R(nb_ns+2,nb_ns+3) = 0.5d0*beta(1)
    R(nb_ns+3,nb_ns+3) = 0.5d0*beta(1)
    
    do i = 2,nb_temp
        R(nb_ns+2+i,nb_ns+3) = -beta(1)*em(i)
    end do
    
    ! Columns ns+4..ns+2+nt
    do j = nb_ns+4,nb_ns+2+nb_temp
        R(1:nb_ns,j) = -beta(j-nb_ns-2)*yi
        R(nb_ns+1,j) = -beta(j-nb_ns-2)*v
        R(nb_ns+2,j) = 0.5d0*beta(j-nb_ns-2)
        R(nb_ns+3,j) = 0.5d0*beta(j-nb_ns-2)
        
        do i = 2,nb_temp
            R(nb_ns+2+i,j) = -beta(j-nb_ns-2)*em(i)
        end do
        R(j,j) = R(j,j) + a2
    end do
    
    ! Negative-split eigenvalues
    l1 = 0.5d0*(u - sqrt(u*u+0.1*a))
    l2 = 0.5d0*((u+a) - sqrt((u+a)*(u+a)+0.1*a))
    l3 = 0.5d0*((u-a) - sqrt((u-a)*(u-a)+0.1*a))
    
    ! Compute L*Lambda
    L(:,1:nb_ns+1) = l1*L(:,1:nb_ns+1)
    L(:,nb_ns+2)   = l2*L(:,nb_ns+2)
    L(:,nb_ns+3)   = l3*L(:,nb_ns+3)
    L(:,nb_ns+4:nb_ns+2+nb_temp) = l1*L(:,nb_ns+4:nb_ns+2+nb_temp) 
    
    ! Finally compute A- = L*Lambda^-*R
    jfr = matmul(L, R)

  END SUBROUTINE Pos_Neg_split_Jac_neqNT_Te_1D_SL
!------------------------------------------------------------------------------!
