!------------------------------------------------------------------------------!
!> This module provides subroutines for variable transformation and eigensystem for 1D calorically perfect gas flows.
  MODULE mod_pg_1D

    IMPLICIT NONE
  
    REAL(KIND=8), SAVE :: gamma, gamma_minus1, gamma_plus1, gamma_minus3, R_gas 
   
    ! Subroutine and functions for handling data.
    CONTAINS 
 
      !------------------------------------------------------! 
      !> This subroutine initializes the specific heat ratios and the specific gas constants.
      !! It has to be called once before the use of the other subroutines provided in the present module
      SUBROUTINE initialize_pg_1D(in_gamma, in_R_gas) 

        REAL(KIND=8), INTENT(IN) :: in_gamma, in_R_gas

        ! Specific heat ratio and gas constant
        gamma = in_gamma
        gamma_plus1  = gamma + 1.d0
        gamma_minus1 = gamma - 1.d0
        gamma_minus3 = gamma - 3.d0
         
        R_gas = in_R_gas 

      END SUBROUTINE initialize_pg_1D

      !------------------------------------------------------!
      !> This subroutine computes the post-shock conditions for a normal shock
      !! wave according to the Rankine-Hugonito jump relations
      SUBROUTINE rh_jump (p1, u1, T1, p2, u2, T2)

        REAL(KIND=8) :: c1, M1, M1s

        REAL(KIND=8), INTENT(IN) :: p1, u1, T1
        REAL(KIND=8), INTENT(OUT) :: p2, u2, T2

        c1  = SQRT(gamma*R_gas*T1)
        M1  = u1/c1
        M1s = M1**2  

        p2 = p1*(2.d0*gamma*M1s - gamma + 1.d0)/gamma_plus1
        u2 = u1 - c1*2.d0/gamma_plus1*(M1 - 1.d0/M1)
        T2 = T1*(2.d0*gamma*M1s - gamma + 1.d0)*(gamma - 1.d0 + 2.d0/M1s)/(gamma_plus1**2)     

      END SUBROUTINE rh_jump 
 
      !------------------------------------------------------!
      !> This subroutine computes the set of conservative variable vector from the set of primitive variable vector
      !! for 1D calorically perfect gas flows.
      SUBROUTINE prim_to_cons_1D (prim, cons)
  
        REAL(KIND=8) :: u, v, rho, p, v2 
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prim
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cons
 
        rho = prim(1)
        u   = prim(2)
        p   = prim(3)
        v2  = u*u
    
        cons(1) = rho
        cons(2) = rho*u
        cons(3) = p/gamma_minus1 + 0.5d0*rho*v2
      
      END SUBROUTINE prim_to_cons_1D
  
      !------------------------------------------------------!
      !> This subroutione computes the set of conservative variable vector from the set of primitive variable vector
      !! for 1D calorically perfect gas flows.
      SUBROUTINE cons_to_prim_1D (cons, prim)
  
        REAL(KIND=8) :: rho, rhou, rhov, rhoE, rhov2
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: prim

        rho   = cons(1)
        rhou  = cons(2)
        rhoE  = cons(3)
        rhov2 = rhou*rhou
      
        prim(1) = rho
        prim(2) = rhou/rho
        prim(3) = gamma_minus1*(rhoE - 0.5d0*rhov2/rho)
    
      END SUBROUTINE cons_to_prim_1D
  
      !------------------------------------------------------!
      !> This subroutine computes the normal inviscid flux for the 1D Euler equations
      !! along the nx direction.
      SUBROUTINE inv_flux_1D (nx, prim, f)
  
        REAL(KIND=8) :: rho, u, v, vn, p, v2, h0, tmp 
  
        REAL(KIND=8), INTENT(IN) :: nx
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prim
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: f 
 
        rho = prim(1)
        u   = prim(2)
        p   = prim(3)
  
        vn = u*nx
        v2 = u*u
        h0 = p*gamma/(rho*gamma_minus1) + 0.5d0*v2
  
        tmp  = rho*vn
        f(1) = tmp
        f(2) = tmp*u + p*nx
        f(3) = tmp*h0
  
      END SUBROUTINE inv_flux_1D 

      !----------------------------------------------------!
      !> This subroutine compitues the inviscid flux Jacobians for the 1D Euler equations.
      SUBROUTINE flux_Jacobian_1D (nx, prim, j)

        REAL(KIND=8) :: rho, u, vn, v2, h0, p, tmp
 
        REAL(KIND=8), INTENT(IN) :: nx
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prim
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: j

        rho = prim(1)
        u   = prim(2)
        p   = prim(3)
  
        vn  = u*nx
        v2  = u*u
        tmp = gamma_minus1*v2
        h0  = p*gamma/(rho*gamma_minus1) + 0.5d0*v2

        ! First column
        j(1,1) = 0.d0
        j(2,1) = 0.5d0*gamma_minus3*v2*nx
        j(3,1) = vn*(0.5d0*tmp - h0)

        ! Second column
        j(1,2) = nx
        j(2,2) = - gamma_minus3*vn
        j(3,2) = (h0 - tmp)*nx
        
        ! Third column
        j(1,3) = 0.d0
        j(2,3) = gamma_minus1*nx
        j(3,3) = gamma*vn

      END SUBROUTINE flux_Jacobian_1D

      !----------------------------------------------------!
      !> This subroutines provides the eigensystem of the 1D Euler equations along the nx direction. 
      SUBROUTINE eigensystem_1D (nx, prim, lambda, right, left) 

        REAL(KIND=8) :: rho, u, vn, v2, h0, p, ek, c, cn
        REAL(KIND=8) :: fac1, fac2, fac3, fac_ov_c, p_ov_rho, u_ov_c, uc

        REAL(KIND=8), INTENT(IN) :: nx
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prim
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: lambda
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: right, left 

        rho = prim(1)
        u   = prim(2)
        p   = prim(3)
        p_ov_rho = p/rho  

        vn = u*nx
        v2 = u*u
        ek = 0.5d0*v2
        c  = DSQRT(gamma*p_ov_rho)
        cn = c*nx
        uc = u*cn
        h0 = p_ov_rho*gamma/gamma_minus1 + ek

        ! Eigenvalues 
        lambda(1) = vn
        lambda(2) = vn - c
        lambda(3) = vn + c
       
        ! First right eigenvector 
        right(1,1)  = 1.d0
        right(2,1)  = u
        right(3,1)  = ek

        ! Second right eigenvector
        right(1,2)  = 1.d0
        right(2,2)  = u - cn
        right(3,2)  = h0 - uc

        ! Third right eigenvector 
        right(1,3)  = 1.d0
        right(2,3)  = u + cn
        right(3,3)  = h0 + uc 

        ! Common factors for left eigenvector matrix
        u_ov_c   = u/c
        fac_ov_c = 0.5d0*u_ov_c
        fac1     = gamma_minus1*u_ov_c
        fac2     = 0.5d0/c
        fac3     = gamma_minus1/(c*c)

        ! First left eigenvector
        left(1,1)  = 1.d0 - 0.5d0*fac3*v2
        left(1,2)  = u*fac3
        left(1,3)  = - fac3

        ! Second left eigenvector
        left(2,1)  = fac_ov_c*(nx + 0.5d0*fac1) 
        left(2,2)  = - fac2*(nx + fac1)
        left(2,3)  = 0.5d0*fac3
         
        ! Third left eigenvector
        left(3,1)  = fac_ov_c*(- nx + 0.5d0*fac1) 
        left(3,2)  = fac2*(nx - fac1)
        left(3,3)  = 0.5d0*fac3

      END SUBROUTINE eigensystem_1D 

      !----------------------------------------------------!
      !> This subroutine computes the positive and negative splits of the inviscid 
      !! flux Jacobian along the nx direction for the 1D Euler equations.
      SUBROUTINE flux_Jacobian_split_1D (nx, prim, Aplus, Aminus)

        REAL(KIND=8) :: rho, p, u
        REAL(KIND=8) :: c, c2, ek, h0, nx_ov_c, ov_c, ov_c2, v2, vn
        REAL(KIND=8) :: l1, l2, l3, l1p, l2p, l3p, l1m, l2m, l3m
        REAL(KIND=8) :: eig_diff, eig_sum
        REAL(KIND=8) :: gamma_minus1, gamma_minus1_ov_c, gamma_minus1_ov_c2 

        REAL(KIND=8), INTENT(IN) :: nx
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prim
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: Aplus, Aminus

        rho = prim(1)
        u   = prim(2)
        p   = prim(3)
        c   = DSQRT(gamma*p/rho)
        c2  = c**2
        v2  = u**2
        ek  = 0.5d0*v2
        vn  = u*nx
        h0  = c2/gamma_minus1 + ek
       
        ov_c  = 1.d0/c
        ov_c2 = ov_c/c
        nx_ov_c = nx*ov_c
        gamma_minus1_ov_c  = gamma_minus1*ov_c
        gamma_minus1_ov_c2 = gamma_minus1_ov_c*ov_c

        ! Eigenvalues 
        l1 = vn 
        l2 = vn - c
        l3 = vn + c

        ! Positive and negative eigenvalues 
        l1p = 0.5d0*(l1 + ABS(l1))
        l2p = 0.5d0*(l2 + ABS(l2))
        l3p = 0.5d0*(l3 + ABS(l3))

        l1m = 0.5d0*(l1 - ABS(l1))
        l2m = 0.5d0*(l2 - ABS(l2))
        l3m = 0.5d0*(l3 - ABS(l3))

        ! A^+ matrix (positive split Jacobian)
        ! Useful factors
        eig_sum  = 0.5d0*(l3p + l2p)
        eig_diff = 0.5d0*(l3p - l2p)      

        ! First column 
        Aplus(1,1) = l1p + (eig_sum - l1p)*gamma_minus1_ov_c2*ek - eig_diff*vn*ov_c
        Aplus(2,1) = l1p*u + (eig_sum - l1p)*gamma_minus1_ov_c2*ek*u - eig_sum*nx*vn + & 
                   & eig_diff*(gamma_minus1_ov_c*ek*nx - u*vn*ov_c)
        Aplus(3,1) = l1p*ek + (eig_sum*h0 - l1p*ek)*0.5d0*gamma_minus1_ov_c2*v2 &
                   & - eig_sum*vn**2 + eig_diff*(0.5d0*gamma_minus1_ov_c*v2*vn - h0*vn*ov_c) 

        ! Second column 
        Aplus(1,2) = (l1p - eig_sum)*gamma_minus1_ov_c2*u + eig_diff*nx_ov_c
        Aplus(2,2) = l1p*gamma_minus1_ov_c2*v2 + eig_sum*(nx**2 - gamma_minus1_ov_c2*v2) + & 
                   & eig_diff*(2.d0 - gamma)*u*nx_ov_c
        Aplus(3,2) = l1p*gamma_minus1_ov_c2*ek*u + eig_sum*(vn*nx - gamma_minus1_ov_c2*u*h0) + &
                   & eig_diff*(h0*nx*ov_c - gamma_minus1_ov_c*u*vn)

        ! Third column 
        Aplus(1,3) = (eig_sum - l1p)*gamma_minus1_ov_c2
        Aplus(2,3) = (eig_sum - l1p)*gamma_minus1_ov_c2*u + eig_diff*gamma_minus1*nx_ov_c
        Aplus(3,3) = - l1p*gamma_minus1_ov_c2*ek + eig_sum*gamma_minus1_ov_c2*h0 + & 
                   &  eig_diff*gamma_minus1_ov_c*vn

        ! A^+-matrix (negative split Jacobian)
        ! Useful factors 
        eig_sum  = 0.5d0*(l3m + l2m)
        eig_diff = 0.5d0*(l3m - l2m)  

        ! First column 
        Aminus(1,1) = l1m + (eig_sum - l1m)*gamma_minus1_ov_c2*ek - eig_diff*vn*ov_c
        Aminus(2,1) = l1m*u + (eig_sum - l1m)*gamma_minus1_ov_c2*ek*u - eig_sum*nx*vn + & 
                   & eig_diff*(gamma_minus1_ov_c*ek*nx - u*vn*ov_c)
        Aminus(3,1) = l1m*ek + (eig_sum*h0 - l1m*ek)*0.5d0*gamma_minus1_ov_c2*v2 &
                   & - eig_sum*vn**2 + eig_diff*(0.5d0*gamma_minus1_ov_c*v2*vn - h0*vn*ov_c)

        ! Second column 
        Aminus(1,2) = (l1m - eig_sum)*gamma_minus1_ov_c2*u + eig_diff*nx_ov_c
        Aminus(2,2) = l1m*gamma_minus1_ov_c2*v2 + eig_sum*(nx**2 - gamma_minus1_ov_c2*v2) + & 
                   & eig_diff*(2.d0 - gamma)*u*nx_ov_c
        Aminus(3,2) = l1m*gamma_minus1_ov_c2*ek*u + eig_sum*(vn*nx - gamma_minus1_ov_c2*u*h0) + &
                   & eig_diff*(h0*nx*ov_c - gamma_minus1_ov_c*u*vn)

        ! Third column
        Aminus(1,3) = (eig_sum - l1m)*gamma_minus1_ov_c2
        Aminus(2,3) = (eig_sum - l1m)*gamma_minus1_ov_c2*u + eig_diff*gamma_minus1*nx_ov_c
        Aminus(3,3) = - l1m*gamma_minus1_ov_c2*ek + eig_sum*gamma_minus1_ov_c2*h0 + & 
                   &  eig_diff*gamma_minus1_ov_c*vn

      END SUBROUTINE flux_Jacobian_split_1D

  END MODULE mod_pg_1D
!------------------------------------------------------------------------------!
