!------------------------------------------------------------------------------!
!> This module provides implementation for subroutine used when computing 1D stagnation line nonequilibrium flows.  
  MODULE mod_neq_1D_SL

    USE mod_general_data,     ONLY: nb_ns, nb_dim, nb_temp, nb_int_temp, nb_te, nb_eq, pos_em,  & 
                                  & pos_u, pos_v, pos_T, pos_Tk, pos_Te, pos_rhou, pos_rhov,    & 
                                  & pos_rhoE, pos_rhoek, pos_rhose, gamma_e, gamma_e_m1,        & 
                                  & ov_gamma_e_m1, Ri
  
    IMPLICIT NONE
  
    ! Subroutines and functions for handling data.
    CONTAINS 
  
      !------------------------------------------------------!
      !> This subroutine computes the conservative variable vector from the primitive variable vector 
      !! for 1D nonequilibrium flows.
      SUBROUTINE prim_to_cons_neq_1D_SL (prim, cons)
   
        USE mod_neq_function_pointer,    ONLY: library_get_energy_densities
  
        INTEGER :: i
        REAL(KIND=8) :: tmp
        REAL(KIND=8) :: rho, ek, u, v
   
        REAL(KIND=8), DIMENSION(nb_temp) :: temp, rho_e  
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prim
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: cons
  
        ! Species densities
        rho = 0.d0   
        DO i = 1,nb_ns 
           tmp = prim(i)
           cons(i) = tmp
           rho     = rho + tmp
        ENDDO
  
        ! Radial and circumferential velocity components
        u = prim(pos_u)
        v = prim(pos_v) 

        ! Kinetic energy per unit mass
        ek = 0.5d0*u**2
  
        ! Radial momentum density
        cons(pos_rhou) = rho*u
  
        ! Circumferential momentum density
        cons(pos_rhov) = rho*v

        ! Temperature vector
        DO i = 1,nb_temp
           temp(i) = prim(pos_T + i - 1)
        ENDDO
   
        ! Total energy per unit volume (total and needed components)
        CALL library_get_energy_densities (prim(1:nb_ns), temp, rho_e)
    
        ! Remaining components of conservative variable vector
        DO i = 1,nb_temp
           cons(pos_rhoE + i - 1) = rho_e(i)
        ENDDO
  
        ! Adding the kinetic energy contribution (rhoE)
        cons(pos_rhoE) = cons(pos_rhoE) + rho*ek 
  
      END SUBROUTINE prim_to_cons_neq_1D_SL
  
      !------------------------------------------------------!
      !> This subroutine computes the primitive variable vector from the conservative variable vector 
      !! for 1D nonequilibrium flows.
      SUBROUTINE cons_to_prim_neq_1D_SL (cons, prim)
  
        USE mod_neq_function_pointer,    ONLY: library_get_temperatures
  
        INTEGER :: i
        REAL(KIND=8) :: rho, rho_ek, u, tmp
        REAL(KIND=8), DIMENSION(nb_temp) :: temp, rho_eint
  
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: prim
  
        ! Species densities 
        rho = 0.d0
        DO i = 1,nb_ns 
           tmp = cons(i)
           prim(i) = tmp
           rho     = rho + tmp
        ENDDO
  
        ! Velocity and kinetic energy density
        u      = cons(pos_rhou)/rho
        rho_ek = 0.5d0*rho*u*u
  
        ! Velocity
        prim(pos_u) = u

        ! Energy densities
        rho_eint(1) = cons(pos_rhoE) - rho_ek 
        DO i = 1,nb_int_temp + nb_te
           rho_eint(i + 1) = cons(pos_rhoE + i)
        ENDDO
             
        ! Computing the temperatures
        CALL library_get_temperatures (prim(1:nb_ns), rho_eint, temp)   

        ! Temperatures 
        DO i = 1,nb_temp
           prim(pos_T + i - 1) = temp(i)
        ENDDO
  
      END SUBROUTINE cons_to_prim_neq_1D_SL
  
      !----------------------------------------------------! 
      !> This subroutine computes the inviscid flux for the 1D stagnation line Euler equations
      !! for nonequilibrium flows 
      SUBROUTINE inviscid_flux_neq_1D_SL (nx, prim, f)

        USE mod_neq_function_pointer,    ONLY: library_get_pressure

        INTEGER :: i
        REAL(KIND=8) :: u, vn, p
        REAL(KIND=8), DIMENSION(nb_eq) :: cons

        REAL(KIND=8), INTENT(IN) :: nx
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prim
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: f

        ! Velocity 
        u  = prim(pos_u)
        vn = u*nx

        ! Mixture static pressure 
        CALL library_get_pressure(prim(1:nb_ns), prim(pos_T:nb_eq), p)

        ! Compute the energy conservative variables 
        CALL prim_to_cons_neq_1D_SL (prim, cons) 
        
        DO i = 1,nb_eq 
           f(i) = cons(i)*vn
        ENDDO

        ! Adding the pressure contribution in the global momentum and global energy equation components
        f(pos_rhou) = f(pos_rhou) + p*nx
        f(pos_rhoE) = f(pos_rhoE) + p*vn 

      END SUBROUTINE inviscid_flux_neq_1D_SL

      !----------------------------------------------------!
      !> This subroutine computes the inviscid flux Jacobian for the 1D stagnation line Euler equations
      !! for nonequilibrium flows (1T case)
      SUBROUTINE flux_Jacobian_neq1T_1D_SL (nx, prim, jf) 

        USE mod_neq_function_pointer,    ONLY: library_get_thermodynamic_data, library_get_density, & 
                                             & library_get_mass_fractions

        INTEGER :: i, j
        REAL(KIND=8) :: eps, tmp    
        REAL(KIND=8) :: c, ek, gamma, h0vn, h0, p, rho, T, u, v, uvn, vn, vvn, vcn
        REAL(KIND=8) :: alpha, ek_alpha, vn_alpha, nx_alpha, nx_ovrho, ov_rho, vn_ovrho
        REAL(KIND=8), DIMENSION(nb_ns) :: ei, yi, yivn, epsi
        REAL(KIND=8), DIMENSION(nb_temp) :: temp, beta

        REAL(KIND=8), INTENT(IN) :: nx 
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prim
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jf 

        ! Physical data
        ! Mixture density and species mass fractions
        CALL library_get_density (prim(1:nb_ns), rho) 

        ! Species mass fractions
        CALL library_get_mass_fractions (rho, prim(1:nb_ns), yi)

        ! Velocity components, kinetic energy per unit mass and temperature 
        u  = prim(pos_u)
        v  = prim(pos_v)
        ek = 0.5d0*u**2
        T  = prim(pos_T)

        ! Temperature
        DO i = 1,nb_temp
           temp(i) = prim(pos_T + i - 1)
        ENDDO 

        ! Compute thermodynamic quantities needed for the flux Jacobian
        CALL library_get_thermodynamic_data (rho, prim(1:nb_ns), temp, c, gamma, p, alpha, beta, ei)

        ! Common factors
        vn       = u*nx
        uvn      = u*vn
        vvn      = v*vn
        ov_rho   = 1.d0/rho
        vn_alpha = vn*alpha
        ek_alpha = ek*alpha
        nx_alpha = nx*alpha
        nx_ovrho = nx*ov_rho
        vn_ovrho = vn*ov_rho 

        ! Total enthalpy per unit mass
        h0 = 0.d0
        DO i = 1,nb_ns 
           h0 = h0 + yi(i)*ei(i)   
        ENDDO
        h0   = h0 + ek + p*ov_rho
        h0vn = h0*vn 

        DO i = 1,nb_ns
           yivn(i) = yi(i)*vn 
           epsi(i) = Ri(i)*T - alpha*ei(i) + ek_alpha
        ENDDO
 
        ! Fill matrix 
        ! Columns j,  j =1,..,nb_ns
        DO j = 1,nb_ns 

           DO i = 1,j - 1
              jf(i,j) = - yivn(i)
           ENDDO

           jf(j,j) = vn - yivn(j)

           DO i = j + 1,nb_ns 
              jf(i,j) = - yivn(i)
           ENDDO

           eps = epsi(j)
           jf(pos_rhou,j) = eps*nx - uvn
           jf(pos_rhov,j) = - vvn
           jf(pos_rhoE,j) = eps*vn - h0vn

        ENDDO 

        ! Column nb_ns + 1
        DO i = 1,nb_ns 
           jf(i,pos_rhou) = prim(i)*nx_ovrho
        ENDDO
 
        jf(pos_rhou,pos_rhou) = 2.d0*vn - vn_alpha
        jf(pos_rhov,pos_rhou) = v*nx
        jf(pos_rhoE,pos_rhou) = h0*nx   - alpha*uvn

        ! Column nb_ns + 2
        DO i = 1,nb_ns 
           jf(i,pos_rhov) = 0.d0 
        ENDDO

        jf(pos_rhou,pos_rhov) = 0.d0
        jf(pos_rhov,pos_rhov) = vn
        jf(pos_rhoE,pos_rhov) = 0.d0

        ! Column nb_ns + 3
        DO i = 1,nb_ns 
           jf(i,pos_rhoE) = 0.d0
        ENDDO

        jf(pos_rhou,pos_rhoE) = nx_alpha
        jf(pos_rhov,pos_rhoE) = 0.d0
        jf(pos_rhoE,pos_rhoE) = vn + vn_alpha

      END SUBROUTINE flux_Jacobian_neq1T_1D_SL   

      !----------------------------------------------------!
      !> This subroutine provides the eigensystem for the 1D stagnation line Euler equations for
      !! nonequilibrium flows (1T case)
      SUBROUTINE eigensystem_neq1T_1D_SL (nx, prim, lambda, right, left)    

        USE mod_neq_function_pointer,    ONLY: library_get_thermodynamic_data, library_get_density, & 
                                             & library_get_mass_fractions

        INTEGER :: i, j
        REAL(KIND=8) :: tmp1, tmp2, tmp3
        REAL(KIND=8) :: alpha, eps, gamma
        REAL(KIND=8) :: c, cn, ek, h0, ov_c2, ov_rho, rho, p, T, u, uc, v, vn, vcn, v_ov_rho
        REAL(KIND=8), DIMENSION(nb_ns) :: ei, yi, yivn, epsi, sigmai
        REAL(KIND=8), DIMENSION(nb_temp) :: temp, beta

        REAL(KIND=8), INTENT(IN) :: nx
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prim
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: lambda
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: right, left

        ! Physical data
        ! Mixture density and species mass fractions
        CALL library_get_density (prim(1:nb_ns), rho) 

        ! Species mass fractions
        CALL library_get_mass_fractions (rho, prim(1:nb_ns), yi)

        ! Velocity components, kinetic energy per unit mass and temperature 
        u  = prim(pos_u)
        v  = prim(pos_v)
        ek = 0.5d0*u**2
        T  = prim(pos_T)

        ! Temperature
        DO i = 1,nb_temp
           temp(i) = prim(pos_T + i - 1)
        ENDDO 

        ! Compute thermodynamic quantities needed for the flux Jacobian
        CALL library_get_thermodynamic_data (rho, prim(1:nb_ns), temp, c, gamma, p, alpha, beta, ei)

        ! Total enthalpy per unit mass
        ov_rho = 1.d0/rho
        h0 = 0.d0
        DO i = 1,nb_ns 
           h0 = h0 + yi(i)*ei(i)   
        ENDDO
        h0 = h0 + ek + p*ov_rho

        ! Common factors
        vn  = u*nx
        uc  = u*c
        cn  = c*nx
        vcn = uc*nx
        v_ov_rho = v*ov_rho

        DO i = 1,nb_ns 
          tmp1   = Ri(i)*T
          tmp2   = ei(i)
          epsi(i)   = tmp1 - alpha*tmp2 
          sigmai(i) = tmp2 - tmp1/alpha 
        ENDDO
        
        epsi   = epsi + alpha*ek 
        sigmai = sigmai + ek

        ! Eigenvalues 
        DO i = 1,nb_ns 
           lambda(i) = vn
        ENDDO
        lambda(pos_rhou) = vn - c 
        lambda(pos_rhov) = vn
        lambda(pos_rhoE) = vn + c
   
        ! Right eigenvector matrix 
        ! Column j,  j = 1,..,nb_ns
        DO j = 1,nb_ns 

           DO i = 1,j - 1
              right(i,j) = 0.d0
           ENDDO

           right(j,j) = 1.d0

           DO i = j + 1,nb_ns
              right(i,j) = 0.d0
           ENDDO

           right(pos_rhou,j) = u
           right(pos_rhov,j) = v
           right(pos_rhoE,j) = sigmai(j)

        ENDDO

        ! Column nb_ns + 1
        DO i = 1,nb_ns 
           right(i,pos_rhou) = yi(i)
        ENDDO
 
        right(pos_rhou,pos_rhou) = u - cn
        right(pos_rhov,pos_rhou) = v
        right(pos_rhoE,pos_rhou) = h0 - vcn

        ! Column nb_ns + 2
        DO i = 1,nb_ns 
           right(i,pos_rhov) = 0.d0
        ENDDO

        right(pos_rhou,pos_rhov) = 0.d0
        right(pos_rhov,pos_rhov) = rho
        right(pos_rhoE,pos_rhov) = 0.d0

        ! Column nb_ns + 3
        DO i = 1,nb_ns 
           right(i,pos_rhoE) = yi(i)
        ENDDO
 
        right(pos_rhou,pos_rhoE) = u + cn
        right(pos_rhov,pos_rhoE) = v
        right(pos_rhoE,pos_rhoE) = h0 + vcn

        ! Left eigenvector matrix
        ! Column j,  j = 1,..,nb_ns
        ov_c2 = 1.d0/c**2
        tmp1  = 0.5d0*vcn*ov_c2
        DO j = 1,nb_ns 

           tmp2 = epsi(j)*ov_c2
           tmp3 = 0.5d0*tmp2

           DO i = 1,j - 1
              left(i,j) = - yi(i)*tmp2
           ENDDO

           left(j,j) = 1.d0 - yi(j)*tmp2

           DO i = j + 1,nb_ns
              left(i,j) = - yi(i)*tmp2
           ENDDO

           left(pos_rhou,j) = tmp3 + tmp1
           left(pos_rhov,j) = - v_ov_rho 
           left(pos_rhoE,j) = tmp3 - tmp1 

        ENDDO

        ! Column nb_ns + 1	
        tmp1 = alpha*u*ov_c2
        tmp2 = 0.5d0*tmp1
        tmp3 = 0.5d0*nx/c

        DO i = 1,nb_ns 
           left(i,pos_rhou) = yi(i)*tmp1
        ENDDO
 
        left(pos_rhou,pos_rhou) = - (tmp3 + tmp2)
        left(pos_rhov,pos_rhou) = 0.d0
        left(pos_rhoE,pos_rhou) = tmp3 - tmp2

        ! Column nb_ns + 2
        DO i = 1,nb_ns 
           left(i,pos_rhov) = 0.d0
        ENDDO

        left(pos_rhou,pos_rhov) = 0.d0
        left(pos_rhov,pos_rhov) = ov_rho
        left(pos_rhoE,pos_rhov) = 0.d0

        ! Column nb_ns + 3
        tmp1 = ov_c2*alpha
        tmp2 = 0.5d0*tmp1

        DO i = 1,nb_ns 
           left(i,pos_rhoE) = - yi(i)*tmp1
        ENDDO
 
        left(pos_rhou,pos_rhoE) = tmp2
        left(pos_rhov,pos_rhoE) = 0.d0
        left(pos_rhoE,pos_rhoE) = tmp2

      END SUBROUTINE eigensystem_neq1T_1D_SL

      !----------------------------------------------------!
      !> This subroutine computes the positive and negative splits of the inviscid 
      !! flux Jacobian of the 1D stagnation line Euler equations for nonequilibrium flows (1T case).
      SUBROUTINE flux_Jacobian_split_neq1T_1D_SL (nx, prim, Aplus, Aminus)

        USE mod_neq_function_pointer,    ONLY: library_get_thermodynamic_data, library_get_density, & 
                                             & library_get_mass_fractions

        INTEGER :: i, j
        REAL(KIND=8) :: tmp1, tmp2
        REAL(KIND=8) :: alpha, eps, gamma, nx2
        REAL(KIND=8) :: c, cn, ek, h0, ov_rho, rho, u, v, u_vn, v_vn, vn, vn2, T, p
        REAL(KIND=8) :: alpha_ov_c2, ov_c, ov_c2, yi_epsj_ov_c2, yi_vn_ov_c, epsj_ov_c2, ov_alpha
        REAL(KIND=8) :: l1, l2, l3, l1p, l2p, l3p, l1m, l2m, l3m, l1m_u, l1p_u, l1p_ov_gamma_m1, l1m_ov_gamma_m1
        REAL(KIND=8) :: eig_diff, eig_sum, eig_diff_vn_ov_c, eig_sum_m_l1p, eig_sum_m_l1m,                     & 
                      & u_eig_sum_m_l1p, v_eig_sum_m_l1p, h0_eig_sum_m_l1p, u_eig_sum_m_l1m, v_eig_sum_m_l1m,  & 
                      & h0_eig_sum_m_l1m, eig_sum_vn_nx, eig_diff_ov_c, eig_diff_nx_ov_c, eig_sum_vn2,         & 
                      & u_eig_sum_m_l1p_alpha_ov_c2, u_eig_sum_m_l1m_alpha_ov_c2,                              & 
                      & v_eig_sum_m_l1p_alpha_ov_c2, v_eig_sum_m_l1m_alpha_ov_c2
        REAL(KIND=8), DIMENSION(nb_ns) :: ei, yi, epsi, sigmai
        REAL(KIND=8), DIMENSION(nb_temp) :: temp, beta

        REAL(KIND=8), INTENT(IN) :: nx
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prim
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: Aplus, Aminus

        ! Initialization 
        Aplus  = 0.d0
        Aminus = 0.d0

        ! Physical data
        ! Mixture density and species mass fractions
        CALL library_get_density (prim(1:nb_ns), rho) 

        ! Species mass fractions
        CALL library_get_mass_fractions (rho, prim(1:nb_ns), yi)

        ! Radial and circumferential velocity components 
        u  = prim(pos_u)
        v  = prim(pos_v)

        ! Kinetic energy per unit mass
        ek = 0.5d0*u**2

        ! Temperatute 
        T  = prim(pos_T)
        DO i = 1,nb_temp
           temp(i) = T
        ENDDO 

        ! Compute thermodynamic quantities needed for the flux Jacobian
        CALL library_get_thermodynamic_data (rho, prim(1:nb_ns), temp, c, gamma, p, alpha, beta, ei)

        ! Total enthalpy per unit mass
        ov_rho = 1.d0/rho
        h0 = 0.d0
        DO i = 1,nb_ns 
           h0 = h0 + yi(i)*ei(i)   
        ENDDO
        h0 = h0 + ek + p*ov_rho

        ! Eigenvalues 
        nx2 = nx**2
        vn  = u*nx
        l1  = vn 
        l2  = vn - c
        l3  = vn + c

        ! Positive and negative eigenvalues 
        l1p = 0.5d0*(l1 + ABS(l1))
        l2p = 0.5d0*(l2 + ABS(l2))
        l3p = 0.5d0*(l3 + ABS(l3))

        l1m = 0.5d0*(l1 - ABS(l1))
        l2m = 0.5d0*(l2 - ABS(l2))
        l3m = 0.5d0*(l3 - ABS(l3))

        ! A^+ matrix (positive split Jacobian)
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
          tmp1   = Ri(i)*T
          tmp2   = ei(i)
          epsi(i)   = tmp1 - alpha*tmp2 
          sigmai(i) = tmp2 - tmp1*ov_alpha 
        ENDDO
        
        epsi   = epsi + alpha*ek 
        sigmai = sigmai + ek

        ! Column i,  i = 1,..,nb_ns 
        DO j = 1,nb_ns 

           epsj_ov_c2 = epsi(j)*ov_c2

           DO i = 1,j - 1
              Aplus(i,j) = yi(i)*(eig_sum_m_l1p*epsj_ov_c2 - eig_diff_vn_ov_c)  
           ENDDO

           Aplus(j,j) = l1p + yi(j)*(eig_sum_m_l1p*epsj_ov_c2 - eig_diff_vn_ov_c) 

           DO i = j + 1,nb_ns 
              Aplus(i,j) = yi(i)*(eig_sum_m_l1p*epsj_ov_c2 - eig_diff_vn_ov_c)
           ENDDO

           Aplus(pos_rhou,j) = epsj_ov_c2*u_eig_sum_m_l1p  + l1p_u - eig_sum_vn_nx + eig_diff_ov_c*(nx*epsi(j) - u_vn)
           Aplus(pos_rhov,j) = epsj_ov_c2*v_eig_sum_m_l1p - eig_diff_ov_c*v_vn
           Aplus(pos_rhoE,j) = epsj_ov_c2*h0_eig_sum_m_l1p + l1p*sigmai(j) + epsi(j)*l1p_ov_gamma_m1 & 
                             & - eig_sum_vn2 + eig_diff_vn_ov_c*(epsi(j) - h0)

        ENDDO

        ! Column nb_ns + 1 
        tmp1 = eig_diff_nx_ov_c  - u_eig_sum_m_l1p_alpha_ov_c2
        DO i = 1,nb_ns 
           Aplus(i,pos_rhou) = tmp1*yi(i)
        ENDDO
  
        Aplus(pos_rhou,pos_rhou) = - u_eig_sum_m_l1p_alpha_ov_c2*u  + eig_sum*nx2 + eig_diff_vn_ov_c*(1.d0 - alpha)
        Aplus(pos_rhov,pos_rhou) = - u_eig_sum_m_l1p_alpha_ov_c2*v  + eig_diff_nx_ov_c*v
        Aplus(pos_rhoE,pos_rhou) = - u_eig_sum_m_l1p_alpha_ov_c2*h0 - l1p_u + eig_sum*vn*nx + eig_diff_ov_c*(h0*nx - alpha*u_vn)

        ! Column nb_bs + 2
        DO i = 1,nb_ns 
           Aplus(i,pos_rhov) = 0.d0
        ENDDO

        Aplus(pos_rhou,pos_rhov) = 0.d0
        Aplus(pos_rhov,pos_rhov) = l1p
        Aplus(pos_rhoE,pos_rhov) = 0.d0

        ! Column nb_ns + 3
        tmp1 = eig_sum_m_l1p*alpha_ov_c2
        DO i = 1,nb_ns 
           Aplus(i,pos_rhoE) = tmp1*yi(i)
        ENDDO

        Aplus(pos_rhou,pos_rhoE) = u_eig_sum_m_l1p_alpha_ov_c2  + alpha*eig_diff_nx_ov_c 
        Aplus(pos_rhov,pos_rhoE) = v_eig_sum_m_l1p_alpha_ov_c2
        Aplus(pos_rhoE,pos_rhoE) = h0_eig_sum_m_l1p*alpha_ov_c2 + l1p + eig_diff_vn_ov_c*alpha

        ! A^- matrix (negative split Jacobian)
        ! Common factors
        eig_sum  = 0.5d0*(l3m + l2m)
        eig_diff = 0.5d0*(l3m - l2m) 

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

        ! Column j,  j = 1,..,nb_ns 
        DO j = 1,nb_ns 

           epsj_ov_c2 = epsi(j)*ov_c2

           DO i = 1,j - 1
              Aminus(i,j) = yi(i)*(eig_sum_m_l1m*epsj_ov_c2 - eig_diff_vn_ov_c)  
           ENDDO

           Aminus(j,j) = l1m + yi(j)*(eig_sum_m_l1m*epsj_ov_c2 - eig_diff_vn_ov_c) 

           DO i = j + 1,nb_ns 
              Aminus(i,j) = yi(i)*(eig_sum_m_l1m*epsj_ov_c2 - eig_diff_vn_ov_c)
           ENDDO

           Aminus(pos_rhou,j) = epsj_ov_c2*u_eig_sum_m_l1m  + l1m_u - eig_sum_vn_nx + eig_diff_ov_c*(nx*epsi(j) - u_vn)
           Aminus(pos_rhov,j) = epsj_ov_c2*v_eig_sum_m_l1m - eig_diff_ov_c*v_vn 
           Aminus(pos_rhoE,j) = epsj_ov_c2*h0_eig_sum_m_l1m + l1m*sigmai(j) + epsi(j)*l1m_ov_gamma_m1 & 
                              & - eig_sum_vn2 + eig_diff_vn_ov_c*(epsi(j) - h0)

        ENDDO

        ! Column nb_ns + 1 
        tmp1 = eig_diff_nx_ov_c - u_eig_sum_m_l1m_alpha_ov_c2
        DO i = 1,nb_ns 
           Aminus(i,pos_rhou) = tmp1*yi(i) 
        ENDDO
  
        Aminus(pos_rhou,pos_rhou) = - u_eig_sum_m_l1m_alpha_ov_c2*u  + eig_sum*nx2 + eig_diff_vn_ov_c*(1.d0 - alpha)
        Aminus(pos_rhov,pos_rhou) = - u_eig_sum_m_l1m_alpha_ov_c2*v  + eig_diff_nx_ov_c*v
        Aminus(pos_rhoE,pos_rhou) = - u_eig_sum_m_l1m_alpha_ov_c2*h0 - l1m_u + eig_sum*vn*nx + eig_diff_ov_c*(h0*nx - alpha*u_vn)

        ! Column nb_ns + 2
        DO i = 1,nb_ns 
           Aminus(i,pos_rhov) = 0.d0
        ENDDO

        Aminus(pos_rhou,pos_rhov) = 0.d0
        Aminus(pos_rhov,pos_rhov) = l1m
        Aminus(pos_rhoE,pos_rhov) = 0.d0

        ! Column nb_ns + 3
        tmp1 = eig_sum_m_l1m*alpha_ov_c2
        DO i = 1,nb_ns 
           Aminus(i,pos_rhoE) = tmp1*yi(i)
        ENDDO

        Aminus(pos_rhou,pos_rhoE) = u_eig_sum_m_l1m_alpha_ov_c2  + alpha*eig_diff_nx_ov_c
        Aminus(pos_rhov,pos_rhoE) = v_eig_sum_m_l1m_alpha_ov_c2
        Aminus(pos_rhoE,pos_rhoE) = h0_eig_sum_m_l1m*alpha_ov_c2 + l1m + eig_diff_vn_ov_c*alpha 

      END SUBROUTINE flux_Jacobian_split_neq1T_1D_SL

  END MODULE mod_neq_1D_SL
!------------------------------------------------------------------------------!
