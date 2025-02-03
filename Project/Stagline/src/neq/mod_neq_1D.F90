!------------------------------------------------------------------------------!
!> This module provides implementation for subroutine used when computing 1D nonequilibrium flows.  
  MODULE mod_neq_1D

    USE mod_general_data,     ONLY: nb_ns, nb_dim, nb_temp, nb_int_temp, nb_te, nb_eq, pos_em,   & 
                                  & pos_u, pos_T, pos_Tk, pos_Te, pos_rhou, pos_rhoE, pos_rhoek, & 
                                  & pos_rhose, gamma_e, gamma_e_m1, ov_gamma_e_m1, Ri
  
    IMPLICIT NONE
  
    ! Subroutines and functions for handling data.
    CONTAINS 
  
      !------------------------------------------------------!
      !> This subroutine computes the conservative variable vector from the primitive variable vector 
      !! for 1D nonequilibrium flows.
      SUBROUTINE prim_to_cons_neq_1D (prim, cons)
   
        USE mod_neq_function_pointer,    ONLY: library_get_energy_densities
  
        INTEGER :: i
        REAL(KIND=8) :: rho, ek, u, tmp
   
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
  
        ! Velocity and kinetic energy per unit mass
        u  = prim(pos_u)
        ek = 0.5d0*u**2
  
        ! Momentum density
        cons(pos_rhou) = rho*u
  
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
  
      END SUBROUTINE prim_to_cons_neq_1D
  
      !------------------------------------------------------!
      !> This subroutine computes the primitive variable vector from the conservative variable vector 
      !! for 1D nonequilibrium flows.
      SUBROUTINE cons_to_prim_neq_1D (cons, prim)
  
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
  
      END SUBROUTINE cons_to_prim_neq_1D
  
      !----------------------------------------------------! 
      !> This subroutine computes the inviscid flux for the 1D Euler equations
      !! for nonequilibrium flows 
      SUBROUTINE inviscid_flux_neq_1D (nx, prim, f)

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
        CALL library_get_pressure(prim(1:nb_ns), prim(nb_ns + 2:nb_eq), p)

        ! Compute the energy conservative variables 
        CALL prim_to_cons_neq_1D (prim, cons) 
        
        DO i = 1,nb_eq 
           f(i) = cons(i)*vn
        ENDDO

        ! Adding the pressure contribution in the global momentum and global energy equation components
        f(pos_rhou) = f(pos_rhou) + p*nx
        f(pos_rhoE) = f(pos_rhoE) + p*vn 

      END SUBROUTINE inviscid_flux_neq_1D

      !----------------------------------------------------!
      !> This subroutine computes the inviscid flux Jacobian for the 1D Euler equations
      !! for nonequilibrium flows (1T case)
      SUBROUTINE flux_Jacobian_neq1T_1D (nx, prim, jf) 

        USE mod_neq_function_pointer,    ONLY: library_get_thermodynamic_data, library_get_density, & 
                                             & library_get_mass_fractions

        INTEGER :: i, j
        REAL(KIND=8) :: eps, tmp    
        REAL(KIND=8) :: c, ek, gamma, h0vn, h0, p, rho, T, u, uvn, vn, vcn
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

        ! Velocity, kinetic energy per unit mass and temperature 
        u  = prim(pos_u)
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
        ! Columns 1,nb_ns
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
           jf(pos_rhoE,j) = eps*vn - h0vn

        ENDDO 

        ! Column nb_ns + 1
        DO i = 1,nb_ns 
           jf(i,pos_rhou) = prim(i)*nx_ovrho
        ENDDO
 
        jf(pos_rhou,pos_rhou) = 2.d0*vn - vn_alpha
        jf(pos_rhoE,pos_rhou) = h0*nx   - alpha*uvn

        ! Column nb_ns + 2
        DO i = 1,nb_ns 
           jf(i,pos_rhoE) = 0.d0
        ENDDO

        jf(pos_rhou,pos_rhoE) = nx_alpha
        jf(pos_rhoE,pos_rhoE) = vn + vn_alpha

      END SUBROUTINE flux_Jacobian_neq1T_1D   

      !----------------------------------------------------!
      !> This subroutine computes the inviscid flux Jacobian for the 1D Euler equations
      !! for nonequilibrium flows (NT case)
      SUBROUTINE flux_Jacobian_neqNT_1D (nx, prim, jf) 

        USE mod_neq_function_pointer,    ONLY: library_get_thermodynamic_data, library_get_density, & 
                                             & library_get_mass_fractions, library_get_energy_densities 

        INTEGER :: i, j, k
        REAL(KIND=8) :: eps, tmp1, tmp2    
        REAL(KIND=8) :: c, ek, gamma, h0vn, h0, p, rho, T, u, uvn, vn, vcn
        REAL(KIND=8) :: alpha, ek_alpha, vn_alpha, nx_alpha, nx_ovrho, ov_rho, vn_ovrho
        REAL(KIND=8), DIMENSION(nb_ns) :: yi, yivn, epsi
        REAL(KIND=8), DIMENSION(nb_ns + nb_ns*nb_int_temp) :: ei
        REAL(KIND=8), DIMENSION(nb_temp) :: temp, beta, rho_eint

        REAL(KIND=8), INTENT(IN) :: nx
        REAL(KIND=8), INTENT(IN), DIMENSION(:) :: prim
        REAL(KIND=8), INTENT(OUT), DIMENSION(:,:) :: jf

        ! Physical data
        ! Mixture density and species mass fractions
        CALL library_get_density (prim(1:nb_ns), rho) 

        ! Species mass fractions
        CALL library_get_mass_fractions (rho, prim(1:nb_ns), yi)

        ! Velocity, kinetic energy per unit mass and temperature 
        u  = prim(pos_u)
        ek = 0.5d0*u**2
        T  = prim(pos_T)

        ! Temperature
        DO i = 1,nb_temp
           temp(i) = prim(pos_T + i - 1)
        ENDDO 

        ! Compute thermodynamic quantities needed for the flux Jacobian
        CALL library_get_thermodynamic_data (rho, prim(1:nb_ns), temp, c, gamma, p, alpha, beta, ei)

        ! Compute energy densities
        CALL library_get_energy_densities (prim(1:nb_ns), temp, rho_eint)

        ! Common factors
        vn       = u*nx
        uvn      = u*vn
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
           tmp1    = Ri(i)*T - alpha*ei(i) + ek_alpha
           tmp2    = 0.d0
           DO k = 1,nb_int_temp
              tmp2 = tmp2 - ei(nb_ns + (i - 1)*nb_int_temp + k)
           ENDDO   
           epsi(i) = tmp1 - alpha*tmp2
        ENDDO

        ! Column i, i = 1,..,nb_ns
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
           jf(pos_rhoE,j) = eps*vn - h0vn

           DO i = 1,nb_int_temp
              jf(pos_rhoE + i,j) = - rho_eint(i + 1)*vn_ovrho
           ENDDO

        ENDDO 

        ! Column nb_ns + 1
        DO i = 1,nb_ns 
           jf(i,pos_rhou) = prim(i)*nx_ovrho
        ENDDO

        jf(pos_rhou,pos_rhou) = 2.d0*vn - vn_alpha
        jf(pos_rhoE,pos_rhou) = h0*nx   - alpha*uvn

        DO i = 1,nb_int_temp
           jf(pos_rhoE + i,pos_rhou) = rho_eint(i + 1)*nx_ovrho
        ENDDO

        ! Column nb_ns + 2
        DO i = 1,nb_ns 
           jf(i,pos_rhoE) = 0.d0
        ENDDO

        jf(pos_rhou,pos_rhoE) = nx_alpha
        jf(pos_rhoE,pos_rhoE) = vn + vn_alpha 

        DO i = 1,nb_int_temp
           jf(pos_rhoE + i,pos_rhoE) = 0.d0
        ENDDO
 
        ! Column nb_ns 2 + k,  k = 1,..,nb_int_temp
        DO j = pos_rhoE + 1,pos_rhoE + nb_int_temp

           DO i = 1,nb_ns 
              jf(i,j) = 0.d0
           ENDDO

           jf(pos_rhou,j) = - nx_alpha
           jf(pos_rhoE,j) = - vn_alpha

           DO i = pos_rhoE + 1,j - 1
              jf(i,j) = 0.d0
           ENDDO

           jf(j,j) = vn

           DO i = j + 1,pos_rhoE + nb_int_temp
             jf(i,j) = 0.d0
           ENDDO

        ENDDO

      END SUBROUTINE flux_Jacobian_neqNT_1D

      !----------------------------------------------------!
      !> This subroutine computes the inviscid flux Jacobian for the 1D Euler equations
      !! for nonequilibrium flows (NT - Te case)
      SUBROUTINE flux_Jacobian_neqNT_Te_1D (nx, prim, jf) 

        USE mod_neq_function_pointer,    ONLY: library_get_thermodynamic_data, library_get_density, & 
                                             & library_get_mass_fractions, library_get_energy_densities 

        INTEGER :: i, j, k
        REAL(KIND=8) :: eps, tmp1, tmp2    
        REAL(KIND=8) :: c, ek, fe, gamma, h0vn, h0, p, pe, rho, T, Te, u, uvn, vn, vcn
        REAL(KIND=8) :: alpha, ek_alpha, ge_m_g, ge_m_g_pe_ov_rho, vn_alpha, nx_alpha, nx_ovrho, & 
                      & ov_rho, vn_ovrho
        REAL(KIND=8), DIMENSION(nb_ns) :: yi, yivn, epsi
        REAL(KIND=8), DIMENSION(nb_ns + nb_ns*nb_int_temp) :: ei
        REAL(KIND=8), DIMENSION(nb_temp) :: temp, beta, rho_eint
       
        REAL(KIND=8), INTENT(IN) :: nx
        REAL(KIND=8), INTENT(IN), DIMENSION(:) :: prim
        REAL(KIND=8), INTENT(OUT), DIMENSION(:,:) :: jf

        ! Physical data
        ! Mixture density and species mass fractions
        CALL library_get_density (prim(1:nb_ns), rho) 
       
        ! Species mass fractions
        CALL library_get_mass_fractions (rho, prim(1:nb_ns), yi)

        ! Velocity, kinetic energy per unit mass and translational temperatures
        ! (heavy-particles and free electrons) 
        u  = prim(pos_u)
        ek = 0.5d0*u**2
        T  = prim(pos_T)
        Te = prim(pos_Te)

        ! Temperature
        DO i = 1,nb_temp
           temp(i) = prim(pos_T + i - 1)
        ENDDO 

        ! Compute thermodynamic quantities needed for the flux Jacobian
        CALL library_get_thermodynamic_data (rho, prim(1:nb_ns), temp, c, gamma, p, alpha, beta, ei)

        ! Compute energy densities
        CALL library_get_energy_densities (prim(1:nb_ns), temp, rho_eint)

        ! Free electron pressure 
        pe = prim(pos_em)*Ri(pos_em)*Te

        ! Common factors
        vn       = u*nx
        uvn      = u*vn
        ov_rho   = 1.d0/rho
        vn_alpha = vn*alpha
        ek_alpha = ek*alpha
        nx_alpha = nx*alpha
        nx_ovrho = nx*ov_rho
        vn_ovrho = vn*ov_rho 
        ge_m_g   = gamma_e - gamma
        ge_m_g_pe_ov_rho = ge_m_g*pe*ov_rho
        fe = ge_m_g*ov_gamma_e_m1*(rho**gamma_e_m1)

        ! Total enthalpy per unit mass
        h0 = 0.d0
        DO i = 1,nb_ns 
           h0 = h0 + yi(i)*ei(i)   
        ENDDO
        h0   = h0 + ek + p*ov_rho
        h0vn = h0*vn 

        ! Free electrons
        yivn(pos_em) = yi(pos_em)*vn
        epsi(pos_em) = ek_alpha + ge_m_g_pe_ov_rho

        ! Heavy particles 
        DO i = pos_em + 1,nb_ns 
           yivn(i) = yi(i)*vn
           tmp1    = Ri(i)*T - alpha*ei(i) + ek_alpha + ge_m_g_pe_ov_rho
           tmp2    = 0.d0
           DO k = 1,nb_int_temp
              tmp2 = tmp2 - ei(nb_ns + (i - 1)*nb_int_temp + k)
           ENDDO   
           epsi(i) = tmp1 - alpha*tmp2
        ENDDO

        ! Column j,  j = 1,..,nb_ns
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
           jf(pos_rhoE,j) = eps*vn - h0vn

           DO i = 1,nb_int_temp + nb_te
              jf(pos_rhoE + i,j) = - rho_eint(i + 1)*vn_ovrho
           ENDDO

        ENDDO 

        ! Column nb_ns + 1
        DO i = 1,nb_ns 
           jf(i,pos_rhou) = prim(i)*nx_ovrho
        ENDDO

        jf(pos_rhou,pos_rhou) = 2.d0*vn - vn_alpha
        jf(pos_rhoE,pos_rhou) = h0*nx   - alpha*uvn

        DO i = 1,nb_int_temp + nb_te 
           jf(pos_rhoE + i,pos_rhou) = rho_eint(i + 1)*nx_ovrho
        ENDDO

        ! Column nb_ns + 2
        DO i = 1,nb_ns 
           jf(i,pos_rhoE) = 0.d0
        ENDDO

        jf(pos_rhou,pos_rhoE) = nx_alpha
        jf(pos_rhoE,pos_rhoE) = vn + vn_alpha 

        DO i = 1,nb_int_temp + nb_te 
           jf(pos_rhoE + i,pos_rhoE) = 0.d0
        ENDDO

        ! Column nb_ns + 2 + k,  k = 1,..,nb_int_temp
        DO j = pos_rhoE + 1,pos_rhoE + nb_int_temp

           DO i = 1,nb_ns 
              jf(i,j) = 0.d0
           ENDDO

           jf(pos_rhou,j) = - nx_alpha
           jf(pos_rhoE,j) = - vn_alpha

           DO i = pos_rhoE + 1,j - 1
              jf(i,j) = 0.d0
           ENDDO

           jf(j,j) = vn

           DO i = j + 1,pos_rhoE + nb_int_temp + nb_te
              jf(i,j) = 0.d0
           ENDDO

        ENDDO

        ! Column nb_ns + 2 + nb_int_temp + nb_te
        DO i = 1,nb_eq
           jf(i,pos_rhose) = 0.d0
        ENDDO

        jf(pos_rhou,pos_rhose) = nx*fe
        jf(pos_rhoE,pos_rhose) = vn*fe

        DO k = 1,nb_int_temp
           jf(pos_rhoek + k - 1,pos_rhose) = 0.d0
        ENDDO 

        jf(pos_rhose,pos_rhose) = vn

      END SUBROUTINE flux_Jacobian_neqNT_Te_1D

      !----------------------------------------------------!
      !> This subroutine provides the eigensystem for the 1D Euler equations for
      !! nonequilibrium flows (1T case)
      SUBROUTINE eigensystem_neq1T_1D (nx, prim, lambda, right, left)    

        USE mod_neq_function_pointer,    ONLY: library_get_thermodynamic_data, library_get_density, & 
                                             & library_get_mass_fractions

        INTEGER :: i, j
        REAL(KIND=8) :: tmp1, tmp2, tmp3
        REAL(KIND=8) :: alpha, eps, gamma
        REAL(KIND=8) :: c, cn, ek, h0, ov_c2, ov_rho, rho, p, T, u, uc, vn, vcn
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

        ! Velocity, kinetic energy per unit mass and temperature 
        u  = prim(pos_u)
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
           right(pos_rhoE,j) = sigmai(j)

        ENDDO

        ! Column nb_ns + 1
        DO i = 1,nb_ns 
           right(i,pos_rhou) = yi(i)
        ENDDO
 
        right(pos_rhou,pos_rhou) = u - cn
        right(pos_rhoE,pos_rhou) = h0 - vcn

        ! Column nb_ns + 2
        DO i = 1,nb_ns 
           right(i,pos_rhoE) = yi(i)
        ENDDO
 
        right(pos_rhou,pos_rhoE) = u + cn
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
        left(pos_rhoE,pos_rhou) = tmp3 - tmp2

        ! Column nb_ns + 2
        tmp1 = ov_c2*alpha
        tmp2 = 0.5d0*tmp1

        DO i = 1,nb_ns 
           left(i,pos_rhoE) = - yi(i)*tmp1
        ENDDO
 
        left(pos_rhou,pos_rhoE) = tmp2
        left(pos_rhoE,pos_rhoE) = tmp2

      END SUBROUTINE eigensystem_neq1T_1D    

      !----------------------------------------------------!
      !> This subroutine provides the eigensystem for the 1D Euler equations for
      !! nonequilibrium flows (NT case)
      SUBROUTINE eigensystem_neqNT_1D (nx, prim, lambda, right, left)    

        USE mod_neq_function_pointer,    ONLY: library_get_thermodynamic_data, library_get_density, & 
                                             & library_get_mass_fractions, library_get_energy_densities 

        INTEGER :: i, j, k, kp, pos_k
        REAL(KIND=8) :: tmp1, tmp2, tmp3
        REAL(KIND=8) :: alpha, eps, gamma
        REAL(KIND=8) :: betak, c, cn, ek, eit, eistar, h0, ov_c2, ov_rho, ov_alpha, rho, p, T, u, uc, & 
                      & vn, vcn, RiT
        REAL(KIND=8), DIMENSION(nb_ns) :: yi, yivn, epsi, sigmai
        REAL(KIND=8), DIMENSION(nb_temp) :: temp, beta, rho_eint
        REAL(KIND=8), DIMENSION(nb_int_temp) :: eintk, ov_betak
        REAL(KIND=8), DIMENSION(nb_ns + nb_ns*nb_int_temp) :: ei

        REAL(KIND=8), INTENT(IN) :: nx
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prim
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: lambda
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: right, left
        
        ! Physical data
        ! Mixture density and species mass fractions
        CALL library_get_density (prim(1:nb_ns), rho) 

        ! Species mass fractions
        CALL library_get_mass_fractions (rho, prim(1:nb_ns), yi)

        ! Velocity, kinetic energy per unit mass and temperature 
        u  = prim(pos_u)
        ek = 0.5d0*u**2
        T  = prim(pos_T)

        ! Temperature
        DO i = 1,nb_temp
           temp(i) = prim(pos_T + i - 1)
        ENDDO 

        ! Compute thermodynamic quantities needed for the flux Jacobian
        CALL library_get_thermodynamic_data (rho, prim(1:nb_ns), temp, c, gamma, p, alpha, beta, ei)

        ! Compute energy densities
        CALL library_get_energy_densities (prim(1:nb_ns), temp, rho_eint)

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
        ov_alpha = 1.d0/alpha

        DO i = 1,nb_ns 
           RiT   = Ri(i)*T
           eit   = ei(i)   
           eistar = 0.d0
           DO k = 1,nb_int_temp
              eistar = eistar + ei(nb_ns + (i - 1)*nb_int_temp + k)
           ENDDO   
           epsi(i)   = RiT - alpha*(eit - eistar)
           sigmai(i) = eit - RiT*ov_alpha
        ENDDO

        epsi   = epsi + alpha*ek 
        sigmai = sigmai + ek

        DO k = 1,nb_int_temp
           eintk(k)    = rho_eint(k + 1)*ov_rho
           ov_betak(k) = 1.d0/beta(k + 1)
        ENDDO

        ! Eigenvalues 
        DO i = 1,nb_ns 
           lambda(i) = vn
        ENDDO

        lambda(pos_rhou) = vn - c 
        lambda(pos_rhoE) = vn + c

        DO k = 1,nb_int_temp
           lambda(pos_rhoek + k - 1) = vn
        ENDDO

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
           right(pos_rhoE,j) = sigmai(j)

           DO k = 1,nb_int_temp
              right(pos_rhoek + k - 1,j) = ei(nb_ns + nb_int_temp*(j - 1) + k)
           ENDDO

        ENDDO

        ! Column nb_ns + 1
        DO i = 1,nb_ns 
           right(i,pos_rhou) = yi(i)
        ENDDO
 
        right(pos_rhou,pos_rhou) = u - cn
        right(pos_rhoE,pos_rhou) = h0 - vcn

        DO k = 1,nb_int_temp
           right(pos_rhoek + k - 1,pos_rhou) = eintk(k)
        ENDDO

        ! Column nb_ns + 2
        DO i = 1,nb_ns 
           right(i,pos_rhoE) = yi(i)
        ENDDO
 
        right(pos_rhou,pos_rhoE) = u + cn
        right(pos_rhoE,pos_rhoE) = h0 + vcn

        DO k = 1,nb_int_temp
           right(pos_rhoek + k - 1,pos_rhoE) = eintk(k)
        ENDDO

        ! Column nb_ns + 2 + k,  k = 1,..,nb_int_temp 
        DO k = 1,nb_int_temp

           pos_k = pos_rhoek + k - 1
           betak = beta(k + 1)

           DO i = 1,nb_ns 
              right(i,pos_k) = 0.d0
           ENDDO

           right(pos_rhou,pos_k) = 0.d0 
           right(pos_rhoE,pos_k) = betak

           DO kp = 1,k - 1
              right(pos_rhoek + kp - 1,pos_k) = 0.d0
           ENDDO

           right(pos_k,pos_k) = betak

           DO kp = k + 1,nb_int_temp
              right(pos_rhoek + kp - 1,pos_k) = 0.d0
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
              left(i,j) = - yi(i)*tmp2
           ENDDO

           left(j,j) = 1.d0 - yi(j)*tmp2

           DO i = j + 1,nb_ns
              left(i,j) = - yi(i)*tmp2
           ENDDO

           left(pos_rhou,j) = tmp3 + tmp1 
           left(pos_rhoE,j) = tmp3 - tmp1 

           DO k = 1,nb_int_temp
              left(pos_rhoek + k - 1,j) = - ei(nb_ns + nb_int_temp*(j - 1) + k)*ov_betak(k)
           ENDDO

        ENDDO

        ! Column nb_ns + 1	
        tmp1 = alpha*u*ov_c2
        tmp2 = 0.5d0*tmp1
        tmp3 = 0.5d0*nx/c

        DO i = 1,nb_ns 
           left(i,pos_rhou) = yi(i)*tmp1
        ENDDO
 
        left(pos_rhou,pos_rhou) = - (tmp3 + tmp2)
        left(pos_rhoE,pos_rhou) = tmp3 - tmp2

        DO k = 1,nb_int_temp
           left(pos_rhoek + k - 1,pos_rhou) = 0.d0
        ENDDO

        ! Column nb_ns + 2
        tmp1 = ov_c2*alpha
        tmp2 = 0.5d0*tmp1

        DO i = 1,nb_ns 
           left(i,pos_rhoE) = - yi(i)*tmp1
        ENDDO
 
        left(pos_rhou,pos_rhoE) = tmp2
        left(pos_rhoE,pos_rhoE) = tmp2

        DO k = 1,nb_int_temp
           left(pos_rhoek + k - 1,pos_rhoE) = 0.d0
        ENDDO 

        ! Column nb_ns + 2 + k,  k = 1,..,nb_int_temp 
        DO k = 1,nb_int_temp

           pos_k = pos_rhoek + k - 1

           DO i = 1,nb_ns 
              left(i,pos_k) = yi(i)*tmp1
           ENDDO

           left(pos_rhou,pos_k) = - tmp2
           left(pos_rhoE,pos_k) = - tmp2

           DO kp = 1,k - 1
              left(pos_rhoek + kp - 1,pos_k) = 0.d0  
           ENDDO

           left(pos_k,pos_k) = ov_betak(k)
 
           DO kp = k + 1,nb_int_temp
              left(pos_rhoek + kp - 1,pos_k) = 0.d0
           ENDDO

        ENDDO

      END SUBROUTINE eigensystem_neqNT_1D

      !----------------------------------------------------!
      !> This subroutine provides the eigensystem for the 1D Euler equations for
      !! nonequilibrium flows (NT - Te case)
      SUBROUTINE eigensystem_neqNT_Te_1D (nx, prim, lambda, right, left)    

        USE mod_neq_function_pointer,    ONLY: library_get_thermodynamic_data, library_get_density, & 
                                             & library_get_mass_fractions, library_get_energy_densities 

        INTEGER :: i, j, k, kp, pos_k
        REAL(KIND=8) :: tmp1, tmp2, tmp3
        REAL(KIND=8) :: alpha, eps, gamma, ge_m_g, ge_m1_pe_ov_rho, ge_m_g_pe_ov_rho, rho_pow_ge,   & 
                      & rho_pow_ge_m1, pe_gamma_e_ov_rhoc2, pe_ge_m1_ov_rho_pow_ge,                 & 
                      & alpha_pe_gamma_e_ov_rhoc2
        REAL(KIND=8) :: betak, c, cn, ek, eit, eistar, fe, fe_ov_c2, h0, ov_c2, ov_rho, ov_alpha,   & 
                      & rho, p, pe, se, T, Te, u, uc, vn, vcn, RiT
        REAL(KIND=8), DIMENSION(nb_ns) :: yi, epsi, sigmai
        REAL(KIND=8), DIMENSION(nb_temp) :: temp, beta, rho_eint
        REAL(KIND=8), DIMENSION(MAX(1,nb_int_temp)) :: eintk, ov_betak
        REAL(KIND=8), DIMENSION(nb_ns + nb_ns*nb_int_temp) :: ei

        REAL(KIND=8), INTENT(IN) :: nx
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prim
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: lambda
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: right, left
                
        ! Physical data
        ! Mixture density and species mass fractions
        CALL library_get_density (prim(1:nb_ns), rho) 

        ! Species mass fractions
        CALL library_get_mass_fractions (rho, prim(1:nb_ns), yi)

        ! Velocity, kinetic energy per unit mass and translational temperatures
        ! (heavy-particle and free electrons) 
        u  = prim(pos_u)
        ek = 0.5d0*u**2
        T  = prim(pos_T)
        Te = prim(pos_Te)

        ! Temperature
        DO i = 1,nb_temp
           temp(i) = prim(pos_T + i - 1)
        ENDDO 

        ! Compute thermodynamic quantities needed for the flux Jacobian
        CALL library_get_thermodynamic_data (rho, prim(1:nb_ns), temp, c, gamma, p, alpha, beta, ei)

        ! Compute energy densities
        CALL library_get_energy_densities (prim(1:nb_ns), temp, rho_eint)

        ! Total enthalpy per unit mass
        ov_rho = 1.d0/rho
        h0 = 0.d0
        DO i = 1,nb_ns 
           h0 = h0 + yi(i)*ei(i)   
        ENDDO
        h0 = h0 + ek + p*ov_rho

        ! Free electron pressure and specific pseudo-entropy
        se = rho_eint(nb_temp)*ov_rho
        pe = prim(pos_em)*Ri(pos_em)*Te

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
        epsi(pos_em)   = 0.d0  
        sigmai(pos_em) = 0.d0 

        ! Heavy-particles 
        DO i = pos_em + 1,nb_ns 
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

        DO k = 1,nb_int_temp
           eintk(k)    = rho_eint(k + 1)*ov_rho
           ov_betak(k) = 1.d0/beta(k + 1)
        ENDDO

        ! Eigenvalues 
        DO i = 1,nb_ns 
           lambda(i) = vn
        ENDDO

        lambda(pos_rhou) = vn - c 
        lambda(pos_rhoE) = vn + c

        DO k = 1,nb_int_temp
           lambda(pos_rhoek + k - 1) = vn
        ENDDO

        lambda(pos_rhose) = vn

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
           right(pos_rhoE,j) = sigmai(j)

           DO k = 1,nb_int_temp
              right(pos_rhoek + k - 1,j) = ei(nb_ns + nb_int_temp*(j - 1) + k)
           ENDDO

           right(pos_rhose,j) = - pe_ge_m1_ov_rho_pow_ge

        ENDDO

        ! Column nb_ns + 1
        DO i = 1,nb_ns 
           right(i,pos_rhou) = yi(i)
        ENDDO
 
        right(pos_rhou,pos_rhou) = u - cn
        right(pos_rhoE,pos_rhou) = h0 - vcn

        DO k = 1,nb_int_temp
           right(pos_rhoek + k - 1,pos_rhou) = eintk(k)
        ENDDO

        right(pos_rhose,pos_rhou) = se

        ! Column nb_ns + 2
        DO i = 1,nb_ns 
           right(i,pos_rhoE) = yi(i)
        ENDDO
 
        right(pos_rhou,pos_rhoE) = u + cn
        right(pos_rhoE,pos_rhoE) = h0 + vcn

        DO k = 1,nb_int_temp
           right(pos_rhoek + k - 1,pos_rhoE) = eintk(k)
        ENDDO

        right(pos_rhose,pos_rhoE) = se

        ! Column nb_ns + 2 + k,  k = 1,..,nb_int_temp 
        DO k = 1,nb_int_temp

           pos_k = pos_rhoek + k - 1
           betak = beta(k + 1)

           DO i = 1,nb_ns 
              right(i,pos_k) = 0.d0
           ENDDO

           right(pos_rhou,pos_k) = 0.d0 
           right(pos_rhoE,pos_k) = betak

           DO kp = 1,k - 1
              right(pos_rhoek + kp - 1,pos_k) = 0.d0
           ENDDO

           right(pos_k,pos_k) = betak

           DO kp = k + 1,nb_int_temp
              right(pos_rhoek + kp - 1,pos_k) = 0.d0
           ENDDO

           right(pos_rhose,pos_k) = 0.d0

        ENDDO

        ! Column nb_ns + 2 + nb_int_temp + nb_te
        DO i = 1,nb_ns 
           right(i,pos_rhose) = 0.d0
        ENDDO

        right(pos_rhou,pos_rhose) = 0.d0
        right(pos_rhoE,pos_rhose) = - ge_m_g*ov_gamma_e_m1/(gamma - 1.d0)

        DO k = 1,nb_int_temp
           right(pos_rhoek + k - 1,pos_rhose) = 0.d0
        ENDDO

        right(pos_rhose,pos_rhose) = 1.d0/rho_pow_ge_m1

        ! Left eigenvector matrix 
        ! Column j,  j = 1,..,nb_ns
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
           left(pos_rhoE,j) = tmp3 - tmp1 

           DO k = 1,nb_int_temp
              left(pos_rhoek + k - 1,j) = - ei(nb_ns + nb_int_temp*(j - 1) + k)*ov_betak(k)
           ENDDO

           left(pos_rhose,j) = ge_m1_pe_ov_rho - epsi(j)*pe_gamma_e_ov_rhoc2

        ENDDO

        ! Column nb_ns + 1	
        tmp1 = alpha*u*ov_c2
        tmp2 = 0.5d0*tmp1
        tmp3 = 0.5d0*nx/c

        DO i = 1,nb_ns 
           left(i,pos_rhou) = yi(i)*tmp1
        ENDDO
 
        left(pos_rhou,pos_rhou) = - (tmp3 + tmp2)
        left(pos_rhoE,pos_rhou) = tmp3 - tmp2

        DO k = 1,nb_int_temp
           left(pos_rhoek + k - 1,pos_rhou) = 0.d0
        ENDDO

        left(pos_rhose,pos_rhou) = alpha_pe_gamma_e_ov_rhoc2*u

        ! Column nb_ns + 2
        tmp1 = ov_c2*alpha
        tmp2 = 0.5d0*tmp1

        DO i = 1,nb_ns 
           left(i,pos_rhoE) = - yi(i)*tmp1
        ENDDO
 
        left(pos_rhou,pos_rhoE) = tmp2
        left(pos_rhoE,pos_rhoE) = tmp2

        DO k = 1,nb_int_temp
           left(pos_rhoek + k - 1,pos_rhoE) = 0.d0
        ENDDO 

        left(pos_rhose,pos_rhoE) = - alpha_pe_gamma_e_ov_rhoc2 

        ! Column nb_ns + 2 + k,  k = 1,..,nb_int_temp 
        DO k = 1,nb_int_temp

           pos_k = pos_rhoek + k - 1

           DO i = 1,nb_ns 
              left(i,pos_k) = yi(i)*tmp1
           ENDDO

           left(pos_rhou,pos_k) = - tmp2
           left(pos_rhoE,pos_k) = - tmp2

           DO kp = 1,k - 1
              left(pos_rhoek + kp - 1,pos_k) = 0.d0  
           ENDDO

           left(pos_k,pos_k) = ov_betak(k)
 
           DO kp = k + 1,nb_int_temp
              left(pos_rhoek + kp - 1,pos_k) = 0.d0
           ENDDO

           left(pos_rhose,pos_k) = alpha_pe_gamma_e_ov_rhoc2

        ENDDO

        ! Column nb_ns + 2 + nb_int_temp + nb_te 
        DO i = 1,nb_ns
           left(i,pos_rhose) = - yi(i)*fe_ov_c2
        ENDDO

        tmp1 = 0.5d0*fe_ov_c2
        left(pos_rhou,pos_rhose) = tmp1
        left(pos_rhoE,pos_rhose) = tmp1

        DO k = 1,nb_int_temp
           left(pos_rhoek + k - 1,pos_rhose) = 0.d0
        ENDDO

        left(pos_rhose,pos_rhose) = rho_pow_ge_m1*(1.d0 - pe_gamma_e_ov_rhoc2*ge_m_g*ov_gamma_e_m1)

      END SUBROUTINE eigensystem_neqNT_Te_1D

      !----------------------------------------------------!
      !> This subroutine computes the positive and negative splits of the inviscid 
      !! flux Jacobian of the 1D Euler equations for nonequilibrium flows (1T case).
      SUBROUTINE flux_Jacobian_split_neq1T_1D (nx, prim, Aplus, Aminus)

        USE mod_neq_function_pointer,    ONLY: library_get_thermodynamic_data, library_get_density, & 
                                             & library_get_mass_fractions

        INTEGER :: i, j
        REAL(KIND=8) :: tmp1, tmp2
        REAL(KIND=8) :: alpha, eps, gamma, nx2
        REAL(KIND=8) :: c, cn, ek, h0, ov_rho, rho, u, u_vn, vn, vn2, T, p
        REAL(KIND=8) :: alpha_ov_c2, ov_c, ov_c2, yi_epsj_ov_c2, yi_vn_ov_c, epsj_ov_c2, ov_alpha
        REAL(KIND=8) :: l1, l2, l3, l1p, l2p, l3p, l1m, l2m, l3m, l1m_u, l1p_u, l1p_ov_gamma_m1, l1m_ov_gamma_m1
        REAL(KIND=8) :: eig_diff, eig_sum, eig_diff_vn_ov_c, eig_sum_m_l1p, eig_sum_m_l1m,    & 
                      & u_eig_sum_m_l1p, h0_eig_sum_m_l1p, u_eig_sum_m_l1m, h0_eig_sum_m_l1m, & 
                      & eig_sum_vn_nx, eig_diff_ov_c, eig_diff_nx_ov_c, eig_sum_vn2,          & 
                      & u_eig_sum_m_l1p_alpha_ov_c2, u_eig_sum_m_l1m_alpha_ov_c2
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

        ! Velocity, kinetic energy per unit mass and temperature 
        u  = prim(pos_u)
        ek = 0.5d0*u**2
        T  = prim(pos_T)

        ! Temperature
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
        h0_eig_sum_m_l1p = h0*eig_sum_m_l1p
        u_eig_sum_m_l1p_alpha_ov_c2 = u_eig_sum_m_l1p*alpha_ov_c2

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

        ! Column j,  j = 1,..,nb_ns 
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
           Aplus(pos_rhoE,j) = epsj_ov_c2*h0_eig_sum_m_l1p + l1p*sigmai(j) + epsi(j)*l1p_ov_gamma_m1 & 
                             & - eig_sum_vn2 + eig_diff_vn_ov_c*(epsi(j) - h0)

        ENDDO

        ! Column nb_ns + 1 
        tmp1 = eig_diff_nx_ov_c  - u_eig_sum_m_l1p_alpha_ov_c2
        DO i = 1,nb_ns 
           Aplus(i,pos_rhou) = tmp1*yi(i)
        ENDDO
  
        Aplus(pos_rhou,pos_rhou) = - u_eig_sum_m_l1p_alpha_ov_c2*u  + eig_sum*nx2 + eig_diff_vn_ov_c*(1.d0 - alpha)
        Aplus(pos_rhoE,pos_rhou) = - u_eig_sum_m_l1p_alpha_ov_c2*h0 - l1p_u + eig_sum*vn*nx + eig_diff_ov_c*(h0*nx - alpha*u_vn)

        ! Column nb_ns + 2
        tmp1 = eig_sum_m_l1p*alpha_ov_c2
        DO i = 1,nb_ns 
           Aplus(i,pos_rhoE) = tmp1*yi(i)
        ENDDO

        Aplus(pos_rhou,pos_rhoE) = u_eig_sum_m_l1p_alpha_ov_c2  + alpha*eig_diff_nx_ov_c 
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
        h0_eig_sum_m_l1m = h0*eig_sum_m_l1m
        u_eig_sum_m_l1m_alpha_ov_c2 = u_eig_sum_m_l1m*alpha_ov_c2

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
           Aminus(pos_rhoE,j) = epsj_ov_c2*h0_eig_sum_m_l1m + l1m*sigmai(j) + epsi(j)*l1m_ov_gamma_m1 & 
                              & - eig_sum_vn2 + eig_diff_vn_ov_c*(epsi(j) - h0)

        ENDDO

        ! Column nb_ns + 1 
        tmp1 = eig_diff_nx_ov_c - u_eig_sum_m_l1m_alpha_ov_c2
        DO i = 1,nb_ns 
           Aminus(i,pos_rhou) = tmp1*yi(i) 
        ENDDO
  
        Aminus(pos_rhou,pos_rhou) = - u_eig_sum_m_l1m_alpha_ov_c2*u  + eig_sum*nx2 + eig_diff_vn_ov_c*(1.d0 - alpha)
        Aminus(pos_rhoE,pos_rhou) = - u_eig_sum_m_l1m_alpha_ov_c2*h0 - l1m_u + eig_sum*vn*nx + eig_diff_ov_c*(h0*nx - alpha*u_vn)

        ! Column nb_ns + 2
        tmp1 = eig_sum_m_l1m*alpha_ov_c2
        DO i = 1,nb_ns 
           Aminus(i,pos_rhoE) = tmp1*yi(i)
        ENDDO

        Aminus(pos_rhou,pos_rhoE) = u_eig_sum_m_l1m_alpha_ov_c2  + alpha*eig_diff_nx_ov_c
        Aminus(pos_rhoE,pos_rhoE) = h0_eig_sum_m_l1m*alpha_ov_c2 + l1m + eig_diff_vn_ov_c*alpha 

      END SUBROUTINE flux_Jacobian_split_neq1T_1D

      !----------------------------------------------------!
      !> This subroutine computes the positive and negative splits of the inviscid 
      !! flux Jacobian of the 1D Euler equations for nonequilibrium flows (NT case).
      SUBROUTINE flux_Jacobian_split_neqNT_1D (nx, prim, Aplus, Aminus)

        USE mod_neq_function_pointer,    ONLY: library_get_thermodynamic_data, library_get_density,     & 
                                             & library_get_mass_fractions, library_get_energy_densities

        INTEGER :: i, j, k, kp
        INTEGER :: pos_k
        REAL(KIND=8) :: tmp1, tmp2, tmp3
        REAL(KIND=8) :: alpha, eps, gamma, nx2
        REAL(KIND=8) :: c, cn, ek, h0, ov_rho, rho, u, u_vn, vn, vn2, T, p
        REAL(KIND=8) :: alpha_ov_c2, alpha_u_ov_c2, ov_c, ov_c2, yi_epsj_ov_c2, yi_vn_ov_c, epsj_ov_c2, ov_alpha
        REAL(KIND=8) :: l1, l2, l3, l1p, l2p, l3p, l1m, l2m, l3m, l1m_u, l1p_u, l1p_ov_gamma_m1, l1m_ov_gamma_m1
        REAL(KIND=8) :: eig_diff, eig_sum, eig_diff_vn_ov_c, eig_diff_vn_ov_c_alpha, eig_sum_m_l1p,             & 
                      & eig_sum_m_l1m, u_eig_sum_m_l1p, h0_eig_sum_m_l1p, u_eig_sum_m_l1m, h0_eig_sum_m_l1m,    & 
                      & eig_sum_vn_nx, eig_diff_ov_c, eig_diff_nx_ov_c, eig_diff_nx_ov_c_alpha, eig_sum_vn2,    & 
                      & eig_sum_m_l1p_alpha_ov_c2, eig_sum_m_l1m_alpha_ov_c2, u_eig_sum_m_l1p_alpha_ov_c2,      & 
                      & u_eig_sum_m_l1m_alpha_ov_c2, h0_eig_sum_m_l1p_alpha_ov_c2, h0_eig_sum_m_l1m_alpha_ov_c2
        REAL(KIND=8), DIMENSION(nb_ns) :: yi, epsi, sigmai, eistar
        REAL(KIND=8), DIMENSION(nb_ns*(nb_int_temp + 1)) :: ei
        REAL(KIND=8), DIMENSION(nb_temp) :: temp, beta, rho_eint
        REAL(KIND=8), DIMENSION(nb_int_temp) :: eintk

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

        ! Velocity, kinetic energy per unit mass and temperature 
        u  = prim(pos_u)
        ek = 0.5d0*u**2
        T  = prim(pos_T)

        ! Temperature
        DO i = 1,nb_temp
           temp(i) = prim(pos_u + i)
        ENDDO 

        ! Compute thermodynamic quantities needed for the flux Jacobian
        CALL library_get_thermodynamic_data (rho, prim(1:nb_ns), temp, c, gamma, p, alpha, beta, ei)

        ! Compute energy densities
        CALL library_get_energy_densities (prim(1:nb_ns), temp, rho_eint) 

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
        vn2   = vn**2 
        ov_c  = 1.d0/c
        ov_c2 = ov_c/c
        ov_alpha    = 1.d0/alpha
        alpha_ov_c2 = alpha*ov_c2
        alpha_u_ov_c2 = alpha_ov_c2*u

        l1p_u = l1p*u

        eig_sum_m_l1p = eig_sum - l1p
        eig_sum_vn_nx = eig_sum*vn*nx
        eig_sum_vn2   = eig_sum*vn2

        u_eig_sum_m_l1p  = u*eig_sum_m_l1p
        h0_eig_sum_m_l1p = h0*eig_sum_m_l1p
        h0_eig_sum_m_l1p_alpha_ov_c2 = h0_eig_sum_m_l1p*alpha_ov_c2
        eig_sum_m_l1p_alpha_ov_c2    = eig_sum_m_l1p*alpha_ov_c2
        u_eig_sum_m_l1p_alpha_ov_c2  = u_eig_sum_m_l1p*alpha_ov_c2

        eig_diff_ov_c    = eig_diff*ov_c
        eig_diff_nx_ov_c = eig_diff_ov_c*nx
        eig_diff_nx_ov_c_alpha = eig_diff_nx_ov_c*alpha
        eig_diff_vn_ov_c       = eig_diff_ov_c*vn
        eig_diff_vn_ov_c_alpha = eig_diff_vn_ov_c*alpha

        l1p_ov_gamma_m1  = l1p/(gamma - 1.d0)

        ! Energy related data
        DO i = 1,nb_ns 
          tmp1  = Ri(i)*T
          tmp2  = ei(i)
          tmp3  = 0.d0
          pos_k = nb_ns + nb_int_temp*(i - 1)
          DO k = 1,nb_int_temp
             tmp3     = tmp3 + ei(pos_k + k)
             eintk(k) = rho_eint(k + 1)*ov_rho
          ENDDO
          epsi(i)   = tmp1 - alpha*(tmp2 - tmp3) 
          sigmai(i) = tmp2 - tmp1*ov_alpha 
          eistar(i) = tmp3
        ENDDO
        
        epsi   = epsi + alpha*ek 
        sigmai = sigmai + ek

        ! Column j,  j = 1,..,nb_ns 
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
           Aplus(pos_rhoE,j) = epsj_ov_c2*h0_eig_sum_m_l1p + l1p*(sigmai(j) - eistar(j)) + epsi(j)*l1p_ov_gamma_m1 & 
                             & - eig_sum_vn2 + eig_diff_vn_ov_c*(epsi(j) - h0)

           DO k = 1,nb_int_temp
              Aplus(pos_rhoek + k - 1,j) = eintk(k)*(epsj_ov_c2*eig_sum_m_l1p - eig_diff_vn_ov_c)  
           ENDDO

        ENDDO

        ! Column nb_ns + 1 
        tmp1 = eig_diff_nx_ov_c  - u_eig_sum_m_l1p_alpha_ov_c2
        DO i = 1,nb_ns 
           Aplus(i,pos_rhou) = tmp1*yi(i)
        ENDDO
  
        Aplus(pos_rhou,pos_rhou) = - u_eig_sum_m_l1p_alpha_ov_c2*u  + eig_sum*nx2 + eig_diff_vn_ov_c*(1.d0 - alpha)
        Aplus(pos_rhoE,pos_rhou) = - u_eig_sum_m_l1p_alpha_ov_c2*h0 - l1p_u + eig_sum*vn*nx + eig_diff_ov_c*(h0*nx - alpha*u_vn)

        DO k = 1,nb_int_temp
           Aplus(pos_rhoek + k - 1,pos_rhou) = eintk(k)*(eig_diff_nx_ov_c - eig_sum_m_l1p*alpha_u_ov_c2)
        ENDDO

        ! Column nb_ns + 2
        tmp1 = eig_sum_m_l1p*alpha_ov_c2
        DO i = 1,nb_ns 
           Aplus(i,pos_rhoE) = tmp1*yi(i)
        ENDDO

        Aplus(pos_rhou,pos_rhoE) = u_eig_sum_m_l1p_alpha_ov_c2  + eig_diff_nx_ov_c_alpha 
        Aplus(pos_rhoE,pos_rhoE) = h0_eig_sum_m_l1p_alpha_ov_c2 + l1p + eig_diff_vn_ov_c_alpha

        DO k = 1,nb_int_temp
           Aplus(pos_rhoek + k - 1,pos_rhoE) = eintk(k)*eig_sum_m_l1p_alpha_ov_c2
        ENDDO

        ! Column nb_ns + 2 + k,  k = 1,..,nb_int_temp
        DO k = 1,nb_int_temp

           pos_k = pos_rhoek + k - 1

           DO i = 1,nb_ns
              Aplus(i,pos_k) = - yi(i)*eig_sum_m_l1p_alpha_ov_c2
           ENDDO

           Aplus(pos_rhou,pos_k) = - u_eig_sum_m_l1p_alpha_ov_c2 - eig_diff_nx_ov_c_alpha
           Aplus(pos_rhoE,pos_k) = - h0_eig_sum_m_l1p_alpha_ov_c2 - eig_diff_vn_ov_c_alpha  

           DO kp = 1,k - 1
              Aplus(pos_rhoek + kp - 1,pos_k) = - eintk(kp)*eig_sum_m_l1p_alpha_ov_c2
           ENDDO

           Aplus(pos_k,pos_k) = l1p - eintk(k)*eig_sum_m_l1p_alpha_ov_c2

           DO kp = k + 1,nb_int_temp
              Aplus(pos_rhoek + kp - 1,pos_k) = - eintk(kp)*eig_sum_m_l1p_alpha_ov_c2
           ENDDO

        ENDDO

        ! A^+ matrix (positive split Jacobian)
        ! Common factors 
        eig_sum  = 0.5d0*(l3m + l2m)
        eig_diff = 0.5d0*(l3m - l2m) 

        l1m_u = l1m*u

        eig_sum_m_l1m = eig_sum - l1m
        eig_sum_vn_nx = eig_sum*vn*nx
        eig_sum_vn2   = eig_sum*vn2

        u_eig_sum_m_l1m  = u*eig_sum_m_l1m
        h0_eig_sum_m_l1m = h0*eig_sum_m_l1m
        h0_eig_sum_m_l1m_alpha_ov_c2 = h0_eig_sum_m_l1m*alpha_ov_c2
        eig_sum_m_l1m_alpha_ov_c2    = eig_sum_m_l1m*alpha_ov_c2
        u_eig_sum_m_l1m_alpha_ov_c2  = u_eig_sum_m_l1m*alpha_ov_c2

        eig_diff_ov_c    = eig_diff*ov_c
        eig_diff_nx_ov_c = eig_diff_ov_c*nx
        eig_diff_nx_ov_c_alpha = eig_diff_nx_ov_c*alpha
        eig_diff_vn_ov_c       = eig_diff_ov_c*vn
        eig_diff_vn_ov_c_alpha = eig_diff_vn_ov_c*alpha

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
           Aminus(pos_rhoE,j) = epsj_ov_c2*h0_eig_sum_m_l1m + l1m*(sigmai(j) - eistar(j)) + epsi(j)*l1m_ov_gamma_m1 & 
                             & - eig_sum_vn2 + eig_diff_vn_ov_c*(epsi(j) - h0)

           DO k = 1,nb_int_temp
              Aminus(pos_rhoek + k - 1,j) = eintk(k)*(epsj_ov_c2*eig_sum_m_l1m - eig_diff_vn_ov_c)  
           ENDDO

        ENDDO

        ! Column nb_ns + 1 
        tmp1 = eig_diff_nx_ov_c  - u_eig_sum_m_l1m_alpha_ov_c2
        DO i = 1,nb_ns 
           Aminus(i,pos_rhou) = tmp1*yi(i)
        ENDDO
  
        Aminus(pos_rhou,pos_rhou) = - u_eig_sum_m_l1m_alpha_ov_c2*u  + eig_sum*nx2 + eig_diff_vn_ov_c*(1.d0 - alpha)
        Aminus(pos_rhoE,pos_rhou) = - u_eig_sum_m_l1m_alpha_ov_c2*h0 - l1m_u + eig_sum*vn*nx + eig_diff_ov_c*(h0*nx - alpha*u_vn)

        DO k = 1,nb_int_temp
           Aminus(pos_rhoek + k - 1,pos_rhou) = eintk(k)*(eig_diff_nx_ov_c - eig_sum_m_l1m*alpha_u_ov_c2)
        ENDDO

        ! Column nb_ns + 2
        tmp1 = eig_sum_m_l1m*alpha_ov_c2
        DO i = 1,nb_ns 
           Aminus(i,pos_rhoE) = tmp1*yi(i)
        ENDDO

        Aminus(pos_rhou,pos_rhoE) = u_eig_sum_m_l1m_alpha_ov_c2  + eig_diff_nx_ov_c_alpha 
        Aminus(pos_rhoE,pos_rhoE) = h0_eig_sum_m_l1m_alpha_ov_c2 + l1m + eig_diff_vn_ov_c_alpha

        DO k = 1,nb_int_temp
           Aminus(pos_rhoek + k - 1,pos_rhoE) = eintk(k)*eig_sum_m_l1m_alpha_ov_c2
        ENDDO

        ! Column nb_ns + 2 + k,  k = 1,..,nb_int_temp
        DO k = 1,nb_int_temp

           pos_k = pos_rhoek + k - 1

           DO i = 1,nb_ns
              Aminus(i,pos_k) = - yi(i)*eig_sum_m_l1m_alpha_ov_c2
           ENDDO

           Aminus(pos_rhou,pos_k) = - u_eig_sum_m_l1m_alpha_ov_c2 - eig_diff_nx_ov_c_alpha
           Aminus(pos_rhoE,pos_k) = - h0_eig_sum_m_l1m_alpha_ov_c2 - eig_diff_vn_ov_c_alpha  

           DO kp = 1,k - 1
              Aminus(pos_rhoek + kp - 1,pos_k) = - eintk(kp)*eig_sum_m_l1m_alpha_ov_c2
           ENDDO

           Aminus(pos_k,pos_k) = l1m - eintk(k)*eig_sum_m_l1m_alpha_ov_c2

           DO kp = k + 1,nb_int_temp
              Aminus(pos_rhoek + kp - 1,pos_k) = - eintk(kp)*eig_sum_m_l1m_alpha_ov_c2
           ENDDO

        ENDDO        
 
      END SUBROUTINE flux_Jacobian_split_neqNT_1D 

      !----------------------------------------------------!
      !> This subroutine computes the positive and negative splits of the inviscid 
      !! flux Jacobian of the 1D Euler equations for nonequilibrium flows (NT - Te case).
      SUBROUTINE flux_Jacobian_split_neqNT_Te_1D (nx, prim, Aplus, Aminus)

        USE mod_neq_function_pointer,    ONLY: library_get_thermodynamic_data, library_get_density, & 
                                             & library_get_mass_fractions, library_get_energy_densities

        INTEGER :: i, j, k, kp
        INTEGER :: pos_k
        REAL(KIND=8) :: tmp1, tmp2, tmp3
        REAL(KIND=8) :: alpha, eps, gamma, ge_m_g, ge_m_g_ov_ge_m1, rho_pow_ge_m1, nx2
        REAL(KIND=8) :: c, cn, ek, h0, ov_rho, rho, u, u_vn, vn, vn2, T, Te, p, pe, se
        REAL(KIND=8) :: alpha_ov_c2, alpha_u_ov_c2, ov_c, ov_c2, yi_epsj_ov_c2, yi_vn_ov_c, epsj_ov_c2, ov_alpha
        REAL(KIND=8) :: l1, l2, l3, l1p, l2p, l3p, l1m, l2m, l3m, l1m_u, l1p_u, l1p_ov_gamma_m1, l1m_ov_gamma_m1
        REAL(KIND=8) :: eig_diff, eig_sum, eig_diff_vn_ov_c, eig_diff_vn_ov_c_alpha, eig_sum_m_l1p,              & 
                      & eig_sum_m_l1m, u_eig_sum_m_l1p, h0_eig_sum_m_l1p, u_eig_sum_m_l1m, h0_eig_sum_m_l1m,     & 
                      & eig_sum_vn_nx, eig_diff_ov_c, eig_diff_nx_ov_c, eig_diff_vn_ov_c_se,                     & 
                      & eig_diff_nx_ov_c_alpha, eig_sum_vn2, eig_sum_m_l1p_alpha_ov_c2,                          & 
                      & eig_sum_m_l1m_alpha_ov_c2, u_eig_sum_m_l1p_alpha_ov_c2, u_eig_sum_m_l1m_alpha_ov_c2,     & 
                      & h0_eig_sum_m_l1p_alpha_ov_c2, h0_eig_sum_m_l1m_alpha_ov_c2,                              &
                      & eig_sum_m_l1p_ov_c2, eig_sum_m_l1m_ov_c2, eig_sum_m_l1p_se, eig_sum_m_l1m_se,            & 
                      & eig_sum_m_l1p_alpha_ov_c2_se, eig_sum_m_l1m_alpha_ov_c2_se
        REAL(KIND=8), DIMENSION(nb_ns) :: yi, epsi, sigmai, eistar
        REAL(KIND=8), DIMENSION(nb_ns*(nb_int_temp + 1)) :: ei
        REAL(KIND=8), DIMENSION(nb_temp) :: temp, beta, rho_eint
        REAL(KIND=8), DIMENSION(MAX(1,nb_int_temp)) :: eintk

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

        ! Velocity, kinetic energy per unit mass and translationa temperatures
        ! (heavy particles and free electrons) 
        u  = prim(pos_u)
        ek = 0.5d0*u**2
        T  = prim(pos_T)
        Te = prim(pos_Te)

        ! Temperature
        DO i = 1,nb_temp
           temp(i) = prim(pos_u + i)
        ENDDO 

        ! Compute thermodynamic quantities needed for the flux Jacobian
        CALL library_get_thermodynamic_data (rho, prim(1:nb_ns), temp, c, gamma, p, alpha, beta, ei)

        ! Compute energy densities
        CALL library_get_energy_densities (prim(1:nb_ns), temp, rho_eint) 

        ! Total enthalpy per unit mass
        ov_rho = 1.d0/rho
        h0 = 0.d0
        DO i = 1,nb_ns 
           h0 = h0 + yi(i)*ei(i)   
        ENDDO
        h0 = h0 + ek + p*ov_rho

        ! Free electron pressure and specific pseudo-entropy
        se = rho_eint(nb_temp)*ov_rho
        pe = prim(pos_em)*Ri(pos_em)*Te

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
        vn2   = vn**2 
        ov_c  = 1.d0/c
        ov_c2 = ov_c/c
        ov_alpha    = 1.d0/alpha
        alpha_ov_c2 = alpha*ov_c2
        alpha_u_ov_c2 = alpha_ov_c2*u

        l1p_u = l1p*u

        eig_sum_m_l1p = eig_sum - l1p
        eig_sum_vn_nx = eig_sum*vn*nx
        eig_sum_vn2   = eig_sum*vn2

        u_eig_sum_m_l1p  = u*eig_sum_m_l1p
        eig_sum_m_l1p_ov_c2 = eig_sum_m_l1p*ov_c2
        eig_sum_m_l1p_se    = se*eig_sum_m_l1p
        h0_eig_sum_m_l1p    = h0*eig_sum_m_l1p
        h0_eig_sum_m_l1p_alpha_ov_c2 = h0_eig_sum_m_l1p*alpha_ov_c2
        eig_sum_m_l1p_alpha_ov_c2    = eig_sum_m_l1p*alpha_ov_c2
        eig_sum_m_l1p_alpha_ov_c2_se = eig_sum_m_l1p_alpha_ov_c2*se
        u_eig_sum_m_l1p_alpha_ov_c2  = u_eig_sum_m_l1p*alpha_ov_c2

        eig_diff_ov_c    = eig_diff*ov_c
        eig_diff_nx_ov_c = eig_diff_ov_c*nx
        eig_diff_nx_ov_c_alpha = eig_diff_nx_ov_c*alpha
        eig_diff_vn_ov_c       = eig_diff_ov_c*vn
        eig_diff_vn_ov_c_se    = eig_diff_vn_ov_c*se
        eig_diff_vn_ov_c_alpha = eig_diff_vn_ov_c*alpha

        l1p_ov_gamma_m1  = l1p/(gamma - 1.d0)       

        ge_m_g = gamma_e - gamma 
        ge_m_g_ov_ge_m1 = ge_m_g*ov_gamma_e_m1
        rho_pow_ge_m1   = rho**gamma_e_m1 

        ! Energy related data
        ! Free electrons 
        epsi(pos_em)   = 0.d0
        sigmai(pos_em) = 0.d0

        ! Heavy-particles 
        DO i = pos_em + 1,nb_ns 
          tmp1  = Ri(i)*T
          tmp2  = ei(i)
          tmp3  = 0.d0
          pos_k = nb_ns + nb_int_temp*(i - 1)
          DO k = 1,nb_int_temp
             tmp3     = tmp3 + ei(pos_k + k)
             eintk(k) = rho_eint(k + 1)*ov_rho
          ENDDO
          epsi(i)   = tmp1 - alpha*(tmp2 - tmp3) 
          sigmai(i) = tmp2 - tmp1*ov_alpha 
          eistar(i) = tmp3
        ENDDO
        
        epsi   = epsi + alpha*ek + ge_m_g*pe*ov_rho
        sigmai = sigmai + ek

        ! Column j,  j = 1,..,nb_ns 
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
           Aplus(pos_rhoE,j) = epsj_ov_c2*h0_eig_sum_m_l1p + l1p*(sigmai(j) - eistar(j)) + epsi(j)*l1p_ov_gamma_m1 & 
                             & - eig_sum_vn2 + eig_diff_vn_ov_c*(epsi(j) - h0)

           DO k = 1,nb_int_temp
              Aplus(pos_rhoek + k - 1,j) = eintk(k)*(epsj_ov_c2*eig_sum_m_l1p - eig_diff_vn_ov_c)  
           ENDDO

           Aplus(pos_rhose,j) = epsj_ov_c2*eig_sum_m_l1p_se - eig_diff_vn_ov_c_se

        ENDDO

        ! Column nb_ns + 1 
        tmp1 = eig_diff_nx_ov_c  - u_eig_sum_m_l1p_alpha_ov_c2
        DO i = 1,nb_ns 
           Aplus(i,pos_rhou) = tmp1*yi(i)
        ENDDO
  
        Aplus(pos_rhou,pos_rhou) = - u_eig_sum_m_l1p_alpha_ov_c2*u  + eig_sum*nx2 + eig_diff_vn_ov_c*(1.d0 - alpha)
        Aplus(pos_rhoE,pos_rhou) = - u_eig_sum_m_l1p_alpha_ov_c2*h0 - l1p_u + eig_sum*vn*nx + eig_diff_ov_c*(h0*nx - alpha*u_vn)

        DO k = 1,nb_int_temp
           Aplus(pos_rhoek + k - 1,pos_rhou) = eintk(k)*(eig_diff_nx_ov_c - eig_sum_m_l1p*alpha_u_ov_c2)
        ENDDO

        Aplus(pos_rhose,pos_rhou) = - u_eig_sum_m_l1p_alpha_ov_c2*se + eig_diff_nx_ov_c*se

        ! Column nb_ns + 2
        tmp1 = eig_sum_m_l1p*alpha_ov_c2
        DO i = 1,nb_ns 
           Aplus(i,pos_rhoE) = tmp1*yi(i)
        ENDDO

        Aplus(pos_rhou,pos_rhoE) = u_eig_sum_m_l1p_alpha_ov_c2  + eig_diff_nx_ov_c_alpha 
        Aplus(pos_rhoE,pos_rhoE) = h0_eig_sum_m_l1p_alpha_ov_c2 + l1p + eig_diff_vn_ov_c_alpha

        DO k = 1,nb_int_temp
           Aplus(pos_rhoek + k - 1,pos_rhoE) = eintk(k)*eig_sum_m_l1p_alpha_ov_c2
        ENDDO

        Aplus(pos_rhose,pos_rhoE) = eig_sum_m_l1p_alpha_ov_c2_se

        ! Column nb_ns + 2 + k,  k = 1,..,nb_int_temp
        DO k = 1,nb_int_temp

           pos_k = pos_rhoek + k - 1

           DO i = 1,nb_ns
              Aplus(i,pos_k) = - yi(i)*eig_sum_m_l1p_alpha_ov_c2
           ENDDO

           Aplus(pos_rhou,pos_k) = - u_eig_sum_m_l1p_alpha_ov_c2 - eig_diff_nx_ov_c_alpha
           Aplus(pos_rhoE,pos_k) = - h0_eig_sum_m_l1p_alpha_ov_c2 - eig_diff_vn_ov_c_alpha  

           DO kp = 1,k - 1
              Aplus(pos_rhoek + kp - 1,pos_k) = - eintk(kp)*eig_sum_m_l1p_alpha_ov_c2
           ENDDO

           Aplus(pos_k,pos_k) = l1p - eintk(k)*eig_sum_m_l1p_alpha_ov_c2

           DO kp = k + 1,nb_int_temp
              Aplus(pos_rhoek + kp - 1,pos_k) = - eintk(kp)*eig_sum_m_l1p_alpha_ov_c2
           ENDDO

           Aplus(pos_rhose,pos_k) = - eig_sum_m_l1p_alpha_ov_c2_se

        ENDDO

        ! Column nb_ns + 2 + nb_int_temp + nb_te
        tmp1 = ge_m_g_ov_ge_m1*rho_pow_ge_m1
        tmp2 = tmp1*eig_sum_m_l1p_ov_c2
        DO i = 1,nb_ns 
           Aplus(i,pos_rhose) = yi(i)*tmp2
        ENDDO 
        
        Aplus(pos_rhou,pos_rhose) = (u_eig_sum_m_l1p*ov_c2 + eig_diff_nx_ov_c)*tmp1
        Aplus(pos_rhoE,pos_rhose) = (h0_eig_sum_m_l1p*ov_c2 + eig_diff_vn_ov_c)*tmp1  

        DO k = 1,nb_int_temp  
           Aplus(pos_rhoek + k - 1,pos_rhose) = eintk(k)*tmp2
        ENDDO

        Aplus(pos_rhose,pos_rhose) = l1p + se*tmp2

        ! A^- matrix (negative split Jacobian)
        ! Common factors 
        eig_sum  = 0.5d0*(l3m + l2m)
        eig_diff = 0.5d0*(l3m - l2m)  

        l1m_u = l1m*u

        eig_sum_m_l1m = eig_sum - l1m
        eig_sum_vn_nx = eig_sum*vn*nx
        eig_sum_vn2   = eig_sum*vn2

        u_eig_sum_m_l1m  = u*eig_sum_m_l1m
        eig_sum_m_l1m_ov_c2 = eig_sum_m_l1m*ov_c2
        eig_sum_m_l1m_se    = se*eig_sum_m_l1m
        h0_eig_sum_m_l1m    = h0*eig_sum_m_l1m
        h0_eig_sum_m_l1m_alpha_ov_c2 = h0_eig_sum_m_l1m*alpha_ov_c2
        eig_sum_m_l1m_alpha_ov_c2    = eig_sum_m_l1m*alpha_ov_c2
        eig_sum_m_l1m_alpha_ov_c2_se = eig_sum_m_l1m_alpha_ov_c2*se
        u_eig_sum_m_l1m_alpha_ov_c2  = u_eig_sum_m_l1m*alpha_ov_c2

        eig_diff_ov_c    = eig_diff*ov_c
        eig_diff_nx_ov_c = eig_diff_ov_c*nx
        eig_diff_nx_ov_c_alpha = eig_diff_nx_ov_c*alpha
        eig_diff_vn_ov_c       = eig_diff_ov_c*vn
        eig_diff_vn_ov_c_se    = eig_diff_vn_ov_c*se
        eig_diff_vn_ov_c_alpha = eig_diff_vn_ov_c*alpha

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
           Aminus(pos_rhoE,j) = epsj_ov_c2*h0_eig_sum_m_l1m + l1m*(sigmai(j) - eistar(j)) + epsi(j)*l1m_ov_gamma_m1 & 
                             & - eig_sum_vn2 + eig_diff_vn_ov_c*(epsi(j) - h0)

           DO k = 1,nb_int_temp
              Aminus(pos_rhoek + k - 1,j) = eintk(k)*(epsj_ov_c2*eig_sum_m_l1m - eig_diff_vn_ov_c)  
           ENDDO

           Aminus(pos_rhose,j) = epsj_ov_c2*eig_sum_m_l1m_se - eig_diff_vn_ov_c_se

        ENDDO

        ! Column nb_ns + 1 
        tmp1 = eig_diff_nx_ov_c  - u_eig_sum_m_l1m_alpha_ov_c2
        DO i = 1,nb_ns 
           Aminus(i,pos_rhou) = tmp1*yi(i)
        ENDDO
  
        Aminus(pos_rhou,pos_rhou) = - u_eig_sum_m_l1m_alpha_ov_c2*u  + eig_sum*nx2 + eig_diff_vn_ov_c*(1.d0 - alpha)
        Aminus(pos_rhoE,pos_rhou) = - u_eig_sum_m_l1m_alpha_ov_c2*h0 - l1m_u + eig_sum*vn*nx + eig_diff_ov_c*(h0*nx - alpha*u_vn)

        DO k = 1,nb_int_temp
           Aminus(pos_rhoek + k - 1,pos_rhou) = eintk(k)*(eig_diff_nx_ov_c - eig_sum_m_l1m*alpha_u_ov_c2)
        ENDDO

        Aminus(pos_rhose,pos_rhou) = - u_eig_sum_m_l1m_alpha_ov_c2*se + eig_diff_nx_ov_c*se

        ! Column nb_ns + 2
        tmp1 = eig_sum_m_l1m*alpha_ov_c2
        DO i = 1,nb_ns 
           Aminus(i,pos_rhoE) = tmp1*yi(i)
        ENDDO

        Aminus(pos_rhou,pos_rhoE) = u_eig_sum_m_l1m_alpha_ov_c2  + eig_diff_nx_ov_c_alpha 
        Aminus(pos_rhoE,pos_rhoE) = h0_eig_sum_m_l1m_alpha_ov_c2 + l1m + eig_diff_vn_ov_c_alpha

        DO k = 1,nb_int_temp
           Aminus(pos_rhoek + k - 1,pos_rhoE) = eintk(k)*eig_sum_m_l1m_alpha_ov_c2
        ENDDO

        Aminus(pos_rhose,pos_rhoE) = eig_sum_m_l1m_alpha_ov_c2_se

        ! Column nb_ns + 2 + k,  k = 1,..,nb_int_temp
        DO k = 1,nb_int_temp

           pos_k = pos_rhoek + k - 1

           DO i = 1,nb_ns
              Aminus(i,pos_k) = - yi(i)*eig_sum_m_l1m_alpha_ov_c2
           ENDDO

           Aminus(pos_rhou,pos_k) = - u_eig_sum_m_l1m_alpha_ov_c2 - eig_diff_nx_ov_c_alpha
           Aminus(pos_rhoE,pos_k) = - h0_eig_sum_m_l1m_alpha_ov_c2 - eig_diff_vn_ov_c_alpha  

           DO kp = 1,k - 1
              Aminus(pos_rhoek + kp - 1,pos_k) = - eintk(kp)*eig_sum_m_l1m_alpha_ov_c2
           ENDDO

           Aminus(pos_k,pos_k) = l1m - eintk(k)*eig_sum_m_l1m_alpha_ov_c2

           DO kp = k + 1,nb_int_temp
              Aminus(pos_rhoek + kp - 1,pos_k) = - eintk(kp)*eig_sum_m_l1m_alpha_ov_c2
           ENDDO

           Aminus(pos_rhose,pos_k) = - eig_sum_m_l1m_alpha_ov_c2_se

        ENDDO

        ! Column nb_ns + 2 + nb_int_temp + nb_te
        tmp1 = ge_m_g_ov_ge_m1*rho_pow_ge_m1
        tmp2 = tmp1*eig_sum_m_l1m_ov_c2
        DO i = 1,nb_ns 
           Aminus(i,pos_rhose) = yi(i)*tmp2
        ENDDO 
        
        Aminus(pos_rhou,pos_rhose) = (u_eig_sum_m_l1m*ov_c2 + eig_diff_nx_ov_c)*tmp1
        Aminus(pos_rhoE,pos_rhose) = (h0_eig_sum_m_l1m*ov_c2 + eig_diff_vn_ov_c)*tmp1  

        DO k = 1,nb_int_temp  
           Aminus(pos_rhoek + k - 1,pos_rhose) = eintk(k)*tmp2
        ENDDO

        Aminus(pos_rhose,pos_rhose) = l1m + se*tmp2

     END SUBROUTINE flux_Jacobian_split_neqNT_Te_1D

  END MODULE mod_neq_1D
!------------------------------------------------------------------------------!
