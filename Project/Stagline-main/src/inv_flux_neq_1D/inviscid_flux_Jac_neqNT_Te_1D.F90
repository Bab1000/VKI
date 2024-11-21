!------------------------------------------------------------------------------!
!> This subroutine computes the inviscid flux Jacobian for the 1D Euler equations for nonequilibrium flows.
!! In this case the flow is characterized by N temperatures plus a separated temperature Te for free electrons.
  SUBROUTINE inviscid_flux_Jac_neqNT_Te_1D (nx, cons, phys_data, a)

    USE mod_general_data,             ONLY: nb_eq, nb_ns, nb_int_temp, nb_te, nb_temp, pos_em, pos_rhou, pos_rhoE,   & 
                                          & pos_rhoek, pos_rhose, pos_u_cell, pos_h0_cell, pos_c_cell, pos_rho_cell, & 
                                          & pos_alpha_cell, pos_gamma_cell, pos_T_cell, pos_ek_cell, pos_ei_cell,    & 
                                          & pos_eki_cell, gamma_e, gamma_e_m1, ov_gamma_e_m1, Ri, epsi, yi, ei
   
    IMPLICIT NONE

    INTEGER :: i, j, k
    REAL(KIND=8) :: eps, tmp1, tmp2    
    REAL(KIND=8) :: c, ek, fe, gamma, h0vn, h0, p, pe, rho, T, Te, u, uvn, vn, vcn
    REAL(KIND=8) :: alpha, ek_alpha, ge_m_g, ge_m_g_pe_ov_rho, vn_alpha, nx_alpha, nx_ovrho, ov_rho, vn_ovrho
    REAL(KIND=8), DIMENSION(nb_eq) :: prim

    REAL(KIND=8), INTENT(IN) :: nx                       !> normal to the cell interface
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons       !> conservative variables
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data  !> physical properties
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: a       !> inviscid flux Jacobian

    ! Physical data
    T     = phys_data(pos_T_cell)
    Te    = phys_data(pos_T_cell + nb_temp - 1)
    u     = phys_data(pos_u_cell)
    ek    = phys_data(pos_ek_cell)
    h0    = phys_data(pos_h0_cell)
    rho   = phys_data(pos_rho_cell)
    alpha = phys_data(pos_alpha_cell)
    gamma = phys_data(pos_gamma_cell)

    ! Free electron pressure 
    pe = cons(pos_em)*Ri(pos_em)*Te

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
    fe   = ge_m_g*ov_gamma_e_m1*(rho**gamma_e_m1)
    h0vn = h0*vn 

    ! Free electrons
    yi(pos_em)   = cons(pos_em)*vn_ovrho
    epsi(pos_em) = ek_alpha + ge_m_g_pe_ov_rho

    ! Heavy particles 
    DO i = pos_em + 1,nb_ns 
       yi(i) = cons(i)*vn_ovrho
       tmp1  = Ri(i)*T - alpha*ei(i) + ek_alpha + ge_m_g_pe_ov_rho
       tmp2  = 0.d0
       DO k = 1,nb_int_temp
          tmp2 = tmp2 - ei(nb_ns + (i - 1)*nb_int_temp + k)
       ENDDO   
       epsi(i) = tmp1 - alpha*tmp2
    ENDDO

    ! Column j,  j = 1,..,nb_ns
    DO j = 1,nb_ns 

       DO i = 1,j - 1
          a(i,j) = - yi(i)
       ENDDO

       a(j,j) = vn - yi(j)

       DO i = j + 1,nb_ns 
          a(i,j) = - yi(i)
       ENDDO

       eps = epsi(j)
       a(pos_rhou,j) = eps*nx - uvn
       a(pos_rhoE,j) = eps*vn - h0vn

       DO i = 1,nb_int_temp + nb_te
          a(pos_rhoE + i,j) = - cons(pos_rhoE + i)*vn_ovrho
       ENDDO

    ENDDO

    ! Column nb_ns + 1
    DO i = 1,nb_ns 
       a(i,pos_rhou) = cons(i)*nx_ovrho
    ENDDO

    a(pos_rhou,pos_rhou) = 2.d0*vn - vn_alpha
    a(pos_rhoE,pos_rhou) = h0*nx   - alpha*uvn

    DO i = 1,nb_int_temp + nb_te
       a(pos_rhoE + i,pos_rhou) = cons(pos_rhoE + i)*nx_ovrho
    ENDDO

    ! Column nb_ns + 2
    DO i = 1,nb_ns 
       a(i,pos_rhoE) = 0.d0
    ENDDO

    a(pos_rhou,pos_rhoE) = nx_alpha
    a(pos_rhoE,pos_rhoE) = vn + vn_alpha 

    DO i = 1,nb_int_temp + nb_te
       a(pos_rhoE + i,pos_rhoE) = 0.d0
    ENDDO

    ! Column nb_ns + 2 + k,  k = 1,..,nb_int_temp
    DO j = pos_rhoE + 1,pos_rhoE + nb_int_temp

       DO i = 1,nb_ns 
          a(i,j) = 0.d0
       ENDDO

       a(pos_rhou,j) = - nx_alpha
       a(pos_rhoE,j) = - vn_alpha

       DO i = pos_rhoE + 1,j - 1
          a(i,j) = 0.d0
       ENDDO

       a(j,j) = vn

       DO i = j + 1,pos_rhoE + nb_int_temp + nb_te
          a(i,j) = 0.d0
       ENDDO

    ENDDO      

    ! Column nb_ns + 2 + nb_int_temp + nb_te
    DO i = 1,nb_eq
       a(i,pos_rhose) = 0.d0
    ENDDO

    a(pos_rhou,pos_rhose) = nx*fe
    a(pos_rhoE,pos_rhose) = vn*fe

    DO k = 1,nb_int_temp
       a(pos_rhoek + k - 1,pos_rhose) = 0.d0
    ENDDO 

    a(pos_rhose,pos_rhose) = vn

  END SUBROUTINE inviscid_flux_Jac_neqNT_Te_1D
!------------------------------------------------------------------------------!
