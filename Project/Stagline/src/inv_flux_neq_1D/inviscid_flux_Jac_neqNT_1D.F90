!------------------------------------------------------------------------------!
!> This subroutine computes the inviscid flux Jacobian for the 1D Euler equations for nonequilibrium flows.
!! In this case the flow is characterized by N temperatures. However, no separated temperature exist for the free electrons.
  SUBROUTINE inviscid_flux_Jac_neqNT_1D (nx, cons, phys_data, a)

    USE mod_general_data,             ONLY: nb_eq, nb_ns, nb_int_temp, pos_u_cell, pos_h0_cell, pos_c_cell,     & 
                                          & pos_rho_cell, pos_alpha_cell, pos_T_cell, pos_ek_cell, pos_ei_cell, & 
                                          & pos_eki_cell, pos_rhou, pos_rhoE, Ri, epsi, yi

    IMPLICIT NONE

    INTEGER :: i, j, k
    REAL(KIND=8) :: eps
    REAL(KIND=8) :: tmp1, tmp2
    REAL(KIND=8) :: ek, h0vn, h0, rho, T, u, uvn, vn, vcn
    REAL(KIND=8) :: alpha, ek_alpha, vn_alpha, nx_alpha, nx_ovrho, ov_rho, vn_ovrho 

    REAL(KIND=8), INTENT(IN) :: nx                       !> normal to the cell interface
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons       !> conservative variables
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data  !> physical properties
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: a       !> inviscid flux Jacobian

    ! Physical data
    T     = phys_data(pos_T_cell)
    u     = phys_data(pos_u_cell)
    ek    = phys_data(pos_ek_cell)
    h0    = phys_data(pos_h0_cell)
    rho   = phys_data(pos_rho_cell)
    alpha = phys_data(pos_alpha_cell)

    ! Common factors
    vn       = u*nx
    uvn      = u*vn
    h0vn     = h0*vn
    ov_rho   = 1.d0/rho
    vn_alpha = vn*alpha
    ek_alpha = ek*alpha
    nx_alpha = nx*alpha
    nx_ovrho = nx*ov_rho
    vn_ovrho = vn*ov_rho

    DO i = 1,nb_ns 
       yi(i) = cons(i)*vn_ovrho
       tmp1  = Ri(i)*T - alpha*phys_data(pos_ei_cell + i - 1) + ek_alpha
       tmp2  = 0.d0
       DO k = 1,nb_int_temp
          tmp2 = tmp2 - phys_data(pos_eki_cell + (i - 1)*nb_int_temp + k - 1)
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

       DO i = 1,nb_int_temp
          a(pos_rhoE + i,j) = - cons(pos_rhoE + i)*vn_ovrho
       ENDDO

    ENDDO 

    ! Column nb_ns + 1
    DO i = 1,nb_ns 
       a(i,pos_rhou) = cons(i)*nx_ovrho
    ENDDO

    a(pos_rhou,pos_rhou) = 2.d0*vn - vn_alpha
    a(pos_rhoE,pos_rhou) = h0*nx   - alpha*uvn

    DO i = 1,nb_int_temp
       a(pos_rhoE + i,pos_rhou) = cons(pos_rhoE + i)*nx_ovrho
    ENDDO

    ! Column nb_ns + 2
    DO i = 1,nb_ns 
       a(i,pos_rhoE) = 0.d0
    ENDDO

    a(pos_rhou,pos_rhoE) = nx_alpha
    a(pos_rhoE,pos_rhoE) = vn + vn_alpha 

    DO i = 1,nb_int_temp
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

       DO i = j + 1,pos_rhoE + nb_int_temp
          a(i,j) = 0.d0
       ENDDO

    ENDDO

  END SUBROUTINE inviscid_flux_Jac_neqNT_1D
!------------------------------------------------------------------------------!
