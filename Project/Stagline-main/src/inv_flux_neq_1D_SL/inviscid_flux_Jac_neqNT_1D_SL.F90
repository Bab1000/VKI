!------------------------------------------------------------------------------!
!> This subroutine computes the inviscid flux Jacobian for the 1D Euler stagnation line equations for nonequilibrium flows.
!! In this case the flow is characterized by N temperatures. However, no separated temperature exist for the free electrons.
  SUBROUTINE inviscid_flux_Jac_neqNT_1D_SL (nx, cons, phys_data, a)

    USE mod_general_data,             ONLY: nb_eq, nb_ns, pos_u_cell,pos_v_cell, pos_h0_cell, pos_c_cell,  &
                                          & pos_rho_cell, pos_alpha_cell, pos_T_cell, pos_ek_cell,         &
                                          & pos_ei_cell, pos_rhou, pos_rhov, pos_rhoE, Ri, epsi, yi,       &
                                          & pos_rhoek, pos_eki_cell, nb_int_temp, eintk
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(KIND=8) :: eps, tmp1
    REAL(KIND=8) :: ek, h0vn, h0, rho, T, u, v, uvn, vn, vvn, vcn
    REAL(KIND=8) :: alpha, ek_alpha, vn_alpha, nx_alpha, nx_ovrho, ov_rho,vn_ovrho
 
    REAL(KIND=8), INTENT(IN) :: nx                       !> normal to the cell interface
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons       !> conservative variables
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data  !> physical properties
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: a       !> inviscid flux Jacobian

    a = 0.d0

    ! Physical data
    T     = phys_data(pos_T_cell)
    u     = phys_data(pos_u_cell)
    v     = phys_data(pos_v_cell)
    ek    = phys_data(pos_ek_cell)
    h0    = phys_data(pos_h0_cell)
    rho   = phys_data(pos_rho_cell)
    alpha = phys_data(pos_alpha_cell)

    ! Common factors
    vn       = u*nx
    uvn      = u*vn
    vvn      = v*vn
    h0vn     = h0*vn
    ov_rho   = 1.d0/rho
    vn_alpha = vn*alpha
    ek_alpha = ek*alpha
    nx_alpha = nx*alpha
    nx_ovrho = nx*ov_rho
    vn_ovrho = vn*ov_rho

    ! Remarks on NT model
    ! alpha = (gamma-1) where cv is related to the translational modes
    ! phys_data(pos_ei_cell) contains the energy density (total witout kinetic)
    ! of species i 
    DO i = 1,nb_ns
       yi(i)   = cons(i)*vn_ovrho
       epsi(i) = Ri(i)*T - alpha*phys_data(pos_ei_cell + i - 1) + ek_alpha
    ENDDO

    DO i = 1, nb_int_temp
      eintk(i) = 0.d0
      DO j = 1, nb_ns
        tmp1 = phys_data(pos_eki_cell + j - 1 +nb_ns*(i-1))
        epsi(j) = epsi(j) - alpha*tmp1
        eintk(i) = eintk(i) + tmp1
      ENDDO
    ENDDO 
    ! Inviscid flux Jacobian 
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
       a(pos_rhov,j) = - vvn
       a(pos_rhoE,j) = eps*vn - h0vn
       DO i = 1, nb_int_temp
          a(pos_rhoE+i, j) = -vn*eintk(i)
       ENDDO

    ENDDO

    ! Column nb_ns + 1
    DO i = 1,nb_ns
       a(i,pos_rhou) = cons(i)*nx_ovrho
    ENDDO

    a(pos_rhou,pos_rhou) = 2.d0*vn - vn_alpha
    a(pos_rhov,pos_rhou) = v*nx
    a(pos_rhoE,pos_rhou) = h0*nx   - alpha*uvn
    DO i = 1, nb_int_temp
       a(pos_rhoE+i, pos_rhou) = eintk(i)
    ENDDO

    ! Column nb_ns + 2
    DO i = 1,nb_ns
       a(i,pos_rhov) = 0.d0
    ENDDO

    a(pos_rhou,pos_rhov) = 0.d0
    a(pos_rhov,pos_rhov) = vn
    a(pos_rhoE,pos_rhov) = 0.d0
    DO i = 1, nb_int_temp
       a(pos_rhoE+i, pos_rhov) = 0.d0
    ENDDO

    ! Column nb_ns + 3
    DO i = 1,nb_ns
       a(i,pos_rhoE) = 0.d0
    ENDDO

    a(pos_rhou,pos_rhoE) = nx_alpha
    a(pos_rhov,pos_rhoE) = 0.d0
    a(pos_rhoE,pos_rhoE) = vn + vn_alpha


    DO i = 1, nb_int_temp
       a(pos_rhoE+i, pos_rhoE) = 0.d0
    ENDDO

    ! Column nb_ns+4, ..., nb_eq
    DO j = 1, nb_int_temp

      DO i = 1, nb_ns
         a(i, pos_rhoE+j) = 0.d0
      ENDDO

      a(pos_rhou,pos_rhoE+j) = - nx_alpha
      a(pos_rhov,pos_rhoE+j) = 0.d0
      a(pos_rhoE,pos_rhoE+j) = - vn_alpha

      DO i = 1, j-1
         a(pos_rhoE+i, pos_rhoE+j) = 0.d0
      ENDDO
      a(pos_rhoE+j, pos_rhoE+j) = vn
      DO i = j+1, nb_int_temp
         a(pos_rhoE+i, pos_rhoE+j) = 0.d0
      ENDDO
    ENDDO

  END SUBROUTINE inviscid_flux_Jac_neqNT_1D_SL
!------------------------------------------------------------------------------!
