!------------------------------------------------------------------------------!
!> This subroutine computes the numerical flux according to the Steger-Warming flux vector splitting 
!! for 1D nonequilibrium gas flows. In this case the gas flow is characterized by NT temperature 
!! plus a separated temperature Te for free electrons. 
  SUBROUTINE Steger_Warming_neqNT_Te_1D (nx, left_data, right_data, u_left, u_right, f) 

    USE mod_general_data,             ONLY: nb_ns, nb_int_temp, nb_te, nb_temp, nb_eq, pos_u_cell, & 
                                          & pos_rho_cell, pos_c_cell, pos_h0_cell, pos_gamma_cell, & 
                                          & pos_pres_cell, pos_T_cell, pos_em, pos_rhou, pos_rhoE, & 
                                          & pos_rhoek, pos_rhose, gamma_e, Ri
    USE mod_numerics_data,            ONLY: fl, fr
    USE mod_neq_function_pointer,     ONLY: library_get_mass_fractions 


    IMPLICIT NONE

    INTEGER :: i, k 
    REAL(KIND=8) :: gamma, gamma_minus1, ge_m_g 
    REAL(KIND=8) :: c, c2, h0, u, vn, rho, p, pe, Te, ov_c2
    REAL(KIND=8) :: fac1, fac2, fac3, fac4, fac5
    REAL(KIND=8) :: p_fac, pe_fac, eig_diff, eig_sum
    REAL(KIND=8) :: l1, l2, l3, l1p, l2p, l3p, l1m, l2m, l3m

    REAL(KIND=8), INTENT(IN) :: nx                        !< normal to the cell interface
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: f          !< numerical flux

    ! Left state (positive split flux F^+)
    ! Physical data
    c     = left_data(pos_c_cell)
    gamma = left_data(pos_gamma_cell)
    h0    = left_data(pos_h0_cell)
    rho   = left_data(pos_rho_cell)
    u     = left_data(pos_u_cell)
    p     = left_data(pos_pres_cell)
    Te    = left_data(pos_T_cell + nb_temp - 1)

    ! Free electron pressure 
    pe = u_left(pos_em)*Ri(pos_em)*Te

    ! Eigenvalues
    vn = u*nx 
    l1 = vn 
    l2 = vn - c 
    l3 = vn + c

    l1p = 0.5d0*(l1 + ABS(l1))
    l2p = 0.5d0*(l2 + ABS(l2))
    l3p = 0.5d0*(l3 + ABS(l3))

    ! Common factors
    c2 = c**2
    ov_c2    = 0.5d0/c2
    eig_sum  = l3p + l2p
    eig_diff = l3p - l2p 

    ge_m_g = gamma_e - gamma
    gamma_minus1 = gamma - 1.d0

    pe_fac = pe*ge_m_g     
    p_fac  = p + pe_fac

    fac1 = (2.d0*gamma_minus1*l1p + eig_sum)*p
    fac2 = eig_sum*pe_fac
    fac3 = eig_diff*p_fac
    fac4 = ov_c2*(fac1 + fac2)/rho
    fac5 = ov_c2*fac3*c

    ! Species continuity equations
    DO i = 1,nb_ns 
       fl(i) = fac4*u_left(i)
    ENDDO

    ! Global momentum equation
    fl(pos_rhou) = fac4*u_left(pos_rhou) + fac5*nx

    ! Global energy equation
    fl(pos_rhoE) = ov_c2*(2.d0*(h0*gamma_minus1 - c2)*l1p*p + h0*eig_sum*p_fac) + fac5*vn  

    ! Internal energy equations
    DO k = 1,nb_int_temp
       fl(pos_rhoek + k - 1) =  fac4*u_left(pos_rhoek + k - 1)
    ENDDO

    ! Free electron pseudo-entropy equation
    fl(pos_rhose) = fac4*u_left(pos_rhose)

    ! Right state (negative split flux F^-)
    ! Physical data
    c     = right_data(pos_c_cell)
    gamma = right_data(pos_gamma_cell)
    h0    = right_data(pos_h0_cell)
    rho   = right_data(pos_rho_cell)
    u     = right_data(pos_u_cell)
    p     = right_data(pos_pres_cell)
    Te    = right_data(pos_T_cell + nb_temp - 1)

    ! Free electron pressure 
    pe = u_right(pos_em)*Ri(pos_em)*Te

    ! Eigenvalues
    vn = u*nx 
    l1 = vn 
    l2 = vn - c 
    l3 = vn + c

    l1m = 0.5d0*(l1 - ABS(l1))
    l2m = 0.5d0*(l2 - ABS(l2))
    l3m = 0.5d0*(l3 - ABS(l3))

    ! Common factors
    c2 = c**2
    ov_c2    = 0.5d0/c2
    eig_sum  = l3m + l2m
    eig_diff = l3m - l2m 

    ge_m_g = gamma_e - gamma
    gamma_minus1 = gamma - 1.d0
   
    pe_fac = pe*ge_m_g     
    p_fac  = p + pe_fac
 
    fac1 = (2.d0*gamma_minus1*l1m + eig_sum)*p
    fac2 = eig_sum*pe_fac
    fac3 = eig_diff*p_fac
    fac4 = ov_c2*(fac1 + fac2)/rho
    fac5 = ov_c2*fac3*c

    ! Species continuity equations
    DO i = 1,nb_ns 
       fr(i) = fac4*u_right(i)
    ENDDO

    ! Global momentum equation
    fr(pos_rhou) = fac4*u_right(pos_rhou) + fac5*nx

    ! Global energy equation
    fr(pos_rhoE) = ov_c2*(2.d0*(h0*gamma_minus1 - c2)*l1m*p + h0*eig_sum*p_fac) + fac5*vn  

     ! Internal energy equations
    DO k = 1,nb_int_temp
       fr(pos_rhoek + k - 1) =  fac4*u_right(pos_rhoek + k - 1)
    ENDDO

    ! Free electron pseudo-entropy equation
    fr(pos_rhose) = fac4*u_right(pos_rhose)  

    ! Sum of positive and negative split fluxes (F = F^+ + F^-)
    f = fl + fr

  END SUBROUTINE Steger_Warming_neqNT_Te_1D
!------------------------------------------------------------------------------!
