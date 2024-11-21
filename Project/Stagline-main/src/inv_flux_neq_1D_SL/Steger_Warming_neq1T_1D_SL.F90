!------------------------------------------------------------------------------!
!> This subroutine computes the numerical flux according to the Steger-Warming flux vector splitting 
!! for 1D stagnation line nonequilibrium gas flows. In this case the flow is described by means of a single temperature.  
  SUBROUTINE Steger_Warming_neq1T_1D_SL (nx, vol_l, vol_r, left_data, right_data, u_left, u_right, f) 

    USE mod_general_data,             ONLY: nb_ns, pos_u_cell, pos_v_cell,  pos_rho_cell, pos_c_cell, pos_h0_cell, & 
                                          & pos_gamma_cell, pos_rhou, pos_rhov, pos_rhoE, yi, nb_eq
    USE mod_numerics_data,            ONLY: fl, fr
    USE mod_neq_function_pointer,     ONLY: library_get_mass_fractions 

    IMPLICIT NONE

    INTEGER :: i 
    REAL(KIND=8) :: gamma, gamma_minus1 
    REAL(KIND=8) :: c, h0, u, v, vn, rho
    REAL(KIND=8) :: fac, fac_l, fac_r, eig_diff, eig_sum
    REAL(KIND=8) :: l1, l2, l3, l1p, l2p, l3p, l1m, l2m, l3m

    REAL(KIND=8), INTENT(IN) :: nx                        !< normal to the cell interface
    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
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
    v     = left_data(pos_v_cell)

    ! Species mass fractions
    CALL library_get_mass_fractions (rho, u_left(1:nb_ns), yi)

    ! Eigenvalues
    vn = u*nx 
    l1 = vn 
    l2 = vn - c 
    l3 = vn + c

    l1p = 0.5d0*(l1 + ABS(l1))
    l2p = 0.5d0*(l2 + ABS(l2))
    l3p = 0.5d0*(l3 + ABS(l3))

    ! Common factors
    eig_sum  = l3p + l2p
    eig_diff = l3p - l2p 

    gamma_minus1 = gamma - 1.d0

    fac   = 2.d0*gamma_minus1*l1p + eig_sum
    fac_l = rho/gamma

    ! Species continuity equations
    DO i = 1,nb_ns 
       fl(i) = yi(i)*fac
    ENDDO

    ! Radial momentum equation
    fl(pos_rhou) = u*fac + eig_diff*c*nx

    ! Circumferential momentum equation
    fl(pos_rhov) = v*fac

    ! Global energy equation
    fl(pos_rhoE) = 2.d0*l1p*(h0*gamma_minus1 - c**2) + eig_sum*h0 + eig_diff*vn*c

    ! Right state (negative split flux F^-)
    ! Physical data
    c     = right_data(pos_c_cell)
    gamma = right_data(pos_gamma_cell)
    h0    = right_data(pos_h0_cell)
    rho   = right_data(pos_rho_cell)
    u     = right_data(pos_u_cell)
    v     = right_data(pos_v_cell)

    ! Species mass fractions
    CALL library_get_mass_fractions (rho, u_right(1:nb_ns), yi)

    ! Eigenvalues
    vn = u*nx 
    l1 = vn 
    l2 = vn - c 
    l3 = vn + c

    l1m = 0.5d0*(l1 - ABS(l1))
    l2m = 0.5d0*(l2 - ABS(l2))
    l3m = 0.5d0*(l3 - ABS(l3))

    ! Common factors
    eig_sum  = l3m + l2m
    eig_diff = l3m - l2m 

    gamma_minus1 = gamma - 1.d0

    fac   = 2.d0*gamma_minus1*l1m + eig_sum
    fac_r = rho/gamma

    ! Species continuity equations
    DO i = 1,nb_ns 
       fr(i) = yi(i)*fac
    ENDDO

    ! Radial momentum equation
    fr(pos_rhou) = u*fac + eig_diff*c*nx

    ! Circumreferential momentum equation 
    fr(pos_rhov) = v*fac

    ! Global energy equation
    fr(pos_rhoE) = 2.d0*l1m*(h0*gamma_minus1 - c**2) + eig_sum*h0 + eig_diff*vn*c

    ! Sum of positive and negative split fluxes (F = F^+ + F^-)
    f = 0.5d0*(fac_l*fl + fac_r*fr)
   
  END SUBROUTINE Steger_Warming_neq1T_1D_SL 
!------------------------------------------------------------------------------!
