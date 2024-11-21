!------------------------------------------------------------------------------!
!> This subroutine computes the convective flux by means of the AUSM scheme 
!! for 1D stagnation line chemically non equilibrium flows.
  SUBROUTINE ausmP_up_as_neq_1D_SL (nx, vol_l, vol_r, left_data, right_data, u_left, u_right, f)
    USE mod_general_data,                 ONLY: nb_eq, nb_ns, pos_u_cell, pos_v_cell, pos_c_cell, pos_pres_cell, & 
                                              & pos_h0_cell, pos_mu_cell, pos_rho, pos_rhou, pos_rhov, pos_rho_cell,   & 
                                              & pos_rhoE, pos_gamma_cell,MachIn

    IMPLICIT NONE

    REAL(KIND=8) :: rhol, ul, vl, pl, cl, h0l, rhor, ur, vr, pr, cr, h0r, unl, unr, gammal, gammar
    REAL(KIND=8) :: a_interface, gamma_interface
    REAL(KIND=8) :: Mnl, Mnr, Mnl_abs, Mnr_abs, M_bar_2, M_o_2, M_o, f_a, M_inf_2
    REAL(KIND=8) :: m_dot_interface, M_interface, MpL, MmR, Mu, P_interface, PpL, PmR, Pu
    REAL(KIND=8) :: alpha_coeff, beta_coeff, sigma_coeff, K_p, K_u
    REAL(KIND=8), DIMENSION(nb_eq) :: Psil, Psir
    INTEGER      :: i_ns

    REAL(KIND=8) :: Ma_function_minus_4, Ma_function_plus_4, P_function_minus_5, P_function_plus_5
    REAL(KIND=8) :: Ma_p_diffusion, P_u_diffusion, a_interface_function

    REAL(KIND=8), INTENT(IN) :: nx                        !< normal to the cell interface
    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: f          !< numerical flux

    ! Data of left and right states
    rhol   = left_data(pos_rho_cell)
    ul     = left_data(pos_u_cell)
    vl     = left_data(pos_v_cell)
    pl     = left_data(pos_pres_cell)
    cl     = left_data(pos_c_cell)
    h0l    = left_data(pos_h0_cell)
    gammal = left_data(pos_gamma_cell)

    rhor   = right_data(pos_rho_cell)
    ur     = right_data(pos_u_cell)
    vr     = right_data(pos_v_cell)
    pr     = right_data(pos_pres_cell)
    cr     = right_data(pos_c_cell)
    h0r    = right_data(pos_h0_cell)
    gammar = right_data(pos_gamma_cell)

    ! GET M_inf_2
    M_inf_2 = MachIn
    !M_inf_2 = 1.d0

    unl = ul*nx
    unr = ur*nx

    ! Parameters for AUSM+up
    beta_coeff  = 0.125d0 
    sigma_coeff = 1.d0
    K_p         = 0.25d0
    K_u         = 0.75d0

    gamma_interface = 0.5D0 * (gammal + gammar)

    ! Normal Mach numbers of left and right states
    a_interface = a_interface_function(unl, unr, h0l, h0r, gamma_interface)
    !a_interface = (cr + cl)/2.d0

    Mnl = unl / a_interface
    Mnr = unr / a_interface

    Mnl_abs = ABS(Mnl)
    Mnr_abs = ABS(Mnr)

    ! Defining the f_a function
    M_bar_2 = (ul * ul + ur * ur) / (2.d0 * a_interface * a_interface)

    M_o_2   = MIN( 1.d0, MAX( M_bar_2, M_inf_2 ) )
    M_o     = SQRT(M_o_2)
    f_a     = M_o * ( 2 - M_o )
    alpha_coeff = 0.1875d0 * ( -4.d0 + 5.d0 * f_a * f_a)

    MpL = Ma_function_plus_4(Mnl, beta_coeff)
    PpL = P_function_plus_5(Mnl, alpha_coeff)
    MmR = Ma_function_minus_4(MnR, beta_coeff)
    PmR = P_function_minus_5(MnR, alpha_coeff)

    ! Computing the pressure diffusion term
    !Mu =  Ma_p_diffusion(pl, pr, rhol, rhor, M_bar_2, a_interface, K_p,sigma_coeff) / f_a
    Mu = 0.d0

    ! Computing the velocity diffusion term
    Pu = P_u_diffusion( rhol, rhor, unl, unr, MnL, MnR, a_interface, K_u, alpha_coeff ) * f_a

    ! Computing the interface Mach
    M_interface = MpL + MmR + Mu

    ! Computing the m_dot_interface
    m_dot_interface = a_interface * M_interface
    IF (M_interface > 0) THEN
        m_dot_interface = m_dot_interface * rhol
    ELSE
        m_dot_interface = m_dot_interface * rhor
    ENDIF

    ! Computing the P_interface
    P_interface = PpL * pl + PmR * pr + Pu

    ! Computing Psil and Psir
    DO i_ns = 1, nb_ns
      Psil(i_ns) = u_left(i_ns) / rhol
    ENDDO
    Psil(pos_rhou) = ul 
    Psil(pos_rhov) = vl
    Psil(pos_rhoE) = h0l

    DO i_ns = 1, nb_ns
      Psir(i_ns) = u_right(i_ns) / rhor
    ENDDO
    Psir(pos_rhou) = ur 
    Psir(pos_rhov) = vr
    Psir(pos_rhoE) = h0r

    ! Numerical flux 
    f = 0.5d0* m_dot_interface * (Psil + Psir) - 0.5d0*ABS(m_dot_interface)*(Psir - Psil)

    ! Adding the pressure contribution in the radial momentum equation
    f(pos_rhou) = f(pos_rhou) + P_interface*nx

  END SUBROUTINE ausmP_up_as_neq_1D_SL

  SUBROUTINE ausmP_up_as_cd_neq_1D_SL (nx, vol_l, vol_r, left_data, right_data, u_left, u_right, f)
    USE mod_general_data,                 ONLY: nb_eq, nb_ns, pos_u_cell, pos_v_cell, pos_c_cell, pos_pres_cell, & 
                                              & pos_h0_cell, pos_mu_cell, pos_rho, pos_rhou, pos_rhov, pos_rho_cell,   & 
                                              & pos_rhoE, pos_gamma_cell

    IMPLICIT NONE

    REAL(KIND=8) :: rhol, ul, vl, pl, cl, h0l, rhor, ur, vr, pr, cr, h0r, unl, unr, gammal, gammar
    REAL(KIND=8) :: a_interface, gamma_interface
    REAL(KIND=8) :: Mnl, Mnr, Mnl_abs, Mnr_abs, M_bar_2, M_o_2, M_o, f_a, M_inf_2
    REAL(KIND=8) :: m_dot_interface, M_interface, MpL, MmR, Mu, P_interface, PpL, PmR, Pu
    REAL(KIND=8) :: alpha_coeff, beta_coeff, sigma_coeff, K_p, K_u
    REAL(KIND=8), DIMENSION(nb_eq) :: Psil, Psir
    INTEGER      :: i_ns

    REAL(KIND=8) :: c, dx, dp, mp, m_visc, p, rho, u, Vp, mul, mur, u_co

    REAL(KIND=8) :: Ma_function_minus_4, Ma_function_plus_4, P_function_minus_5, P_function_plus_5
    REAL(KIND=8) :: Ma_p_diffusion, P_u_diffusion, a_interface_function
    REAL(KIND=8) :: viscous_flow_preconditioning_mass_flow


    REAL(KIND=8), INTENT(IN) :: nx                        !< normal to the cell interface
    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: f          !< numerical flux

    ! Data of left and right states
    rhol   = left_data(pos_rho_cell)
    ul     = left_data(pos_u_cell)
    vl     = left_data(pos_v_cell)
    pl     = left_data(pos_pres_cell)
    cl     = left_data(pos_c_cell)
    h0l    = left_data(pos_h0_cell)
    gammal = left_data(pos_gamma_cell)
    mul    = left_data(pos_mu_cell)

    rhor   = right_data(pos_rho_cell)
    ur     = right_data(pos_u_cell)
    vr     = right_data(pos_v_cell)
    pr     = right_data(pos_pres_cell)
    cr     = right_data(pos_c_cell)
    h0r    = right_data(pos_h0_cell)
    gammar = right_data(pos_gamma_cell)
    mur  = right_data(pos_mu_cell)

    ! GET M_inf_2
    !M_inf_2 = 0.000000101102d0
    M_inf_2 = 1.d-2

    unl = ul*nx
    unr = ur*nx

    ! Parameters for AUSM+up
    beta_coeff  = 0.125d0 
    sigma_coeff = 1.d0
    K_p         = 0.25d0
    K_u         = 0.75d0

    gamma_interface = 0.5D0 * (gammal + gammar)

    ! Normal Mach numbers of left and right states
    a_interface = a_interface_function(unl, unr, h0l, h0r, gamma_interface)
    !a_interface = (cr + cl)/2.d0

    Mnl = unl / a_interface
    Mnr = unr / a_interface

    Mnl_abs = ABS(Mnl)
    Mnr_abs = ABS(Mnr)

    ! Defining the f_a function
    M_bar_2 = (ul * ul + ur * ur) / (2.d0 * a_interface * a_interface)

    M_o_2   = MIN( 1.d0, MAX( M_bar_2, M_inf_2 ) )
    M_o     = SQRT(M_o_2)
    f_a     = M_o * ( 2 - M_o )

    alpha_coeff = 0.1875d0 * ( -4.d0 + 5.d0 * f_a * f_a)

    MpL = Ma_function_plus_4(Mnl, beta_coeff)
    PpL = P_function_plus_5(Mnl, alpha_coeff)
    MmR = Ma_function_minus_4(MnR, beta_coeff)
    PmR = P_function_minus_5(MnR, alpha_coeff)

    ! Computing the pressure diffusion term
    Mu =  Ma_p_diffusion(pl, pr, rhol, rhor, M_bar_2, a_interface, K_p,sigma_coeff) ! / f_a 

    ! Computing the velocity diffusion term
    Pu = P_u_diffusion( rhol, rhor, unl, unr, MnL, MnR, a_interface, K_u, alpha_coeff ) * f_a  

    ! Computing the interface Mach
    M_interface = MpL + MmR + Mu

    ! Computing the m_dot_interface
    m_dot_interface = a_interface * M_interface
    IF (M_interface > 0) THEN
        m_dot_interface = m_dot_interface * rhol
    ELSE
        m_dot_interface = m_dot_interface * rhor
    ENDIF

    ! Computing the P_interface
    P_interface = PpL * pl + PmR * pr + Pu

    ! Computing Psil and Psir
    DO i_ns = 1, nb_ns
      Psil(i_ns) = u_left(i_ns) / rhol
    ENDDO
    Psil(pos_rhou) = ul 
    Psil(pos_rhov) = vl
    Psil(pos_rhoE) = h0l

    DO i_ns = 1, nb_ns
      Psir(i_ns) = u_right(i_ns) / rhor
    ENDDO
    Psir(pos_rhou) = ur 
    Psir(pos_rhov) = vr
    Psir(pos_rhoE) = h0r

    ! Compute pre-conditioning velocity
    rho = 0.5d0*(rhol + rhor)
    m_visc = 0.5d0*(mul + mur)
    dx = 0.5d0*(vol_l + vol_r)
    u  = 0.5d0*(ul + ur)
    u_co = 1.d0
  
    ! Interface pressure diffusion term (added for pressure velocity coupling at low speeds) 
    mp = 0.d0!viscous_flow_preconditioning_mass_flow(pl, pr, u, u_co, m_visc, rho , dx) * nx
 
    ! Numerical flux 
    f = 0.5d0 * ( m_dot_interface + mp ) * (Psil + Psir) - 0.5d0*ABS(m_dot_interface )*(Psir - Psil)

    ! Adding the pressure contribution in the radial momentum equation
    f(pos_rhou) = f(pos_rhou) + P_interface*nx

  END SUBROUTINE ausmP_up_as_cd_neq_1D_SL

!------------------------------------------------------------------------------!
! Implementation of the necessary functions for computing the ausm schemes
!------------------------------------------------------------------------------!

FUNCTION Ma_function_plus_1 (Ma)
IMPLICIT NONE
REAL(KIND=8), INTENT(IN) :: Ma
REAL(KIND=8)             :: Ma_function_plus_1

Ma_function_plus_1 = 0.5D0 * (Ma + abs(Ma))
END FUNCTION Ma_function_plus_1

!------------------------------------------------------------------------------!

FUNCTION Ma_function_minus_1 (Ma)
IMPLICIT NONE
REAL(KIND=8), INTENT(IN) :: Ma
REAL(KIND=8)             :: Ma_function_minus_1

Ma_function_minus_1 = 0.5D0 * (Ma - abs(Ma))
END FUNCTION Ma_function_minus_1

!------------------------------------------------------------------------------!

FUNCTION Ma_function_plus_2 (Ma)
IMPLICIT NONE
REAL(KIND=8), INTENT(IN) :: Ma
REAL(KIND=8)             :: Ma_function_plus_2

Ma_function_plus_2 = 0.25D0 * (Ma + 1.D0)**2.D0
END FUNCTION Ma_function_plus_2

!------------------------------------------------------------------------------!

FUNCTION Ma_function_minus_2 (Ma)
IMPLICIT NONE
REAL(KIND=8), INTENT(IN) :: Ma
REAL(KIND=8)             :: Ma_function_minus_2

Ma_function_minus_2 = -0.25D0 * (Ma - 1.D0)**2.D0
END FUNCTION Ma_function_minus_2

!------------------------------------------------------------------------------!

FUNCTION Ma_function_plus_4 (Ma, beta_coeff)
IMPLICIT NONE
REAL(KIND=8), INTENT(IN) :: Ma
REAL(KIND=8), INTENT(IN) :: beta_coeff
REAL(KIND=8)             :: Ma_function_plus_1, Ma_function_plus_2, Ma_function_minus_2
REAL(KIND=8)             :: Ma_function_plus_4

IF ( abs(Ma) >= 1.D0 ) THEN
    Ma_function_plus_4 = Ma_function_plus_1(Ma)
ELSE 
    Ma_function_plus_4 = Ma_function_plus_2(Ma) * (1.D0 - 16.D0 * beta_coeff * Ma_function_minus_2(Ma))
END IF
END FUNCTION Ma_function_plus_4

!------------------------------------------------------------------------------!

FUNCTION Ma_function_minus_4 (Ma, beta_coeff)
IMPLICIT NONE
REAL(KIND=8), INTENT(IN) :: Ma
REAL(KIND=8), INTENT(IN) :: beta_coeff
REAL(KIND=8)             :: Ma_function_minus_1, Ma_function_plus_2, Ma_function_minus_2
REAL(KIND=8)             :: Ma_function_minus_4

IF ( abs(Ma) >= 1.D0 ) THEN
    Ma_function_minus_4 = Ma_function_minus_1(Ma)
ELSE 
    Ma_function_minus_4 = Ma_function_minus_2(Ma) * (1.D0 + 16.D0 * beta_coeff * Ma_function_plus_2(Ma))
END IF
END FUNCTION Ma_function_minus_4

!------------------------------------------------------------------------------!

FUNCTION Ma_p_diffusion (pL, pR, rhoL, rhoR, M_bar_2, alpha_interface, K_p, sigma_coeff)
IMPLICIT NONE
REAL(KIND=8), INTENT(IN) :: pL, pR, rhoL, rhoR
REAL(KIND=8), INTENT(IN) :: M_bar_2, alpha_interface
REAL(KIND=8), INTENT(IN) :: K_p, sigma_coeff
REAL(KIND=8)             :: Ma_p_diffusion

Ma_p_diffusion = - K_p * MAX(1.D0 - sigma_coeff * M_bar_2, 0.D0) * (pR - pL) / ((rhoR + rhoL) / 2.D0 * alpha_interface**2.D0 )

END FUNCTION Ma_p_diffusion

!------------------------------------------------------------------------------!

FUNCTION P_function_plus_5 (Ma, alpha_coeff)
IMPLICIT NONE
REAL(KIND=8), INTENT(IN) :: Ma
REAL(KIND=8), INTENT(IN) :: alpha_coeff
REAL(KIND=8)             :: Ma_function_plus_1, Ma_function_plus_2, Ma_function_minus_2
REAL(KIND=8)             :: P_function_plus_5

IF ( abs(Ma) >= 1.D0 ) THEN
    P_function_plus_5 = Ma_function_plus_1(Ma) / Ma
ELSE 
    P_function_plus_5 = Ma_function_plus_2(Ma) * ((2.D0 - Ma ) - 16.D0 * alpha_coeff * Ma * Ma_function_minus_2(Ma) )
END IF
END FUNCTION P_function_plus_5

!------------------------------------------------------------------------------!

FUNCTION P_function_minus_5 (Ma, alpha_coeff)
IMPLICIT NONE
REAL(KIND=8), INTENT(IN) :: Ma
REAL(KIND=8), INTENT(IN) :: alpha_coeff
REAL(KIND=8)             :: Ma_function_minus_1, Ma_function_plus_2, Ma_function_minus_2
REAL(KIND=8)             :: P_function_minus_5

IF ( abs(Ma) >= 1.D0 ) THEN
    P_function_minus_5 = Ma_function_minus_1(Ma) / Ma
ELSE 
    P_function_minus_5 = Ma_function_minus_2(Ma) * (( -2.D0 - Ma ) + 16.D0 * alpha_coeff * Ma * Ma_function_plus_2(Ma) )
END IF
END FUNCTION P_function_minus_5

!------------------------------------------------------------------------------!

FUNCTION P_u_diffusion ( rhoL, rhoR, uL, uR, MaL, MaR, alpha_interface, K_u, alpha_coeff )
IMPLICIT NONE
REAL(KIND=8), INTENT(IN) :: rhoL, rhoR, uL, uR, MaL, MaR
REAL(KIND=8), INTENT(IN) :: alpha_interface
REAL(KIND=8), INTENT(IN) :: K_u, alpha_coeff
REAL(KIND=8)             :: P_function_plus_5, P_function_minus_5
REAL(KIND=8)             :: P_u_diffusion

P_u_diffusion = - K_u * P_function_plus_5(MaL, alpha_coeff) * P_function_minus_5(MaR, alpha_coeff) &
                                                                           * (rhoL + rhoR) * alpha_interface * (uR - uL)

END FUNCTION P_u_diffusion

!------------------------------------------------------------------------------!

FUNCTION a_interface_function(unl, unr, h0l, h0r, gamma_interface)
IMPLICIT NONE
REAL(KIND=8), INTENT(IN) :: unl, unr, h0l, h0r
REAL(KIND=8), INTENT(IN) :: gamma_interface
REAL(KIND=8)             :: a_star_2, a_star, a_hat_L, a_hat_R
REAL(KIND=8)             :: a_interface_function

a_star_2 = 2.d0 * (gamma_interface - 1.d0) / (gamma_interface + 1.d0) * (h0l + h0r) / 2.d0
a_star = SQRT(a_star_2)
a_hat_L = a_star_2 / MAX(a_star, unl) ! Check sign
a_hat_R = a_star_2 / MAX(a_star, -unr)

a_interface_function = MIN(a_hat_L, a_hat_R)

END FUNCTION a_interface_function

!------------------------------------------------------------------------------!

FUNCTION weighting_function(Ma_interface, Ma_limit)
IMPLICIT NONE
REAL(KIND=8), INTENT(IN) :: Ma_interface
REAL(KIND=8), INTENT(IN) :: Ma_limit
REAL(KIND=8)             :: Ma_ratio
REAL(KIND=8)             :: weighting_function

Ma_ratio = Ma_interface / Ma_limit

IF( Ma_ratio < 1.d0 ) THEN
    weighting_function = -6.d0 * Ma_ratio**5.d0 + 15.d0 * Ma_ratio**4.d0 - 10.d0 * Ma_ratio**3.d0 + 1.d0
ELSE
    weighting_function = 0.d0
ENDIF

END FUNCTION weighting_function

!------------------------------------------------------------------------------!

FUNCTION viscous_flow_preconditioning_mass_flow(p_L, p_R, u_local, u_cutoff, mu_interface, rho_interface, l_ref)
IMPLICIT NONE
REAL(KIND=8), INTENT(IN) :: p_L, p_R
REAL(KIND=8), INTENT(IN) :: u_local, u_cutoff
REAL(KIND=8), INTENT(IN) :: mu_interface, rho_interface, l_ref
REAL(KIND=8)             :: u_viscous, inv_velocity
REAL(KIND=8)             :: viscous_flow_preconditioning_mass_flow

u_viscous = mu_interface / (rho_interface * l_ref)

inv_velocity = 1.D0 / MAX( ABS(u_local), u_viscous, u_cutoff )

viscous_flow_preconditioning_mass_flow = - (p_R - p_L) * inv_velocity

END FUNCTION
