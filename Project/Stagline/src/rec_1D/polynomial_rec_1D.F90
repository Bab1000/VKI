!------------------------------------------------------------------------------!
!> This subroutine applies a 1D polynomial reconstruction for 1D calorically perfect gas flows.
  SUBROUTINE poly_rec_1D (ull, ul, ur, urr, prop_ll, prop_l, prop_r, prop_rr, & 
                       &  u_left, u_right, prop_left, prop_right)

    USE mod_general_data,       ONLY: pos_u_cell, pos_pres_cell, pos_rho_cell, nb_eq
    USE mod_function_pointer,   ONLY: limiter, get_cons_phys_from_prim

    IMPLICIT NONE

    REAL(KIND=8), PARAMETER :: eps = 1.d-20
    REAL(KIND=8), PARAMETER :: fac = 0.5d0
    REAL(KIND=8), PARAMETER :: omega_l = fac, omega_r = - fac
    REAL(KIND=8) :: rl_rho, rr_rho, rl_v, rr_v, rl_p, rr_p
    REAL(KIND=8) :: lim_l_rho, lim_r_rho, lim_l_v, lim_r_v, lim_l_p, lim_r_p
    REAL(KIND=8) :: rholl, rhol, rhor, rhorr, vll, vl, vr, vrr, pll, pl, pr, prr
    REAL(KIND=8) :: drho, dv, dp, drhol, drhor, dvl, dvr, dpl, dpr
    REAL(KIND=8), DIMENSION(nb_eq) :: prim_left, prim_right

    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: ull         !< conservative variables of left-left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: ul          !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: ur          !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: urr         !< conservative variables of right-right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: prop_ll     !< physical properties of left-left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: prop_l      !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: prop_r      !< physical properties of right state     
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: prop_rr     !< physical properties of right-right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_left      !< conservative variables of reconstructed left state 
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_right     !< conservative variables of reconstructed right state 
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: prop_left   !< physical properties of reconstructed left state 
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: prop_right  !< physical properties of reconstructed right state

    ! Density, velocity and pressure of all states
    rholl = ull(1)
    rhol  = ul(1)
    rhor  = ur(1)
    rhorr = urr(1)

    vll = prop_ll(pos_u_cell)
    vl  = prop_l(pos_u_cell)
    vr  = prop_r(pos_u_cell)
    vrr = prop_rr(pos_u_cell)

    pll = prop_ll(pos_pres_cell)
    pl  = prop_l(pos_pres_cell)
    pr  = prop_r(pos_pres_cell)
    prr = prop_rr(pos_pres_cell) 

    ! Density, velocity and pressure jumps
    drho  = rhor  - rhol + eps
    drhol = rhol  - rholl + eps
    drhor = rhorr - rhor + eps

    dv  = vr  - vl + eps
    dvl = vl  - vll + eps
    dvr = vrr - vr + eps

    dp  = pr  - pl + eps
    dpl = pl  - pll + eps
    dpr = prr - pr + eps

    ! Ratios of consecutive differences (density)
    rl_rho = drho/drhol
    rr_rho = drho/drhor 

    ! Ratios of consecutive differences (velocity)
    rl_v = dv/dvl
    rr_v = dv/dvr 

    ! Ratios of consecutive differences (pressure)
    rl_p = dp/dpl
    rr_p = dp/dpr

    ! Primitive variable reconstruction
    ! Density
    lim_l_rho = limiter(rl_rho)
    lim_r_rho = limiter(rr_rho)

    prim_left(1)  = rhol + omega_l*lim_l_rho*drhol
    prim_right(1) = rhor + omega_r*lim_r_rho*drhor

    ! Velocity
    lim_l_v = limiter(rl_v)
    lim_r_v = limiter(rr_v)
     
    prim_left(2)  = vl + omega_l*lim_l_v*dvl
    prim_right(2) = vr + omega_r*lim_r_v*dvr
     
    ! Pressure 
    lim_l_p = limiter(rl_p)
    lim_r_p = limiter(rr_p)

    prim_left(3)  = pl + omega_l*lim_l_p*dpl
    prim_right(3) = pr + omega_r*lim_r_p*dpr

    ! Get the conservative variables and the physical properties from the reconstructed primitive variables 
    ! of left and and right states
    CALL get_cons_phys_from_prim(prim_left, u_left, prop_left)
    CALL get_cons_phys_from_prim(prim_right, u_right, prop_right)

  END SUBROUTINE poly_rec_1D 
!------------------------------------------------------------------------------!
