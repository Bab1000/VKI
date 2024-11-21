!------------------------------------------------------------------------------!
!> This subroutine applies a 1D polynomial reconstruction for 1D stagnation line calorically perfect gas flows.
!! The reconstruction is performed on the set of variables V = [rho  u  v  T]^T.
  SUBROUTINE poly_rec_1D_SL (vol_ll, vol_l, vol_r, vol_rr, ull, ul, ur, urr, prop_ll, prop_l, prop_r, prop_rr, & 
                           & u_left, u_right, prop_left, prop_right)

    USE mod_general_data,       ONLY: pos_u_cell, pos_v_cell, pos_pres_cell, pos_rho_cell, pos_T_cell, & 
                                    & pos_rho, pos_u, pos_p, nb_eq, nb_dim, R_gas
    USE mod_function_pointer,   ONLY: limiter, get_cons_phys_from_prim

    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8), PARAMETER :: eps = 1.d-20
    REAL(KIND=8), PARAMETER :: fac = 0.5d0
    REAL(KIND=8), PARAMETER :: omega_l = fac, omega_r = - fac
    REAL(KIND=8) :: rl_rho, rr_rho, rl_v, rr_v, rl_T, rr_T
    REAL(KIND=8) :: lim_l_rho, lim_r_rho, lim_l_v, lim_r_v, lim_l_T, lim_r_T
    REAL(KIND=8) :: rholl, rhol, rhor, rhorr, vll, vl, vr, vrr, Tll, Tl, Tr, Trr
    REAL(KIND=8) :: drho, dv, dp, dT, drhol, drhor, dul, dvl, dvr, dTl, dTr
    REAL(KIND=8), DIMENSION(nb_eq) :: prim_left, prim_right

    REAL(KIND=8),               INTENT(IN)  :: vol_l 
    REAL(KIND=8),               INTENT(IN)  :: vol_r
    REAL(KIND=8),               INTENT(IN)  :: vol_ll
    REAL(KIND=8),               INTENT(IN)  :: vol_rr
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

    ! Density 
    rholl = ull(pos_rho)
    rhol  = ul(pos_rho)
    rhor  = ur(pos_rho)
    rhorr = urr(pos_rho)

    ! Jumps
    drho  = rhor  - rhol + eps
    drhol = rhol  - rholl + eps
    drhor = rhorr - rhor + eps

    ! Ratios of consecutive differences
    rl_rho = drho/drhol
    rr_rho = drho/drhor 

    ! Gradient limiting
    lim_l_rho = limiter(rl_rho)
    lim_r_rho = limiter(rr_rho)

    ! Reconstruction
    prim_left(pos_rho)  = rhol + omega_l*lim_l_rho*drhol
    prim_right(pos_rho) = rhor + omega_r*lim_r_rho*drhor

    ! Radial and circumferential velocity
    DO i = 1,nb_dim

       vll = prop_ll(pos_u_cell + i - 1)
       vl  = prop_l(pos_u_cell + i - 1)
       vr  = prop_r(pos_u_cell + i - 1)
       vrr = prop_rr(pos_u_cell + i - 1)

       ! Jumps
       dv  = vr  - vl + eps
       dvl = vl  - vll + eps
       dvr = vrr - vr + eps  

       ! Ratios of consecutive differences
       rl_v = dv/dvl
       rr_v = dv/dvr  

       ! Gradient limiting
       lim_l_v = limiter(rl_v)
       lim_r_v = limiter(rr_v)

       ! Reconstruction
       prim_left(pos_u + i - 1)  = vl + omega_l*lim_l_v*dvl
       prim_right(pos_u + i - 1) = vr + omega_r*lim_r_v*dvr

    ENDDO

    ! Temperature 
    Tll = prop_ll(pos_T_cell)
    Tl  = prop_l(pos_T_cell)
    Tr  = prop_r(pos_T_cell)
    Trr = prop_rr(pos_T_cell) 

    ! Jumps
    dT  = Tr  - Tl + eps
    dTl = Tl  - Tll + eps
    dTr = Trr - Tr + eps

    ! Ratios of consecutive differences 
    rl_T = dT/dTl
    rr_T = dT/dTr

    ! Gradient limiting
    lim_l_T = limiter(rl_T)
    lim_r_t = limiter(rr_T)     

    ! Reconstruction
    Tl = Tl + omega_l*lim_l_T*dTl
    Tr = Tr + omega_r*lim_r_T*dtr

    ! Reconstructed pressure 
    prim_left(pos_p)  = prim_left(pos_rho)*Tl*R_gas
    prim_right(pos_p) = prim_right(pos_rho)*Tr*R_gas

    ! Get the conservative variables and the physical properties from the reconstructed primitive variables 
    ! of left and and right states
    CALL get_cons_phys_from_prim(prim_left, u_left, prop_left)
    CALL get_cons_phys_from_prim(prim_right, u_right, prop_right)

  END SUBROUTINE poly_rec_1D_SL

!------------------------------------------------------------------------------!
!> This subroutine applies a 1D polynomial reconstruction for 1D stagnation line calorically perfect gas flows.
!! The reconstruction is performed on the set of variables V = [rho  u  v  T]^T, accounting for metrics (A. Turchi)
  SUBROUTINE poly_rec_1D_SL_metr (vol_ll, vol_l, vol_r, vol_rr, ull, ul, ur, urr, prop_ll, prop_l, prop_r, prop_rr, & 
                           & u_left, u_right, prop_left, prop_right)

    USE mod_general_data,       ONLY: pos_u_cell, pos_v_cell, pos_pres_cell, pos_rho_cell, pos_T_cell, & 
                                    & pos_rho, pos_u, pos_p, nb_eq, nb_dim, R_gas
    USE mod_function_pointer,   ONLY: limiter, get_cons_phys_from_prim

    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8), PARAMETER :: eps = 1.d-20
    REAL(KIND=8) :: rl_rho, rr_rho, rl_v, rr_v, rl_T, rr_T
    REAL(KIND=8) :: lim_l_rho, lim_r_rho, lim_l_v, lim_r_v, lim_l_T, lim_r_T
    REAL(KIND=8) :: rholl, rhol, rhor, rhorr, vll, vl, vr, vrr, Tll, Tl, Tr, Trr
    REAL(KIND=8) :: drho, dv, dp, dT, drhol, drhor, dul, dvl, dvr, dTl, dTr
    REAL(KIND=8) :: vol_sum_l, vol_sum, vol_sum_r 
    REAL(KIND=8) :: omega_l 
    REAL(KIND=8) :: omega_r 
    REAL(KIND=8), DIMENSION(nb_eq) :: prim_left, prim_right

    REAL(KIND=8),               INTENT(IN)  :: vol_l 
    REAL(KIND=8),               INTENT(IN)  :: vol_r
    REAL(KIND=8),               INTENT(IN)  :: vol_ll
    REAL(KIND=8),               INTENT(IN)  :: vol_rr
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

    vol_sum_l = vol_ll + vol_l
    vol_sum   = vol_l + vol_r
    vol_sum_r = vol_r + vol_rr

    omega_l =   vol_l / vol_sum_l 
    omega_r = - vol_r / vol_sum_r


    ! Density 
    rholl = ull(pos_rho)
    rhol  = ul(pos_rho)
    rhor  = ur(pos_rho)
    rhorr = urr(pos_rho)

    ! Jumps
    drho  = rhor  - rhol + eps
    drhol = rhol  - rholl + eps
    drhor = rhorr - rhor + eps
           
    ! Ratio of consecutive differences (species densities)
    rl_rho = (drho * vol_sum_l )  / (drhol * vol_sum)
    rr_rho = (drho * vol_sum_r )  / (drhor * vol_sum)

    ! Gradient limiting
    lim_l_rho = limiter(rl_rho)
    lim_r_rho = limiter(rr_rho)

    ! Reconstruction
    prim_left(pos_rho)  = rhol + omega_l*lim_l_rho*drhol
    prim_right(pos_rho) = rhor + omega_r*lim_r_rho*drhor

    ! Radial and circumferential velocity
    DO i = 1,nb_dim

       vll = prop_ll(pos_u_cell + i - 1)
       vl  = prop_l(pos_u_cell + i - 1)
       vr  = prop_r(pos_u_cell + i - 1)
       vrr = prop_rr(pos_u_cell + i - 1)

       ! Jumps
       dv  = vr  - vl + eps
       dvl = vl  - vll + eps
       dvr = vrr - vr + eps  

       ! Ratios of consecutive differences
       rl_v = (dv * vol_sum_l) / (dvl * vol_sum)
       rr_v = (dv * vol_sum_r) / (dvr * vol_sum) 

       ! Gradient limiting
       lim_l_v = limiter(rl_v)
       lim_r_v = limiter(rr_v)

       ! Reconstruction
       prim_left(pos_u + i - 1)  = vl + omega_l*lim_l_v*dvl
       prim_right(pos_u + i - 1) = vr + omega_r*lim_r_v*dvr

    ENDDO

    ! Temperature 
    Tll = prop_ll(pos_T_cell)
    Tl  = prop_l(pos_T_cell)
    Tr  = prop_r(pos_T_cell)
    Trr = prop_rr(pos_T_cell) 

    ! Jumps
    dT  = Tr  - Tl + eps
    dTl = Tl  - Tll + eps
    dTr = Trr - Tr + eps

    ! Ratios of consecutive differences 
    rl_T = (dT * vol_sum_l )  / (dTl * vol_sum)
    rr_T = (dT * vol_sum_r )  / (dTr * vol_sum)

    ! Gradient limiting
    lim_l_T = limiter(rl_T)
    lim_r_t = limiter(rr_T)     

    ! Reconstruction
    Tl = Tl + omega_l*lim_l_T*dTl
    Tr = Tr + omega_r*lim_r_T*dtr

    ! Reconstructed pressure 
    prim_left(pos_p)  = prim_left(pos_rho)*Tl*R_gas
    prim_right(pos_p) = prim_right(pos_rho)*Tr*R_gas

    ! Get the conservative variables and the physical properties from the reconstructed primitive variables 
    ! of left and and right states
    CALL get_cons_phys_from_prim(prim_left, u_left, prop_left)
    CALL get_cons_phys_from_prim(prim_right, u_right, prop_right)

  END SUBROUTINE poly_rec_1D_SL_metr
