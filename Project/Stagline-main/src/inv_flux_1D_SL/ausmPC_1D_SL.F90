!------------------------------------------------------------------------------!
!> This subroutine computes the convective flux by means of the preconditioned AUSM scheme 
!! for 1D stagnation line calorically perfect gase flows.
  SUBROUTINE ausmPC_1D_SL (nx, vol_l, vol_r, left_data, right_data, u_left, u_right, f)

    USE mod_general_data,                 ONLY: nb_eq, pos_u_cell, pos_v_cell, pos_c_cell, pos_pres_cell, & 
                                              & pos_h0_cell, pos_mu_cell, pos_rho, pos_rhou, pos_rhov,    & 
                                              & pos_rhoE
    USE mod_function_pointer,             ONLY: get_prec_vel_1D_SL

    IMPLICIT NONE

    REAL(KIND=8) :: c, dx, dp, mdot, mp, mu, p, rho, u, Vp
    REAL(KIND=8) :: cl, cr, h0l, h0r, mul, mur, pl, pr, rhol, rhor, ul, ur, unl, unr, vl, vr
    REAL(KIND=8), DIMENSION(nb_eq) :: Phil, Phir

    REAL(KIND=8), INTENT(IN) :: nx                        !< normal to the cell interface
    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: f          !< numerical flux

    ! Data of left and right states
    rhol = u_left(pos_rho)
    ul   = left_data(pos_u_cell)
    vl   = left_data(pos_v_cell)
    pl   = left_data(pos_pres_cell)
    cl   = left_data(pos_c_cell)
    h0l  = left_data(pos_h0_cell)
    mul  = left_data(pos_mu_cell)

    rhor = u_right(pos_rho)
    ur   = right_data(pos_u_cell)
    vr   = right_data(pos_v_cell)
    pr   = right_data(pos_pres_cell)
    cr   = right_data(pos_c_cell)
    h0r  = right_data(pos_h0_cell)
    mur  = right_data(pos_mu_cell)
  
    unl = ul*nx
    unr = ur*nx

    ! Interface density, pressure, speed of sound and viscosity
    rho = 0.5d0*(rhol + rhor)
    p   = 0.5d0*(pl + pr)
    c   = 0.5d0*(cl + cr)
    mu  = 0.5d0*(mul + mur)

    ! Interface mass flux 
    mdot = 0.5d0*rho*(unl + unr)

    ! Compute pre-conditioning velocity
    dx = 0.5d0*(vol_l + vol_r)
    dp = ABS(pr - pr)
    u  = MAX(1.d0,0.5d0*ABS(ul + ur))
    CALL get_prec_vel_1D_SL (dx, rho, mu, u, dp, c, Vp)
  
    ! Interface pressure diffusion term (added for pressure velocity coupling at low speeds) 
    mp = - (pr - pl)*nx/Vp

    ! Left and right advection vectors 
    Phil(pos_rho)  = 1.d0
    Phil(pos_rhou) = ul 
    Phil(pos_rhov) = vl
    Phil(pos_rhoE) = h0l

    Phir(pos_rho)  = 1.d0
    Phir(pos_rhou) = ur 
    Phir(pos_rhov) = vr
    Phir(pos_rhoE) = h0r

    ! Numerical flux 
    f = 0.5d0*(mdot + mp)*(Phil + Phir) - 0.5d0*ABS(mdot)*(Phir - Phil)

    ! Adding the pressure contribution in the radial momentum equation
    f(pos_rhou) = f(pos_rhou) + p*nx

  END SUBROUTINE ausmPC_1D_SL
!------------------------------------------------------------------------------!
