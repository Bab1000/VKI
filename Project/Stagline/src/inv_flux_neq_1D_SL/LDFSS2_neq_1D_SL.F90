!------------------------------------------------------------------------------!
!> This subroutine computes the convective flux by means of the preconditioned AUSM scheme 
!! for 1D stagnation line calorically perfect gase flows.
  SUBROUTINE LDFSS2_neq_1D_SL (nx, vol_l, vol_r, left_data, right_data, u_left, u_right, f)

    USE mod_general_data,             ONLY: nb_eq, nb_ns, nb_temp, nb_int_temp, nb_te, pos_u_cell, pos_v_cell,     & 
                                          & pos_h0_cell, pos_pres_cell, pos_gamma_cell, pos_c_cell, pos_rho_cell,  & 
                                          & pos_rhou, pos_rhov, pos_rhoE, pos_T_cell, pos_mu_cell
    USE mod_function_pointer,             ONLY: get_prec_vel_1D_SL
 

    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8) :: rhol, rhor, h0l, h0r, pl, pr, ul, ur, vl, vr, al, ar, dx
    REAL(KIND=8) :: ahalf, ml, mr, betal, betar, mhalf, mhalfp, mhalfm
    REAL(KIND=8) :: alphap, alpham, cvlp, cvlm, cp, cm, plp, prm, dlp, drm
    REAL(KIND=8), DIMENSION(nb_eq) :: phil, phir
    REAL(KIND=8), PARAMETER :: delta = 1.0d0

    REAL(KIND=8), INTENT(IN) :: nx                        !< normal to the cell interface
    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: f          !< numerical flux


    ! Data of left and right states
    rhol = left_data(pos_rho_cell)
    rhor = right_data(pos_rho_cell)
    pl   = left_data(pos_pres_cell)   ! pressure
    pr   = right_data(pos_pres_cell)
    ul   = left_data(pos_u_cell)      ! radial velocoity
    ur   = right_data(pos_u_cell)
    vl   = left_data(pos_v_cell)      ! tangential velocity
    vr   = right_data(pos_v_cell)
    al   = left_data(pos_c_cell)      ! sound speed
    ar   = right_data(pos_c_cell)
    h0l  = left_data(pos_h0_cell)
    h0r  = right_data(pos_h0_cell)
    
    ahalf = 0.5d0 * (al + ar)
    ml    = ul/ahalf
    mr    = ur/ahalf
    
    betal = -max(0.0d0, 1.0d0-int(abs(ml)))
    betar = -max(0.0d0, 1.0d0-int(abs(mr)))
    mhalf = 0.25d0*betal*betar*(sqrt(0.5d0*(ml*ml+mr*mr))-1.0d0)**2.0d0
    
    mhalfp = mhalf*max(0.0d0, 1.0d0-(pl-pr)/(pl+pr)-2.0d0*abs(pl-pr)/pl)
    mhalfm = mhalf*max(0.0d0, 1.0d0+(pl-pr)/(pl+pr)-2.0d0*abs(pl-pr)/pr)
    
    alphap = 0.5d0*(1.0d0+sign(1.0d0,ml))
    alpham = 0.5d0*(1.0d0-sign(1.0d0,mr))
    
    cvlp = alphap*(1.0d0+betal)*ml-0.25d0*betal*(ml+1.0d0)**2.0d0
    cvlm = alpham*(1.0d0+betar)*mr+0.25d0*betar*(mr-1.0d0)**2.0d0
    
    cp = cvlp - mhalfp
    cm = cvlm + mhalfm
    
    plp = 0.25d0*(2.0d0-ml)*(ml+1.0d0)**2.0d0
    prm = 0.25d0*(2.0d0+mr)*(mr-1.0d0)**2.0d0
    
    dlp = alphap*(1.0d0+betal)-betal*plp
    drm = alpham*(1.0d0+betar)-betar*prm
    
    phil = u_left;  phil(pos_rhoE) = phil(pos_rhoE)+pl
    phir = u_right; phir(pos_rhoE) = phir(pos_rhoE)+pr
    
    dx = 0.5d0*(vol_l+vol_r)
    f = ahalf*(cp*phil + cm*phir)
    f(pos_rhou) = f(pos_rhou) + nx*(dlp*pl + drm*pr)
    
    !write(*,*) ">>>>", u_left, "<<<<"
    !write(*,*) ">>>>", u_right, "<<<<"
    !write(*,*) ">>>>", f, "<<<<"
    !stop


  END SUBROUTINE LDFSS2_neq_1D_SL
!------------------------------------------------------------------------------!
