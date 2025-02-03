!------------------------------------------------------------------------------!
!> This subroutine computes the convective flux by means of the flux-vector-splitting proposed by van Leer
!  for 1D stagnation line calorically perfect gas flow.
  SUBROUTINE van_Leer_1D_SL (nx, vol_l, vol_r, left_data, right_data, u_left, u_right, f)
  
    USE mod_general_data,                ONLY: nb_eq, pos_u_cell, pos_v_cell, pos_c_cell, pos_h0_cell,   & 
                                             & pos_pres_cell, pos_rho_cell, pos_rho, pos_rhou, pos_rhov, & 
                                             & pos_rhoE, gamma
    USE mod_numerics_data,               ONLY: fl, fr

    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8) :: tmp1, tmp2, tmp3, tmp4, tmp5
    REAL(KIND=8) :: gamma_minus1, gamma_plus1
    REAL(KIND=8) :: ml, mr, cl, cr, h0l, h0r, Mnl, Mnr, pr, pl, rhol, rhor, ul, ur, vl, vr, unl, unr
  
    REAL(KIND=8), INTENT(IN) :: nx                        !< normal to the cell interface
    REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
    REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: f          !< numerical flux 

    ! Common factors
    gamma_plus1  = gamma + 1.d0 
    gamma_minus1 = gamma - 1.d0 
  
    ! Data of left and right states
    ul   = left_data(pos_u_cell)
    vl   = left_data(pos_v_cell)
    cl   = left_data(pos_c_cell)
    rhol = left_data(pos_rho_cell)
    
    ur   = right_data(pos_u_cell)    
    vr   = right_data(pos_v_cell)
    cr   = right_data(pos_c_cell)
    rhor = right_data(pos_rho_cell)
    
    unl = ul*nx
    unr = ur*nx
    
    ! Left and right state normal Mach numbers 
    Mnl = unl/cl
    Mnr = unr/cr
 
    ! Left state
    IF (Mnl.GE.1.d0) THEN
  
       ! Supersonic flow 
       pl  = left_data(pos_pres_cell)
       h0l = left_data(pos_h0_cell)  
      
       ! Mass flow
       ml = rhol*unl
 
       ! Inviscid flux
       fl(pos_rho)  = ml 
       fl(pos_rhou) = ml*ul + pl*nx
       fl(pos_rhov) = ml*vl
       fl(pos_rhoE) = ml*h0l

    ELSE 
  
       IF (Mnl.GT.-1.d0) THEN
  
          ! Subsonic flow
          tmp1 = Mnl + 1.d0
          ml   = 0.25d0*rhol*cl*tmp1*tmp1
  
          tmp1 = cl/gamma
          tmp2 = 2.d0 - Mnl
          tmp3 = gamma_minus1*ul + 2.d0*cl
          tmp4 = 2.d0*gamma_minus1*gamma_plus1 
          tmp5 = tmp1*tmp2

          fl(pos_rho)  = 1.d0
          fl(pos_rhou) = ul + tmp5
          fl(pos_rhov) = vl
          fl(pos_rhoE) = tmp3*tmp3/tmp4

          fl = ml*fl

       ELSE

          fl = 0.d0 

       ENDIF
  
    ENDIF
  
    ! Right state 
    IF (Mnr.LE.-1.d0) THEN
  
       ! Supersonic flow 
       pr  = right_data(pos_pres_cell)
       h0r = right_data(pos_h0_cell)  
 
       ! Mass flow
       mr = rhor*unr  

       ! Inviscid flux
       fr(pos_rho)  = mr 
       fr(pos_rhou) = mr*ur + pr*nx
       fr(pos_rhov) = mr*vr
       fr(pos_rhoE) = mr*h0r

    ELSE 
  
       IF (Mnr.LT.1.d0) THEN
  
          ! Subsonic flow
          tmp1 = Mnr - 1.d0
          mr   = - 0.25d0*rhor*cr*tmp1*tmp1
 
          tmp1 = cr/gamma
          tmp2 = - 2.d0 - Mnr
          tmp3 = gamma_minus1*ur - 2.d0*cr
          tmp4 = 2.d0*gamma_minus1*gamma_plus1 
          tmp5 = tmp1*tmp2

          fr(pos_rho)  = 1.d0
          fr(pos_rhou) = ur + tmp5
          fr(pos_rhov) = vr
          fr(pos_rhoE) = tmp3*tmp3/tmp4

          fr = mr*fr
 
       ELSE
 
          fr = 0.d0 

       ENDIF
  
    ENDIF
  
    ! Adding left and right contributes
    f = fl + fr
 
  END SUBROUTINE van_Leer_1D_SL
!------------------------------------------------------------------------------!
