!------------------------------------------------------------------------------!
!> This subroutine computes the convective flux by means of the flux-vector-splitting 
!! proposed by van Leer for 1D stagnation line nonequilibrium flows.
  SUBROUTINE van_leer_neq_1D_SL (nx, vol_l, vol_r, left_data, right_data, u_left, u_right, f)
 
    USE mod_general_data,             ONLY: nb_eq, nb_ns, nb_temp, nb_int_temp, nb_te, pos_u_cell, pos_v_cell,     & 
                                          & pos_h0_cell, pos_pres_cell, pos_gamma_cell, pos_c_cell, pos_rho_cell,  & 
                                          & pos_rhou, pos_rhov, pos_rhoE, pos_T_cell
    USE mod_numerics_data,            ONLY: fl, fr
     
    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8) :: tmp1, tmp2, tmp3
    REAL(KIND=8) :: cl, cr, gl, gr, h0l, h0r, Mnl, Mnr, pr, pl, rhol, rhor, ul, ur, vl, vr, unl, unr
    REAL(KIND=8) :: fml, fmr, fmom_l, fmom_r, fene_l, fene_r
    
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
    cl   = left_data(pos_c_cell)
    ul   = left_data(pos_u_cell)
    vl   = left_data(pos_v_cell)
    
    rhor = right_data(pos_rho_cell)
    cr   = right_data(pos_c_cell)
    ur   = right_data(pos_u_cell)
    vr   = right_data(pos_v_cell)    

    unl = ul*nx
    unr = ur*nx
  
    ! Normal Mach numbers of left and right states
    Mnl = unl/cl
    Mnr = unr/cr
   
    ! Left state
    IF (Mnl.GE.1.d0) THEN
  
       ! Supersonic flow
       DO i = 1,nb_eq
          fl(i) = u_left(i)*unl 
       ENDDO
  
       ! Adding the pressure contribution to the global momentum and energy equations
       pl = left_data(pos_pres_cell)
       fl(pos_rhou) = fl(pos_rhou) + pl*nx
       fl(pos_rhoE) = fl(pos_rhoE) + pl*unl
  
    ELSE  
  
      IF (Mnl.GT.-1.d0) THEN
                
         gl  = left_data(pos_gamma_cell)
         h0l = left_data(pos_h0_cell)

         ! Subsonic flow
         tmp1   = Mnl + 1.d0
         fml    = 0.25d0*cl*tmp1*tmp1
         tmp2   = rhol*fml
         fmom_l = (-unl + 2.d0*cl)*nx/gl + ul
         tmp3   = unl - cl
         fene_l = h0l - tmp3*tmp3/(gl + 1.d0)
  
         DO i = 1,nb_ns 
            fl(i) = u_left(i)*fml
         ENDDO

         fl(pos_rhou) = tmp2*fmom_l 
         fl(pos_rhov) = tmp2*vl 
         fl(pos_rhoE) = tmp2*fene_l 
  
         DO i = 1,nb_int_temp + nb_te
            fl(pos_rhoE + i) = u_left(pos_rhoE + i)*fml
         ENDDO
  
       ELSE 
 
         fl = 0.d0
 
       ENDIF
  
    ENDIF
  
    ! Right state
    IF (Mnr.LE.-1d0) THEN
  
       ! Supersonic flow
       DO i = 1,nb_eq
          fr(i) = u_right(i)*unr 
       ENDDO
  
       ! Adding the pressure contribution to the global momentum and energy equations
       pr = right_data(pos_pres_cell)
       fr(pos_rhou) = fr(pos_rhou) + pr*nx
       fr(pos_rhoE) = fr(pos_rhoE) + pr*unr
  
    ELSE

      IF (Mnr.LT.1.d0) THEN
                 
         gr  = right_data(pos_gamma_cell)
         h0r = right_data(pos_h0_cell)

         ! Subsonic flow
         tmp1   = Mnr - 1.d0
         fmr    = - 0.25d0*cr*tmp1*tmp1
         tmp2   = rhor*fmr
         fmom_r = (-unr - 2.d0*cr)*nx/gr + ur
         tmp3   = unr + cr
         fene_r = h0r - tmp3*tmp3/(gr + 1.d0)
  
         DO i = 1,nb_ns 
            fr(i) = u_right(i)*fmr
         ENDDO

         fr(pos_rhou) = tmp2*fmom_r 
         fr(pos_rhov) = tmp2*vr
         fr(pos_rhoE) = tmp2*fene_r
  
         DO i = 1,nb_int_temp + nb_te 
            fr(pos_rhoE + i) = u_right(pos_rhoE + i)*fmr
         ENDDO 
  
      ELSE 

         fr = 0.d0

      ENDIF
  
    ENDIF
  
    ! Adding left and right contributes
    f = fl + fr

  END SUBROUTINE van_Leer_neq_1D_SL
!------------------------------------------------------------------------------!
