!------------------------------------------------------------------------------!
!> This subroutine computes the convective flux by means of the flux-vector-splitting proposed by van Leer 
!! for 1D calorically perfect gas flows.
  SUBROUTINE van_Leer_1D (nx, left_data, right_data, u_left, u_right, f)
  
    USE mod_general_data,                ONLY: nb_eq, pos_u_cell, pos_c_cell, pos_h0_cell, pos_pres_cell, &
                                             & pos_rho_cell, gamma
    USE mod_numerics_data,               ONLY: fl, fr

    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8) :: tmp
    REAL(KIND=8) :: gamma_minus1, gamma_plus1
    REAL(KIND=8) :: ml, mr
    REAL(KIND=8) :: cl, cr, h0l, h0r, Mnl, Mnr, pr, pl, rhol, rhor, ul, ur, unl, unr
 
    REAL(KIND=8), INTENT(IN) :: nx                        !< normal to the cell interface
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
    cl   = left_data(pos_c_cell)
    rhol = left_data(pos_rho_cell)
    
    ur   = right_data(pos_u_cell)    
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
       pl = left_data(pos_pres_cell)
  
       DO i = 1,nb_eq 
          fl(i) = u_left(i)*unl  
       ENDDO
       
       ! Adding pressure contributions to global momentum and energy equations
       fl(2) = fl(2) + pl*nx 
       fl(3) = fl(3) + pl*unl

    ELSE 
  
       IF (Mnl.GT.-1.d0) THEN
  
          ! Subsonic flow
          h0l = left_data(pos_h0_cell)  
          tmp = Mnl + 1.d0
          ml  = 0.25d0*rhol*cl*tmp*tmp
  
          tmp = unl - cl  
          fl(1) = ml
          fl(2) = ml*((-unl + 2.d0*cl)*nx/gamma + ul) 
          fl(3) = ml*(h0l - tmp*tmp/gamma_plus1)
  
       ELSE

          fl = 0.d0 

       ENDIF
  
    ENDIF
  
    ! Right state 
    IF (Mnr.LE.-1.d0) THEN
  
       ! Supersonic flow 
       pr = right_data(pos_pres_cell)
  
       DO i = 1,nb_eq 
          fr(i) = u_right(i)*unr  
       ENDDO
       
       ! Adding pressure contributions to global momentum and energy equations
       fr(2) = fr(2) + pr*nx 
       fr(3) = fr(3) + pr*unr

    ELSE 
  
       IF (Mnr.LT.1.d0) THEN
  
          ! Subsonic flow
          h0r = right_data(pos_h0_cell)
          tmp = Mnr - 1.d0
          mr  = - 0.25d0*rhor*cr*tmp*tmp
  
          tmp = unr + cr
  
          fr(1) = mr
          fr(2) = mr*((- unr - 2.d0*cr)*nx/gamma + ur)
          fr(3) = mr*(h0r - tmp*tmp/gamma_plus1)
  
       ELSE
 
          fr = 0.d0 

       ENDIF
  
    ENDIF
  
    ! Adding left and right contributes
    f = fl + fr
  
    END SUBROUTINE van_Leer_1D
!------------------------------------------------------------------------------!
