!------------------------------------------------------------------------------!
!> This subroutine computes the Roe's averaged state for 1D calorically perfect gas flows.
  SUBROUTINE roe_avgState_1D (left_data, right_data, roe_data)

    USE mod_general_data,     ONLY: nb_eq, pos_u_cell, pos_pres_cell, pos_h0_cell, pos_c_cell, & 
                                  & pos_ek_cell, pos_rho_cell, gamma


    IMPLICIT NONE

    REAL(KIND=8) :: gamma_minus1
    REAL(KIND=8) :: tmp
    REAL(KIND=8) :: a, b, r
    REAL(KIND=8) :: c, ek, h0, u, rho, p
    REAL(KIND=8) :: h0l, h0r, rhol, rhor, ul, ur
    
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: roe_data   !< physical properties of Roe's averaged state

    ! Common factor 
    gamma_minus1 = gamma - 1.d0

    ! Velocity, total enthalpy and density of left and right states
    ul   = left_data(pos_u_cell)
    h0l  = left_data(pos_h0_cell)
    rhol = left_data(pos_rho_cell)
   
    ur   = right_data(pos_u_cell)
    h0r  = right_data(pos_h0_cell)
    rhor = right_data(pos_rho_cell)

    ! Roe's averaged state
    r = DSQRT(rhor/rhol)
    a = 1.d0/(1.d0 + r)
    b = a*r

    rho = DSQRT(rhol*rhor)
    u   = a*ul + b*ur
    ek  = 0.5d0*u**2
    h0  = a*h0l + b*h0r
    tmp = gamma_minus1*(h0 - ek)
    c   = DSQRT(tmp)
    p   = rho*tmp/gamma

    roe_data(pos_u_cell)    = u
    roe_data(pos_h0_cell)   = h0
    roe_data(pos_ek_cell)   = ek
    roe_data(pos_c_cell)    = c
    roe_data(pos_pres_cell) = p
    roe_data(pos_rho_cell)  = rho 

  END SUBROUTINE roe_avgState_1D
!------------------------------------------------------------------------------!
