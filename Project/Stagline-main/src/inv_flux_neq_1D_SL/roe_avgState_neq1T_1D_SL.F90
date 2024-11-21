!------------------------------------------------------------------------------!
!> This subroutine computes the Roe's averaged state for 1D stagnation line nonequilibrium flows. 
!! In this case only one temperature is used.
  SUBROUTINE roe_avgState_neq1T_1D_SL (u_left, u_right, left_data, right_data, u_Roe, data_Roe)

    USE mod_general_data,             ONLY: nb_ns, nb_eq, nb_temp, pos_u_cell, pos_v_cell, pos_h0_cell, & 
                                          & pos_c_cell, pos_pres_cell, pos_rho_cell, pos_T_cell,        &
                                          & pos_ei_cell, pos_gamma_cell, pos_ek_cell, pos_alpha_cell,   & 
                                          & pos_beta_cell, pos_rhou, pos_rhov, pos_rhoE, yi, rhoi, ei,  & 
                                          & beta, temp
    USE mod_neq_function_pointer,     ONLY: library_get_thermodynamic_data

    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8) :: a, b, r
    REAL(KIND=8) :: tmp1, tmp2
    REAL(KIND=8) :: h0l, h0r, rhol, rhor, Tl, Tr, ul, ur, vl, vr 
    REAL(KIND=8) :: c, c2, cn, ek, h0, rho, rhou, rhov, rhoE, p, T, u, v
    REAL(KIND=8) :: alpha, gamma

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !> conservative variabbles of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !> conservative variabbles of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !> physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !> physical properties of right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_Roe      !> conservative variables of Roe's averaged state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: data_Roe   !> physical properties of Roe's averaged state

    ! Data of left and right states
    Tl   = left_data(pos_T_cell)
    ul   = left_data(pos_u_cell)
    vl   = left_data(pos_v_cell)
    h0l  = left_data(pos_h0_cell)
    rhol = left_data(pos_rho_cell)

    Tr   = right_data(pos_T_cell)
    ur   = right_data(pos_u_cell)
    vr   = right_data(pos_v_cell)
    h0r  = right_data(pos_h0_cell)
    rhor = right_data(pos_rho_cell)

    ! Roe's averaged state
    r = DSQRT(rhor/rhol)
    a = 1.d0/(1.d0 + r)
    b = a*r
    
    ! Average density, total enthalpy, temperature, velocity and kinetic energy
    ! The formula used for the temperature is rigorous only in the case of constant specific heats 
    rho = DSQRT(rhol*rhor) 
    h0  = a*h0l + b*h0r
    u   = a*ul + b*ur
    v   = a*vl + b*vr
    T   = a*Tl + b*Tr
    ek  = 0.5d0*u**2
    
    ! Average species mass fractions and densities
    tmp1 = a/rhol
    tmp2 = b/rhor 
    DO i = 1,nb_ns 
       yi(i) = u_left(i)*tmp1 + u_right(i)*tmp2
    ENDDO
    rhoi = yi*rho

    ! Thermodynamic data of Roe' averaged state
    temp(1) = T
    CALL library_get_thermodynamic_data (rho, rhoi, temp, c, gamma, p, alpha, beta, ei) 

    ! Total momentum components and energy density 
    rhou = rho*u
    rhov = rho*v
    rhoE = rho*h0 - p
 
    ! Roe's averaged state conservative variable vector 
    DO i = 1,nb_ns 
       u_Roe(i) = rhoi(i)
    ENDDO
    u_Roe(pos_rhou) = rhou 
    u_Roe(pos_rhov) = rhov
    u_Roe(pos_rhoE) = rhoE

    ! Roes' averaged state physical data vector 
    data_Roe(pos_T_cell)     = T  
    data_Roe(pos_u_cell)     = u
    data_Roe(pos_v_cell)     = v
    data_Roe(pos_h0_cell)    = h0
    data_Roe(pos_c_cell)     = c
    data_Roe(pos_gamma_cell) = gamma
    data_Roe(pos_pres_cell)  = p
    data_Roe(pos_rho_cell)   = rho
    data_Roe(pos_ek_cell)    = ek
    data_Roe(pos_alpha_cell) = alpha 
    data_Roe(pos_beta_cell)  = beta(1)

    DO i = 1,nb_ns 
       data_Roe(pos_ei_cell + i - 1) = ei(i)  
    ENDDO

  END SUBROUTINE roe_avgState_neq1T_1D_SL 
!------------------------------------------------------------------------------!
