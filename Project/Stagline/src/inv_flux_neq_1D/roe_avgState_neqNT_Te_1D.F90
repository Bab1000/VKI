!------------------------------------------------------------------------------!
!> This subroutine computes the Roe's averaged state for a 1D nonequilibrium flow. 
!! In this case the flow is characterized by N temperatures plus a separated temperature Te for free electrons.
  SUBROUTINE roe_avgState_neqNT_Te_1D (u_left, u_right, left_data, right_data, u_Roe, data_Roe)

    USE mod_general_data,             ONLY: nb_ns, nb_eq, nb_temp, nb_int_temp, pos_u_cell, pos_h0_cell,  & 
                                          & pos_c_cell, pos_pres_cell, pos_rho_cell, pos_T_cell,          & 
                                          & pos_ei_cell, pos_gamma_cell, pos_ek_cell, pos_alpha_cell,     & 
                                          & pos_beta_cell, pos_em, pos_rhou, pos_rhoE, pos_rhoek,         & 
                                          & pos_rhose, pos_Te, gamma_e_m1, yi, rhoi, ei, Ri, beta, temp
    USE mod_neq_function_pointer,     ONLY: library_get_thermodynamic_data

    IMPLICIT NONE

    INTEGER :: i, k
    REAL(KIND=8) :: a, b, r
    REAL(KIND=8) :: tmp1, tmp2
    REAL(KIND=8) :: h0l, h0r, rhol, rhor, Tl, Tr, Sel, Ser, ul, ur 
    REAL(KIND=8) :: c, c2, cn, ek, h0, rho, rhou, rhoE, p, T, Te, u, Se
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
    h0l  = left_data(pos_h0_cell)
    rhol = left_data(pos_rho_cell)
    Sel  = u_left(nb_eq)

    Tr   = right_data(pos_T_cell)
    ur   = right_data(pos_u_cell)
    h0r  = right_data(pos_h0_cell)
    rhor = right_data(pos_rho_cell)
    Ser  = u_right(nb_eq)

    ! Roe's averaged state
    r = DSQRT(rhor/rhol)
    a = 1.d0/(1.d0 + r)
    b = a*r

    ! Average density, total enthalpy, temperatures, velocity, kinetic energy
    ! and pseudo-entropy. The formula used for the temperatures is rigorous only in the case of constant specific heats
    rho = DSQRT(rhol*rhor) 
    h0  = a*h0l  + b*h0r
    u   = a*ul   + b*ur
    DO i = 1,1 + nb_int_temp
       temp(i) = a*left_data(pos_T_cell + i - 1) + b*right_data(pos_T_cell + i - 1)  
    ENDDO
    Se  = a*Sel  + b*Ser
    ek  = 0.5d0*u**2

    ! Average species mass fractions and densities
    tmp1 = a/rhol
    tmp2 = b/rhor 
    DO i = 1,nb_ns 
       yi(i) = u_left(i)*tmp1 + u_right(i)*tmp2
    ENDDO
    rhoi = yi*rho

    ! Free electron temperature computed from the free electron pseudo-entropy
    Te = Se*(rho**gamma_e_m1)/(Ri(pos_em)*rhoi(pos_em))
    temp(pos_Te) = Te

    ! Thermodynamic data of Roe' averaged state
    !temp(1) = T
    !temp(2) = Te
    CALL library_get_thermodynamic_data (rho, rhoi, temp, c, gamma, p, alpha, beta, ei) 

    ! Total momentum and energy density 
    rhou = rho*u
    rhoE = rho*h0 - p
 
    ! Roe's averaged state conservative variable vector 
    DO i = 1,nb_ns 
       u_Roe(i) = rhoi(i)
    ENDDO
    u_Roe(pos_rhou)  = rhou 
    u_Roe(pos_rhoE)  = rhoE

    DO k = 1,nb_int_temp
       u_Roe(pos_rhoek + k - 1) = a*u_left(pos_rhoek+ k - 1) + b*u_right(pos_rhoek+ k - 1) 
    ENDDO

    u_Roe(pos_rhose) = Se

    ! Roes' averaged state physical data vector 
    DO i = 1,nb_temp 
       data_Roe(pos_T_cell + i - 1) = temp(i)
    ENDDO  
    data_Roe(pos_u_cell)     = u
    data_Roe(pos_h0_cell)    = h0
    data_Roe(pos_c_cell)     = c
    data_Roe(pos_gamma_cell) = gamma
    data_Roe(pos_pres_cell)  = p
    data_Roe(pos_rho_cell)   = rho
    data_Roe(pos_ek_cell)    = ek
    data_Roe(pos_alpha_cell) = alpha 

    DO i = 1,nb_int_temp + 1
       data_Roe(pos_beta_cell + i - 1) = beta(i)
    ENDDO

    DO i = 1,nb_ns + nb_ns*nb_int_temp
       data_Roe(pos_ei_cell + i - 1) = ei(i)  
    ENDDO

  END SUBROUTINE roe_avgState_neqNT_Te_1D 
!------------------------------------------------------------------------------!
