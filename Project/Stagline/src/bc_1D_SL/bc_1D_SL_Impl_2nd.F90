!------------------------------------------------------------------------------!
!> This subroutine applies a supersonic inlet boundary condition for 1D stagnation line 
!! calorically perfect gas flows.
  SUBROUTINE sup_in_1D_SL_Impl_2nd (id, phys_data1, u_phys1, phys_data2, u_phys2, & 
                                 &  ghost_data1, u_ghost1, ghost_data2, u_ghost2, &
                                 &  ju_ghost)

    USE mod_general_data,          ONLY: pos_u_cell, pos_v_cell, pos_ek_cell, pos_pres_cell, pos_rho_cell, nb_eq 
    USE mod_domain_boundary,       ONLY: boundary, boundary_data, get_boundary_inlet 
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim

    IMPLICIT NONE

    REAL(KIND=8) :: ov_3
    REAL(KIND=8) :: rho_in, u_in, v_in, p_in
    REAL(KIND=8) :: rhoG, rhoP, uG, uP, vG, vP, ekG, ekP
    REAL(KIND=8) :: rhoG_ov_rhoP, uP_rhoG_ov_rhoP
    REAL(KIND=8), DIMENSION(nb_eq) :: inlet_data, prim
    TYPE(boundary) :: bound

    INTEGER, INTENT(IN) :: id                               !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys1       !< conservative variables of the first physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys2       !< conservative variables of the second physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data1    !< physical properties of the first physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data2    !< physical properties of the second physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost1     !< conservative variables of the first ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost2     !< conservative variables of the second ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data1  !< physical properties of the first ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data2  !< physical properties of the second ghost cell
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: ju_ghost   !< ghost state conservative variable Jacobian

    ! Inlet values of density, velocity components and pressure
    bound      = boundary_data(id)
    inlet_data = get_boundary_inlet(nb_eq, bound)

    rho_in = inlet_data(1)
    u_in   = inlet_data(2)
    v_in   = inlet_data(3)
    p_in   = inlet_data(4)

    ! Common factor
    ov_3 = 1.d0/3.d0

    ! Ghost state 1
    ! Primitive variables
    prim(1) = (4.d0*rho_in - phys_data2(pos_rho_cell))*ov_3
    prim(2) = (4.d0*u_in   - phys_data2(pos_u_cell))*ov_3
    prim(3) = (4.d0*v_in   - phys_data2(pos_v_cell))*ov_3
    prim(4) = (4.d0*p_in   - phys_data2(pos_pres_cell))*ov_3

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost1, ghost_data1)
 
    ! Ghost state 2
    ! Primitive variables
    prim(1) = 2.0*rho_in - phys_data2(pos_rho_cell)
    prim(2) = 2.d0*u_in  - phys_data2(pos_u_cell)
    prim(3) = 2.d0*v_in  - phys_data2(pos_v_cell)
    prim(4) = 2.d0*p_in  - phys_data2(pos_pres_cell)

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost2, ghost_data2)    
 
    ! Evaluation of the matrix dUg/dUp for implicti boundary condition 
    ! (the Jacobian matrix is computed as for a first order scheme for sake of robustness)
    ! Density, velocity components and kinetic energy per unit mass of physical state 1
    rhoP = phys_data1(pos_rho_cell)
    uP   = phys_data1(pos_u_cell)
    vP   = phys_data1(pos_v_cell)
    ekP  = phys_data1(pos_ek_cell)
    
    ! Density, velocity components and kinetic energy per unit mass of ghost state 1
    rhoG = ghost_data1(pos_rho_cell)
    uG   = ghost_data1(pos_u_cell)
    vG   = ghost_data1(pos_v_cell)
    ekG  = ghost_data1(pos_ek_cell) 
    rhoG_ov_rhoP = rhoG/rhoP
    uP_rhoG_ov_rhoP = uP*rhoG_ov_rhoP    

    ! First column 
    ju_ghost(1,1) = - 1.d0
    ju_ghost(2,1) = - uG + uP_rhoG_ov_rhoP
    ju_ghost(3,1) = - vG + vP*rhoG_ov_rhoP
    ju_ghost(4,1) = - ekG - ekP + uG*uP_rhoG_ov_rhoP

    ! Second column
    ju_ghost(1,2) = 0.d0
    ju_ghost(2,2) = - rhoG_ov_rhoP
    ju_ghost(3,2) = 0.d0
    ju_ghost(4,2) = uP - uG*rhoG_ov_rhoP

    ! Third column 
    ju_ghost(1,3) = 0.d0
    ju_ghost(2,3) = 0.d0
    ju_ghost(3,3) = - rhoG_ov_rhoP
    ju_ghost(4,3) = 0.d0

    ! Fourth column 
    ju_ghost(1,4) = 0.d0
    ju_ghost(2,4) = 0.d0
    ju_ghost(3,4) = 0.d0
    ju_ghost(4,4) = - 1.d0

  END SUBROUTINE sup_in_1D_SL_Impl_2nd
!------------------------------------------------------------------------------!
!> This subroutine applies a slip wall boundary condition for 1D stagnation line 
!! calorically perfect gas flows. 
  SUBROUTINE slip_wall_1D_SL_Impl_2nd (id, phys_data1, u_phys1, phys_data2, u_phys2, & 
                                    &  ghost_data1, u_ghost1, ghost_data2, u_ghost2, & 
                                    &  ju_ghost)

    USE mod_general_data,          ONLY: pos_u_cell, pos_v_cell, pos_ek_cell, pos_pres_cell, pos_rho_cell, nb_eq
    USE mod_domain_boundary,       ONLY: boundary_data, get_boundary_rhoin, get_boundary_pin
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim 
 
    IMPLICIT NONE

    REAL(KIND=8) :: rhoG, rhoP, uG, uP, vG, vP, ekG, ekP
    REAL(KIND=8) :: rhoP_ov_rhoG, uP_rhoP_ov_rhoG
    REAL(KIND=8), DIMENSION(nb_eq) :: prim

    INTEGER, INTENT(IN) :: id                               !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys1       !< conservative variables of the first physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys2       !< conservative variables of the second physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data1    !< physical properties of the first physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data2    !< physical properties of the second physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost1     !< conservative variables of the first ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost2     !< conservative variables of the second ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data1  !< physical properties of the first ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data2  !< physical properties of the second ghost cell
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: ju_ghost   !< ghost state conservative variable Jacobian

    ! Density, velocity components and kinetic energy of physical state 1
    rhoP = phys_data1(pos_rho_cell)
    uP   = phys_data1(pos_u_cell)
    vP   = phys_data1(pos_v_cell)
    ekP  = phys_data1(pos_ek_cell)

    ! Ghost state 1
    ! Primitive variables 
    prim(1) = rhoP
    prim(2) = - uP
    prim(3) = vP
    prim(4) = (4.d0*phys_data1(pos_pres_cell) - phys_data2(pos_pres_cell))/3.d0 
   
    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost1, ghost_data1)
 
    ! Ghost state 2
    u_ghost2    = u_ghost1
    ghost_data2 = ghost_data1
 
    ! Evaluation of the matrix dUg/dUp for implicti boundary condition
    ! (the Jacobian matrix is computed as for a first order scheme for sake of robustness) 
    ! Density, velocity components and kinetic energy per unit mass of ghost state 1
    rhoG = ghost_data1(pos_rho_cell)
    uG   = ghost_data1(pos_u_cell)
    vG   = ghost_data1(pos_v_cell)
    ekG  = ghost_data1(pos_ek_cell) 
    rhoP_ov_rhoG = rhoG/rhoP
    uP_rhoP_ov_rhoG = uP*rhoP_ov_rhoG

    ! First column 
    ju_ghost(1,1) = 1.d0
    ju_ghost(2,1) = uG + uP_rhoP_ov_rhoG
    ju_ghost(3,1) = vG - vP*rhoP_ov_rhoG
    ju_ghost(4,1) = ekG + ekP + uG*uP_rhoP_ov_rhoG

    ! Second column
    ju_ghost(1,2) = 0.d0
    ju_ghost(2,2) = - rhoP_ov_rhoG
    ju_ghost(3,2) = 0.d0
    ju_ghost(4,2) = - uP - uG*rhoP_ov_rhoG

    ! Third column 
    ju_ghost(1,3) = 0.d0
    ju_ghost(2,3) = 0.d0
    ju_ghost(3,3) = rhoP_ov_rhoG
    ju_ghost(4,3) = 0.d0

    ! Fourth column 
    ju_ghost(1,4) = 0.d0
    ju_ghost(2,4) = 0.d0
    ju_ghost(3,4) = 0.d0
    ju_ghost(4,4) = 1.d0    

  END SUBROUTINE slip_wall_1D_SL_Impl_2nd
!------------------------------------------------------------------------------!
!> This subroutine applies an isothermal slip wall boundary condition for 1D stagnation line 
!! calorically perfect gas flows. 
  SUBROUTINE no_slip_iso_Twall_1D_SL_Impl_2nd (id, phys_data1, u_phys1, phys_data2, u_phys2, & 
                                            &  ghost_data1, u_ghost1, ghost_data2, u_ghost2, & 
                                            &  ju_ghost)
    
    USE mod_general_data,          ONLY: pos_u_cell, pos_v_cell, pos_pres_cell, pos_T_cell, pos_rho_cell, pos_ek_cell, & 
                                       & nb_eq, R_gas, gamma
    USE mod_domain_boundary,       ONLY: boundary, boundary_data, get_boundary_Twall
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim 
 
    IMPLICIT NONE

    REAL(KIND=8), PARAMETER :: Tmin = 400.d0
    REAL(KIND=8) :: u1, v1, p1, p2, T2, Twall
    REAL(KIND=8) :: pG1, TG1, TG2
    REAL(KIND=8) :: rhoP, uP, vP, p, TP
    REAL(KIND=8) :: rhoG, pG, uG, vG, TG
    REAL(KIND=8) :: rhoG_ov_rhoP, ekP, EG, ov_TG
    REAL(KIND=8) :: fac1, fac2, fac3, fac4, fac5
    REAL(KIND=8) :: gamma_minus1 
    REAL(KIND=8), DIMENSION(nb_eq) :: prim
    TYPE(boundary) :: bound

    INTEGER, INTENT(IN) :: id                               !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys1       !< conservative variables of the first physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys2       !< conservative variables of the second physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data1    !< physical properties of the first physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data2    !< physical properties of the second physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost1     !< conservative variables of the first ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost2     !< conservative variables of the second ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data1  !< physical properties of the first ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data2  !< physical properties of the second ghost cell
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: ju_ghost   !< ghost state conservative variable Jacobian

    ! Wall temperature
    bound = boundary_data(id)
    Twall = get_boundary_Twall(bound)

    ! Pressure of physical states 1 and 2
    p1 = phys_data1(pos_pres_cell)
    p2 = phys_data2(pos_pres_cell)

    ! Temperature of physical state 2
    T2 = phys_data2(pos_T_cell)

    ! Velocity components of physical state 1
    u1 = phys_data1(pos_u_cell)
    v1 = phys_data1(pos_v_cell)

    ! Ghost state 1
    ! Primitive variables
    pG1 = (4.d0*p1 - p2)/3.d0
    TG1 = MAX((4.d0*Twall - T2)/3.d0,Tmin)

    prim(1) = pG1/(R_gas*TG1)
    prim(2) = - u1
    prim(3) = - v1
    prim(4) = pG1

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost1, ghost_data1)
 
    ! Ghost state 2
    ! Primitive variables
    TG2 = MAX(2.d0*Twall - T2,Tmin)

    prim(1) = pG1/(R_gas*TG2)
    prim(2) = - u1
    prim(3) = - v1
    prim(4) = pG1

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost2, ghost_data2) 

    ! Evaluation of the matrix dUg/dUp for implicit boundary condition
    ! (the Jacobian matrix is computed as for a first order scheme for sake of robustness)
    ! Common factors
    gamma_minus1 = gamma - 1.d0 
    TP    = phys_data1(pos_T_cell) 
    ekP   = phys_data1(pos_ek_cell)
    rhoP  = phys_data1(pos_rho_cell)
    rhoG  = ghost_data1(pos_rho_cell)
    uG    = ghost_data1(pos_u_cell)
    vG    = ghost_data1(pos_v_cell)
    uP    = u1
    vP    = v1
    TG    = TG1
    EG    = u_ghost1(4)/rhoG
    ov_TG = 1.d0/TG 
    fac1  = (1.d0 + TP*ov_TG)*ov_TG 
    fac2  = gamma_minus1/R_gas*ekP - TP 
    fac3  = (TP*ov_TG + fac2*fac1)
    fac4  = fac1*gamma_minus1/R_gas
    fac5  = fac4*uP
    rhoG_ov_rhoP = rhoG/rhoP 

    ! First column 
    ju_ghost(1,1) = fac3
    ju_ghost(2,1) = uG*fac3 + rhoG_ov_rhoP*uP
    ju_ghost(3,1) = vG*fac3 + rhoG_ov_rhoP*vP
    ju_ghost(4,1) = EG*fac3 + rhoG_ov_rhoP*(uG*uP - fac2*R_gas/gamma_minus1)

    ! Second column 
    ju_ghost(1,2) = - fac5
    ju_ghost(2,2) = - fac5*uG - rhoG_ov_rhoP
    ju_ghost(3,2) = - fac5*vG
    ju_ghost(4,2) = - fac5*EG + rhoG_ov_rhoP*(uP - UG)    

    ! Third column 
    ju_ghost(1,3) = 0.d0
    ju_ghost(2,3) = 0.d0
    ju_ghost(3,3) = - rhoG_ov_rhoP
    ju_ghost(4,3) = 0.d0

    ! Fourth column
    ju_ghost(1,4) = fac4
    ju_ghost(2,4) = fac4*uG
    ju_ghost(3,4) = fac4*vG
    ju_ghost(4,4) = - 1.d0 + fac4*EG
    
  END SUBROUTINE no_slip_iso_Twall_1D_SL_Impl_2nd
!------------------------------------------------------------------------------!
