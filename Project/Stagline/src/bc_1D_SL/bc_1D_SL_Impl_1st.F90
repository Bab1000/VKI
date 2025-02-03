!------------------------------------------------------------------------------!
!> This subroutine applies a supersonic inlet boundary condition for 1D stagnation line 
!! calorically perfect gas flows.
  SUBROUTINE sup_in_1D_SL_Impl (id, phys_data, u_phys, ghost_data, u_ghost, ju_ghost)

    USE mod_general_data,          ONLY: pos_u_cell, pos_v_cell, pos_pres_cell, pos_rho_cell, pos_ek_cell, nb_eq 
    USE mod_domain_boundary,       ONLY: boundary, boundary_data, get_boundary_inlet 
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim

    IMPLICIT NONE

    REAL(KIND=8) :: rho_in, u_in, v_in, p_in
    REAL(KIND=8) :: rhoG, rhoP, uG, uP, vG, vP, ekG, ekP
    REAL(KIND=8) :: rhoG_ov_rhoP, uP_rhoG_ov_rhoP
    REAL(KIND=8), DIMENSION(nb_eq) :: inlet_data, prim
    TYPE(boundary) :: bound

    INTEGER, INTENT(IN) :: id                              !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys       !< conservative variables of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data    !< physical properties of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost     !< conservative variables of ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data  !< physical properties of ghost cell
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: ju_ghost  !< ghost state conservative variable Jacobian

    ! Inlet values of density, velocity components and pressure
    bound      = boundary_data(id)
    inlet_data = get_boundary_inlet(nb_eq, bound)

    rho_in = inlet_data(1)
    u_in   = inlet_data(2)
    v_in   = inlet_data(3)
    p_in   = inlet_data(4)

    ! Density, velocity components and kinetic energy per unit mass of physical state
    rhoP = phys_data(pos_rho_cell)
    uP   = phys_data(pos_u_cell)
    vP   = phys_data(pos_v_cell)
    ekP  = phys_data(pos_ek_cell)

    ! Primitive variables (linear extrapolation used for density, velocity and pressure)
    prim(1) = 2.d0*rho_in - rhoP
    prim(2) = 2.d0*u_in   - uP
    prim(3) = 2.d0*v_in   - vP
    prim(4) = 2.d0*p_in   - phys_data(pos_pres_cell) 

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost, ghost_data)

    ! Evaluation of the matrix dUg/dUp for implicti boundary condition
    ! Density, velocity components and kinetic energy per unit mass of ghost state
    rhoG = ghost_data(pos_rho_cell)
    uG   = ghost_data(pos_u_cell)
    vG   = ghost_data(pos_v_cell)
    ekG  = ghost_data(pos_ek_cell) 
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

  END SUBROUTINE sup_in_1D_SL_Impl
!------------------------------------------------------------------------------!
!> This subroutine applies a slip wall boundary condition for 1D stagnation line 
!! calorically perfect gas flows. 
  SUBROUTINE slip_wall_1D_SL_Impl (id, phys_data, u_phys, ghost_data, u_ghost, ju_ghost)

    USE mod_general_data,          ONLY: pos_u_cell, pos_v_cell, pos_pres_cell, pos_ek_cell, pos_rho_cell, nb_eq
    USE mod_domain_boundary,       ONLY: boundary_data, get_boundary_rhoin, get_boundary_pin
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim 
 
    IMPLICIT NONE

    REAL(KIND=8) :: rhoG, rhoP, uG, uP, vG, vP, ekG, ekP
    REAL(KIND=8) :: rhoP_ov_rhoG, uP_rhoP_ov_rhoG
    REAL(KIND=8), DIMENSION(nb_eq) :: prim

    INTEGER, INTENT(IN) :: id                              !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys       !< conservative variables of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data    !< physical properties of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost     !< conservative variables of ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data  !< physical properties of ghost cell
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: ju_ghost  !< ghost state conservative variable Jacobian

    ! Density, velocity components and kinetic energy per unit mass of physical state 
    rhoP = phys_data(pos_rho_cell)
    uP   = phys_data(pos_u_cell)
    vP   = phys_data(pos_v_cell)
    ekP  = phys_data(pos_ek_cell)

    ! Ghost state
    ! Primitive variables
    prim(1) = rhoP
    prim(2) = - uP
    prim(3) = vP
    prim(4) = phys_data(pos_pres_cell)

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost, ghost_data)

    ! Evaluation of the matrix dUg/dUp for implicti boundary condition
    ! Density, velocity components and kinetic energy per unit mass of ghost state
    rhoG = ghost_data(pos_rho_cell)
    uG   = ghost_data(pos_u_cell)
    vG   = ghost_data(pos_v_cell)
    ekG  = ghost_data(pos_ek_cell) 
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

  END SUBROUTINE slip_wall_1D_SL_Impl
!------------------------------------------------------------------------------!
!> This subroutine applies a no slip wall boundary condition for 1D stagnation line 
!! calorically perfect gas flows. 
  SUBROUTINE no_slip_iso_Twall_1D_SL_Impl (id, phys_data, u_phys, ghost_data, u_ghost, ju_ghost)

    USE mod_general_data,          ONLY: pos_u_cell, pos_v_cell, pos_pres_cell, pos_T_cell, pos_rho_cell, pos_ek_cell, & 
                                       & nb_eq, R_gas, gamma
    USE mod_domain_boundary,       ONLY: boundary, boundary_data, get_boundary_Twall
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim 
 
    IMPLICIT NONE

    REAL(KIND=8), PARAMETER :: Tmin = 400.d0
    REAL(KIND=8) :: Twall
    REAL(KIND=8) :: rhoP, uP, vP, p, TP
    REAL(KIND=8) :: rhoG, pG, uG, vG, TG
    REAL(KIND=8) :: rhoG_ov_rhoP, ekP, EG, ov_TG
    REAL(KIND=8) :: fac1, fac2, fac3, fac4, fac5
    REAL(KIND=8) :: gamma_minus1 
    REAL(KIND=8), DIMENSION(nb_eq) :: prim
    TYPE(boundary) :: bound

    INTEGER, INTENT(IN) :: id                              !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys       !< conservative variables of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data    !< physical properties of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost     !< conservative variables of ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data  !< physical properties of ghost cell
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: ju_ghost  !< ghost state conservative variable Jacobian

    ! Wall temperature
    bound = boundary_data(id)
    Twall = get_boundary_Twall(bound)

    ! Pressure, density, kinetic energy per unit mass, velocity components and temperature of physical state  
    p    = phys_data(pos_pres_cell)
    rhoP = phys_data(pos_rho_cell)
    uP   = phys_data(pos_u_cell)
    vP   = phys_data(pos_v_cell)
    ekP  = phys_data(pos_ek_cell)
    TP   = phys_data(pos_T_cell)

    ! Ghost state
    ! Temperature
    TG = MAX(2.d0*Twall - TP,Tmin)

    ! Primitive variables
    rhoG = p/(R_gas*TG)
    uG   = - uP
    vG   = - vP
    pG   = p

    prim(1) = rhoG
    prim(2) = uG
    prim(3) = vG
    prim(4) = pG
    
    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost, ghost_data)
        
    ! Evaluation of the matrix dUg/dUp for implicit boundary condition
    ! Common factors
    gamma_minus1 = gamma - 1.d0
    EG    = u_ghost(4)/rhoG
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

  END SUBROUTINE no_slip_iso_Twall_1D_SL_Impl 
!------------------------------------------------------------------------------!
