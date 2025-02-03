!------------------------------------------------------------------------------!
!> This subroutine applies a supersonic inlet boundary condition for 1D stagnation line 
!! calorically perfect gas flows.
  SUBROUTINE sup_in_1D_SL_Expl (id, phys_data, u_phys, ghost_data, u_ghost)

    USE mod_general_data,          ONLY: pos_u_cell, pos_v_cell, pos_pres_cell, pos_rho_cell, nb_eq 
    USE mod_domain_boundary,       ONLY: boundary, boundary_data, get_boundary_inlet 
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim

    IMPLICIT NONE

    REAL(KIND=8) :: rho_in, u_in, v_in, p_in
    REAL(KIND=8), DIMENSION(nb_eq) :: inlet_data, prim
    TYPE(boundary) :: bound

    INTEGER, INTENT(IN) :: id                              !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys       !< conservative variables of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data    !< physical properties of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost     !< conservative variables of ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data  !< physical properties of ghost cell

    ! Inlet values of density, velocity components and pressure
    bound      = boundary_data(id)
    inlet_data = get_boundary_inlet(nb_eq, bound)

    rho_in = inlet_data(1)
    u_in   = inlet_data(2)
    v_in   = inlet_data(3)
    p_in   = inlet_data(4)

    ! Primitive variables (linear extrapolation used for density, velocity and pressure)
    prim(1) = 2.d0*rho_in - phys_data(pos_rho_cell)
    prim(2) = 2.d0*u_in   - phys_data(pos_u_cell)
    prim(3) = 2.d0*v_in   - phys_data(pos_v_cell)
    prim(4) = 2.d0*p_in   - phys_data(pos_pres_cell) 

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost, ghost_data)

  END SUBROUTINE sup_in_1D_SL_Expl
!------------------------------------------------------------------------------!
!> This subroutine applies a slip wall boundary condition for 1D stagnation line 
!! calorically perfect gas flows. 
  SUBROUTINE slip_wall_1D_SL_Expl (id, phys_data, u_phys, ghost_data, u_ghost)

    USE mod_general_data,          ONLY: pos_u_cell, pos_v_cell, pos_pres_cell, pos_rho_cell, nb_eq
    USE mod_domain_boundary,       ONLY: boundary_data, get_boundary_rhoin, get_boundary_pin
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim 
 
    IMPLICIT NONE

    REAL(KIND=8) :: rho, u, v, p
    REAL(KIND=8), DIMENSION(nb_eq) :: prim
 
    INTEGER, INTENT(IN) :: id                              !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys       !< conservative variables of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data    !< physical properties of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost     !< conservative variables of ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data  !< physical properties of ghost cell

    ! Density, velocity components and pressure of physical state 
    rho = phys_data(pos_rho_cell)
    u   = phys_data(pos_u_cell)
    v   = phys_data(pos_v_cell)
    p   = phys_data(pos_pres_cell)

    ! Ghost state
    ! Primitive variables
    prim(1) = rho
    prim(2) = - u
    prim(3) = v
    prim(4) = p 

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost, ghost_data)
 
  END SUBROUTINE slip_wall_1D_SL_Expl
!------------------------------------------------------------------------------!
!> This subroutine applies an isothermal no slip wall boundary condition for 1D stagnation line 
!! calorically perfect gas flows. 
  SUBROUTINE no_slip_iso_Twall_1D_SL_Expl (id, phys_data, u_phys, ghost_data, u_ghost)

    USE mod_general_data,          ONLY: pos_u_cell, pos_v_cell, pos_pres_cell, pos_T_cell, nb_eq, R_gas
    USE mod_domain_boundary,       ONLY: boundary, boundary_data, get_boundary_Twall
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim 
 
    IMPLICIT NONE

    REAL(KIND=8) :: u, v, p, T, Tghost, Twall
    REAL(KIND=8), DIMENSION(nb_eq) :: prim
    TYPE(boundary) :: bound

    INTEGER, INTENT(IN) :: id                              !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys       !< conservative variables of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data    !< physical properties of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost     !< conservative variables of ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data  !< physical properties of ghost cell

    ! Wall temperature
    bound = boundary_data(id)
    Twall = get_boundary_Twall(bound)
    
    ! Pressure, velocity components and temperature of physical state  
    p = phys_data(pos_pres_cell)
    u = phys_data(pos_u_cell)
    v = phys_data(pos_v_cell)
    T = phys_data(pos_T_cell)

    ! Ghost state
    ! Temperature
    Tghost = MAX(2.d0*Twall - T,0.9d0*Twall)

    ! Primitive variables
    prim(1) = p/(R_gas*Tghost)
    prim(2) = - u
    prim(3) = - v
    prim(4) = p

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost, ghost_data)

  END SUBROUTINE no_slip_iso_Twall_1D_SL_Expl
!------------------------------------------------------------------------------!
