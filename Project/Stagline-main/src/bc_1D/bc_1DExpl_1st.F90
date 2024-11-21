!------------------------------------------------------------------------------!
!> This subroutine applies a supersonic inlet boundary condition for 1D calorically perfct gas flows.
  SUBROUTINE sup_in_1DExpl (id, phys_data, u_phys, ghost_data, u_ghost)

    USE mod_general_data,          ONLY: pos_u_cell, pos_pres_cell, pos_rho_cell, nb_eq 
    USE mod_domain_boundary,       ONLY: boundary, boundary_data, get_boundary_inlet 
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim

    IMPLICIT NONE

    REAL(KIND=8) :: rho_in, u_in, p_in
    REAL(KIND=8), DIMENSION(nb_eq) :: inlet_data, prim
    TYPE(boundary) :: bound

    INTEGER, INTENT(IN) :: id                              !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys       !< conservative variables of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data    !< physical properties of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost     !< conservative variables of ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data  !< physical properties of ghost cell

    ! Inlet values of density, velocity and pressure
    bound      = boundary_data(id)
    inlet_data = get_boundary_inlet(nb_eq, bound)

    rho_in = inlet_data(1)
    u_in   = inlet_data(2)
    p_in   = inlet_data(3)

    ! Primitive variables (linear extrapolation used for density, velocity and pressure)
    prim(1) = 2.d0*rho_in - phys_data(pos_rho_cell)
    prim(2) = 2.d0*u_in   - phys_data(pos_u_cell)
    prim(3) = 2.d0*p_in   - phys_data(pos_pres_cell) 

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost, ghost_data)

  END SUBROUTINE sup_in_1DExpl
!------------------------------------------------------------------------------!
!> This subroutine applies a supersonic outlet boundary condition.
  SUBROUTINE sup_out_Expl (id, phys_data, u_phys, ghost_data, u_ghost)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: id                              !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys       !< conservative variables of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data    !< physical properties of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost     !< conservative variables of ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data  !< physical properties of ghost cell
 
    ! Ghost state conservative variables and physical properties
    u_ghost    = u_phys
    ghost_data = phys_data 

  END SUBROUTINE sup_out_Expl 
!------------------------------------------------------------------------------!
!> This subroutine applies a subsonic inlet boundary condition for 1D calorically perfect gas flows. 
  SUBROUTINE sub_in_rhoin_pin_1DExpl  (id, phys_data, u_phys, ghost_data, u_ghost)

    USE mod_general_data,          ONLY: pos_u_cell, pos_pres_cell, pos_rho_cell, nb_eq
    USE mod_domain_boundary,       ONLY: boundary_data, get_boundary_rhoin, get_boundary_pin
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim 
 
    IMPLICIT NONE

    REAL(KIND=8) :: rho_in, p_in
    REAL(KIND=8), DIMENSION(nb_eq) :: prim

    INTEGER, INTENT(IN) :: id                              !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys       !< conservative variables of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data    !< physical properties of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost     !< conservative variables of ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data  !< physical properties of ghost cell

    ! Inlet values of density and pressure
    rho_in = get_boundary_rhoin (boundary_data (id))
    p_in   = get_boundary_pin (boundary_data (id))  
    
    ! Primitive variables (linear extrapolation used for density and pressure)
    prim(1) = 2.d0*rho_in - phys_data(pos_rho_cell)
    prim(2) = phys_data(pos_u_cell)
    prim(3) = 2.d0*p_in   - phys_data(pos_pres_cell) 

    ! Compute the ghost state conservative variables and physical properties from primitive variables
    CALL get_cons_phys_from_prim(prim, u_ghost, ghost_data)

  END SUBROUTINE sub_in_rhoin_pin_1DExpl
!------------------------------------------------------------------------------!
!> This subroutine applies a subsonic outlet boundary condition for 1D calorically perfect gas flows. 
  SUBROUTINE sub_out_pout_1DExpl (id, phys_data, u_phys, ghost_data, u_ghost)

    USE mod_general_data,          ONLY: pos_u_cell, pos_pres_cell, pos_rho_cell, nb_eq
    USE mod_domain_boundary,       ONLY: boundary_data, get_boundary_pout
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim

    IMPLICIT NONE

    REAL(KIND=8) :: pout 
    REAL(KIND=8), DIMENSION(nb_eq) :: prim

    INTEGER, INTENT(IN) :: id                              !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys       !< conservative variables of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data    !< physical properties of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost     !< conservative variables of ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data  !< physical properties of ghost cell

    ! Imposed outlet pressure 
    pout = get_boundary_pout (boundary_data (id))

    ! Primitive variables (linear extrapolation used for static pressure)
    prim(1) = phys_data(pos_rho_cell) 
    prim(2) = phys_data(pos_u_cell)
    prim(3) = 2.d0*pout - phys_data(pos_pres_cell)

    ! Compute the ghost state conservative variables and physical properties from primitive variables
    CALL get_cons_phys_from_prim(prim, u_ghost, ghost_data)

  END SUBROUTINE sub_out_pout_1DExpl
!------------------------------------------------------------------------------!
