!------------------------------------------------------------------------------!
!> This subroutine applies a supersonic outlet boundary condition for 1D calorically perfect gas flows.
  SUBROUTINE sup_out_1DImpl (id, phys_data, u_phys, ghost_data, u_ghost, ju_ghost)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: id                              !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys       !< conservative variables of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data    !< physical properties of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost     !< conservative variables of ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data  !< physical properties of ghost cell
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: ju_ghost  !< ghost state conservative variable Jacobian    

    ! Ghost state conservative variables and physical properties
    u_ghost    = u_phys
    ghost_data = phys_data 

    ! Evaluation of dUg/dUp matrix for implicit treatment of boundary condition  
    ! First column 
    ju_ghost(1,1) = 1.d0
    ju_ghost(2,1) = 0.d0
    ju_ghost(3,1) = 0.d0

    ! Second column 
    ju_ghost(1,2) = 0.d0
    ju_ghost(2,2) = 1.d0
    ju_ghost(3,2) = 0.d0

    ! Third column 
    ju_ghost(1,3) = 0.d0
    ju_ghost(2,3) = 0.d0
    ju_ghost(3,3) = 1.d0 

  END SUBROUTINE sup_out_1DImpl
!------------------------------------------------------------------------------!
!> This subroutine applies a supersonic inlet boundary condition for 1D calorically perfect gas flows.
  SUBROUTINE sup_in_1DImpl (id, phys_data, u_phys, ghost_data, u_ghost, ju_ghost)

    USE mod_general_data,          ONLY: nb_eq, pos_u_cell, pos_pres_cell, pos_c_cell, pos_rho_cell
    USE mod_domain_boundary,       ONLY: boundary_data, boundary, get_boundary_inlet
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim

    IMPLICIT NONE

    REAL(KIND=8) :: rho_in, u_in, p_in
    REAL(KIND=8) :: rho_phys, vel_phys, p_phys
    REAL(KIND=8) :: tmp1, tmp2, tmp3
    REAL(KIND=8), DIMENSION(nb_eq) :: inlet_data, prim
    TYPE(boundary) :: bound
 
    INTEGER, INTENT(IN) :: id                              !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys       !< conservative variables of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data    !< physical properties of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost     !< conservative variables of ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data  !< physical properties of ghost cell
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: ju_ghost  !< ghost state conservative variable Jacobian

    ! Inlet values of density, velocity and pressure
    bound      = boundary_data(id)
    inlet_data = get_boundary_inlet(nb_eq, bound)

    rho_in = inlet_data(1)
    u_in   = inlet_data(2)
    p_in   = inlet_data(3)

    rho_phys = phys_data(pos_rho_cell)
    vel_phys = phys_data(pos_u_cell)
    p_phys   = phys_data(pos_pres_cell) 

    ! Primitive variables (linear extrapolation is used for all variables)
    prim(1) = 2.d0*rho_in - rho_phys
    prim(2) = 2.d0*u_in   - vel_phys
    prim(3) = 2.d0*p_in   - p_phys 

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost, ghost_data)
    
    ! Evaluation of dUg/dUp matrix for implicit boundary condition 
    ! Common factors
    tmp1 = 2.d0*u_in
    tmp2 = vel_phys/rho_phys
    tmp3 = rho_in*u_in    

    ! First column 
    ju_ghost(1,1) = - 1.d0    
    ju_ghost(2,1) = - tmp1 + 2.d0*rho_in*tmp2    
    ju_ghost(3,1) = 4.d0*tmp3*tmp2 - tmp1*u_in + rho_in*tmp2*u_in

    ! Second column 
    ju_ghost(1,2) = 0.d0    
    ju_ghost(2,2) = - 2.d0*rho_in/rho_phys + 1.d0    
    ju_ghost(3,2) = - 4.d0*tmp3/rho_phys + tmp1 + 2.d0*rho_in*tmp2    

    ! Third column 
    ju_ghost(1,3) = 0.d0    
    ju_ghost(2,3) = 0.d0    
    ju_ghost(3,3) = - 1.d0    

  END SUBROUTINE sup_in_1DImpl 
!------------------------------------------------------------------------------!
!> This subroutine applies a subsonic inlet boundary condition for 1D calorically perfect gas flows. 
  SUBROUTINE sub_in_rhoin_pin_1DImpl (id, phys_data, u_phys, ghost_data, u_ghost, ju_ghost)

    USE mod_general_data,          ONLY: pos_u_cell, pos_pres_cell, pos_rho_cell, nb_eq
    USE mod_domain_boundary,       ONLY: boundary_data, get_boundary_rhoin, get_boundary_pin
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim

    IMPLICIT NONE
 
    REAL(KIND=8) :: rho_in, p_in
    REAL(KIND=8) :: rho, rho_phys, vel_phys, p_phys
    REAL(KIND=8) :: tmp1, tmp2
    REAL(KIND=8), DIMENSION(nb_eq) :: prim

    INTEGER, INTENT(IN) :: id                              !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys       !< conservative variables of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data    !< physical properties of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost     !< conservative variables of ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data  !< physical properties of ghost cell
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: ju_ghost  !< ghost state conservative variable Jacobian

    ! Inlet values of density and pressure
    rho_in = get_boundary_rhoin (boundary_data (id))
    p_in   = get_boundary_pin (boundary_data (id))  
    
    rho_phys = phys_data(pos_rho_cell)
    vel_phys = phys_data(pos_u_cell)
    p_phys   = phys_data(pos_pres_cell)

    ! Ghost state primitive variables (linear extrapolation used for density and pressure)
    prim(1) = 2.d0*rho_in - rho_phys 
    prim(2) = vel_phys 
    prim(3) = 2.d0*p_in - p_phys  

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost, ghost_data)

    ! Evaluation of dUg/dUp matrix for implicit boundary condition 
    ! Useful quantities
    rho  = prim(1)
    tmp1 = rho/rho_phys
    tmp2 = vel_phys*(1.d0 + tmp1)

    ! First column 
    ju_ghost(1,1) = - 1.d0 
    ju_ghost(2,1) = - tmp2
    ju_ghost(3,1) = - vel_phys*tmp2 

    ! Second column
    ju_ghost(1,2) = 0.d0
    ju_ghost(2,2) = tmp1
    ju_ghost(3,2) = tmp2

    ! Third column 
    ju_ghost(1,3) = 0.d0  
    ju_ghost(2,3) = 0.d0  
    ju_ghost(3,3) = - 1.d0 
    
  END SUBROUTINE sub_in_rhoin_pin_1DImpl
!------------------------------------------------------------------------------!
!> This subroutine applies a subsonic outlet boundary conditions for 1D calorically perfect gas flows. 
  SUBROUTINE sub_out_pout_1DImpl (id, phys_data, u_phys, ghost_data, u_ghost, ju_ghost)

    USE mod_general_data,          ONLY: pos_u_cell, pos_pres_cell, pos_ek_cell, pos_rho_cell, nb_eq
    USE mod_domain_boundary,       ONLY: boundary_data, get_boundary_pout
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim

    IMPLICIT NONE

    REAL(KIND=8) :: pout
    REAL(KIND=8) :: rho_phys, p_phys, vel_phys
    REAL(KIND=8) :: u, ek
    REAL(KIND=8), DIMENSION(nb_eq) :: prim

    INTEGER, INTENT(IN) :: id                              !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys       !< conservative variables of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data    !< physical properties of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost     !< conservative variables of ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data  !< physical properties of ghost cell
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: ju_ghost  !< ghost state conservative variable Jacobian

    ! Imposed outlet pressure 
    pout = get_boundary_pout (boundary_data (id))

    ! Density, velocity and pressure of physical state
    rho_phys = phys_data(pos_rho_cell)
    vel_phys = phys_data(pos_u_cell)
    p_phys   = phys_data(pos_pres_cell)

    ! Primitive variabels (linear extrapolation used for static pressure)
    prim(1) = rho_phys 
    prim(2) = vel_phys
    prim(3) = 2.d0*pout - p_phys

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost, ghost_data)

    ! Evaluation of dUg/dUp matrix for implicit boundary condition 
    u  = ghost_data(pos_u_cell)
    ek = ghost_data(pos_ek_cell)

    ! First column 
    ju_ghost(1,1) = 1.d0
    ju_ghost(2,1) = 0.d0
    ju_ghost(3,1) = - 2.d0*ek

    ! Second column 
    ju_ghost(1,2) = 0.d0
    ju_ghost(2,2) = 1.d0
    ju_ghost(3,2) = -2.d0*u 

    ! Third column
    ju_ghost(1,3) = 0.d0
    ju_ghost(2,3) = 0.d0
    ju_ghost(3,3) = -1.d0

  END  SUBROUTINE sub_out_pout_1DImpl
!------------------------------------------------------------------------------~!
