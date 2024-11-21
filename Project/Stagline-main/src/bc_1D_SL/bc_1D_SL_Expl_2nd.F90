!------------------------------------------------------------------------------!
!> This subroutine applies a supersonic inlet boundary condition for 1D stagnation line nonequilibrium flows.
  SUBROUTINE sup_in_1D_SL_Expl_2nd (id, phys_data1, u_phys1, phys_data2, u_phys2, & 
                                 &  ghost_data1, u_ghost1, ghost_data2, u_ghost2)

    USE mod_general_data,          ONLY: pos_u_cell, pos_v_cell, pos_pres_cell, pos_rho_cell, nb_eq 
    USE mod_domain_boundary,       ONLY: boundary, boundary_data, get_boundary_inlet 
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim

    IMPLICIT NONE

    REAL(KIND=8) :: ov_3
    REAL(KIND=8) :: rho_in, u_in, v_in, p_in
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

  END SUBROUTINE sup_in_1D_SL_Expl_2nd
!------------------------------------------------------------------------------!
!> This subroutine applies a slip wall boundary condition for 1D stagnation line nonequilibrium flows. 
  SUBROUTINE slip_wall_1D_SL_Expl_2nd (id, phys_data1, u_phys1, phys_data2, u_phys2,   & 
                                    &  ghost_data1, u_ghost1, ghost_data2, u_ghost2)

    USE mod_general_data,          ONLY: pos_u_cell, pos_v_cell, pos_pres_cell, pos_rho_cell, nb_eq
    USE mod_domain_boundary,       ONLY: boundary_data, get_boundary_rhoin, get_boundary_pin
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim 
 
    IMPLICIT NONE

    REAL(KIND=8) :: rho1, u1, v1, p1, p2
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
 
    ! Density, velocity components and pressure of physical state 
    rho1 = phys_data1(pos_rho_cell)
    u1   = phys_data1(pos_u_cell)
    v1   = phys_data1(pos_v_cell)
    p1   = phys_data1(pos_pres_cell)
    p2   = phys_data2(pos_pres_cell)

    ! Ghost state 1
    ! Primitive variables 
    prim(1) = rho1
    prim(2) = - u1
    prim(3) = v1
    prim(4) = (4.d0*p1 - p2)/3.d0 
   
    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost1, ghost_data1)
 
    ! Ghost state 2
    u_ghost2    = u_ghost1
    ghost_data2 = ghost_data1

  END SUBROUTINE slip_wall_1D_SL_Expl_2nd
!------------------------------------------------------------------------------!
!> This subroutine applies an isothrmal no slip wall boundary condition for 1D stagnation line 
!! calorically perfect gas flows. 
  SUBROUTINE no_slip_iso_Twall_1D_SL_Expl_2nd (id, phys_data1, u_phys1, phys_data2, u_phys2,   & 
                                           &  ghost_data1, u_ghost1, ghost_data2, u_ghost2)

    USE mod_general_data,          ONLY: pos_u_cell, pos_v_cell, pos_pres_cell, pos_rho_cell, pos_T_cell, nb_eq, R_gas
    USE mod_domain_boundary,       ONLY: boundary, boundary_data, get_boundary_Twall
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim 
 
    IMPLICIT NONE

    REAL(KIND=8) :: u1, v1, p1, T1, T2, Twall
    REAL(KIND=8) :: TG1, TG2
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

    ! Wall temperature
    bound = boundary_data(id)
    Twall = get_boundary_Twall(bound)

    ! Pressure of physical state 1
    p1 = phys_data1(pos_pres_cell)

    ! Temperature of physical states 1 and 2
    T1 = phys_data1(pos_T_cell)
    T2 = phys_data2(pos_T_cell)

    ! Velocity components of physical state 1
    u1 = phys_data1(pos_u_cell)
    v1 = phys_data1(pos_v_cell)

    ! Ghost state 1
    ! Primitive variables
    TG1 = MAX(Twall - 0.5d0*(T2 - T1),0.9d0*Twall)

    prim(1) = p1/(R_gas*TG1)
    prim(2) = - u1
    prim(3) = - v1
    prim(4) = p1

    ! Compute the ghost state conservative variables and physical properties from primitive variables
    CALL get_cons_phys_from_prim(prim, u_ghost1, ghost_data1)

    ! Ghost state 2
    u_ghost2 = u_ghost1
    ghost_data2 = ghost_data1

  END SUBROUTINE no_slip_iso_Twall_1D_SL_Expl_2nd
!------------------------------------------------------------------------------!
