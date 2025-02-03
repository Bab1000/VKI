!------------------------------------------------------------------------------!
!> This subroutine applies a supersonic inlet boundary condition in case of use
!! of two ghost cells for 1D calorically perfect gas flows. A linear extrapolation is used for conservative variables. 
!! Physical data are then computed from the conservative variables. The first ghost cells is that closer to the domain boundary.
  SUBROUTINE sup_in_1DImpl_2nd (id, phys_data1, u_phys1, phys_data2, u_phys2, & 
                              &  ghost_data1, u_ghost1, ghost_data2, u_ghost2, &
                              &  ju_ghost)

    USE mod_general_data,          ONLY: pos_u_cell, pos_pres_cell, pos_rho_cell, nb_eq
    USE mod_domain_boundary,       ONLY: boundary_data, get_boundary_inlet, boundary
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim

    IMPLICIT NONE

    REAL(KIND=8) :: rho_in, u_in, p_in
    REAL(KIND=8) :: rho1, rho2, u1, u2, p1, p2
    REAL(KIND=8) :: tmp1, tmp2, tmp3
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

    ! Inlet values of density, velocity and pressure
    bound      = boundary_data(id)
    inlet_data = get_boundary_inlet(nb_eq, bound)
    
    rho_in = inlet_data(1)
    u_in   = inlet_data(2)
    p_in   = inlet_data(3) 

    ! Density, velocity and pressure of physical states
    rho1  = phys_data1(pos_rho_cell)
    u1    = phys_data1(pos_u_cell)
    p1    = phys_data1(pos_pres_cell)

    rho2  = phys_data2(pos_rho_cell)
    u2    = phys_data2(pos_u_cell)
    p2    = phys_data2(pos_pres_cell)

    ! Ghost state 1
    ! Primitive variables (linear extrapolation for all variables)
    prim(1) = 2.d0*rho_in - rho1
    prim(2) = 2.d0*u_in   - u1
    prim(3) = 2.d0*p_in   - p1 

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost1, ghost_data1)

    ! Ghost state 2 
    ! Primitive variables (linear extrapolation for all variables)
    prim(1) = 2.d0*rho_in - rho2
    prim(2) = 2.d0*u_in   - u2
    prim(3) = 2.d0*p_in   - p2 

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost2, ghost_data2)

    ! Evaluation of dUg/dUp matrix for implicit of boundary condition 
    ! (the Jacobian matrix is computed as for a first order scheme for sake of robustness)
    ! Common factors
    tmp1 = 2.d0*u_in
    tmp2 = u1/rho1
    tmp3 = rho_in*u_in    

    ! First column 
    ju_ghost(1,1) = - 1.d0    
    ju_ghost(2,1) = - tmp1 + 2.d0*rho_in*tmp2    
    ju_ghost(3,1) = 4.d0*tmp3*tmp2 - tmp1*u_in + rho_in*tmp2*u_in

    ! Second column 
    ju_ghost(1,2) = 0.d0    
    ju_ghost(2,2) = - 2.d0*rho_in/rho1 + 1.d0    
    ju_ghost(3,2) = - 4.d0*tmp3/rho1 + tmp1 + 2.d0*rho_in*tmp2    

    ! Third column 
    ju_ghost(1,3) = 0.d0    
    ju_ghost(2,3) = 0.d0    
    ju_ghost(3,3) = - 1.d0

  END  SUBROUTINE sup_in_1DImpl_2nd
!------------------------------------------------------------------------------!
!> This subroutine applies a supersonic outlet boundary condition in case of use
!! of two ghost cells for 1D calorically perfect gas flows. A linear extrapolation is used for conservative variables. 
!! Physical data are then computed from the conservative variables. The first ghost cells is that closer to the domain boundary.
  SUBROUTINE sub_out_pout_1DImpl_2nd (id, phys_data1, u_phys1, phys_data2, u_phys2, & 
                              &  ghost_data1, u_ghost1, ghost_data2, u_ghost2, &
                              &  ju_ghost)

    USE mod_general_data,          ONLY: pos_u_cell, pos_pres_cell, pos_rho_cell, nb_eq
    USE mod_domain_boundary,       ONLY: boundary_data, get_boundary_pout
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim

    IMPLICIT NONE

    REAL(KIND=8) :: p_out
    REAL(KIND=8) :: rho1, rho2, u1, u2, p1, p2
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

    ! Outlet value of pressure
    p_out = get_boundary_pout (boundary_data (id))
    
    ! Density, velocity and pressure of physical states
    rho1  = phys_data1(pos_rho_cell)
    u1    = phys_data1(pos_u_cell)
    p1    = phys_data1(pos_pres_cell)

    rho2  = phys_data2(pos_rho_cell)
    u2    = phys_data2(pos_u_cell)
    p2    = phys_data2(pos_pres_cell)

    ! Ghost cell 1
    ! Primitive variables (linear extrapolation is used for pressure)
    prim(1) = 2.d0*rho1  - rho2
    prim(2) = 2.d0*u1    - u2
    prim(3) = 2.d0*p_out - p1

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost1, ghost_data1)

    ! Ghost cell 2
    ! Primitive variables (linear extrapolation is used for pressure)
    prim(1) = 3.d0*rho1  - 2.d0*rho2
    prim(2) = 3.d0*u1    - 2.d0*u2
    prim(3) = 2.d0*p_out - p2

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost2, ghost_data2)

    ! Evaluation of dUg/dUp matrix for implicit treatment 
    ! (the Jacobian matrix is computed as for a first order scheme for sake of robustness)
    ! First column 
    ju_ghost(1,1) = 1.d0
    ju_ghost(2,1) = 0.d0
    ju_ghost(3,1) = - u1*u1

    ! Second column 
    ju_ghost(1,2) = 0.d0
    ju_ghost(2,2) = 1.d0
    ju_ghost(3,2) = -2.d0*u1 

    ! Third column
    ju_ghost(1,3) = 0.d0
    ju_ghost(2,3) = 0.d0
    ju_ghost(3,3) = -1.d0 

  END SUBROUTINE sub_out_pout_1DImpl_2nd
!------------------------------------------------------------------------------!
!> This subroutine applies a supersonic outlet boundary condition in case of use
!! of two ghost cells for 1D calorically perfect gas flows. A linear extrapolation is used for conservative variables. 
!! Physical data are then computed from the conservative variables. The first ghost cells is that closer to the domain boundary.
  SUBROUTINE sup_out_1DImpl_2nd (id, phys_data1, u_phys1, phys_data2, u_phys2, & 
                              &  ghost_data1, u_ghost1, ghost_data2, u_ghost2, &
                              &  ju_ghost)

    USE mod_general_data,          ONLY: pos_u_cell, pos_pres_cell, pos_rho_cell, nb_eq
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim

    IMPLICIT NONE

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

    ! Ghost state 1
    ! Primitive variables (linear extrapolation is used)
    prim(1) = 2.d0*phys_data1(pos_rho_cell) - phys_data2(pos_rho_cell)
    prim(2) = phys_data1(pos_u_cell) 
    prim(3) = 2.d0*phys_data1(pos_pres_cell) - phys_data2(pos_pres_cell)

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost1, ghost_data1)

    ! Ghost state 2
    ! Primitive variables (linear extrapolation is used)
    prim(1) = 2.d0*ghost_data1(pos_rho_cell) - phys_data1(pos_rho_cell)
    prim(2) = phys_data1(pos_u_cell) 
    prim(3) = 2.d0*ghost_data1(pos_pres_cell) - phys_data1(pos_pres_cell)

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost2, ghost_data2)

    ! Evaluation of dUg/dUp matrix for implicit boundary condition 
    ! (the Jacobian matrix is computed as for a first order scheme for sake of robustness)
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

  END SUBROUTINE sup_out_1DImpl_2nd 
!------------------------------------------------------------------------------!
!> This subroutine applies a subsonic inlet boundary condition for 1D calorically perfect gas flows
!! in case of use of two ghost cells. A parabolic extrapolation is used for primitive variables that are being 
!! imposed (static density and pressure) while a linear extrapolation is used for the velocity. 
!! Conservative variables are then computed from primitive variables. The first ghost cells is that close to the domain boundary.
  SUBROUTINE sub_in_rhoin_pin_1DImpl_2nd (id, phys_data1, u_phys1, phys_data2, u_phys2,  & 
                                       &  ghost_data1, u_ghost1, ghost_data2, u_ghost2,  & 
                                       &  ju_ghost)

    USE mod_general_data,          ONLY: pos_u_cell, pos_pres_cell, pos_rho_cell, nb_eq
    USE mod_domain_boundary,       ONLY: boundary_data, get_boundary_rhoin, get_boundary_pin
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim

    IMPLICIT NONE

    REAL(KIND=8) :: tmp1, tmp2
    REAL(KIND=8) :: rho_in, p_in
    REAL(KIND=8) :: rho1, rho2, u1, u2, p1, p2
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

    ! Inlet values of density and pressure
    rho_in = get_boundary_rhoin (boundary_data (id))
    p_in   = get_boundary_pin (boundary_data (id))  
    
    ! Density, velocity and pressure of physical states
    rho1 = phys_data1(pos_rho_cell)
    u1   = phys_data1(pos_u_cell)
    p1   = phys_data1(pos_pres_cell)

    rho2 = phys_data2(pos_rho_cell)
    u2   = phys_data2(pos_u_cell)
    p2   = phys_data2(pos_pres_cell)

    ! Ghost state 1
    ! Primitive variables (linear extrapolation is used for pressure)
    prim(1) = 2.d0*rho_in - rho1
    prim(2) = 2.d0*u1 - u2
    prim(3) = 2.d0*p_in - p1

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost1, ghost_data1) 

    ! Ghost cell 2
    prim(1) = 2.d0*rho_in - rho2
    prim(2) = 3.d0*u1 - 2.d0*u2
    prim(3) = 2.d0*p_in - p2

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost2, ghost_data2)

    ! Evaluation of dUg/dUp matrix for implicit boundary condition 
    ! (the Jacobian matrix is computed as for a first order scheme for sake of robustness) 
    ! Useful quantities
    tmp1 = 2.d0*u1/rho1
    tmp2 = rho_in*tmp1

    ! First column 
    ju_ghost(1,1) = - 1.d0 
    ju_ghost(2,1) = - tmp2
    ju_ghost(3,1) = - u1*tmp2 

    ! Second column
    ju_ghost(1,2) = 0.d0
    ju_ghost(2,2) = 2.d0*rho_in/rho1 - 1.d0
    ju_ghost(3,2) = tmp2

    ! Third column 
    ju_ghost(1,3) = 0.d0  
    ju_ghost(2,3) = 0.d0  
    ju_ghost(3,3) = - 1.d0  

  END SUBROUTINE sub_in_rhoin_pin_1DImpl_2nd 
!------------------------------------------------------------------------------!
