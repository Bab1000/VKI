!------------------------------------------------------------------------------!
!> This subroutine applies a supersonic inlet boundary conditions for 1D stagnation line nonequilibrium flows.
  SUBROUTINE sup_in_neq_1D_SL_Impl(id, phys_data, u_phys, ghost_data, u_ghost, ju_ghost)

    USE mod_general_data,          ONLY: nb_ns, nb_temp, nb_eq, pos_u, pos_v, pos_T, pos_u_cell, & 
                                       & pos_v_cell, pos_T_cell                                         
    USE mod_domain_boundary,       ONLY: boundary, boundary_data, get_boundary_inlet 
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim
 
    IMPLICIT NONE

    INTEGER :: i, j
    REAL(KIND=8) :: rho_in, u_in, v_in, p_in
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

    ! Primitive variables (linear extrapolation used for density, velocity and temperatures)
    ! Species densities
    DO i = 1,nb_ns 
       prim(i) = 2.d0*inlet_data(i) - u_phys(i)
    ENDDO

    ! Velocity components 
    prim(pos_u) = 2.d0*inlet_data(pos_u) - phys_data(pos_u_cell)
    prim(pos_v) = 2.d0*inlet_data(pos_v) - phys_data(pos_v_cell)

    ! Temperature(s)
    DO i = 1,nb_temp
       prim(pos_T + i - 1) = 2.d0*inlet_data(pos_T + i - 1) - phys_data(pos_T_cell + i - 1)
    ENDDO

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost, ghost_data)

    !! APPROXIMATION !!
    DO j = 1,nb_eq 

       DO i = 1,j - 1
          ju_ghost(i,j) = 0.d0 
       ENDDO

       ju_ghost(i,i) = - 1.d0

       DO i = j + 1,nb_eq 
          ju_ghost(i,j) = 0.d0
       ENDDO

    ENDDO

    !PRINT*
    !WRITE(*,'(A)')'In "sup_in_neq_1D_SL_Impl", boundary conditions Jacobian to be implemented!'
    !PRINT*
    !STOP

  END SUBROUTINE sup_in_neq_1D_SL_Impl 
!------------------------------------------------------------------------------!
