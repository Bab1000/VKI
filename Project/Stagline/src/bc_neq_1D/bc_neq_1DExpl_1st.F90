!------------------------------------------------------------------------------!
!> This subroutine applies a supersonic inlet boundary conditions for 1D nonequilibrium flows.
  SUBROUTINE sup_in_neq_1DExpl(id, phys_data, u_phys, ghost_data, u_ghost)

    USE mod_general_data,             ONLY: nb_eq, nb_ns, nb_temp, pos_u_cell, pos_T_cell, pos_u 
    USE mod_domain_boundary,          ONLY: boundary_data, get_boundary_inlet, boundary 
    USE mod_function_pointer,         ONLY: get_cons_phys_from_prim

    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8), DIMENSION(nb_eq) :: inlet_data, prim
    TYPE(boundary) :: bound

    INTEGER, INTENT(IN) :: id                              !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys       !< conservative variables of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data    !< physical properties of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost     !< conservative variables of ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data  !< physical properties of ghost cell

    ! Inlet data (species densities, velocity and temperature(s))
    bound      = boundary_data(id)
    inlet_data = get_boundary_inlet (nb_eq, bound)

    ! Ghost state primitive variables (linear extrapolation used for all variables)
    ! Species densities
    DO i = 1,nb_ns 
       prim(i) = 2.d0*inlet_data(i) - u_phys(i)
    ENDDO

    ! Velocity
    prim(pos_u) = 2.d0*inlet_data(pos_u) - phys_data(pos_u_cell) 

    ! Temperature(s)
    DO i = 1,nb_temp
       prim(pos_u + i) = 2.d0*inlet_data(pos_u + i) - phys_data(pos_T_cell + i - 1)
    ENDDO

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost, ghost_data)
    
  END SUBROUTINE sup_in_neq_1DExpl 
!------------------------------------------------------------------------------!
!> This subroutine applies a subsonic inlet boundary condition for 1D nonequilibrium flows.
  SUBROUTINE sub_in_neq_rhoiin_T_1DExpl (id, phys_data, u_phys, ghost_data, u_ghost)

    USE mod_general_data,             ONLY: nb_ns, nb_temp, nb_eq, pos_u, pos_T, pos_u_cell, pos_T_cell
    USE mod_domain_boundary,          ONLY: boundary, boundary_data, get_boundary_rhoi, get_boundary_Tvec
    USE mod_function_pointer,         ONLY: get_cons_phys_from_prim

    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8), DIMENSION(nb_ns) :: rhoi_in
    REAL(KIND=8), DIMENSION(nb_temp) :: temp_in
    REAL(KIND=8), DIMENSION(nb_eq) :: prim
    TYPE(boundary) :: bound

    INTEGER, INTENT(IN) :: id                              !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys       !< conservative variables of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data    !< physical properties of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost     !< conservative variables of ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data  !< physical properties of ghost cell
  
    ! Species densities at inlet
    bound   = boundary_data(id)
    rhoi_in = get_boundary_rhoi (nb_ns, bound)

    ! Temperature vector at inlet
    temp_in = get_boundary_Tvec (nb_temp, bound)
    
    ! Ghost state primivite variables (linear extrapolation used for species densities and temperature(s))
    ! Species densities
    DO i = 1,nb_ns
       prim(i) = 2.d0*rhoi_in(i) - u_phys(i) 
    ENDDO

    ! Velocity 
    prim(pos_u) = phys_data(pos_u_cell)

    ! Temperature(s)
    DO i = 1,nb_temp
       prim(pos_T + i - 1) = 2.d0*temp_in(i) - phys_data(pos_T_cell + i - 1)
    ENDDO

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost, ghost_data)

  END SUBROUTINE sub_in_neq_rhoiin_T_1DExpl
!------------------------------------------------------------------------------!
!> This subroutine applies a subsonic outlet bondary condition for 1D nonequilibrium flows.
!! The physical variable imposed is the static temperature.
  SUBROUTINE sub_out_neq_Tout_1DExpl (id, phys_data, u_phys, ghost_data, u_ghost)

    USE mod_general_data,             ONLY: nb_ns, nb_temp, nb_eq, pos_u_cell, pos_T_cell, pos_u, pos_T 
    USE mod_domain_boundary,          ONLY: boundary, boundary_data, get_boundary_Tout
    USE mod_function_pointer,         ONLY: get_cons_phys_from_prim

    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8) :: Tout
    REAL(KIND=8), DIMENSION(nb_eq) :: prim
    TYPE(boundary) :: bound

    INTEGER, INTENT(IN) :: id                              !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys       !< conservative variables of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data    !< physical properties of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost     !< conservative variables of ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data  !< physical properties of ghost cell

    ! Outlet temperature
    bound = boundary_data(id)
    Tout  = get_boundary_Tout (bound)   
  
    ! Ghost state primitive variables (linear extrapolation used for the temperature)
    ! Species densities
    DO i = 1,nb_ns
       prim(i) = u_phys(i)
    ENDDO

    ! Velocity 
    prim(pos_u) = phys_data(pos_u_cell)

    ! Temperature(s)
    prim(pos_T) = 2.d0*Tout - phys_data(pos_T_cell)
    DO i = 2,nb_temp
       prim(pos_T + i - 1) = phys_data(pos_T_cell + i - 1)
    ENDDO
    
    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost, ghost_data)

  END SUBROUTINE sub_out_neq_Tout_1DExpl
!------------------------------------------------------------------------------!
