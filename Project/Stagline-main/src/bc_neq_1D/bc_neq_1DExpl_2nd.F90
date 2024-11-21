!------------------------------------------------------------------------------!
!> This subroutine applies a supersonic inlet boundary condition in case of use
!! of two ghost cells for 1D nonequilibrium flows. A linear extrapolation is used 
!! for conservative variables. Physical data are then computed from the conservative
!! variable values. The first ghost cells is that close to the domain boundary.
  SUBROUTINE sup_in_neq_1DExpl_2nd (id, phys_data1, u_phys1, phys_data2, u_phys2, & 
                                 &  ghost_data1, u_ghost1, ghost_data2, u_ghost2)

    USE mod_general_data,             ONLY: nb_ns, nb_eq, nb_temp, pos_u_cell, pos_T_cell, pos_u
    USE mod_domain_boundary,          ONLY: boundary, boundary_data, get_boundary_inlet 
    USE mod_function_pointer,         ONLY: get_cons_phys_from_prim

    IMPLICIT NONE

    INTEGER :: i
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

    ! Inlet data
    bound      = boundary_data(id)
    inlet_data = get_boundary_inlet (nb_eq, bound)

    ! Ghost state 1
    ! Primitive variables (linear extrapolation used for all variables)
    ! Species densities
    DO i = 1,nb_ns 
       prim(i) = 2.d0*inlet_data(i) - u_phys1(i)
    ENDDO

    ! Velocity 
    prim(pos_u) = 2.d0*inlet_data(pos_u) - phys_data1(pos_u_cell)

    ! Temperature(s)
    ! Temperatures 
    DO i = 1,nb_temp 
       prim(pos_u + i) = 2.d0*inlet_data(pos_u + i) - phys_data1(pos_T_cell + i - 1)
    ENDDO

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost1, ghost_data1)

    ! Ghost state 2
    ! Primitive variables (linear extrapolation used for all variables)
    ! Species densities
    DO i = 1,nb_ns 
       prim(i) = 2.d0*inlet_data(i) - u_phys2(i)
    ENDDO

    ! Velocity 
    prim(pos_u) = 2.d0*inlet_data(pos_u) - phys_data2(pos_u_cell)

    ! Temperature(s)
    ! Temperatures 
    DO i = 1,nb_temp 
       prim(pos_u + i) = 2.d0*inlet_data(pos_u + i) - phys_data2(pos_T_cell + i - 1)
    ENDDO

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost2, ghost_data2)

  END SUBROUTINE sup_in_neq_1DExpl_2nd 
!------------------------------------------------------------------------------!
!> This subroutine applies a supersonic outlet boundary condition in case of use
!! of two ghost cells for 1D nonequilibrium flows. No extrapolation is used for conservative variables
!! and physical data for sake of robustness. The first ghost cells is that close to the domain boundary.
  SUBROUTINE sup_out_neq_1DExpl_2nd (id, phys_data1, u_phys1, phys_data2, u_phys2, & 
                                  &  ghost_data1, u_ghost1, ghost_data2, u_ghost2)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: id                               !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys1       !< conservative variables of the first physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys2       !< conservative variables of the second physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data1    !< physical properties of the first physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data2    !< physical properties of the second physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost1     !< conservative variables of the first ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost2     !< conservative variables of the second ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data1  !< physical properties of the first ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data2  !< physical properties of the second ghost cell

    ! Ghost cell 1
    u_ghost1 = u_phys1
    ghost_data1 = phys_data1

    ! Ghost state 2
    u_ghost2 = u_ghost1
    ghost_data2 = ghost_data1

  END SUBROUTINE  sup_out_neq_1DExpl_2nd 
!------------------------------------------------------------------------------!
!> This subroutine applies a subsonic inlet boundary condition in case of use
!! of two ghost cells for 1D nonequilibrium flows. A linear extrapolation is used 
!! for conservative variables. Physical data are then computed from the conservative
!! variable values. The first ghost cells is that close to the domain boundary.
  SUBROUTINE sub_in_neq_rhoiin_T_1DExpl_2nd (id, phys_data1, u_phys1, phys_data2, u_phys2, & 
                               &  ghost_data1, u_ghost1, ghost_data2, u_ghost2)

    USE mod_general_data,             ONLY: nb_ns, nb_temp, nb_eq, pos_u_cell, pos_T_cell, pos_u, pos_T
    USE mod_domain_boundary,          ONLY: boundary, boundary_data, get_boundary_rhoi, get_boundary_Tvec
    USE mod_function_pointer,         ONLY: get_cons_phys_from_prim

    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8), DIMENSION(nb_ns) :: rhoi_in
    REAL(KIND=8), DIMENSION(nb_temp) :: temp_in
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

    ! Species densities at inlet
    bound   = boundary_data(id)
    rhoi_in = get_boundary_rhoi (nb_ns, bound)

    ! Temperature vector at inlet
    temp_in = get_boundary_Tvec (nb_temp, bound)

    ! Ghost state 1
    ! Primitive variables
    ! Species densities
    DO i = 1,nb_ns
       prim(i) = 2.d0*rhoi_in(i) - u_phys1(i) 
    ENDDO

    ! Velocity 
    prim(pos_u) = 2.d0*phys_data1(pos_u_cell) - phys_data2(pos_u_cell)

    ! Temperature(s)
    DO i = 1,nb_temp
       prim(pos_T + i - 1) = 2.d0*temp_in(i) - phys_data1(pos_T_cell + i - 1)
    ENDDO
   
    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost1, ghost_data1)

    ! Ghost state 2
    ! Species and mixture densities
    DO i = 1,nb_ns
       prim(i) = 2.d0*rhoi_in(i) - u_phys2(i) 
    ENDDO

    ! Velocity and kinetic energy per unit mass
    prim(pos_u) = 3.d0*phys_data1(pos_u_cell) - 2.d0*phys_data2(pos_u_cell)

    ! Tmperature(s)
    DO i = 1,nb_temp
       prim(pos_T + i - 1) = 2.d0*temp_in(i) - phys_data2(pos_T_cell + i - 1)
    ENDDO

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost2, ghost_data2) 

  END SUBROUTINE sub_in_neq_rhoiin_T_1DExpl_2nd
!------------------------------------------------------------------------------!

