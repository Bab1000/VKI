!------------------------------------------------------------------------------!
!> This subroutine applies a supersonic inlet boundary conditions for 1D stagnation line nonequilibrium flows.
  SUBROUTINE sup_in_neq_1D_SL_Impl_2nd(id, phys_data1, u_phys1, phys_data2, u_phys2, & 
                                    &  ghost_data1, u_ghost1, ghost_data2, u_ghost2, & 
                                    &  ju_ghost)


    USE mod_general_data,          ONLY: nb_ns, nb_temp, nb_eq, pos_u, pos_v, pos_T, pos_u_cell, & 
                                       & pos_v_cell, pos_T_cell                                         
    USE mod_domain_boundary,       ONLY: boundary, boundary_data, get_boundary_inlet 
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim
 
    IMPLICIT NONE

    INTEGER :: i, j
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
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: ju_ghost   !< ghost state conservative variable Jacobian

    ! Inlet values of density, velocity components and pressure
    bound      = boundary_data(id)
    inlet_data = get_boundary_inlet(nb_eq, bound)

    ! Common factor
    ov_3 = 1.d0/3.d0

    ! Ghost state 1
    ! Primitive variables
    ! Species densities
    DO i = 1,nb_ns 
       prim(i) = (4.d0*inlet_data(i) - u_phys2(i))*ov_3
    ENDDO

    ! Velocity components 
    prim(pos_u) = (4.d0*inlet_data(pos_u) - phys_data2(pos_u_cell))*ov_3
    prim(pos_v) = (4.d0*inlet_data(pos_v) - phys_data2(pos_v_cell))*ov_3

    ! Temperature(s)
    DO i = 1,nb_temp
       prim(pos_T + i - 1) = (4.d0*inlet_data(pos_T + i - 1) - phys_data2(pos_T_cell + i - 1))*ov_3
    ENDDO

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost1, ghost_data1)

    ! Ghost state 2 
    ! Primitive variables 
    ! Species densities
    DO i = 1,nb_ns 
       prim(i) = 2.d0*inlet_data(i) - u_phys2(i)
    ENDDO

    ! Velocity components 
    prim(pos_u) = 2.d0*inlet_data(pos_u) - phys_data2(pos_u_cell)
    prim(pos_v) = 2.d0*inlet_data(pos_v) - phys_data2(pos_v_cell)

    ! Temperature(s)
    DO i = 1,nb_temp
       prim(pos_T + i - 1) = 2.d0*inlet_data(pos_T + i - 1) - phys_data2(pos_T_cell + i - 1)
    ENDDO

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost2, ghost_data2)

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
    !WRITE(*,'(A)')'In "sup_in_neq_1D_SL_Impl_2nd", boundary conditions Jacobian to be implemented!'
    !PRINT*
    !STOP

  END SUBROUTINE sup_in_neq_1D_SL_Impl_2nd
!------------------------------------------------------------------------------!
