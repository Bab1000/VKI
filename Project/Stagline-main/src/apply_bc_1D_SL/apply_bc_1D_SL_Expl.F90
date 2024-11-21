!------------------------------------------------------------------------------!
! This computes the ghost states according to the selected boundary condition in case of 
! an explicit time-integration scheme for 1D stagnation line flows.
  SUBROUTINE apply_bc_1D_SL_Expl(bound_id, phys_data1, u_phys1, phys_data2, u_phys2, & 
                               & ghost_data1, u_ghost1, ghost_data2, u_ghost2)

    USE mod_numerics_data,        ONLY: poly_rec
    USE mod_function_pointer,     ONLY: get_ghost_state_Expl_1D_SL_1st, get_ghost_state_Expl_1D_SL_2nd

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: bound_id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys1, phys_data1, u_phys2, phys_data2 
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost1, ghost_data1, u_ghost2, ghost_data2
  
    SELECT CASE(poly_rec)

      ! 1st order accurate solution
      CASE('constant')
        CALL get_ghost_state_Expl_1D_SL_1st(bound_id)%bc_procedure(bound_id, phys_data1, u_phys1, & 
                                                                 & ghost_data1, u_ghost1)

        u_ghost2    = u_ghost1
        ghost_data2 = ghost_data1

      ! 2nd order accurate solution
      CASE('linear')
        CALL get_ghost_state_Expl_1D_SL_2nd(bound_id)%bc_procedure(bound_id, phys_data1, u_phys1, phys_data2,   & 
                                                                 & u_phys2, ghost_data1, u_ghost1, ghost_data2, & 
                                                                 & u_ghost2)

    END SELECT

  END SUBROUTINE apply_bc_1D_SL_Expl
!------------------------------------------------------------------------------!
