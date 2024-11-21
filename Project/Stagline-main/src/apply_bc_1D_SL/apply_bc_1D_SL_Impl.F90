!------------------------------------------------------------------------------!
! This computes the ghost states according to the selected boundary condition in case of 
! an implicit time-integration scheme for 1D stagnation line flows.
  SUBROUTINE apply_bc_1D_SL_Impl(bound_id, phys_data1, u_phys1, phys_data2, u_phys2, & 
                               & ghost_data1, u_ghost1, ghost_data2, u_ghost2, jb)

    USE mod_numerics_data,        ONLY: poly_rec, bc_Jac
    USE mod_function_pointer,     ONLY: get_ghost_state_Expl_1D_SL_1st, get_ghost_state_Expl_1D_SL_2nd, & 
                                      & get_ghost_state_Impl_1D_SL_1st, get_ghost_state_Impl_1D_SL_2nd

    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: bound_id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys1, phys_data1, u_phys2, phys_data2 
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost1, ghost_data1, u_ghost2, ghost_data2
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jb

    SELECT CASE(poly_rec)

      ! 1st order accurate solution
      CASE('constant')
       
        SELECT CASE(bc_Jac(bound_id))
       
          ! Analytical boundary condition Jacobian
          CASE('analytical')
            CALL get_ghost_state_Impl_1D_SL_1st(bound_id)%bc_procedure(bound_id, phys_data1, u_phys1, & 
                                                                     & ghost_data1, u_ghost1, jb)
          
          ! Numerical boundary condition Jacobian
          CASE('numerical')
            CALL get_ghost_state_Expl_1D_SL_1st(bound_id)%bc_procedure(bound_id, phys_data1, u_phys1, & 
                                                                     & ghost_data1, u_ghost1) 
  
            ! Jacobian by means of numerical differentiation
            CALL bc_num_Jac_1st (bound_id, u_phys1, u_ghost1, jb)

        END SELECT

        u_ghost2    = u_ghost1
        ghost_data2 = ghost_data1

      ! 2nd order accurate solution
      CASE('linear')

        SELECT CASE(bc_Jac(bound_id))

          ! Analytical boundary condition Jacobian
          CASE('analytical')
            CALL get_ghost_state_Impl_1D_SL_2nd(bound_id)%bc_procedure(bound_id, phys_data1, u_phys1, phys_data2,   & 
                                                                     & u_phys2, ghost_data1, u_ghost1, ghost_data2, & 
                                                                     & u_ghost2, jb)
          ! Numerical boundary condition Jacobian
          CASE('numerical')
            CALL get_ghost_state_Expl_1D_SL_2nd(bound_id)%bc_procedure(bound_id, phys_data1, u_phys1, phys_data2,   & 
                                                                     & u_phys2, ghost_data1, u_ghost1, ghost_data2, & 
                                                                     & u_ghost2) 


            ! Jacobian by means of numerical differentiation
            CALL bc_num_Jac_2nd (bound_id, u_phys1, u_phys2, phys_data2, u_ghost1, jb)

        END SELECT

    END SELECT

10 FORMAT(100E14.6)
15 FORMAT(A)

    CONTAINS
   
      !----------------------------------------------------!
      ! This subroutine evaluates the boundary conditions Jacobian numerically 
      SUBROUTINE bc_num_Jac_1st (id, u_phys, u_ghost, jac)

        USE mod_general_data,         ONLY: nb_eq, nb_prop, eta
        USE mod_function_pointer,     ONLY: get_phys_from_cons

        INTEGER :: i, j
        REAL(KIND=8) :: du, ov_du, u
        REAL(KIND=8), DIMENSION(nb_eq) :: u_phys_pert, u_ghost_pert
        REAL(KIND=8), DIMENSION(nb_prop) :: phys_data_pert, ghost_data_pert

        INTEGER, INTENT(IN) :: id
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys, u_ghost
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jac 
        
        DO j = 1,nb_eq

           ! Perturb conservative variables 
           DO i = 1,nb_eq 
              u_phys_pert(i) = u_phys(i)
           ENDDO

           u  = u_phys(j)
           du = eta*MAX(ABS(u),1.d-20)*SIGN(1.d0,u)

           u_phys_pert(j) = u + du

           ! Get the physical properties corresponding to the perturbed conservative variables
           CALL get_phys_from_cons(u_phys_pert, phys_data_pert) 

           ! Get the boundary condition to the perturbed conservative variables
           CALL get_ghost_state_Expl_1D_SL_1st(id)%bc_procedure(id, phys_data_pert, u_phys_pert, & 
                                                              & ghost_data_pert, u_ghost_pert)

           ! Numerical differentiation
           ov_du = 1.d0/du
           DO i = 1,nb_eq
              jac(i,j) = (u_ghost_pert(i) - u_ghost(i))*ov_du
           ENDDO

        ENDDO

      END SUBROUTINE bc_num_Jac_1st

      !----------------------------------------------------!
      ! This subroutine evaluates the boundary conditions Jacobian numerically 
      SUBROUTINE bc_num_Jac_2nd (id, u1, u2, data2, ug1, jac)

        USE mod_general_data,         ONLY: nb_eq, nb_prop, eta
        USE mod_function_pointer,     ONLY: get_phys_from_cons

        INTEGER :: i, j
        REAL(KIND=8) :: du, ov_du, u
        REAL(KIND=8), DIMENSION(nb_eq) :: up1, ugp1, ug2
        REAL(KIND=8), DIMENSION(nb_prop) :: datap1, datag1p, datag2

        INTEGER, INTENT(IN) :: id
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u1, u2, data2, ug1
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jac 
       
        DO j = 1,nb_eq

           ! Perturb conservative variables 
           DO i = 1,nb_eq 
              up1(i) = u1(i)
           ENDDO

           u  = u1(j)
           du = eta*MAX(ABS(u),1.d-20)*SIGN(1.d0,u)            
           
           up1(j) = u + du

           ! Get the physical properties corresponding to the perturbed conservative variables
           CALL get_phys_from_cons(up1, datap1) 

           ! Get the boundary condition to the perturbed conservative variables
           CALL get_ghost_state_Expl_1D_SL_2nd(id)%bc_procedure(id, datap1, up1, data2, u2, datag1p,& 
                                                              & ugp1, datag2, ug2)

           ! Numerical differentiation
           ov_du = 1.d0/du
           DO i = 1,nb_eq
              jac(i,j) = (ugp1(i) - ug1(i))*ov_du
           ENDDO

        ENDDO

      END SUBROUTINE bc_num_Jac_2nd

  END SUBROUTINE apply_bc_1D_SL_Impl
!------------------------------------------------------------------------------!
