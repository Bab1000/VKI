!------------------------------------------------------------------------------!
!> This subroutine computes the physical data for all physical cells given the whole solution 
!! conservative variable vector. It is used for calorically perfect gas and nonequilibrium flows.
  SUBROUTINE physical_data()

    USE mod_general_data,          ONLY: nb_eq, nb_cells, nb_prop, start_u_phys, start_prop_phys, & 
                                       & u_vec, phys_prop, u, cell_prop
    USE mod_function_pointer,      ONLY: get_phys_from_cons

    IMPLICIT NONE

    INTEGER :: i, j
    INTEGER :: i1, i2

    ! Update physical properties
    ! Index initialization
    i1 = start_u_phys    
    i2 = start_prop_phys
  
    ! Loop over physical cells 
    DO i = 1,nb_cells

       ! Cell conservative variable vector
       DO j = 1,nb_eq
          u_vec(j) = u(i1 + j)
       ENDDO

       ! Get physical properties from conservative variables
       CALL get_phys_from_cons (u_vec, phys_prop)
       DO j = 1,nb_prop 
          cell_prop(i2 + j) = phys_prop(j)
       ENDDO 

       ! Index update
       i1 = i1 + nb_eq
       i2 = i2 + nb_prop

    ENDDO

  END SUBROUTINE physical_data
!------------------------------------------------------------------------------!
