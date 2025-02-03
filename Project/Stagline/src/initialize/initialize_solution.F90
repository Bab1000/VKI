!------------------------------------------------------------------------------!
!> This subroutine initializes the conservative variables for all cells given the initial field in primitive variables.
  SUBROUTINE initialize_solution()

    USE mod_general_data,         ONLY: nb_tot_cells, nb_eq, model_name, p, u
    USE mod_function_pointer,     ONLY: get_cons_from_prim

    IMPLICIT NONE
  
    INTEGER :: i, j
    REAL(KIND=8), DIMENSION(nb_eq):: cons, prim

    ! Loop over all cells (physical + ghost)
    DO i = 1,nb_tot_cells
       
       DO j = 1,nb_eq 
          prim(j) = p((i - 1)*nb_eq + j)
       ENDDO

       ! Compute cell conservative variables
       CALL get_cons_from_prim (prim, cons)
  
       ! Allocating initial solution for conservative variables
       DO j = 1,nb_eq 
          u((i - 1)*nb_eq + j) = cons(j) 
       ENDDO
  
    ENDDO
    
  END SUBROUTINE initialize_solution 
!------------------------------------------------------------------------------!

