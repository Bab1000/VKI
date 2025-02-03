!------------------------------------------------------------------------------!
!> This subroutine post-processes the solution.
  SUBROUTINE post_process_solution ()

    USE mod_general_data,               ONLY: nb_cells, nb_ns, nb_temp, nb_eq, nb_prop,   & 
                                            & start_u_phys, start_prop_phys, pos_T_cell,  & 
                                            & pos_u_cell, u, cell_prop, xc, flag_post_CR, &
                                            & simulation_id, ios_dir
    USE mod_neq_function_pointer,       ONLY: library_write_fvmcc_solution 

    IMPLICIT NONE

    INTEGER, PARAMETER :: out1 = 1, out2 = 2
    INTEGER :: i, is
    INTEGER :: i1, i2
    REAL(KIND=8) :: x, xold, v, vold
    REAL(KIND=8), DIMENSION(nb_ns) :: rhoi
    REAL(KIND=8), DIMENSION(nb_temp) :: temp

    IF (flag_post_CR.EQV..TRUE.) THEN 

       WRITE(*,5)'solver_fvmcc_F90:: solution post-processing'
       PRINT*

       ! Output files
       IF (ios_dir) THEN
           OPEN(UNIT=out1,FILE='./'//TRIM(simulation_id)//'_post_flowfield.dat',STATUS='unknown')
           OPEN(UNIT=out2,FILE='./'//TRIM(simulation_id)//'_pop.dat',STATUS='unknown')
       ELSE
           OPEN(UNIT=out1,FILE='../output/'//TRIM(simulation_id)//'_post_flowfield.dat',STATUS='unknown')
           OPEN(UNIT=out2,FILE='../output/'//TRIM(simulation_id)//'_pop.dat',STATUS='unknown')
       ENDIF

       ! Index initization
       i1 = start_u_phys 
       i2 = start_prop_phys

       xold = xc(3)
       vold = cell_prop(i2 + pos_u_cell)
 
       DO i = 1,nb_cells

          ! Centroid location
          x = xc(i + 2)

          ! Species densities
          DO is = 1,nb_ns
             rhoi(is) = u(i1 + is)
          ENDDO

          ! Velocity 
          v = cell_prop(i2 + pos_u_cell)

          ! Temperature(s)
          DO is = 1,nb_temp
             temp(is) = cell_prop(i2 + pos_T_cell + is - 1)
          ENDDO
         
          ! Post-processing files 
          CALL library_write_fvmcc_solution (out1, out2, i, xold, x, vold, v, rhoi, temp)

          ! Index update
          i1 = i1 + nb_eq
          i2 = i2 + nb_prop

          ! Variable update
          xold = x
          vold = v 

       ENDDO

       CLOSE(out1)
       CLOSE(out2)

    ENDIF

5 FORMAT(A)

  END SUBROUTINE post_process_solution
!------------------------------------------------------------------------------!
