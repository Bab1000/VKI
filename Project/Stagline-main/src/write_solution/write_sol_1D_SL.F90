!------------------------------------------------------------------------------!
!> This subroutine writes the solution for 1D stagnation line calorically perfect gas simulations. 
!! An additional file (restart.dat) is written in order to enable a re-start from a previous solution.
  SUBROUTINE write_sol_1D_SL ()

    USE mod_general_data,    ONLY: nb_eq, nb_dim, nb_prop, nb_cells, pos_u_cell, pos_pres_cell, pos_c_cell,  &
                                 & pos_h0_cell, pos_T_cell, start_u_phys, start_prop_phys, xc, u, cell_prop, & 
                                 & flag_extra_data, simulation_id, ios_dir

    IMPLICIT NONE

    INTEGER :: i, j
    INTEGER :: i1, i2
    INTEGER, PARAMETER :: out1 = 1
    INTEGER, PARAMETEr :: out2 = 2 
    REAL(KIND=8) :: Mach
    REAL(KIND=8), DIMENSION(nb_eq) :: cons
    REAL(KIND=8), DIMENSION(nb_prop) :: cell_data
    
    IF (ios_dir) THEN
        OPEN(UNIT=out1,FILE='./'//TRIM(simulation_id)//'_flowfield.dat',STATUS='unknown')
        OPEN(UNIT=out2,FILE='./'//TRIM(simulation_id)//'_restart.dat',STATUS='unknown')
    ELSE
        OPEN(UNIT=out1,FILE='../output/'//TRIM(simulation_id)//'_flowfield.dat',STATUS='unknown')
        OPEN(UNIT=out2,FILE='../output/'//TRIM(simulation_id)//'_restart.dat',STATUS='unknown')
    ENDIF

    ! Flowfield file header
    WRITE(out1,10)'# solver_fvmcc_F90:: 1D stagnation line solution file'

    IF (flag_extra_data.EQV..TRUE.) THEN
       WRITE(out1,10)'# Calorically perfect gas flow'
       WRITE(out1,10)'# Primitive variables and extra data'
       WRITE(out1,10)'# xc  rho  u  v  p  T  M  h0'
    ELSE 
       WRITE(out1,10)'# Calorically perfect gas flow'
       WRITE(out1,10)'# Primitive variables'
       WRITE(out1,10)'# xc  rho  u  v  p'
    ENDIF

    ! Index initialization
    i1 = start_u_phys
    i2 = start_prop_phys

    ! Write data to files
    ! Extra data included
    IF (flag_extra_data.EQV..TRUE.) THEN

       DO i = 1,nb_cells
      
          ! Conservative variables 
          DO j = 1,nb_eq
             cons(j) = u(i1 + j)
          ENDDO

          ! Physical properties
          DO j = 1,nb_prop
             cell_data(j) = cell_prop(i2 + j)
          ENDDO

          ! Mach number
          Mach = ABS(cell_data(pos_u_cell))/cell_data(pos_c_cell) 
          WRITE(out1,20)xc(i + 2),cons(1),cell_data(pos_u_cell:pos_u_cell + nb_dim - 1),  & 
                      & cell_data(pos_pres_cell),cell_data(pos_T_cell),Mach,cell_data(pos_h0_cell)
          WRITE(out2,20)cons(1),cell_data(pos_u_cell:pos_u_cell + nb_dim - 1),cell_data(pos_pres_cell)

          ! Index update
          i1 = i1 + nb_eq
          i2 = i2 + nb_prop

       ENDDO

    ! No extra data included
    ELSE 

       DO i = 1,nb_cells
      
          ! Conservative variables 
          DO j = 1,nb_eq
             cons(j) = u(i1 + j)
          ENDDO

          ! Physical properties
          DO j = 1,nb_prop
             cell_data(j) = cell_prop(i2 + j)
          ENDDO

          WRITE(out1,20)xc(i + 2),cons(1),cell_data(pos_u_cell:pos_u_cell + nb_dim - 1),  & 
                      & cell_data(pos_pres_cell)
          WRITE(out2,20)cons(1),cell_data(pos_u_cell:pos_u_cell + nb_dim - 1),cell_data(pos_pres_cell)

          ! Index update
          i1 = i1 + nb_eq
          i2 = i2 + nb_prop

       ENDDO    

    ENDIF

    CLOSE(out1)
    CLOSE(out2)

10  FORMAT(A)
20  FORMAT(10000(E20.10,1X))

  END SUBROUTINE write_sol_1D_SL
!------------------------------------------------------------------------------!
