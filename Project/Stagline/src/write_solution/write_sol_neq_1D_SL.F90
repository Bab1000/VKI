!------------------------------------------------------------------------------!
!> This subroutine writes the solution for 1D stagnation line nonequilibrium flow simulations.
!! An additional file (restart.dat) is written in order to enable a re-start from a previous solution.
  SUBROUTINE write_sol_neq_1D_SL ()

    USE mod_general_data,    ONLY: nb_eq, nb_dim, nb_ns, nb_temp, nb_prop, nb_cells, pos_u_cell,     &
                                 & pos_pres_cell, pos_T_cell, pos_rho_cell, pos_c_cell, pos_h0_cell, & 
                                 & start_u_phys, start_prop_phys, xc, u, cell_prop, flag_extra_data, &
                                 & simulation_id, ios_dir
   USE mod_neq_function_pointer, ONLY: library_get_molar_fractions, library_species_name, library_comp_tol
 
   USE mod_general_ablation,  ONLY: mdot, uwall
   USE mod_surface_properties
 
    IMPLICIT NONE

    INTEGER :: i, j
    INTEGER :: i1, i2
    INTEGER, PARAMETER :: out1 = 1
    INTEGER, PARAMETER :: out2 = 2 
    INTEGER, PARAMETER :: out3 = 3
    REAL(KIND=8) :: Mach, T, xi_i, rhoi
    REAL(KIND=8), DIMENSION(nb_eq) :: cons
    REAL(KIND=8), DIMENSION(nb_prop) :: cell_data
    REAL(KIND=8), DIMENSION(nb_ns) :: xi
    

   CHARACTER(LEN=20),DIMENSION(nb_ns) :: density_filename
   CHARACTER(LEN=20),DIMENSION(nb_ns) :: mole_frac_filename
   CHARACTER(LEN=30):: temp 
   CHARACTER(LEN=1000):: density_ans
   CHARACTER(LEN=1000):: mole_frac_ans
   CHARACTER(LEN=10), DIMENSION(nb_ns):: species_name
 
    IF (ios_dir) THEN
    OPEN(UNIT=out1,FILE='./'//TRIM(simulation_id)//'_flowfield.dat',STATUS='unknown')
    OPEN(UNIT=out2,FILE='./'//TRIM(simulation_id)//'_restart.dat',STATUS='unknown')
    OPEN(UNIT=out3,FILE='./'//TRIM(simulation_id)//'_restartBC.dat',STATUS='unknown')
    ELSE
    OPEN(UNIT=out1,FILE='../output/'//TRIM(simulation_id)//'_flowfield.dat',STATUS='unknown')
    OPEN(UNIT=out2,FILE='../output/'//TRIM(simulation_id)//'_restart.dat',STATUS='unknown')
    OPEN(UNIT=out3,FILE='../output/'//TRIM(simulation_id)//'_restartBC.dat',STATUS='unknown')
    ENDIF


    CALL library_species_name(species_name)

    DO i=1,nb_ns
    WRITE (density_filename(i) ,"('rhoi(',A,') ')")TRIM(species_name(i))
    WRITE (mole_frac_filename(i) ,"('xi(',A,') ')")TRIM(species_name(i))
    ENDDO 
    temp=''
    mole_frac_ans=''
    density_ans= '' 

    DO i=1,int(nb_ns/2) 
    temp = ' "'//TRIM(density_filename(i*2-1))//'"  "'// TRIM(density_filename(i*2))//'"'
    density_ans = TRIM(density_ans)//TRIM(temp)
    temp = ' "'//TRIM(mole_frac_filename(i*2-1))//'" "'// TRIM(mole_frac_filename(i*2))//'"'
    mole_frac_ans = TRIM(mole_frac_ans)//TRIM(temp)
    ENDDO
   
   IF (MOD(nb_ns,2).NE.0) THEN   

   density_ans= TRIM(density_ans)// ' "' //TRIM(density_filename(nb_ns))//'"'
   mole_frac_ans= TRIM(mole_frac_ans)// ' "' //TRIM(mole_frac_filename(nb_ns))//'"'

   ENDIF

            
    ! Flowfield file header
    WRITE(out1,10)'# solver_fvmcc_F90:: 1D stagnation line solution file'

    IF (flag_extra_data.EQV..TRUE.) THEN
       WRITE(out1,10)'# Nonequilibrium flow'
       WRITE(out1,10)'# Primitive variables and extra data'

       IF (nb_temp.EQ.1) THEN 

          WRITE(out1,10)'#"xc"  "u"  "v"  "T"  "p"  "rho"  "h0"  "Mf"'//TRIM(mole_frac_ans)//TRIM(density_ans)
                                                                                                  
       ELSE   
          WRITE(out1,10)'#"xc"  "u"  "v"  "T"  "Tint"  "p"  "rho"  "h0"  "Mf"'//TRIM(mole_frac_ans)//&
                                                                        &TRIM(density_ans)
       ENDIF
      

    !   WRITE(out1,10)'# ZONE T="'//TRIM(simulation_id)//'"'

   ELSE 
       WRITE(out1,10)'# Nonequilibrium flow'
       WRITE(out1,10)'# Primitive variables'
       WRITE(out1,10)'# "xc"  "u"  "v"  "T"  "Tint" '//TRIM(mole_frac_ans)//TRIM(density_ans)
   !    WRITE(out1,10)'# ZONE T="'//TRIM(simulation_id)//'"'

    ENDIF

    ! Index initialization
    i1 = start_u_phys
    i2 = start_prop_phys

    ! Write data to files
    ! Extra data included
    IF (flag_extra_data.EQV..TRUE.) THEN
      WRITE(out3,20) rho_i_surface(1:nb_ns), u_surface, v_surface,T_surface(1:nb_temp)
       DO i = 1,nb_cells
      
          ! Conservative variables 
          DO j = 1,nb_eq
             cons(j) = u(i1 + j)
          ENDDO

          ! Physical properties
          DO j = 1,nb_prop
             cell_data(j) = cell_prop(i2 + j)
          ENDDO
             
           CALL library_get_molar_fractions (cons(1:nb_ns), xi) 
           CALL library_comp_tol(xi)
               
          ! Frozen Mach number
          Mach = ABS(cell_data(pos_u_cell))/cell_data(pos_c_cell)

         DO j = 1,nb_ns
          if (cons(j) < 1.D-18) THEN
           cons(j) = 1.D-18
           endif 
         enddo

         WRITE(out1,20)xc(i + 2),cell_data(pos_u_cell:pos_u_cell + nb_dim - 1),  & 
                       & cell_data(1:nb_temp),cell_data(pos_pres_cell),cell_data(pos_rho_cell),  & 
                      & cell_data(pos_h0_cell),Mach, xi, cons(1:nb_ns)   
          WRITE(out2,20)cons(1:nb_ns),cell_data(pos_u_cell:pos_u_cell + nb_dim - 1),cell_data(1:nb_temp)
          WRITE(out3,20)cons(1:nb_ns),cell_data(pos_u_cell:pos_u_cell + nb_dim - 1),cell_data(1:nb_temp)

              
          ! Index update
          i1 = i1 + nb_eq
          i2 = i2 + nb_prop

       ENDDO

    ! No extra data included
    ELSE
      WRITE(out3,20) rho_i_surface(1:nb_ns),u_surface, v_surface,T_surface(1:nb_temp)
       DO i = 1,nb_cells
      
          ! Conservative variables 
          DO j = 1,nb_eq
             cons(j) = u(i1 + j)
          ENDDO

          ! Physical properties
          DO j = 1,nb_prop
             cell_data(j) = cell_prop(i2 + j)
          ENDDO

          ! Frozen Mach number
          Mach = ABS(cell_data(pos_u_cell))/cell_data(pos_c_cell)

          WRITE(out1,20)xc(i + 2),cell_data(pos_u_cell:pos_u_cell + nb_dim - 1),  & 
                      & cell_data(1:nb_temp), xi, cons(1:nb_ns)
          WRITE(out2,20)cons(1:nb_ns),cell_data(pos_u_cell:pos_u_cell + nb_dim - 1),cell_data(1:nb_temp)
          WRITE(out3,20)cons(1:nb_ns),cell_data(pos_u_cell:pos_u_cell + nb_dim - 1),cell_data(1:nb_temp)

          ! Index update
          i1 = i1 + nb_eq
          i2 = i2 + nb_prop

       ENDDO

    ENDIF

    CLOSE(out1)
    CLOSE(out2)
    CLOSE(out3)

10  FORMAT(A)
20  FORMAT(10000(E20.10,1X))

  END SUBROUTINE write_sol_neq_1D_SL
!------------------------------------------------------------------------------!
