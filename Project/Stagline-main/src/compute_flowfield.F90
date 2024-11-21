!------------------------------------------------------------------------------!
!> This subroutine computes the flowfield by means of method of line approach (separated time and space discretization).
  SUBROUTINE compute_flowfield ()

#include "config.h"

    USE mod_general_data,        ONLY: nb_cells, read_rate, write_rate, inter_cfl, nb_ns, nb_eq, u, u_old, du, adapt_cfl, &
                                 & simulation_id, dres, ios_dir, log_cfl_number, log_cfl_ite, log_cfl, log_cfl_int, &
                                 & cell_prop, nb_prop, pos_T_cell, ser_cfl, gamma_ser, alpha_ser, volumes,          &
                                 & ablation_library_name, nb_temp
    USE mod_numerics_data,      ONLY: nb_stages, time_step, stage_coeff, cfl
    USE mod_function_pointer,   ONLY: evolve_solution
    USE mod_radiation,           ONLY: coupling_period, flag_radiation, compute_rad_source_term, wall_Rad

#ifdef CARBONABLA
    USE mod_ablation_pointers,  ONLY: procedure_ablation_write_sol
#endif


    IMPLICIT NONE

    INTEGER, PARAMETER :: out1 = 15
    INTEGER :: i, j
    INTEGER :: it, k
    INTEGER :: pos
    INTEGER :: cell_number_max_residual, equation_number_max_residual 
    REAL(KIND=8) :: tmp, tmp_res, alpha, res0, max_residual
    REAL (KIND=8) :: max_cfl, min_cfl, pow_max, pow_min
    REAL(KIND=8) :: time
    REAL(KIND=8), DIMENSION(nb_eq) :: res
    REAL(KIND=8), DIMENSION(nb_temp) :: free_stream_energies
    INTEGER :: previous_iter   
    LOGICAL :: stop_cond
    EXTERNAL :: write_sol, physical_data, update_cfl 

    ! Explicit interface for the function "ev_stop_cond"
    INTERFACE
      FUNCTION ev_stop_cond (it, level, res)
        INTEGER, INTENT(IN) :: it
        REAL(KIND=8), INTENT(IN) :: level
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: res
        LOGICAL :: ev_stop_cond
      END FUNCTION ev_stop_cond
    END INTERFACE 

    ! Convergence file
    IF(ios_dir) THEN
        OPEN(UNIT=out1,FILE='./'//TRIM(simulation_id)//'_convergence.dat',STATUS='unknown')
    ELSE
    OPEN(UNIT=out1,FILE='../output/'//TRIM(simulation_id)//'_convergence.dat',STATUS='unknown')
    ENDIF

    WRITE(out1,5)'# Convergence file'
    WRITE(out1,5)'# Iter  CFL  res'
    CALL flush(out1)

    ! Single-stage time-integration scheme 
    IF (nb_stages.EQ.1) THEN

       it   = 0
       res  = 0.d0
       res0 = 0.d0
       time = 0.d0
       log_cfl_int = 1
       previous_iter = 0

       stop_cond = ev_stop_cond (it, time, res)
       
       DO WHILE (stop_cond.EQV..FALSE.)
       
          !CALL shock_sensor()

          IF(flag_radiation .AND. (mod(it,coupling_period).EQ.0) ) THEN
            
            ! Make sure the restart file is up to date
            CALL write_sol()
            ! Update the radition source terms
            CALL compute_rad_source_term()
          ENDIF

          ! Compute solution variation vector (in terms of conserative variables)
          CALL evolve_solution()
          
          !! Get the freestream energies
           !write(*,*) 'Freestream energies:'
           !do i = 1,nb_temp
           !    free_stream_energies(i) = u(nb_cells*(nb_eq)+nb_ns+2+i)
           !    write(*,*) free_stream_energies(i)
           !end do
           
          max_residual = 0.0
          ! Update solution and compute its residual (L2 norm) 
          pos = nb_eq
          res = 0.d0
          DO i = 1,nb_cells
             pos = pos + nb_eq
             DO j = 1,nb_eq
                tmp = du(pos + j)
                tmp_res = dres(pos + j)
                DO WHILE (j <= nb_ns .AND. u(pos+j)-tmp < 0.0d0)
                    tmp = 0.5d0*tmp
                END DO
                
                if (j > nb_ns+3 .and. u(pos+j)-tmp < 0.0d0) then
                    write(*,*) 'Internal energy going negative in cell ', i
                    write(*,*) 'energy density = ', u(pos+j)/sum(u(pos+1:pos+1+nb_ns))
                    write(*,*) 'e(i-1)         = ', u(pos+j-nb_eq)
                    write(*,*) 'e(i+1)         = ', u(pos+j+nb_eq)
                    write(*,*) 'update         = ', -tmp
                end if
                    u(pos + j) = u(pos + j) - tmp
                !end if
                
                IF (tmp_res> max_residual) THEN

                max_residual = tmp_res
                cell_number_max_residual = i+2
                equation_number_max_residual = j

                ENDIF  
                  
                res(j)     = res(j) + tmp_res*tmp_res*volumes(i+2)
             ENDDO
             
               !write(*,*) "updated cell ", i
               !DO j = 1,nb_eq
               !   write(*,*) i, u(pos+j)
              !END DO
          ENDDO

         ! IF(mod(it+1,10).EQ.0)  CALL update_mesh()
         !  CALL update_mesh()

          ! Update physical properties 
          !write(*,*) "calling physical_data() for cell "
          CALL physical_data()
          !DO i = 1,nb_cells
          !   IF (cell_prop(nb_prop*(i+1)+pos_T_cell)<=50.d0) PRINT *, "T below 50 K", i
          !   IF (cell_prop(nb_prop*(i+1)+pos_T_cell+1)<=50.d0) PRINT *, "Tve below 50 K", i
          !ENDDO
          
          ! Residual in log10 scale
          !res = DLOG10(DSQRT(res/nb_cells)) 
          res = DLOG10(DSQRT(res/sum(volumes(3:nb_cells+2)))) 
          
          it   = it + 1
          time = time + time_step

          ! Test stop condition
          stop_cond = ev_stop_cond (it, time, res)
         
          WRITE(out1,20)it,cfl,res(nb_ns + 1:nb_eq) 
          CALL flush(out1)

          ! Solution output file
          IF (MOD(it,write_rate).EQ.0.AND.(write_rate.GT.1)) THEN
             WRITE(*,5)'solver_fvmcc_F90:: Solution output update'
             CALL write_sol()
             WRITE(*,*) "Max residual at position in the equation", DLOG(max_residual), cell_number_max_residual, &
          & equation_number_max_residual
          ENDIF

          IF (ablation_library_name .NE.'none') THEN
             CALL procedure_ablation_write_sol()
          ELSE
          ENDIF

          ! Update the CFL number (history log specified by the user)
          !IF ((log_cfl.EQV..TRUE.).AND.(MOD(it, previous_iter + log_cfl_ite(log_cfl_int)).EQ.0)) THEN
          !   CALL update_cfl()
          !   previous_iter = it
          !ENDIF

          ! Update the CFL number (if specified by the user)
          IF ((inter_cfl.EQV..TRUE.).AND.(MOD(it,read_rate).EQ.0)) THEN
             CALL update_cfl()
          ENDIF
          
          ! Update the CFL based on the SER strategy
          IF (ser_cfl) THEN
              IF (dlog10(sum(10.0d0**res)) < res0 * 1.01) THEN
                  cfl = cfl * 1.05d0
                  res0 = dlog10(sum(10.0d0**res))
              ELSE IF (dlog10(sum(10.0d0**res)) > res0 * 1.2) THEN
                  cfl = cfl * 0.5d0
                  res0 = dlog10(sum(10.0d0**res))
              END IF
    
              
              !cfl = min(cfl, max(1.0, real(it/100)))              
              !cfl = gamma_ser*(res0/sum(10.0D0**res))**alpha_ser
          ENDIF

           if (adapt_cfl) then
              max_cfl = 1.10d0
              min_cfl = min(0.25D0,cfl) 
              tmp_res = 10.d0**res(nb_eq)
            !   pow_max = 1.5d0
            !   pow_min = 1.0d0
                pow_max = 1.5
                pow_min = 1.5

              if (it > 1) then
               if (res(nb_eq) < 5.0) then
                  ! max_cfl = 1.25d0
                  min_cfl = 1.0D0 
               end if

               if (it > 3000) then
                  ! pow_max = 2.0
                  ! pow_min = 0.5
                  pow_max = 2.0
                  pow_min = 1.5
                  min_cfl = min(100.0D0,cfl) 
               endif

                  if (tmp_res <= 1.D0*res0) then
                      cfl = cfl * min((res0 / tmp_res)**pow_max, max_cfl)
                  else if (tmp_res > 1.10D0*res0) then
                       cfl = max(cfl * (res0 / tmp_res)**pow_min, min_cfl)
                     !  cfl = max(cfl/2.0D0, min_cfl)
                  endif
              else
                  res0 = tmp_res
              end if

              res0 = 0.9*res0 + 0.1*tmp_res
          endif  

       ENDDO

    ! Multi-stage time-integration scheme
    ELSE

       it   = 0
       res  = 0.d0 
       time = 0.d0
       log_cfl_int = 1
       previous_iter = 0

       stop_cond = ev_stop_cond (it, time, res)

       DO WHILE (stop_cond.EQV..FALSE.)

          u_old = u

          DO k = 1,nb_stages

             ! Stage coefficient  
             alpha = stage_coeff(k)

             ! Compute solution variation vector (in terms of conserative variables)
             CALL evolve_solution()

             ! Update solution (stage k) 
             pos = nb_eq
             res = 0.d0
             DO i = 1,nb_cells
                pos = pos + nb_eq
                DO j = 1,nb_eq
                   tmp = du(pos + j)*alpha
                   u(pos + j) = u_old(pos + j) - tmp
                ENDDO
             ENDDO
          
             ! Update physical properties (stage k)
             CALL physical_data ()
 
          ENDDO

          ! Solution residual (L2 norm) 
          pos = nb_eq
          res = 0.d0
          DO i = 1,nb_cells
             pos = pos + nb_eq
             DO j = 1,nb_eq
                tmp_res = dres(pos + j)
                res(j) = res(j) + tmp_res*tmp_res
             ENDDO
          ENDDO

          ! Residual in log10 scale
          res = DLOG10(DSQRT(res/nb_cells)) 

          it   = it + 1
          time = time + time_step

          ! Test stop condition
          stop_cond = ev_stop_cond (it, time, res)
         
          WRITE(out1,20)it,cfl,res(nb_ns + 1:nb_eq) 
          CALL flush(out1)

          ! Solution output file
          IF (MOD(it,write_rate).EQ.0.AND.(write_rate.GT.1)) THEN
             WRITE(*,5)'solver_fvmcc_F90:: Solution output update'
             CALL write_sol()
          ENDIF
         
         ! Update the CFL number (history log specified by the user)

         IF ((log_cfl.EQV..TRUE.).AND.(MOD(it, previous_iter + log_cfl_ite(log_cfl_int)).EQ.0)) THEN

         CALL update_cfl()
         previous_iter = it
         ENDIF


          ! Update the CFL number (if specified by the user)
          IF ((inter_cfl.EQV..TRUE.).AND.(MOD(it,read_rate).EQ.0)) THEN
             CALL update_cfl()
          ENDIF

       ENDDO

    ENDIF

    CLOSE(out1)

    PRINT*
    WRITE(*,5)'**********************************'
    PRINT*
    WRITE(*,5)'solver_fvmcc_F90:: Solution is done!'
    PRINT*
    WRITE(*,5)'**********************************'
    PRINT*
    WRITE(*,10)'Time',time,'[s]'
    PRINT*
    WRITE(*,30)'Iter',it,'res',res(nb_ns + 1:nb_eq)
    PRINT*
    
    WRITE(*,5)'**********************************'
    PRINT*
    WRITE(*,5)'solver_fvmcc_F90:: Writing output files'
    PRINT*

5  FORMAT(A)
10 FORMAT(A4,1X,E14.6,1X,A4)
20 FORMAT(I8,1X,F20.5,1X,100F10.5)
30 FORMAT(A4,1X,I6,1X,A3,1X,100F10.5)

   CONTAINS

      SUBROUTINE update_mesh()

         USE mod_general_data,             ONLY: xc, volumes
         USE mod_adapt,                    ONLY: adapt_mesh, interpolate

         IMPLICIT NONE
         REAL(KIND=8), DIMENSION(nb_cells+1) :: x_nodes
         REAL(KIND=8), DIMENSION(nb_cells)   :: prim_T
         REAL(KIND=8), DIMENSION(nb_eq, nb_cells)   :: cons
         INTEGER :: i

         x_nodes(1) = xc(3)-volumes(3)*0.5d0
         DO i = 1, nb_cells
            x_nodes(i+1) = x_nodes(i)+volumes(i+2)
         ENDDO

         CALL physical_data()
         DO i = 1, nb_cells
            prim_T(i) = cell_prop((i+2)*nb_prop+pos_T_cell)
         ENDDO

         ! Computes the new mesh
         CALL adapt_mesh(x_nodes, prim_T)

         DO i = 1,nb_cells
            xc(i + 2)      = 0.5d0*(x_nodes(i) + x_nodes(i + 1))
            volumes(i + 2) = x_nodes(i + 1) - x_nodes(i)
         ENDDO
         volumes(1) = volumes(3)
         volumes(2) = volumes(3)
         xc(2) = xc(3) - volumes(3)
         xc(1) = xc(2) - volumes(3)
         volumes(nb_cells + 3) = volumes(nb_cells + 2)
         volumes(nb_cells + 4) = volumes(nb_cells + 2)
         xc(nb_cells + 3) = xc(nb_cells + 2) + volumes(nb_cells + 2)
         xc(nb_cells + 4) = xc(nb_cells + 3) + volumes(nb_cells + 2)
         ! Writes the mesh
         CALL system('rm mesh.dat')
         OPEN(UNIT=1000,FILE='mesh.dat')
         WRITE(1000,*) nb_cells+1
         DO i = 1,nb_cells+1
            WRITE(1000,*) x_nodes(i)
         ENDDO
         CLOSE(1000)
         
         ! Interpolate conservative variables on the new mesh
         cons = RESHAPE(u(2*nb_eq+1:(nb_cells+2)*nb_eq), (/nb_eq, nb_cells/))
         CALL interpolate(cons)
         u(2*nb_eq+1:(nb_cells+2)*nb_eq) = PACK(cons,.true.)

      END SUBROUTINE

      !SUBROUTINE shock_sensor()
!
      !USE mod_general_data,             ONLY: xc, cell_prop, nb_cells, nb_prop, pos_pres_cell, &
      !                                  & start_prop_phys,shockIndex, shockPosition
!
      !IMPLICIT NONE
      !INTEGER :: i2,i, shockIndex_old, shockIndexStar
      !REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: Pressure, xcenter, gradP
      !REAL(KIND=8), PARAMETER :: alpha_relax = 0.5
!
      !allocate(Pressure(nb_cells))
      !allocate(xcenter(nb_cells))
      !allocate(gradP(nb_cells-1))
     !
      !i2 = start_prop_phys
     !
      !shockIndex_old = shockIndex
!
      !DO i = 1, nb_cells
      !  Pressure(i) = cell_prop(i2 + pos_pres_cell)
      !  xcenter(i) = xc(i+2)
      !  i2 = i2 + nb_prop
      !ENDDO
!
      ! DO i = nb_cells-1, 1, -1
      !  gradP(i) = (Pressure(i+1)-Pressure(i))/(xcenter(i+1)-xcenter(i))
      ! if (gradP(i) < -100.D0) THEN
      !  shockIndexStar = i-1
      !  shockPosition = xcenter(i-1)
      !  exit
      !  ENDIF
      ! ENDDO
      ! shockIndex = int((1-alpha_relax)*shockIndex_old+alpha_relax*shockIndexStar) 
      ! shockIndex = shockIndex-1
      !WRITE(*,*) "Shock position at", shockIndex, xcenter(shockIndex), "relaxed", shockIndexStar, shockIndex_old
      !deallocate(Pressure)
      !deallocate(xcenter)
      !deallocate(gradP)
      !END SUBROUTINE shock_sensor 


   END SUBROUTINE compute_flowfield 
!------------------------------------------------------------------------------!


