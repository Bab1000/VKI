! 1-D VKI stagnation-line code solution interpolator
!
!> @author
!> Alessandro Turchi 
!>
!> @brief
!> This tool interpolates a solution of the VKI stagnation-line code over a new mesh.
!
!>\b DESCRIPTION: \n 
!> This tool interpolates a solution of the VKI stagnation-line code over a new mesh.
!!
!! It takes as input the following files:
!!      * solution_old.dat: the restart.dat file of the simulation on the old mesh writen by the VKI stagnation-line code
!!      * mesh_old.dat:     the old mesh file
!!      * mesh_new.dat:     the new mesh file
!!
!! It generates the following output:
!!      * solution_new.dat: the file to be used as restart file for the simulation with the new mesh
!!      * solution_new.tec: file identical to the solution_new.dat but including the new mesh in the first column. To be used to check the interpolation results.
!!      * solution_new.tec: file identical to the solution_old.dat but including the old mesh in the first column. To be used to compare with the interpolated result (solution_new.tec).

PROGRAM mesh_interp

IMPLICIT NONE

integer, parameter :: prec = selected_real_kind(15, 307)
integer, parameter :: in1 = 11, in2 = 12,  in3 = 13
integer, parameter :: out1 = 14, out2 = 15, out3 = 16
integer :: ios
integer :: i, j, k
integer :: new_nodes_num , old_nodes_num
integer :: species, ntemps
real(kind=prec), allocatable, dimension(:) :: new_node, old_node
real(kind=prec), allocatable, dimension(:) :: new_cell_center, old_cell_center
real(kind=prec), allocatable, dimension(:) :: delta_new, delta_old
real(kind=prec), allocatable, dimension(:,:) :: sol_old, sol_new
real(kind=prec), allocatable, dimension(:,:) :: dsol_ov_dspace
character(len=40) :: output
character(len=3) :: string

        
        OPEN(UNIT=in1,FILE='mesh_new.dat',STATUS='old',IOSTAT=ios)

            IF (ios.NE.0) THEN 
               WRITE(*,5)'Attention: mesh_new.dat file not found'
               PRINT*
               STOP
            ENDIF

        OPEN(UNIT=in2,FILE='mesh_old.dat',STATUS='old',IOSTAT=ios)

            IF (ios.NE.0) THEN 
               WRITE(*,5)'Attention: mesh_old.dat file not found'
               PRINT*
               STOP
            ENDIF

        OPEN(UNIT=in3,FILE='solution_old.dat',STATUS='old',IOSTAT=ios)

            IF (ios.NE.0) THEN 
               WRITE(*,5)'Attention: solution_old.dat file not found'
               PRINT*
               STOP
            ENDIF

        WRITE(*,*)'Insert species number:'         
        READ(*,*)  species

        WRITE(*,*) 'Insert number of temperatures:'
        READ(*,*)  ntemps

        OPEN(UNIT=out1,FILE='solution_new.dat',STATUS='new',IOSTAT=ios)

            IF (ios.NE.0) THEN 
               WRITE(*,5)'Attention: solution_new.dat file already exists'
               PRINT*
               STOP
            ENDIF

        OPEN(UNIT=out2,FILE='solution_new.tec',STATUS='new',IOSTAT=ios)

            IF (ios.NE.0) THEN 
               WRITE(*,5)'Attention: solution_new.tec file already exists'
               PRINT*
               STOP
            ENDIF

        OPEN(UNIT=out3,FILE='solution_old.tec',STATUS='new',IOSTAT=ios)

            IF (ios.NE.0) THEN 
               WRITE(*,5)'Attention: solution_old.tec file already exists'
               PRINT*
               STOP

            ENDIF
        

        READ(in1,*) new_nodes_num

        ALLOCATE(new_node(new_nodes_num))
        ALLOCATE(delta_new(new_nodes_num-1))
        ALLOCATE(new_cell_center(new_nodes_num-1))
        ALLOCATE(sol_new(new_nodes_num-1,species+2+ntemps))

        DO i=1,new_nodes_num
                READ(in1,10) new_node(i)
        ENDDO

        DO i=1,new_nodes_num-1
                delta_new(i) = new_node(i+1)-new_node(i)
        ENDDO
        
        DO i=1,new_nodes_num-1
                new_cell_center(i) = new_node(i)+delta_new(i)/2.d0
        ENDDO


        READ(in2,*) old_nodes_num

        ALLOCATE(old_node(old_nodes_num))
        ALLOCATE(delta_old(old_nodes_num-1))
        ALLOCATE(old_cell_center(old_nodes_num-1))
        ALLOCATE(sol_old(old_nodes_num-1,species+2+ntemps))

        DO i=1, old_nodes_num
                READ(in2,*) old_node(i)
        ENDDO

        DO i=1,old_nodes_num-1
                delta_old(i) = old_node(i+1)-old_node(i)
        ENDDO
        
        DO i=1,old_nodes_num-1
                old_cell_center(i) = old_node(i)+delta_old(i)/2.d0
        ENDDO
        DO i = 1 , old_nodes_num-1
       
                READ(in3,*) (sol_old(i , k),k=1,species+2+ntemps)
        ENDDO

        ALLOCATE(dsol_ov_dspace(old_nodes_num-2,species+2+ntemps))
        DO k = 1, species+2+ntemps
                DO i = 1, old_nodes_num-2
                        dsol_ov_dspace(i,k) = (sol_old(i+1, k) - sol_old(i,k))/(old_cell_center(i+1) - old_cell_center(i))
                ENDDO
        ENDDO

        DO i = 1, new_nodes_num-1 

                DO j = 1, old_nodes_num-2

                        IF (new_cell_center(i).ge.old_cell_center(j).and.new_cell_center(i).lt.old_cell_center(j+1)) THEN
                                DO k = 1, species+2+ntemps
                                        sol_new(i,k) = sol_old(j,k) +  dsol_ov_dspace(j,k) *(new_cell_center(i)-old_cell_center(j))
                                ENDDO

                        ELSEIF (new_cell_center(i).lt.old_cell_center(1)) THEN
                                DO k = 1, species+2+ntemps
                                        sol_new(i,k) = sol_old(1,k) + dsol_ov_dspace(1,k)*(new_cell_center(i)-old_cell_center(1))
                                ENDDO
                        ELSEIF (new_cell_center(i).ge.old_cell_center(old_nodes_num-1)) THEN
                                DO k = 1, species+2+ntemps
                                        sol_new(i,k) = sol_old(old_nodes_num-1,k) + dsol_ov_dspace(old_nodes_num-2,k) &
                                                     & *(new_cell_center(i)-old_cell_center(old_nodes_num-1))
                                ENDDO
                        ELSE
                        ENDIF
                ENDDO
        ENDDO

        DO i = 1, new_nodes_num-1 
                DO k = 1, species
                        IF (sol_new(i,k) < 0.0D0) THEN
                        sol_new(i,k) = 1.0D-16
                        ENDIF
                ENDDO
        ENDDO

        DO i = 1, new_nodes_num-1
                WRITE(out1,8) (sol_new(i,k), k = 1, species+2+ntemps)
                WRITE(out2,9) new_cell_center(i),  (sol_new(i,k), k = 1, species+2+ntemps)
        ENDDO


        DO i = 1, old_nodes_num-1
                WRITE(out3,9) old_cell_center(i), (sol_old(i,k), k = 1, species+2+ntemps)
        ENDDO


DEALLOCATE(new_node,delta_new,new_cell_center,sol_new)
DEALLOCATE(old_node,delta_old,old_cell_center,sol_old)
DEALLOCATE(dsol_ov_dspace)

5 FORMAT(A)
8 FORMAT(100(e15.8,X))
9 FORMAT(100(e15.8,X))
10  FORMAT(1e15.8)
END PROGRAM  mesh_interp

