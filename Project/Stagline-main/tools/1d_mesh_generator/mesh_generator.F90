! 1-D mesh generator for the VKI stagnation-line code
!
!> @author
!> Alessandro Turchi 
!>
!> @brief
!> This tool generates a 1-D mesh in the right format for the VKI stagnation-line code.
!
!>\b DESCRIPTION: \n 
!> This tool generates a 1-D mesh in the right format for the VKI stagnation-line code.
!! It takes as input the followings:
!!    * stagnation point radius (distance from body center): real 
!!    * free-stream boundary radius (distance from body center): real
!!    * number of volumes for the mesh: integer
!!    * stretching function (parabola, hyperbolic tangent and hyperbolic sine are implemented): string
!!    * stretching intensity (~0 for equispaced mesh, >0 to a refinement next to the wall): real
!!         note: the stretching intensity has to be greater than 0 (use 0.000000001 to have an equispaced mesh) 
!!
!! It generates the following output:
!!
!!   * mesh.dat: is the file containing the mesh to be used with the VKI stagnation-line code
!!   * mesh.tec: is a file containing two columns: the cell index and the stencil size. Use this file to check the mesh spacing.

PROGRAM mesh_generator

IMPLICIT NONE

integer, parameter :: prec = selected_real_kind(15, 307)
integer, parameter :: out1 = 11, out2 = 12
integer :: ios
integer :: i
integer :: points, volumes 
real(kind=prec) :: in_radius, out_radius
real(kind=prec) :: delta
real(kind=prec) :: a,b,c
real(kind=prec) :: intens=1.d0 
real(kind=prec), allocatable, dimension(:) :: mesh
real(kind=prec), allocatable, dimension(:) :: stencil
character*20 ::stretch_func
        write(*,*)'Insert stagnation point radius:'
        read(*,*), in_radius

        write(*,*)'Insert outer boundary radius:'
        read(*,*), out_radius

        write(*,*)'Insert number of volumes:'
        read(*,*), volumes

        write(*,*)'Insert stretching function:'
        write(*,*)'(select among "parabola", "hyper_tangent" or "hyper_sine")'
        read(*,*), stretch_func

        points = volumes+1

        allocate(mesh(points))
        allocate(stencil(volumes))
        
        delta=(out_radius-in_radius)/volumes
        do i=1, points
                mesh(i)= (i-1)*delta
        end do

        select case(stretch_func)

                case('parabola')


                        write(*,*)'Insert stretching intensity (0<):'
                        write(*,*)'(higher values mean stronger stretching near the wall)'
                        read(*,*) intens
                        
                        do i=1, points
                                mesh(i)= (i-1)*delta+in_radius
                        end do
                        intens=10.d0**intens
                        intens=1.d0/intens
                        a= (1-intens)*(in_radius-out_radius)/((in_radius*in_radius)-(out_radius*out_radius)+ &
                                (2.d0*in_radius)*(out_radius-in_radius))
                        b= -2.d0*a*in_radius+intens
                        c= in_radius+a*in_radius*in_radius-in_radius*intens
                        mesh = a*mesh*mesh+b*mesh+c

                        
                        write(*,*)'The selected stretching function is a parabola!!!'

                case('hyper_tangent')
                      
                        write(*,*)'Insert stretching intensity (>0):'
                        read(*,*) intens

                        mesh=1.d0+dtanh((mesh/(out_radius-in_radius)-1.d0)*intens)/dtanh(intens)            
                        mesh=mesh*(out_radius-in_radius)+in_radius

                        write(*,*)'The selected stretching function is a hyperbolic tangent!!!'

                case('hyper_sine')
                      
                        write(*,*)'Insert stretching intensity (>0):'
                        read(*,*) intens

                        mesh=dsinh(mesh/(out_radius-in_radius)*intens)/dsinh(intens)            
                        mesh=mesh*(out_radius-in_radius)+in_radius

                        write(*,*)'The selected stretching function is a hyperbolic sine!!!'

                case default
                        
                        write(*,*)'ATTENTION: unknown stretching function!!! The mesh file will not be generated!!'
                        STOP        
        end select


        OPEN(UNIT=out1,FILE='mesh.dat',STATUS='unknown',IOSTAT=ios)

        write(out1,*) points
        write(out1,10) (mesh(i),i=1,points)     


        IF (ios.NE.0) THEN
                PRINT*
                WRITE(*,20)'ATTENTION: something went wrong with the mesh.dat file writing!!!'          
                PRINT*
                STOP
        ENDIF
        close(out1)

        OPEN(UNIT=out2,FILE='mesh.tec',STATUS='unknown',IOSTAT=ios)


        do i=1, volumes
                stencil(i)= mesh(i+1) - mesh(i) 
        end do

        write(out2,30) (real(i), stencil(i),i=1,volumes)     


        IF (ios.NE.0) THEN
                PRINT*
                WRITE(*,20)'ATTENTION: something went wrong with the mesh.tec file writing!!!'          
                PRINT*
                STOP
        ENDIF
        close(out2)

deallocate(mesh)
10  FORMAT(1e15.8)
20 FORMAT(A)
30  FORMAT(2e15.8)
END PROGRAM mesh_generator


