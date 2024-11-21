
module mod_adapt
    implicit none
    
    integer :: ns
    integer :: nt
    integer :: nx
    
    character(len=100) :: mesh_file
    character(len=100) :: restart_file
    
    real(kind=8), allocatable, dimension(:,:) :: mesh
    real(kind=8), allocatable, dimension(:,:) :: new_mesh
    real(kind=8), allocatable, dimension(:,:) :: prim
    real(kind=8), allocatable, dimension(:,:) :: new_prim
    
    real(kind=8), allocatable, dimension(:) :: dtdx1
    real(kind=8), allocatable, dimension(:) :: dtdx2
    real(kind=8), allocatable, dimension(:) :: cdf

    contains
    
    !! Get the parameters from the input files
    subroutine initialize()
        implicit none        
        integer :: i
        character(len=100) ns_str, nt_str
        
        ! First read in the current grid and temperature vectors
        !write(*,*)'Insert grid file name:'
        !read(*,*), mesh_file
        !write(*,*) 'Insert the restart file:'
        !read(*,*) restart_file
        !write(*,*) 'Insert number of species:'
        !read(*,*) ns
        !write(*,*) 'Insert number of temperatures:'
        !read(*,*) nt
        
        i = command_argument_count()
        if (i .ne. 4) then
            write(*,*) 'usage: adapt mesh restart NS NT'
            stop
        end if
        
        call get_command_argument(1, mesh_file)
        call get_command_argument(2, restart_file)
        call get_command_argument(3, ns_str); read(ns_str,*) ns
        call get_command_argument(4, nt_str); read(nt_str,*) nt
        
        open(10,file=mesh_file)
        read(10,*) nx
        
        ! Allocate matrices
        allocate(mesh(nx, 2))
        allocate(new_mesh(nx, 2))
        allocate(prim(ns+nt+2,nx-1)) ! will be transposed after reading
        allocate(new_prim(nx-1,ns+nt+2))
        allocate(dtdx1(nx-1))
        allocate(dtdx2(nx-1))
        allocate(cdf(nx-1))
        
        ! Read in the mesh
        do i = 1,nx
            read(10,*) mesh(i,1)
        end do
        close(10)
        
        ! Read in the primitive values
        open(10,file=restart_file)
        read(10,*) prim
        prim = transpose(prim)
        
        ! Compute the cell centers
        do i = 1,nx-1
            mesh(i,2) = 0.5d0*(mesh(i,1)+mesh(i+1,1))
        end do
    end subroutine initialize
    
    subroutine evaluate_derivatives()
        implicit none
        integer :: i
        real(kind=8) :: a, b, x1, x2, x3, y1, y2, y3
        
        ! Derivatives at first point
        x1 = mesh(1,2); y1 = prim(1,ns+3)
        x2 = mesh(2,2); y2 = prim(2,ns+3)
        x3 = mesh(3,2); y3 = prim(3,ns+3)
        
        a = ((y2-y1)*(x1-x3)+(y3-y1)*(x2-x1))/((x1-x3)*(x2*x2-x1*x1)+(x2-x1)*(x3*x3-x1*x1))
        b = ((y2-y1)-a*(x2*x2-x1*x1))/(x2-x1)
    
        dtdx1(1) = 2.0d0*a*x1 + b
        dtdx2(1) = 2.0d0*a
        
        ! Compute derivatives away from boundaries
        do i = 2,nx-2
            x1 = mesh(i-1,2); y1 = prim(i-1,ns+3)
            x2 = mesh(i,2);   y2 = prim(i,ns+3)
            x3 = mesh(i+1,2); y3 = prim(i+1,ns+3)
            
            a = ((y2-y1)*(x1-x3)+(y3-y1)*(x2-x1))/((x1-x3)*(x2*x2-x1*x1)+(x2-x1)*(x3*x3-x1*x1))
            b = ((y2-y1)-a*(x2*x2-x1*x1))/(x2-x1)
            !c = y1-a*x1*x1-b*x1
        
            dtdx1(i) = 2.0d0*a*x2 + b
            dtdx2(i) = 2.0d0*a
        end do
        
        ! Derivatives at last point
        x1 = mesh(nx-3,2); y1 = prim(nx-3,ns+3)
        x2 = mesh(nx-2,2); y2 = prim(nx-2,ns+3)
        x3 = mesh(nx-1,2); y3 = prim(nx-1,ns+3)
        
        a = ((y2-y1)*(x1-x3)+(y3-y1)*(x2-x1))/((x1-x3)*(x2*x2-x1*x1)+(x2-x1)*(x3*x3-x1*x1))
        b = ((y2-y1)-a*(x2*x2-x1*x1))/(x2-x1)
    
        dtdx1(nx-1) = 2.0d0*a*x3 + b
        dtdx2(nx-1) = 2.0d0*a        
    end subroutine
    
    subroutine compute_cdf()
        implicit none
        integer :: i
        real(kind=8) :: f1, f2, m1, m2, mcdf
        real(kind=8), parameter :: alpha = 0.5
        real(kind=8), parameter :: ratio = 5.0e-4
        
        m1 = maxval(abs(dtdx1))
        m2 = maxval(abs(dtdx2))
        
        cdf(1) = 0.0
        do i = 2,nx-1
            f1 = (1.0d0-alpha)*abs(dtdx1(i-1)/m1) + alpha*abs(dtdx2(i-1)/m2)
            f2 = (1.0d0-alpha)*abs(dtdx1(i)/m1) + alpha*abs(dtdx2(i)/m2)
            cdf(i) = 0.5d0*(f1+f2)
        end do
        
        mcdf = maxval(cdf)
        do i = 2,nx-1
            cdf(i) = cdf(i-1) + (cdf(i)+ratio*mcdf)*(mesh(i,2)-mesh(i-1,2))
        end do
        
        ! Rescale
        cdf = cdf / cdf(nx-1)
    end subroutine
    
    subroutine adapt_mesh()
        implicit none
        integer :: i, j
        real(kind=8) :: x1, x2, y1, y2, x
        real(kind=8), dimension(ns+nt+2) :: v1, v2
        
        ! Compute the derivatives and CDF
        call evaluate_derivatives()
        call compute_cdf()
        
        ! Interpolate CDF to grid cell
        j = 2
        do i = 2,nx-1
            x = dble(i-1)/dble(nx-1)
            do while (cdf(j) < x .and. j < nx-1)
                j = j+1
            end do
            
            x1 = cdf(j-1); y1 = mesh(j-1,2)
            x2 = cdf(j);   y2 = mesh(j,2)
            
            new_mesh(i,1) = (x-x1)/(x1-x2)*(y1-y2)+y1
        end do
        
        ! Boundary points are unchanged
        new_mesh(1,1) = mesh(1,1)
        new_mesh(nx,1) = mesh(nx,1)
        
        ! Now compute the cell centers
        do i = 1,nx-1
            new_mesh(i,2) = 0.5d0*(new_mesh(i,1)+new_mesh(i+1,1))
        end do
        
        ! Interpolate the primative variables onto the new mesh (cell centers)
        j = 2
        do i = 1,nx-1
            x = new_mesh(i,2)
            do while (mesh(j,2) < x .and. j < nx-1)
                j = j+1
            end do
            
            x1 = mesh(j-1,2); v1 = prim(j-1,:)
            x2 = mesh(j,2);   v2 = prim(j,:)
        
            new_prim(i,:) = (x-x1)/(x1-x2)*(v1-v2)+v1
        end do
    end subroutine

end module mod_adapt

program adapt
    use mod_adapt
    implicit none
    integer :: i, j
    
    call initialize()
    call adapt_mesh()
    
    ! Output the results
    open(10,file=mesh_file)
    write(10,*) nx
    do i = 1,nx
        write(10,*) new_mesh(i,1)
    end do
    close(10)
    
    open(10,file=trim(mesh_file)//'.old')
    write(10,*) nx
    do i = 1,nx
        write(10,*) mesh(i,1)
    end do
    close(10)
    
    open(10,file=restart_file)
    do i = 1,nx-1
        write(10,*) (new_prim(i,j), j=1,ns+nt+2)
    end do
    close(10)
    
    open(10,file=trim(restart_file)//'.old')
    do i = 1,nx-1
        write(10,*) (prim(i,j), j=1,ns+nt+2)
    end do
    close(10)
    
end program adapt


