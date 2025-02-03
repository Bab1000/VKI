module mod_adapt

    implicit none
    
    integer :: nx
    real(kind=8), parameter :: alpha = 0.5
    real(kind=8), parameter :: ratio = 0.01
    
    real(kind=8), allocatable, dimension(:,:) :: mesh
    real(kind=8), allocatable, dimension(:,:) :: new_mesh
    real(kind=8), allocatable, dimension(:) :: prim 
    real(kind=8), allocatable, dimension(:) :: dtdx1
    real(kind=8), allocatable, dimension(:) :: dtdx2
    real(kind=8), allocatable, dimension(:) :: cdf

    contains
    
    subroutine evaluate_derivatives()
        implicit none
        integer :: i
        real(kind=8) :: a, b, x1, x2, x3, y1, y2, y3
        
        ! Derivatives at first point
        x1 = mesh(1,2); y1 = prim(1)
        x2 = mesh(2,2); y2 = prim(2)
        x3 = mesh(3,2); y3 = prim(3)
        
        a = ((y2-y1)*(x1-x3)+(y3-y1)*(x2-x1))/((x1-x3)*(x2*x2-x1*x1)+(x2-x1)*(x3*x3-x1*x1))
        b = ((y2-y1)-a*(x2*x2-x1*x1))/(x2-x1)
    
        dtdx1(1) = 2.0d0*a*x1 + b
        dtdx2(1) = 2.0d0*a
        
        ! Compute derivatives away from boundaries
        do i = 2,nx-2
            x1 = mesh(i-1,2); y1 = prim(i-1)
            x2 = mesh(i,2);   y2 = prim(i)
            x3 = mesh(i+1,2); y3 = prim(i+1)
            
            a = ((y2-y1)*(x1-x3)+(y3-y1)*(x2-x1))/((x1-x3)*(x2*x2-x1*x1)+(x2-x1)*(x3*x3-x1*x1))
            b = ((y2-y1)-a*(x2*x2-x1*x1))/(x2-x1)
            !c = y1-a*x1*x1-b*x1
        
            dtdx1(i) = 2.0d0*a*x2 + b
            dtdx2(i) = 2.0d0*a
        end do
        
        ! Derivatives at last point
        x1 = mesh(nx-3,2); y1 = prim(nx-3)
        x2 = mesh(nx-2,2); y2 = prim(nx-2)
        x3 = mesh(nx-1,2); y3 = prim(nx-1)
        
        a = ((y2-y1)*(x1-x3)+(y3-y1)*(x2-x1))/((x1-x3)*(x2*x2-x1*x1)+(x2-x1)*(x3*x3-x1*x1))
        b = ((y2-y1)-a*(x2*x2-x1*x1))/(x2-x1)
    
        dtdx1(nx-1) = 2.0d0*a*x3 + b
        dtdx2(nx-1) = 2.0d0*a        
    end subroutine
    
    subroutine compute_cdf()
        implicit none
        integer :: i
        real(kind=8) :: f1, f2, m1, m2, mcdf
        
        m1 = maxval(abs(dtdx1))
        m2 = maxval(abs(dtdx2))
        
        cdf(1) = 0.0
        do i = 2,nx-1
            f1 = (1.0d0-alpha)*abs(dtdx1(i-1)/m1) + alpha*abs(dtdx2(i-1)/m2)
            f2 = (1.0d0-alpha)*abs(dtdx1(i)/m1) + alpha*abs(dtdx2(i)/m2)
            cdf(i) = 0.5d0*(f1+f2)
        end do
        
        mcdf = maxval(cdf)
        cdf(1) = (cdf(2)*0.5d0+ratio*10*mcdf)*(mesh(2,1)-mesh(1,1))
        !do i = 2,10
        !    cdf(i) = cdf(i-1) + (cdf(i)+ratio*(11.d0-i)*mcdf)*(mesh(i,1)-mesh(i-1,1))
        !end do
        do i = 2,nx-1
            cdf(i) = cdf(i-1) + (cdf(i))*(mesh(i,1)-mesh(i-1,1)) + &
                     & ratio*mcdf*(1.0d0-(0.5d0*(mesh(i,1)+mesh(i-1,1))-mesh(1,1))/ &
                     & (mesh(nx-1,1)-mesh(1,1))+0.2d0)*(mesh(i,1)-mesh(i-1,1))
        end do
 
        ! Rescale
        cdf = cdf / cdf(nx-1)
   
    end subroutine
    
    subroutine adapt_mesh(mesh_inout, prim_in)
        implicit none
        real(kind=8), dimension(:), intent(inout) :: mesh_inout
        real(kind=8), dimension(:), intent(in) :: prim_in
        integer :: i, j
        real(kind=8) :: x1, x2, y1, y2, x
       
        nx = size(mesh_inout)
        ! Allocate matrices
        if(allocated(mesh)) deallocate(mesh)
        if(allocated(new_mesh)) deallocate(new_mesh)
        if(allocated(prim)) deallocate(prim)
        if(allocated(dtdx1)) deallocate(dtdx1)
        if(allocated(dtdx2)) deallocate(dtdx2)
        if(allocated(cdf)) deallocate(cdf)
        allocate(mesh(nx, 2))
        allocate(new_mesh(nx, 2))
        allocate(prim(nx-1))
        allocate(dtdx1(nx-1))
        allocate(dtdx2(nx-1))
        allocate(cdf(nx-1))
        
        ! Mesh nodes
        mesh(1:nx,1) = mesh_inout(:)
        ! Mesh cell centers
        do i = 1,nx-1
            mesh(i,2) = 0.5d0*(mesh(i,1)+mesh(i+1,1))
        end do

        ! field data to adapt with
        prim = prim_in
        
        ! Compute the derivatives and CDF
        call evaluate_derivatives()
        call compute_cdf()

        ! Interpolate CDF to grid cell
        j = 2
        do i = 2,nx-1
            x = dble(i-1)/dble(nx-1)
            
            do while (cdf(j) < x .and. j < nx-1)
                j = j+1  !PROBLEM HERE
            end do
    
            x1 = cdf(j-1); y1 = mesh(j-1,2)
            x2 = cdf(j);   y2 = mesh(j,2)
            new_mesh(i,1) = (x-x1)/(x1-x2)*(y1-y2)+y1
        end do
  
        ! Boundary points are unchanged
        new_mesh(1,1) = mesh(1,1)
        new_mesh(nx,1) = mesh(nx,1)

        do i = 1, nx-1
            new_mesh(i,2) = 0.5d0*(new_mesh(i,1)+new_mesh(i+1,1))
        enddo

        
        mesh_inout(:) = new_mesh(:,1)
    end subroutine

    subroutine interpolate(field)
        implicit none

        real(kind=8), dimension(:,:), intent(inout) :: field

        real(kind=8), dimension(size(field,1),size(field,2)) :: new_field
        real(kind=8), dimension(size(field,1)) :: v1, v2
        real(kind=8) :: x1, x2, x
        integer :: i,j

        ! Interpolate the field onto the new mesh (cell centers)
        j = 2
        do i = 1,nx-1
            x = new_mesh(i,2)
            do while (mesh(j,2) < x .and. j < nx-1)
                j = j+1
            end do

            x1 = mesh(j-1,2); v1 = field(:,j-1)
            x2 = mesh(j,2);   v2 = field(:,j)

            new_field(:,i) = (x-x1)/(x1-x2)*(v1-v2)+v1
        end do

        field = new_field

    end subroutine

end module mod_adapt
