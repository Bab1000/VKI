!------------------------------------------------------------------------------!
! This subroutine provides the geometrical source term (area variation) for a calorically perfect gas flow. 
  SUBROUTINE source_term_quasi1D (cell_id, cell_data, cons, s)

    USE mod_general_data,            ONLY: geom_source, pos_u_cell, pos_h0_cell
 
    IMPLICIT NONE

    REAL(KIND=8) :: dlnadx, rhou_dlnadx
    REAL(KIND=8) :: rhou, u, h0

    INTEGER, INTENT(IN) :: cell_id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cell_data, cons
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s

    ! Logarithmic derivative of area variation
    dlnadx = geom_source (cell_id)
 
    rhou = cons(2)
    u    = cell_data(pos_u_cell) 
    h0   = cell_data(pos_h0_cell)

    ! Useful quantity
    rhou_dlnadx = rhou*dlnadx

    ! Source term 
    s(1) = - rhou_dlnadx 
    s(2) = - rhou_dlnadx*u 
    s(3) = - rhou_dlnadx*h0 

  END SUBROUTINE source_term_quasi1D
!------------------------------------------------------------------------------!
! This subroutine provides the geometrical source term (area variation) and the related Jacobian for a calorically perfect gas flow. 
  SUBROUTINE source_term_quasi1D_Jac (source_id, cell_id, cell_data, cons, s, js)

    USE mod_general_data,            ONLY: gamma, geom_source, pos_u_cell, pos_h0_cell
 
    IMPLICIT NONE

    REAL(KIND=8) :: gamma_minus1, gamma_minus1_v2
    REAL(KIND=8) :: dlnadx, udlnadx, rhou_dlnadx
    REAL(KIND=8) :: h0, rhou, u, v2

    INTEGER, INTENT(IN) :: source_id, cell_id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cell_data, cons
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: js

    ! Logarithmic derivative of nozzle area 
    dlnadx = geom_source (cell_id)

    rhou = cons(2)
    u    = cell_data(pos_u_cell) 
    h0   = cell_data(pos_h0_cell)

    ! Useful quantities 
    v2   = u*u
    udlnadx     = u*dlnadx
    rhou_dlnadx = rhou*dlnadx
    gamma_minus1    = gamma - 1.d0
    gamma_minus1_v2 = gamma_minus1*v2

    ! Source term 
    s(1) = - rhou_dlnadx 
    s(2) = - rhou_dlnadx*u 
    s(3) = - rhou_dlnadx*h0 

    ! Source term Jacobian
    ! First column 
    js(1,1) = 0.d0     
    js(2,1) = v2*dlnadx
    js(3,1) = - (0.5d0*gamma_minus1_v2 - h0)*udlnadx
      
    ! Second column
    js(1,2) = - dlnadx
    js(2,2) = - 2.d0*udlnadx
    js(3,2) = - (h0 - gamma_minus1_v2)*dlnadx           

    ! Third column
    js(1,3) = 0.d0
    js(2,3) = 0.d0
    js(3,3) = - gamma*udlnadx      

  END SUBROUTINE source_term_quasi1D_Jac
!------------------------------------------------------------------------------!
