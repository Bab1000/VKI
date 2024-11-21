!------------------------------------------------------------------------------!
!> This subroutine computes the inviscid flux Jacobian for the 1D stagnation line Euler equations
!! for a calorically perfect gas.
  SUBROUTINE inviscid_flux_Jac_1D_SL (nx, phys_data, a) 

    USE mod_general_data,         ONLY: pos_u_cell, pos_v_cell, pos_h0_cell, pos_c_cell, pos_ek_cell, gamma

    IMPLICIT NONE

    REAL(KIND=8) :: gamma_minus1, gamma_plus1, gamma_minus3, gamma_minus1_v2
    REAL(KIND=8) :: c, ek, h0, u, v, vn

    REAL(KIND=8), INTENT(IN) :: nx                        !< normal to the cell interface
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data   !< physical properties
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: a        !< inviscid flux Jacobian

    ! Velocity components, total enthalpy and speed of sound
    u  = phys_data(pos_u_cell)
    v  = phys_data(pos_v_cell)
    c  = phys_data(pos_c_cell)
    ek = phys_data(pos_ek_cell)
    h0 = phys_data(pos_h0_cell)  
    vn = u*nx 
  
    ! Common factors
    gamma_plus1     = gamma + 1.d0 
    gamma_minus1    = gamma - 1.d0 
    gamma_minus3    = gamma - 3.d0  
    gamma_minus1_v2 = 2.d0*ek*gamma_minus1 

    ! Inviscid flux Jacobian 
    ! First column
    a(1,1) = 0.d0
    a(2,1) = 0.5d0*gamma_minus3*u*vn 
    a(3,1) = - v*vn
    a(4,1) = vn*(0.5d0*gamma_minus1_v2 - h0)

    ! Second column
    a(1,2) = nx
    a(2,2) = - gamma_minus3*vn
    a(3,2) = v*nx
    a(4,2) = (h0 - gamma_minus1_v2)*nx
        
    ! Third column
    a(1,3) = 0.d0
    a(2,3) = 0.d0
    a(3,3) = vn
    a(4,3) = 0.d0

    !Fourth column
    a(1,4) = 0.d0
    a(2,4) = gamma_minus1*nx
    a(3,4) = 0.d0
    a(4,4) = gamma*vn 

   END SUBROUTINE inviscid_flux_Jac_1D_SL
!------------------------------------------------------------------------------!
