!-----------------------------------------------------------------------------!
!> This subroutine computes the inviscid source term for 1D stagnation line calorically perfect gas flows (sphere case).
  SUBROUTINE source_term_inv_1D_SL_sph (cell_id, r, prop, cons, s)

    USE mod_general_data,        ONLY: pos_u_cell, pos_v_cell, pos_pres_cell, pos_h0_cell,   & 
                                     & pos_rho_cell, pos_rho, pos_rhou, pos_rhov, pos_rhoE,  & 
                                     & p_inf

    IMPLICIT NONE

    REAL(KIND=8) :: rho, h0, p, u, v
    REAL(KIND=8) :: ov_r, m2_rho_vel_sum_ov_r, m3_rho_vel_sum_ov_r

    INTEGER, INTENT(IN) :: cell_id                  !< cell identifier
    REAL(KIND=8), INTENT(IN) :: r                   !< cell centroid radial location
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prop  !< physical properties
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons  !< conservative variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s    !< source term 
 
    ! Density, velocity components, pressure and specific total enthalpy
    rho = prop(pos_rho_cell)
    u   = prop(pos_u_cell)
    v   = prop(pos_v_cell)
    p   = prop(pos_pres_cell)
    h0  = prop(pos_h0_cell)

    ! Common factors
    ov_r = 1.d0/r
    m2_rho_vel_sum_ov_r = 2.d0*rho*(u + v)*ov_r
    m3_rho_vel_sum_ov_r = 1.5d0*m2_rho_vel_sum_ov_r

    ! Source term
    s(pos_rho)  = - m2_rho_vel_sum_ov_r
    s(pos_rhou) = - m2_rho_vel_sum_ov_r*u
    s(pos_rhov) = - (m3_rho_vel_sum_ov_r*v - 2.d0*(p - p_inf)*ov_r)
    s(pos_rhoE) = - m2_rho_vel_sum_ov_r*h0

  END SUBROUTINE source_term_inv_1D_SL_sph
!-----------------------------------------------------------------------------!
