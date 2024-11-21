!------------------------------------------------------------------------------!
!> This subroutine applies a 1D polynomial reconstruction for 1D nonequilibrium flows.
  SUBROUTINE poly_rec_neq_1D (ull, ul, ur, urr, prop_ll, prop_l, prop_r, prop_rr, & 
                            & u_left, u_right, prop_left, prop_right)

    USE mod_general_data,              ONLY: nb_ns, nb_temp, nb_eq, pos_u_cell, pos_T_cell, pos_u, pos_T
    USE mod_function_pointer,          ONLY: limiter, get_cons_phys_from_prim

    IMPLICIT NONE

    INTEGER :: i, pos 
    REAL(KIND=8), PARAMETER :: eps = 1.d-40
    REAL(KIND=8), PARAMETER :: fac = 0.5d0
    REAL(KIND=8), PARAMETER :: omega_l = fac
    REAL(KIND=8), PARAMETER :: omega_r = - fac
    REAL(KIND=8) :: ratio_l, ratio_r, lim_l, lim_r
    REAL(KIND=8) :: drhol, drhor, delta_rho, dTl, dTr, delta_T, dvl, dvr, delta_v
    REAL(KIND=8) :: densl, densr, rhol, rhor, Tl, Tr, vl, vr
    REAL(KIND=8), DIMENSION(nb_eq) :: prim_left, prim_right

    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: ull         !< conservative variables of left-left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: ul          !< conservative variables of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: ur          !< conservative variables of right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: urr         !< conservative variables of right-right state
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: prop_ll     !< physical properties of left-left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: prop_l      !< physical properties of left state
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: prop_r      !< physical properties of right state     
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: prop_rr     !< physical properties of right-right state
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_left      !< conservative variables of reconstructed left state 
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_right     !< conservative variables of reconstructed right state 
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: prop_left   !< physical properties of reconstructed left state 
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: prop_right  !< physical properties of reconstructed right state

    ! Species density reconstruction 
    DO i = 1,nb_ns 

       densl = ul(i)
       densr = ur(i)

       drhol     = densl  - ull(i) + eps
       drhor     = urr(i) - densr  + eps
       delta_rho = densr - densl + eps
           
       ! Ratio of consecutive differences (species densities)
       ratio_l = delta_rho/drhol
       ratio_r = delta_rho/drhor
       
       ! Slope limiters of left and right states
       lim_l = limiter(ratio_l)
       lim_r = limiter(ratio_r)

       ! Reconstructed species density
       prim_left(i)  = densl + omega_l*lim_l*drhol 
       prim_right(i) = densr + omega_r*lim_r*drhor

    ENDDO

    ! Velocity reconstruction
    vl = prop_l(pos_u_cell)
    vr = prop_r(pos_u_cell)

    dvl     = vl  - prop_ll(pos_u_cell) + eps
    dvr     = prop_rr(pos_u_cell) - vr  + eps
    delta_v = vr - vl + eps

    ! Ratio of consecutive differences (velocity)
    ratio_l = delta_v/dvl
    ratio_r = delta_v/dvr

    ! Slope limiters of left and right states
    lim_l = limiter(ratio_l)
    lim_r = limiter(ratio_r)
     
    ! Reconstructed velocities
    prim_left(pos_u)  = vl + omega_l*lim_l*dvl 
    prim_right(pos_u) = vr + omega_r*lim_r*dvr

    ! Temperature(s)
    DO i = 1,nb_temp

       pos = pos_T_cell + i - 1

       Tl = prop_l(pos)
       Tr = prop_r(pos)

       dTl     = Tl - prop_ll(pos) + eps
       dTr     = prop_rr(pos) - Tr + eps
       delta_T = Tr - Tl + eps

       ! Ratio of consecutive differences (temperatures)
       ratio_l = delta_T/dTl
       ratio_r = delta_T/dTr

       ! Slope limiters of left and right states
       lim_l = limiter(ratio_l)
       lim_r = limiter(ratio_r)

       ! Reconstructed temperature
       prim_left(pos_T + i - 1)  = Tl + omega_l*lim_l*dTl
       prim_right(pos_T + i - 1) = Tr + omega_r*lim_r*dTr 

    ENDDO
   
    ! Get the conservative variables and the physical properties from the reconstructed primitive variables 
    ! of left and and right states
    CALL get_cons_phys_from_prim(prim_left, u_left, prop_left)
    CALL get_cons_phys_from_prim(prim_right, u_right, prop_right)
 
  END SUBROUTINE poly_rec_neq_1D
!------------------------------------------------------------------------------!
