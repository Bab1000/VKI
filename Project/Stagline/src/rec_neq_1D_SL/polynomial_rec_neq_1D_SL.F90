!------------------------------------------------------------------------------!
!> This subroutine applies a 1D polynomial reconstruction for 1D stagnation line nonequilibrium flows.
  SUBROUTINE poly_rec_neq_1D_SL (vol_ll, vol_l, vol_r, vol_rr, ull, ul, ur, urr, prop_ll, prop_l, prop_r, prop_rr, & 
                               & u_left, u_right, prop_left, prop_right)

    USE mod_general_data,              ONLY: nb_ns, nb_temp, nb_eq, nb_dim, pos_u_cell, pos_T_cell, pos_u, pos_T
    USE mod_function_pointer,          ONLY: limiter, get_cons_phys_from_prim

    IMPLICIT NONE

    INTEGER :: i, pos 
    REAL(KIND=8), PARAMETER :: eps = 1.d-40
    REAL(KIND=8), PARAMETER :: fac = 0.5d0
    REAL(KIND=8), PARAMETER :: lambda = 0.2

    REAL(KIND=8), PARAMETER :: omega_l = fac
    REAL(KIND=8), PARAMETER :: omega_r = - fac
    REAL(KIND=8) :: ratio_l, ratio_r, lim_l, lim_r
    REAL(KIND=8) :: drhol, drhor, delta_rho, dTl, dTr, delta_T, dvl, dvr, dv
    REAL(KIND=8) :: densl, densr, rhol, rhor, Tl, Tr, vll, vl, vr, vrr
    REAL(KIND=8), DIMENSION(nb_eq) :: prim_left, prim_right

    REAL(KIND=8),               INTENT(IN)  :: vol_l 
    REAL(KIND=8),               INTENT(IN)  :: vol_r
    REAL(KIND=8),               INTENT(IN)  :: vol_ll
    REAL(KIND=8),               INTENT(IN)  :: vol_rr
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
      ! prim_left(i)  = densl + omega_l*lim_l*drhol 
      ! prim_right(i) = densr + omega_r*lim_r*drhor
       prim_left(i)  = lambda*densl + (1.D0-lambda)*(densl + omega_l*lim_l*drhol)
       prim_right(i) = lambda*densr + (1.D0-lambda)*(densr + omega_r*lim_r*drhor)
       if((abs(drhol) < eps).OR.(abs(delta_rho) < eps)) prim_left(i) = densl
       if((abs(drhor) < eps).OR.(abs(delta_rho) < eps)) prim_right(i) = densr
    ENDDO

    ! Radial and circumferential velocity
    DO i = 1,nb_dim

       pos = pos_u_cell + i - 1

       vll = prop_ll(pos)
       vl  = prop_l(pos)
       vr  = prop_r(pos)
       vrr = prop_rr(pos)

       ! Jumps
       dv  = vr  - vl + eps
       dvl = vl  - vll + eps
       dvr = vrr - vr + eps  

       ! Ratios of consecutive differences
       ratio_l = dv/dvl
       ratio_r = dv/dvr  

       ! Gradient limiting
       lim_l = limiter(ratio_l)
       lim_r = limiter(ratio_r)

       ! Reconstruction
       !prim_left(pos_u + i - 1)  = vl + omega_l*lim_l*dvl
       !prim_right(pos_u + i - 1) = vr + omega_r*lim_r*dvr
       prim_left(pos_u + i - 1)  =  lambda*vl + (1.D0-lambda)*(vl + omega_l*lim_l*dvl)
       prim_right(pos_u + i - 1) =  lambda*vr + (1.D0-lambda)*(vr + omega_r*lim_r*dvr)
       if((abs(dvl) < eps).OR.(abs(dv) < eps)) prim_left(pos_u + i - 1) = vl
       if((abs(dvr) < eps).OR.(abs(dv) < eps)) prim_right(pos_u + i -1) = vr

    ENDDO

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
      ! prim_left(pos_T + i - 1)  = Tl + omega_l*lim_l*dTl
      ! prim_right(pos_T + i - 1) = Tr + omega_r*lim_r*dTr 
       prim_left(pos_T + i - 1)  = lambda*Tl + (1.D0-lambda)*(Tl + omega_l*lim_l*dTl)
       prim_right(pos_T + i - 1) = lambda*Tr +(1.D0-lambda)*(Tr + omega_r*lim_r*dTr)

    ENDDO

    ! Get the conservative variables and the physical properties from the reconstructed primitive variables 
    ! of left and and right states
    CALL get_cons_phys_from_prim(prim_left, u_left, prop_left)
    CALL get_cons_phys_from_prim(prim_right, u_right, prop_right)
      
  END SUBROUTINE poly_rec_neq_1D_SL
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!> This subroutine applies a 1D polynomial reconstruction for 1D stagnation line nonequilibrium flows, accounting for metrics (A. Turchi)
  SUBROUTINE poly_rec_neq_1D_SL_metr (vol_ll, vol_l, vol_r, vol_rr, ull, ul, ur, urr, prop_ll, prop_l, prop_r, prop_rr, & 
                               & u_left, u_right, prop_left, prop_right)

    USE mod_general_data,              ONLY: nb_ns, nb_temp, nb_eq, nb_dim, pos_u_cell, pos_T_cell, pos_u, pos_T
    USE mod_function_pointer,          ONLY: limiter, get_cons_phys_from_prim

    IMPLICIT NONE

    INTEGER :: i, pos, j 
    REAL(KIND=8), PARAMETER :: eps = 1.d-40
    REAL(KIND=8), PARAMETER :: lambda = 0.0

    REAL(KIND=8) :: ratio_l, ratio_r, lim_l, lim_r
    REAL(KIND=8) :: drhol, drhor, delta_rho, dTl, dTr, delta_T, dvl, dvr, dv
    REAL(KIND=8) :: densl, densr, rhol, rhor, Tl, Tr, vll, vl, vr, vrr
    REAL(KIND=8) :: vol_sum_l, vol_sum, vol_sum_r 
    REAL(KIND=8) :: omega_l 
    REAL(KIND=8) :: omega_r 
    REAL(KIND=8), DIMENSION(nb_eq) :: prim_left, prim_right

    REAL(KIND=8),               INTENT(IN)  :: vol_l 
    REAL(KIND=8),               INTENT(IN)  :: vol_r
    REAL(KIND=8),               INTENT(IN)  :: vol_ll
    REAL(KIND=8),               INTENT(IN)  :: vol_rr
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

    vol_sum_l = vol_ll + vol_l
    vol_sum   = vol_l + vol_r
    vol_sum_r = vol_r + vol_rr

    omega_l =   vol_l / vol_sum_l 
    omega_r = - vol_r / vol_sum_r

    ! Species density reconstruction 
    DO i = 1,nb_ns 

       densl = ul(i)
       densr = ur(i)

       drhol     = densl  - ull(i) + eps
       drhor     = urr(i) - densr  + eps
       delta_rho = densr - densl + eps
          
       ! Ratio of consecutive differences (species densities)
       ratio_l = (delta_rho * vol_sum_l )  / (drhol * vol_sum)
       ratio_r = (delta_rho * vol_sum_r )  / (drhor * vol_sum)
       
       ! Slope limiters of left and right states
       lim_l = limiter(ratio_l)
       lim_r = limiter(ratio_r)

       ! Reconstructed species density
       prim_left(i)  = lambda*densl + (1.D0-lambda)*(densl + omega_l*lim_l*drhol)
       prim_right(i) = lambda*densr + (1.D0-lambda)*(densr + omega_r*lim_r*drhor)
       if((abs(drhol) < eps).OR.(abs(delta_rho) < eps)) prim_left(i) = densl
       if((abs(drhor) < eps).OR.(abs(delta_rho) < eps)) prim_right(i) = densr

    ENDDO

    ! Radial and circumferential velocity
    DO i = 1,nb_dim

       pos = pos_u_cell + i - 1

       vll = prop_ll(pos)
       vl  = prop_l(pos)
       vr  = prop_r(pos)
       vrr = prop_rr(pos)

       ! Jumps
       dv  = vr  - vl + eps
       dvl = vl  - vll + eps
       dvr = vrr - vr + eps  

       ! Ratios of consecutive differences
       ratio_l = (dv * vol_sum_l) / (dvl * vol_sum)
       ratio_r = (dv * vol_sum_r) / (dvr * vol_sum) 

       ! Gradient limiting
       lim_l = limiter(ratio_l)
       lim_r = limiter(ratio_r)

       ! Reconstruction
       prim_left(pos_u + i - 1)  =  lambda*vl + (1.D0-lambda)*(vl + omega_l*lim_l*dvl)
       prim_right(pos_u + i - 1) =  lambda*vr + (1.D0-lambda)*(vr + omega_r*lim_r*dvr)
       if((abs(dvl) < eps).OR.(abs(dv) < eps)) prim_left(pos_u + i - 1) = vl
       if((abs(dvr) < eps).OR.(abs(dv) < eps)) prim_right(pos_u + i -1) = vr

    ENDDO

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
       prim_left(pos_T + i - 1)  = lambda*Tl + (1.D0-lambda)*(Tl + omega_l*lim_l*dTl)
       prim_right(pos_T + i - 1) = lambda*Tr +(1.D0-lambda)*(Tr + omega_r*lim_r*dTr) 
       if((abs(dTl) < eps).OR.(abs(delta_T) < eps)) prim_left(pos_T + i - 1) = Tl
       if((abs(dTr) < eps).OR.(abs(delta_T) < eps)) prim_right(pos_T + i - 1) = Tr

    ENDDO
   
    ! Get the conservative variables and the physical properties from the reconstructed primitive variables 
    ! of left and and right states
    CALL get_cons_phys_from_prim(prim_left, u_left, prop_left)
    CALL get_cons_phys_from_prim(prim_right, u_right, prop_right)


  END SUBROUTINE poly_rec_neq_1D_SL_metr
!------------------------------------------------------------------------------!
