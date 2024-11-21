!-----------------------------------------------------------------------------!
!> This subroutine computes the inviscid source term and the related Jacobian for 1D stagnation line 
!! nonequilibrium flows (cylinder case).
  SUBROUTINE source_term_inv_neq_1D_SL_cyl_Jac (cell_id, source_id, r, prop, cons, s, js)

    USE mod_general_data,            ONLY: nb_ns, nb_temp, pos_rho_cell, pos_u_cell, pos_v_cell, pos_h0_cell,   & 
                                         & pos_pres_cell, pos_ei_cell, pos_ek_cell, pos_alpha_cell, pos_T_cell,  & 
                                         & pos_rhou, pos_rhov, pos_rhoE, p_inf, yi, epsi, Ri, nb_eq
    USE mod_neq_function_pointer,    ONLY: library_get_mass_fractions

    IMPLICIT NONE

    INTEGER :: i, j
    REAL(KIND=8) :: alpha, ek, rho, h0, u, v, p, T
    REAL(KIND=8) :: ov_r, u_ov_r, v_ov_r, h0_ov_r
    REAL(KIND=8) :: vel_sum_ov_r, rho_vel_sum_ov_r, vel_sum_ov_r_u, vel_sum_ov_r_v, vel_sum_ov_r_h0 
    REAL(KIND=8), DIMENSION(nb_ns) :: yi_ov_r

    INTEGER, INTENT(IN) :: cell_id                   !< cell identifier
    INTEGER, INTENT(IN) :: source_id                 !< source term identifier
    REAL(KIND=8), INTENT(IN) :: r                    !< cell centroid radial location
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: prop   !< physical properties
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons   !< conservative variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s     !< source term
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: js  !< source term Jacobian  
 
    ! Temperature, density, velocity components, kinetic energy per unit mass, pressure and specific total enthalpy
    T     = prop(pos_T_cell)
    rho   = prop(pos_rho_cell)
    u     = prop(pos_u_cell)
    v     = prop(pos_v_cell)
    p     = prop(pos_pres_cell)
    h0    = prop(pos_h0_cell)
    alpha = prop(pos_alpha_cell)
    ek    = prop(pos_ek_cell)

    ! Common factors
    ov_r = 1.d0/r
    u_ov_r  = u*ov_r
    v_ov_r  = v*ov_r
    h0_ov_r = h0*ov_r
    vel_sum_ov_r   = (u + v)*ov_r
    vel_sum_ov_r_u = vel_sum_ov_r*u
    vel_sum_ov_r_v = vel_sum_ov_r*v
    vel_sum_ov_r_h0 = vel_sum_ov_r*h0
    rho_vel_sum_ov_r = rho*vel_sum_ov_r

    ! Species mass fractions
    CALL library_get_mass_fractions(rho, cons(1:nb_ns), yi)

    DO i = 1,nb_ns 
       yi_ov_r(i) = yi(i)*ov_r
       epsi(i) = Ri(i)*T - alpha*prop(pos_ei_cell + i - 1) 
    ENDDO
    epsi = epsi + alpha*ek

    ! Source term
    ! Species continuity equations
    DO i = 1,nb_ns 
       s(i) = - vel_sum_ov_r*cons(i)
    ENDDO

    ! Radial and circumferential momentum equations
    s(pos_rhou) = - rho_vel_sum_ov_r*u
    s(pos_rhov) = - 2.d0*(rho_vel_sum_ov_r*v - (p - p_inf)*ov_r)

    ! Global energy equation
    s(pos_rhoE) = - rho_vel_sum_ov_r*h0

    ! Internal energy and free-electron pseudo-entropy equations
    DO i = 2,nb_temp
       s(pos_rhoE + i - 1) = - vel_sum_ov_r*cons(pos_rhoE + i - 1)
    ENDDO

    ! Source term Jacobian
    ! Column i,  i = 1,..,nb_ns
    DO j = 1,nb_ns 

       DO i = 1,j - 1
          js(i,j) = vel_sum_ov_r*yi(i)
       ENDDO

       js(j,j) = vel_sum_ov_r*(yi(j) - 1.d0)

       DO i = j + 1,nb_ns 
          js(i,j) = vel_sum_ov_r*yi(i)
       ENDDO

       js(pos_rhou,j) = vel_sum_ov_r_u
       js(pos_rhov,j) = - 2.d0*epsi(j)*ov_r + 2.d0*vel_sum_ov_r_v
       js(pos_rhoE,j) = vel_sum_ov_r_h0 - vel_sum_ov_r*epsi(j)

       DO i = 2,nb_temp
          js(pos_rhoE+i-1,j) = 0.d0
       ENDDO

    ENDDO

    ! Column nb_ns + 1 
    DO i = 1,nb_ns 
       js(i,pos_rhou) = - yi_ov_r(i)
    ENDDO

    js(pos_rhou,pos_rhou) = - (vel_sum_ov_r + u_ov_r)
    js(pos_rhov,pos_rhou) = - 2.d0*(v_ov_r + alpha*u_ov_r)
    js(pos_rhoE,pos_rhou) = alpha*vel_sum_ov_r_u - h0_ov_r

    DO i = 2,nb_temp
       js(pos_rhoE+i-1,pos_rhou) = 0.d0
    ENDDO

    ! Column nb_ns + 2
    DO i = 1,nb_ns 
       js(i,pos_rhov) = - yi_ov_r(i)
    ENDDO

    js(pos_rhou,pos_rhov) = - u_ov_r
    js(pos_rhov,pos_rhov) = - 2.d0*(v_ov_r + vel_sum_ov_r)
    js(pos_rhoE,pos_rhov) = - h0_ov_r

    DO i = 2,nb_temp
       js(pos_rhoE+i-1,pos_rhov) = 0.d0
    ENDDO

    ! Column nb_ns + 3 
    DO i = 1,nb_ns 
       js(i,pos_rhoE) = 0.d0
    ENDDO

    js(pos_rhou,pos_rhoE) = 0.d0
    js(pos_rhov,pos_rhoE) = - 2.d0*alpha*ov_r
    js(pos_rhoE,pos_rhoE) = - vel_sum_ov_r*(alpha + 1)

    DO i = 2,nb_temp
       js(pos_rhoE+i-1,pos_rhoE) = 0.d0
    ENDDO

    ! Column j,  j = nb_ns+4,..,nb_eq 
    DO j = nb_ns+4, nb_eq

       DO i = 1,j - 1
          js(i,j) = 0.d0
       ENDDO

       js(j,j) = - vel_sum_ov_r

       DO i = j + 1,nb_eq
          js(i,j) = 0.d0
       ENDDO

    ENDDO

  END SUBROUTINE source_term_inv_neq_1D_SL_cyl_Jac
!-----------------------------------------------------------------------------!
