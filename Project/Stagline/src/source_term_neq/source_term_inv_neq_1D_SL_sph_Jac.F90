!-----------------------------------------------------------------------------!
!> This subroutine computes the inviscid source term and the related Jacobian for 1D stagnation line 
!! nonequilibrium flows (sphere case).
  SUBROUTINE source_term_inv_neq_1D_SL_sph_Jac (cell_id, source_id, r, prop, cons, s, js)

    USE mod_general_data,            ONLY: nb_ns, nb_temp, pos_rho_cell, pos_u_cell, pos_v_cell, pos_h0_cell,   & 
                                         & pos_pres_cell, pos_ei_cell, pos_ek_cell, pos_alpha_cell, pos_T_cell,  & 
                                         & pos_rhou, pos_rhov, pos_rhoE, p_inf, yi, epsi, Ri, nb_eq,            &
                                         & nb_int_temp, pos_eki_cell, eintk
    USE mod_neq_function_pointer,    ONLY: library_get_mass_fractions

    IMPLICIT NONE

    INTEGER :: i, j
    REAL(KIND=8) :: alpha, ek, rho, h0, u, v, p, T, tmp1
    REAL(KIND=8) :: ov_r, u_ov_r, v_ov_r, m2_h0_ov_r
    REAL(KIND=8) :: m2_vel_sum_ov_r, m3_vel_sum_ov_r, m2_rho_vel_sum_ov_r, m3_rho_vel_sum_ov_r
    REAL(KIND=8) :: m2_vel_sum_ov_r_u, m3_vel_sum_ov_r_v, m2_vel_sum_ov_r_h0
    REAL(KIND=8), DIMENSION(nb_ns) :: m2_yi_ov_r

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
    u_ov_r = u*ov_r
    v_ov_r = v*ov_r
    m2_h0_ov_r = 2.d0*h0*ov_r
    m2_vel_sum_ov_r = 2.d0*(u + v)*ov_r
    m2_vel_sum_ov_r_h0 = m2_vel_sum_ov_r*h0
    m3_vel_sum_ov_r = 1.5d0*m2_vel_sum_ov_r
    m2_vel_sum_ov_r_u = m2_vel_sum_ov_r*u
    m3_vel_sum_ov_r_v = m3_vel_sum_ov_r*v
    m2_rho_vel_sum_ov_r = rho*m2_vel_sum_ov_r
    m3_rho_vel_sum_ov_r = 1.5d0*m2_rho_vel_sum_ov_r

    ! Species mass fractions
    CALL library_get_mass_fractions(rho, cons(1:nb_ns), yi)

    DO i = 1,nb_ns 
       m2_yi_ov_r(i) = 2.d0*yi(i)*ov_r
       epsi(i) = Ri(i)*T - alpha*prop(pos_ei_cell + i - 1) 
    ENDDO
    epsi = epsi + alpha*ek

    DO i = 1, nb_int_temp
      eintk(i) = 0.d0
      DO j = 1, nb_ns
        tmp1 = prop(pos_eki_cell + j - 1 +nb_ns*(i-1))
        epsi(j) = epsi(j) - alpha*tmp1
        eintk(i) = eintk(i) + tmp1
      ENDDO
    ENDDO

    ! Source term
    ! Species continuity equations
    DO i = 1,nb_ns 
       s(i) = - m2_vel_sum_ov_r*cons(i)
    ENDDO

    ! Radial and circumferential momentum equations
    s(pos_rhou) = - m2_rho_vel_sum_ov_r*u
    s(pos_rhov) = - (m3_rho_vel_sum_ov_r*v - 2.d0*(p - p_inf)*ov_r)

    ! Global energy equation
    s(pos_rhoE) = - m2_rho_vel_sum_ov_r*h0

    ! Internal energy and free-electron pseudo-entropy equations
    DO i = 2,nb_temp
       s(pos_rhoE + i - 1) = - m2_vel_sum_ov_r*cons(pos_rhoE + i - 1)
    ENDDO

    ! Source term Jacobian 
    ! Column i,  i = 1,..,nb_ns 
    DO j = 1,nb_ns 

       DO i = 1,j - 1
          js(i,j) = m2_vel_sum_ov_r*yi(i)
       ENDDO

       js(j,j) = m2_vel_sum_ov_r*(yi(j) - 1.d0)

       DO i = j + 1,nb_ns 
          js(i,j) = m2_vel_sum_ov_r*yi(i)
       ENDDO

       js(pos_rhou,j) = m2_vel_sum_ov_r_u
       js(pos_rhov,j) = m3_vel_sum_ov_r_v + 2.d0*ov_r*epsi(j)
       js(pos_rhoE,j) = m2_vel_sum_ov_r_h0 - m2_vel_sum_ov_r*epsi(j)

       DO i = 1,nb_int_temp
          js(pos_rhoE+i,j) = m2_vel_sum_ov_r*eintk(i)
       ENDDO

    ENDDO

    ! Column nb_ns + 1
    DO i = 1,nb_ns 
       js(i,pos_rhou) = - m2_yi_ov_r(i)
    ENDDO

    js(pos_rhou,pos_rhou) = - (m2_vel_sum_ov_r + 2.d0*u_ov_r)
    js(pos_rhov,pos_rhou) = - (3.d0*v_ov_r + 2.d0*alpha*u_ov_r)
    js(pos_rhoE,pos_rhou) = - m2_h0_ov_r  + alpha*m2_vel_sum_ov_r_u 

    DO i = 1,nb_int_temp
       js(pos_rhoE+i,pos_rhou) = - 2.0*ov_r*eintk(i)
    ENDDO

    ! Column nb_ns + 2
    DO i = 1,nb_ns 
       js(i,pos_rhov) = - m2_yi_ov_r(i)
    ENDDO

    js(pos_rhou,pos_rhov) = - 2.d0*u_ov_r
    js(pos_rhov,pos_rhov) = - m3_vel_sum_ov_r - 3.d0*v_ov_r
    js(pos_rhoE,pos_rhov) = - m2_h0_ov_r

    DO i = 1,nb_int_temp
       js(pos_rhoE+i,pos_rhov) = - 2.d0*ov_r*eintk(i)
    ENDDO

    ! Column nb_ns + 3
    DO i = 1,nb_ns 
       js(i,pos_rhoE) = 0.d0
    ENDDO
   
    js(pos_rhou,pos_rhoE) = 0.d0
    js(pos_rhov,pos_rhoE) = 2.d0*alpha*ov_r
    js(pos_rhoE,pos_rhoE) = - m2_vel_sum_ov_r*(1.d0 + alpha)

    DO i = 1,nb_int_temp
       js(pos_rhoE+i,pos_rhoE) = 0.d0
    ENDDO

    ! Column j,  j = nb_ns+4,..,nb_eq 
    DO j = 1, nb_int_temp

       DO i = 1, nb_ns
          js(i, pos_rhoE+j) = 0.d0
       ENDDO

       js(pos_rhou, pos_rhoE+j) = 0.d0
       js(pos_rhov, pos_rhoE+j) = -2.0*alpha*ov_r
       js(pos_rhoE, pos_rhoE+j) = m2_vel_sum_ov_r*alpha

       DO i = 1,j - 1
          js(pos_rhoE+i,pos_rhoE+j) = 0.d0
       ENDDO

       js(pos_rhoE+j,pos_rhoE+j) = - m2_vel_sum_ov_r

       DO i = j + 1, nb_int_temp
          js(pos_rhoE+i,pos_rhoE+j) = 0.d0
       ENDDO

    ENDDO

  END SUBROUTINE source_term_inv_neq_1D_SL_sph_Jac
!-----------------------------------------------------------------------------!
