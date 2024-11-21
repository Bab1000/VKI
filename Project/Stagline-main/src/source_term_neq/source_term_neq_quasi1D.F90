!------------------------------------------------------------------------------!
!> This subroutine provides the area source term for quasi 1D nonequilibrium flows.
  SUBROUTINE source_term_neq_quasi1D (cell_id, cell_data, cons, s)

    USE mod_general_data,            ONLY: nb_eq, nb_ns, nb_int_temp, nb_te, nb_temp, pos_u_cell,  & 
                                         &  pos_h0_cell, pos_rho_cell, pos_rhou, pos_rhoE, geom_source 
 
    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8) :: tmp
    REAL(KIND=8) :: dlnadx, dlnau
    REAL(KIND=8) :: h0, rho, u

    INTEGER, INTENT(IN) :: cell_id                       !< cell identifier
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cell_data  !< physical properties
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons       !< conservative variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s         !< source term

    ! Logarithmic derivative of nozzle area
    dlnadx =  geom_source (cell_id)

    u   = cell_data(pos_u_cell)
    h0  = cell_data(pos_h0_cell)
    rho = cell_data(pos_rho_cell)
 
    ! Common factor 
    dlnau = - dlnadx*u

    ! Source term due to area variation
    DO i = 1,nb_ns
       s(i) = cons(i)*dlnau
    ENDDO
          
    tmp = rho*dlnau  
    s(pos_rhou) = u*tmp
    s(pos_rhoE) = h0*tmp

    ! Internal energy and free electron pseudo-entropy equations 
    DO i = 1,nb_int_temp + nb_te
       s(pos_rhoE + i) = cons(pos_rhoE + i)*dlnau
    ENDDO

  END SUBROUTINE source_term_neq_quasi1D
!------------------------------------------------------------------------------!
! This subroutine provides the area source term and the related jacobian for quasi 1D nonequilibrium flows.
  SUBROUTINE source_term_neq_quasi1D_Jac (source_id, cell_id, cell_data, cons, s, js)

    USE mod_general_data,            ONLY: nb_eq, nb_ns, nb_temp, nb_int_temp, pos_u_cell, pos_h0_cell,  &
                                         & pos_rho_cell, pos_T_cell, pos_ek_cell, pos_alpha_cell,        & 
                                         & pos_ei_cell, pos_ek_cell, pos_eki_cell, pos_beta_cell,        & 
                                         & pos_betak_cell, pos_rhou, pos_rhoE, Ri, epsi, geom_source    

    IMPLICIT NONE

    INTEGER :: i, j, k, kp, pos
    REAL(KIND=8) :: tmp1, tmp2, tmp3, tmp4
    REAL(KIND=8) :: dlnadx, dlnau
    REAL(KIND=8) :: alpha
    REAL(KIND=8) :: ek, h0, rho, T, u, v2

    INTEGER, INTENT(IN) :: cell_id                       !< cell identifier
    INTEGER, INTENT(IN) :: source_id                     !< source term identifier
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cell_data  !< physical properties
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons       !< conservative variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s         !< source term
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: js      !< source term Jacobian

    ! Logarithmic derivative of nozzle area
    dlnadx = geom_source (cell_id)

    T     = cell_data(pos_T_cell)
    u     = cell_data(pos_u_cell)
    h0    = cell_data(pos_h0_cell)
    rho   = cell_data(pos_rho_cell)
    ek    = cell_data(pos_ek_cell)
    alpha = cell_data(pos_alpha_cell)

    ! Common factor 
    dlnau = - dlnadx*u

    ! Source term due to area variation
    DO i = 1,nb_ns
       s(i) = cons(i)*dlnau
    ENDDO
          
    tmp1 = rho*dlnau  
    s(pos_rhou) = u*tmp1
    s(pos_rhoE) = h0*tmp1

    ! Internal energy equations
    DO k = 1,nb_int_temp
       s(pos_rhoE + k) = cons(pos_rhoE + k)*dlnau
    ENDDO

    ! Components of species energy depending on translational temperature only
    DO i = 1,nb_ns 
       tmp1 = cell_data(pos_ei_cell + i - 1)
       DO k = 1,nb_int_temp
          tmp1 = tmp1 - cell_data(pos_eki_cell + (i - 1)*nb_int_temp + k - 1)
       ENDDO
       epsi(i) = - alpha*tmp1 + Ri(i)*T
    ENDDO
    epsi = epsi + alpha*ek

    ! Source term Jacobian due area variation
    ! Useful common factors
    v2   = u*u
    tmp1 = u/rho*dlnadx    
    tmp2 = v2*dlnadx
    tmp3 = dlnadx/rho
    tmp4 = alpha*dlnau 

    ! Columns i,  i = 1,..,nb_ns
    DO j = 1,nb_ns 

       DO i = 1,j - 1 
          js(i,j) = cons(i)*tmp1
       ENDDO

       js(j,j) = dlnau + cons(j)*tmp1

       DO i = j + 1,nb_ns 
          js(i,j) = cons(i)*tmp1
       ENDDO 

       js(pos_rhou,j) = tmp2
       js(pos_rhoE,j) = (epsi(j) - h0)*dlnau

       DO k = 1,nb_int_temp 
          js(pos_rhoE + k,j) = cons(pos_rhoE + k)*tmp1
       ENDDO

    ENDDO 

    ! Column nb_ns + 1
    DO i = 1,nb_ns 
       js(i,pos_rhou) = - cons(i)*tmp3     
    ENDDO

    js(pos_rhou,pos_rhou) = 2.d0*dlnau
    js(pos_rhoE,pos_rhou) = - (h0 - alpha*v2)*dlnadx
 
    DO k = 1,nb_int_temp 
       js(pos_rhoE + k,pos_rhou) = - cons(pos_rhoE + k)*tmp3
    ENDDO

    ! Column nb_ns + 2
    DO i = 1,nb_ns 
       js(i,pos_rhoE) = 0.d0
    ENDDO

    js(pos_rhou,pos_rhoE) = 0.d0
    js(pos_rhoE,pos_rhoE) = (1.d0 + alpha)*dlnau 
 
    DO k = 1,nb_int_temp
       js(pos_rhoE + k,pos_rhoE) = 0.d0
    ENDDO

    ! Column nb_ns + 2 + k,  k = 1,..,nb_int_temp
    DO k = 1,nb_int_temp 

       pos = pos_rhoE + k

       DO i = 1,nb_ns 
          js(i,pos) = 0.d0
       ENDDO
 
       js(pos_rhou,pos) = 0.d0
       js(pos_rhoE,pos) = - tmp4

       DO kp = 1,k - 1
          js(pos_rhoE + kp,pos) = 0.d0  
       ENDDO

       js(pos,pos) = dlnau

       DO kp = k + 1,nb_int_temp
          js(pos_rhoE + kp,pos) = 0.d0
       ENDDO

    ENDDO

  END SUBROUTINE source_term_neq_quasi1D_Jac 
!------------------------------------------------------------------------------!
! This subroutine provides the area source term and the related jacobian for quasi 1D nonequilibrium flows.
! Thermal nonequilibrium between heavy particles and free electrons is accounted for by means 
! of a separate free electron pseudo-entropy equation
  SUBROUTINE source_term_neq_el_quasi1D_Jac (source_id, cell_id, cell_data, cons, s, js)

    USE mod_general_data,            ONLY: nb_eq, nb_ns, nb_temp, nb_int_temp, nb_te,                & 
                                         & pos_u_cell, pos_h0_cell, pos_rho_cell, pos_T_cell,        & 
                                         & pos_ek_cell, pos_alpha_cell, pos_ei_cell, pos_ek_cell,    &                                     
                                         & pos_eki_cell, pos_beta_cell, pos_betak_cell,              &
                                         & pos_gamma_cell, pos_em, pos_rhou, pos_rhoE, pos_rhoek,    & 
                                         & pos_rhose, Ri, epsi, gamma_e, ov_gamma_e_m1, geom_source    

    IMPLICIT NONE

    INTEGER :: i, j, k, kp, pos
    REAL(KIND=8) :: ge_m_g
    REAL(KIND=8) :: tmp1, tmp2, tmp3, tmp4
    REAL(KIND=8) :: dlnadx, dlnau
    REAL(KIND=8) :: alpha, gamma
    REAL(KIND=8) :: ek, h0, rho, T, Te, u, v2, pe

    INTEGER, INTENT(IN) :: cell_id                       !< cell identifier
    INTEGER, INTENT(IN) :: source_id                     !< source term identifier
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cell_data  !< physical properties
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons       !< conservative variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s         !< source term
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: js      !< source term Jacobian
    
    ! Logarithmic derivative of nozzle area
    dlnadx = geom_source (cell_id)

    T     = cell_data(pos_T_cell)
    Te    = cell_data(pos_T_cell + nb_temp - 1)
    u     = cell_data(pos_u_cell)
    h0    = cell_data(pos_h0_cell)
    rho   = cell_data(pos_rho_cell)
    ek    = cell_data(pos_ek_cell)
    alpha = cell_data(pos_alpha_cell)
    gamma = cell_data(pos_gamma_cell)

    ! Common factors 
    dlnau  = - dlnadx*u
    ge_m_g = gamma_e - gamma

    ! Source term due to area variation
    DO i = 1,nb_ns
       s(i) = cons(i)*dlnau
    ENDDO
          
    tmp1 = rho*dlnau  
    s(pos_rhou) = u*tmp1
    s(pos_rhoE) = h0*tmp1

    ! Internal energy conservation equations 
    DO k = 1,nb_int_temp
       s(pos_rhoek + k - 1) = cons(pos_rhoek + k - 1)*dlnau
    ENDDO

    ! Free electron pseudo-entropy equation
    s(pos_rhose) = cons(pos_rhose)*dlnau

    ! Components of species energy depending on translational temperature of heavy particles 
    pe = cons(pos_em)*Ri(pos_em)*Te
    epsi(pos_em) = 0.d0
    DO i = pos_em + 1,nb_ns 
       tmp1 = cell_data(pos_ei_cell + i - 1)
       DO k = 1,nb_int_temp
          tmp1 = tmp1 - cell_data(pos_eki_cell + (i - 1)*nb_int_temp + k - 1)
       ENDDO
       epsi(i) = - alpha*tmp1 + Ri(i)*T
    ENDDO
    epsi = epsi + alpha*ek + pe*ge_m_g/rho

    ! Source term Jacobian due area variation
    ! Useful common factors
    v2   = u*u
    tmp1 = u/rho*dlnadx    
    tmp2 = v2*dlnadx
    tmp3 = dlnadx/rho
    tmp4 = alpha*dlnau

    ! Columns i,  i = 1,..,nb_ns
    DO j = 1,nb_ns 

       DO i = 1,j - 1 
          js(i,j) = cons(i)*tmp1
       ENDDO

       js(j,j) = dlnau + cons(j)*tmp1

       DO i = j + 1,nb_ns 
          js(i,j) = cons(i)*tmp1
       ENDDO 

       js(pos_rhou,j) = tmp2
       js(pos_rhoE,j) = (epsi(j) - h0)*dlnau

       DO k = 1,nb_int_temp
          js(pos_rhoek + k - 1,j) = cons(pos_rhoek + k - 1)*tmp1
       ENDDO

       js(pos_rhose,j) = cons(pos_rhose)*tmp1

    ENDDO 

    ! Column nb_ns + 1
    DO i = 1,nb_ns 
       js(i,pos_rhou) = - cons(i)*tmp3     
    ENDDO

    js(pos_rhou,pos_rhou) = 2.d0*dlnau
    js(pos_rhoE,pos_rhou) = - (h0 - alpha*v2)*dlnadx
 
    DO k = 1,nb_int_temp
       js(pos_rhoek + k - 1,pos_rhou) = - cons(pos_rhoek + k - 1)*tmp3
    ENDDO

    js(pos_rhose,pos_rhou) = - cons(pos_rhose)*tmp3

    ! Column nb_ns + 2
    DO i = 1,nb_ns 
       js(i,pos_rhoE) = 0.d0
    ENDDO

    js(pos_rhou,pos_rhoE) = 0.d0
    js(pos_rhoE,pos_rhoE) = (1.d0 + alpha)*dlnau 
 
    DO k = 1,nb_int_temp
       js(pos_rhoek + k - 1,pos_rhoE) = 0.d0
    ENDDO

    js(pos_rhose,pos_rhoE) = 0.d0

    ! Column nb_ns + 2 + k,  k = 1,..,nb_int_temp
    DO k = 1,nb_int_temp 

       pos = pos_rhoek + k - 1

       DO i = 1,nb_ns 
          js(i,pos) = 0.d0
       ENDDO
 
       js(pos_rhou,pos) = 0.d0
       js(pos_rhoE,pos) = - tmp4

       DO kp = 1,k - 1
          js(pos_rhoek + kp - 1,pos) = 0.d0  
       ENDDO

       js(pos,pos) = dlnau

       DO kp = k + 1,nb_int_temp
          js(pos_rhoek + kp - 1,pos) = 0.d0
       ENDDO

       js(pos_rhose,pos) = 0.d0

    ENDDO

    ! Column nb_ns + 2 + nb_int_temp + nb_te
    DO i = 1,nb_ns 
       js(i,pos_rhose) = 0.d0
    ENDDO
 
    js(pos_rhou,pos_rhose) = 0.d0
    js(pos_rhoE,pos_rhose) = dlnau*ge_m_g*ov_gamma_e_m1*(rho**gamma_e)

    DO k = 1,nb_int_temp
       js(pos_rhoek + k - 1,pos_rhose) = 0.d0  
    ENDDO

    js(pos_rhose,pos_rhose) = dlnau

  END SUBROUTINE source_term_neq_el_quasi1D_Jac 
!------------------------------------------------------------------------------!

