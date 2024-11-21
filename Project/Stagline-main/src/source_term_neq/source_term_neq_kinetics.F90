!------------------------------------------------------------------------------!
!> This subroutine provides the source terms due to collisional and radiative
!! processes for nonequilibrium flows.
  SUBROUTINE source_term_neq_kinetics (cell_id, cell_data, cons, s)

    USE mod_general_data,            ONLY: nb_ns, nb_temp 
    USE mod_neq_function_pointer,    ONLY: library_get_source 

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: cell_id                       !< cell identifier
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cell_data  !< physical properties 
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons       !< conservative variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s         !< source term

    ! Calling the thermodynamic library
    CALL library_get_source (cons(1:nb_ns), cell_data(1:nb_temp), s)

  END SUBROUTINE source_term_neq_kinetics
!------------------------------------------------------------------------------!
!> This subroutine provides the source terms due to collisional and radiative
!! processes for nonequilibrium flows. The source term for the free electron
!! energy is multiplied by a factor in order to obtain the correspondent quantity
!! in the free-electron pseudo-entropy equation. 
  SUBROUTINE source_term_neq_el_kinetics (cell_id, cell_data, cons, s)

    USE mod_general_data,            ONLY: nb_ns, nb_te, nb_temp, nb_eq, pos_rho_cell, gamma_e_m1 
    USE mod_neq_function_pointer,    ONLY: library_get_source 

    IMPLICIT NONE

    REAL(KIND=8) :: rho

    INTEGER, INTENT(IN) :: cell_id                       !< cell identifier
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cell_data  !< physical properties 
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons       !< conservative variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s         !< source term 
 
    ! Calling the thermodynamic library
    CALL library_get_source (cons(1:nb_ns), cell_data(1:nb_temp), s)

    ! Density 
    rho = cell_data(pos_rho_cell)

    ! Source term in the free-electron pseudo-entropy equation
    s(nb_eq) = s(nb_eq)*(gamma_e_m1)/(rho**gamma_e_m1)

  END SUBROUTINE source_term_neq_el_kinetics
!------------------------------------------------------------------------------!
!> This subroutine provides the source term and the related Jacobian due to collisional and radiative
!! processes for nonequilibrium flows.
  SUBROUTINE source_term_neq_kinetics_Jac (source_id, cell_id, cell_data, cons, s, js)

    USE mod_general_data,            ONLY: nb_ns, nb_dim, nb_temp, nb_int_temp, pos_u_cell, pos_ek_cell, pos_beta_cell, & 
                                         & pos_betak_cell, pos_ei_cell, pos_eki_cell, pos_T_cell, pos_u_cell,           & 
                                         & pos_rhou, pos_rhoE, pos_rhoek, rhoi, temp, eiT, ov_betak, vel
    USE mod_numerics_data,           ONLY: js_line
    USE mod_neq_function_pointer,    ONLY: library_get_source_Jac

    IMPLICIT NONE

    INTEGER :: i, j, d, l, k, p
    INTEGER :: dim_js, ns_p1, pos, pos_d, pos_k
    REAL(KIND=8) :: fac, tmp, tmp1, tmp2
    REAL(KIND=8) :: beta, inv_beta, ek

    INTEGER, INTENT(IN) :: source_id                     !< source term identifier
    INTEGER, INTENT(IN) :: cell_id                       !< cell identifier
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cell_data  !< physical properties 
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons       !< conservative variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s         !< source term
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: js      !< source term Jacobian

    ! Useful quantity 
    ns_p1 = nb_ns + 1

    ! Physical data 
    ! Species densities 
    DO i = 1,nb_ns 
       rhoi(i) = cons(i)
    ENDDO
    
    ! Velocity components 
    DO d = 1,nb_dim 
       vel(d) = cell_data(pos_u_cell + d - 1)
    ENDDO
 
    ! Temperatures 
    DO i = 1,nb_temp
       temp(i) = cell_data(pos_T_cell + i - 1)
    ENDDO
    
    ! Energy components in thermal equilibrium with translation of heavy-particles (at temperature T)
    DO i = 1,nb_ns
       tmp1 = cell_data(pos_ei_cell + i - 1) 
       DO k = 1,nb_int_temp 
          tmp1 = tmp1 - cell_data(pos_eki_cell + (i - 1)*nb_int_temp + k - 1)
       ENDDO
       eiT(i) = tmp1
    ENDDO

    ! Kinetic energy per unit mass and beta factors
    ek   = cell_data(pos_ek_cell)
    beta = cell_data(pos_beta_cell)

    inv_beta = 1.d0/beta
    DO k = 1,nb_int_temp
       ov_betak(k) = 1.d0/cell_data(pos_betak_cell + k - 1)
    ENDDO

    ! Source term due to kinetics and related Jacobian (with respect to primitive variables)
    CALL library_get_source_Jac (rhoi, temp, s, js_line)

    ! Number of partial derivatives
    dim_js = nb_ns + nb_temp

    ! Chain rule application for source term Jacobian with respect to
    ! conservative variables (dS_dU = dS_dP*dP_dU)
    ! Column i,  i = 1,..,nb_ns
    DO j = 1,nb_ns 

       ! Species continuity equations
       tmp1 = inv_beta*(ek - eiT(j))
       DO i = 1,nb_ns
          pos  = (i - 1)*dim_js
          tmp2 = js_line(pos + j) + tmp1*js_line(pos + ns_p1)
          DO k = 1,nb_int_temp 
             tmp2 = tmp2 - js_line(pos + ns_p1 + k)*cell_data(pos_eki_cell + (j - 1)*nb_int_temp + k - 1)*ov_betak(k)
          ENDDO 
          js(i,j) = tmp2
       ENDDO

       ! Global momenum equation
       DO d = 1,nb_dim
          js(pos_rhou + d - 1,j) = 0.d0 
       ENDDO
          
       ! Global energy equation
       js(pos_rhoE,j) = 0.d0

       ! Internal energy equations
       DO k = 1,nb_int_temp
          pos  = (nb_ns + k - 1)*dim_js
          tmp2 = js_line(pos + j) + tmp1*js_line(pos + ns_p1)
          DO l = 1,nb_int_temp 
             tmp2 = tmp2 - js_line(pos + ns_p1 + l)*cell_data(pos_eki_cell + (j - 1)*nb_int_temp + l - 1)*ov_betak(l)
          ENDDO 
          js(pos_rhoek + k - 1,j) = tmp2
       ENDDO

    ENDDO

    ! Columns nb_ns + d,  d = 1,..,nb_dim 
    DO d = 1,nb_dim 

       pos_d = pos_rhou + d - 1

       ! Species continuity equations
       tmp1 = - vel(d)*inv_beta
       DO i = 1,nb_ns 
           js(i,pos_d) = tmp1*js_line((i - 1)*dim_js + ns_p1)
       ENDDO

       ! Global momenum equation
       DO p = 1,nb_dim 
          js(pos_rhou + p - 1,pos_d) = 0.d0 
       ENDDO

       ! Global energy equation
       js(pos_rhoE,pos_d) = 0.d0

       ! Internal energy equations
       DO k = 1,nb_int_temp 
          js(pos_rhoek + k - 1,pos_d) = tmp1*js_line((nb_ns + k - 1)*dim_js + ns_p1)
       ENDDO

    ENDDO

    ! Column nb_ns + nb_dim + 1
    ! Species continuity equations 
    DO i = 1,nb_ns 
       js(i,pos_rhoE) = inv_beta*js_line((i - 1)*dim_js + ns_p1) 
    ENDDO

    ! Global momenum equation
    DO d = 1,nb_dim 
       js(pos_rhou + d - 1,pos_rhoE) = 0.d0 
    ENDDO

    ! Global energy equation
    js(pos_rhoE,pos_rhoE) = 0.d0

    ! Internal energy equations 
    DO k = 1,nb_int_temp
       js(pos_rhoek + k - 1,pos_rhoE) = inv_beta*js_line((nb_ns + k - 1)*dim_js + ns_p1)
    ENDDO

    ! Columns nb_ns + nb_dim + 1 + k,  k = 1,..,nb_int_temp
    DO k = 1,nb_int_temp 

       pos_k = pos_rhoek + k - 1

       ! Species continuity equations
       DO i = 1,nb_ns 
          pos = (i - 1)*dim_js
          js(i,pos_k) = - inv_beta*js_line(pos + ns_p1) + js_line(pos + ns_p1 + k)*ov_betak(k)
       ENDDO

       ! Global momenum equation
       DO d = 1,nb_dim 
          js(pos_rhou + d - 1,pos_k) = 0.d0
       ENDDO
 
       ! Global energy equation
       js(pos_rhoE,pos_k) = 0.d0

       ! Internal energy equations 
       DO l = 1,nb_int_temp
          pos = (nb_ns + l - 1)*dim_js 
          js(pos_rhoek + l - 1,pos_k) = - inv_beta*js_line(pos + ns_p1) + js_line(pos + ns_p1 + k)*ov_betak(k)
       ENDDO

    ENDDO

  END SUBROUTINE source_term_neq_kinetics_Jac
!------------------------------------------------------------------------------!
!> This subroutine provides the source term and the related Jacobian due to collisional and radiative
!! processes for nonequilibrium flows.
!! Thermal nonequilibrium between heavy particles and free electrons is accounted for by means 
!! of a separate free electron pseudo-entropy equation
  SUBROUTINE source_term_neq_el_kinetics_Jac (source_id, cell_id, cell_data, cons, s, js)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: source_id                     !< source term identifier
    INTEGER, INTENT(IN) :: cell_id                       !< cell identifier
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cell_data  !< physical properties 
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: cons       !< conservative variables
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: s         !< source term
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: js      !< source term Jacobian

    s = 0.d0
    js = 0.d0

    PRINT*,'to be implemented!'
    STOP

  END SUBROUTINE source_term_neq_el_kinetics_Jac
!------------------------------------------------------------------------------!
