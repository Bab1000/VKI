!------------------------------------------------------------------------------!
!> This subroutine applies a supersonic outlet boundary condition for 1D nonequilibrium flows.
  SUBROUTINE sup_out_neq_1DImpl (id, phys_data, u_phys, ghost_data, u_ghost, ju_ghost)

    USE mod_general_data,          ONLY: nb_eq 

    IMPLICIT NONE

    INTEGER :: i, j

    INTEGER, INTENT(IN) :: id                              !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys       !< conservative variables of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data    !< physical properties of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost     !< conservative variables of ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data  !< physical properties of ghost cell
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: ju_ghost  !< ghost state conservative variable Jacobian

    ! Ghost state conservative variable and physical property vectors
    u_ghost    = u_phys
    ghost_data = phys_data 

    ! Evaluation of dUg/dUp matrix for implicit boundary condition  
    DO j = 1,nb_eq

       DO i = 1,j - 1 
          ju_ghost(i,j) = 0.d0
       ENDDO

       ju_ghost(j,j) = 1.d0

       DO i = j + 1,nb_eq 
          ju_ghost(i,j) = 0.d0
       ENDDO

    ENDDO

  END SUBROUTINE sup_out_neq_1DImpl
!------------------------------------------------------------------------------!
!> This subroutine applies a supersonic inlet boundary conditions for 1D nonequilibrium flows.
  SUBROUTINE sup_in_neq_1DImpl (id, phys_data, u_phys, ghost_data, u_ghost, ju_ghost)

    USE mod_general_data,             ONLY: nb_eq, nb_ns, nb_temp, pos_u_cell, pos_T_cell, pos_u 
    USE mod_domain_boundary,          ONLY: boundary_data, get_boundary_inlet, boundary 
    USE mod_function_pointer,         ONLY: get_cons_phys_from_prim

    IMPLICIT NONE

    INTEGER :: i, j
    REAL(KIND=8), DIMENSION(nb_eq) :: inlet_data, prim
    TYPE(boundary) :: bound

    INTEGER, INTENT(IN) :: id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys, phys_data
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost, ghost_data
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: ju_ghost

    ! Inlet data (species densities, velocity and temperature(s))
    bound      = boundary_data(id)
    inlet_data = get_boundary_inlet (nb_eq, bound)

    ! Ghost state primitive variables (linear extrapolation used for all variables)
    ! Species densities
    DO i = 1,nb_ns 
       prim(i) = 2.d0*inlet_data(i) - u_phys(i)
    ENDDO

    ! Velocity
    prim(pos_u) = 2.d0*inlet_data(pos_u) - phys_data(pos_u_cell) 

    ! Temperature(s)
    DO i = 1,nb_temp
       prim(pos_u + i) = 2.d0*inlet_data(pos_u + i) - phys_data(pos_T_cell + i - 1)
    ENDDO

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost, ghost_data)

    !! APPROXIMATION !!
    DO j = 1,nb_eq 

       DO i = 1,j - 1
          ju_ghost(i,j) = 0.d0 
       ENDDO

       ju_ghost(i,i) = - 1.d0

       DO i = j + 1,nb_eq 
          ju_ghost(i,j) = 0.d0
       ENDDO

    ENDDO

    PRINT*
    WRITE(*,'(A)')'In "sup_in_neq_1DImpl", boundary conditions Jacobian to be implemented!'
    PRINT*
    STOP
 
  END SUBROUTINE sup_in_neq_1DImpl
!------------------------------------------------------------------------------!
!> This subroutine applies a subsonic inlet boundary condition for 1D nonequilibrium flows. 
  SUBROUTINE sub_in_neq_rhoiin_T_1DImpl (id, phys_data, u_phys, ghost_data, u_ghost, ju_ghost)
 
    USE mod_general_data,             ONLY: nb_ns, nb_eq, nb_temp, nb_int_temp, nb_te, pos_h0_cell, pos_u_cell,  & 
                                          & pos_ek_cell, pos_T_cell, pos_pres_cell, pos_rho_cell,                & 
                                          & pos_alpha_cell, pos_beta_cell, pos_betak_cell, pos_ei_cell,          & 
                                          & pos_eki_cell, pos_em, pos_u, pos_T, pos_rhou, pos_rhoE, pos_rhoek,   & 
                                          & pos_rhose, Re, gamma_e, gamma_e_m1, ov_gamma_e_m1, beta, ei, eiT
    USE mod_domain_boundary,          ONLY: boundary, boundary_data, get_boundary_rhoi, get_boundary_Tvec
    USE mod_function_pointer,         ONLY: get_cons_phys_from_prim

    IMPLICIT NONE

    INTEGER :: i, j, k, kp
    INTEGER :: pos, pos_k
    REAL(KIND=8) :: tmp1
    REAL(KIND=8) :: alpha, betaG_ov_betaP, ek, rho, u
    REAL(KIND=8) :: rhoP, rhoG, rho_sum, rho_sum_ov_rhoP_u, rho_sum_ov_rhoP_u2, ov_rhoP 
    REAL(KIND=8) :: rho_emP, rho_emG, rho_emP_ov_rho_emG, rho_emG_ov_rho_emP, rhoG_pow_ge_m1, peP, TeP, TeG
    REAL(KIND=8) :: fac1_em, fac2_em, fac3_em
    REAL(KIND=8), DIMENSION(nb_ns) :: rhoi_in
    REAL(KIND=8), DIMENSION(nb_temp) :: temp_in
    REAL(KIND=8), DIMENSION(MAX(1,nb_int_temp)) :: betakG_ov_betakP
    REAL(KIND=8), DIMENSION(nb_eq) :: prim
    TYPE(boundary) :: bound

    INTEGER, INTENT(IN) :: id                              !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys       !< conservative variables of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data    !< physical properties of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost     !< conservative variables of ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data  !< physical properties of ghost cell
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: ju_ghost  !< ghost state conservative variable Jacobian
   
    ! Species densities at inlet
    bound   = boundary_data(id)
    rhoi_in = get_boundary_rhoi (nb_ns, bound)

    ! Temperature vector at inlet
    temp_in = get_boundary_Tvec (nb_temp, bound)

    ! Ghost state primitive variables (linear extrapolation used for species densities and temperature(s))
    ! Species densities
    DO i = 1,nb_ns
       prim(i) = 2.d0*rhoi_in(i) - u_phys(i) 
    ENDDO

    ! Velocity 
    prim(pos_u) = phys_data(pos_u_cell)

    ! Temperature(s)
    DO i = 1,nb_temp
       prim(pos_T + i - 1) = 2.d0*temp_in(i) - phys_data(pos_T_cell + i - 1)
    ENDDO

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost, ghost_data)

    ! Evaluation of dUg/dUp matrix for implicit boundary condition  
    ! Common factors
    u  = phys_data(pos_u_cell)  
    ek = phys_data(pos_ek_cell) 
    rho  = ghost_data(pos_rho_cell) 
    rhoP = phys_data(pos_rho_cell) 
    ov_rhoP   = 1.d0/rhoP
    rho_sum   = rho + rhoP
    rho_sum_ov_rhoP_u  = u*rho_sum/rhoP
    rho_sum_ov_rhoP_u2 = u*rho_sum_ov_rhoP_u
   
    betaG_ov_betaP = ghost_data(pos_beta_cell)/phys_data(pos_beta_cell)
    DO k = 1,nb_int_temp
       betakG_ov_betakP(k) = ghost_data(pos_betak_cell + k - 1)/phys_data(pos_betak_cell + k - 1) 
    ENDDO
   
    DO i = 1,nb_ns + nb_int_temp*nb_ns 
       ei(i) = ghost_data(pos_ei_cell + i - 1)
    ENDDO

    ! Matrix dUg/dUp entries
    IF (nb_te.EQ.0) THEN

       ! Components of species energy depending on translational temperature of heavy-particles
       DO i = 1,nb_ns 
          tmp1  = phys_data(pos_ei_cell + i - 1)
          DO k = 1,nb_int_temp
             tmp1  = tmp1 - phys_data(pos_eki_cell + (i - 1)*nb_int_temp + k - 1)
          ENDDO
          eiT(i) = tmp1
       ENDDO

       ! Column i,  i = 1,..,nb_ns
       DO j = 1,nb_ns 
 
          DO i = 1,j - 1
             ju_ghost(i,j) = 0.d0
          ENDDO

          ju_ghost(j,j) = - 1.d0

          DO i = j + 1,nb_ns 
             ju_ghost(i,j) = 0.d0
          ENDDO
 
          ju_ghost(pos_rhou,j) = - rho_sum_ov_rhoP_u

          tmp1 = - ei(j) + betaG_ov_betaP*(eiT(j) - ek) + ek - rho_sum_ov_rhoP_u2 
          DO k = 1,nb_int_temp
             tmp1 = tmp1 + phys_data(pos_eki_cell + (k - 1)*nb_int_temp + j - 1)*betakG_ov_betakP(k)
          ENDDO
          ju_ghost(pos_rhoE,j) = tmp1

          DO k = 1,nb_int_temp
             ju_ghost(pos_rhoE + k,j) = phys_data(pos_eki_cell + (k - 1)*nb_int_temp + j - 1)*betakG_ov_betakP(k)  &
                                   &  - ei(nb_ns + (k - 1)*nb_int_temp + j) 
          ENDDO

       ENDDO

       ! Column nb_ns + 1
       DO i = 1,nb_ns 
          ju_ghost(i,pos_rhou) = 0.d0
       ENDDO
  
       ju_ghost(pos_rhou,pos_rhou) = rho*ov_rhoP
       ju_ghost(pos_rhoE,pos_rhou) = betaG_ov_betaP*u + u_ghost(pos_rhou)*ov_rhoP
       
       DO k = 1,nb_int_temp
          ju_ghost(pos_rhoE + k,pos_rhou) = 0.d0 
       ENDDO

       ! Column nb_ns + 2
       DO i = 1,nb_ns 
          ju_ghost(i,pos_rhoE) = 0.d0
       ENDDO
  
       ju_ghost(pos_rhou,pos_rhoE) = 0.d0
       ju_ghost(pos_rhoE,pos_rhoE) = - betaG_ov_betaP

       DO k = 1,nb_int_temp
          ju_ghost(pos_rhoE + k,pos_rhoE) = 0.d0 
       ENDDO

       ! Column nb_ns + 2 + k,  k = 1,..,nb_int_temp
       DO k = 1,nb_int_temp

          tmp1 = betakG_ov_betakP(k) 

          pos_k = pos_rhoek + k - 1

          DO i = 1,nb_ns 
             ju_ghost(i,pos_k) = 0.d0
          ENDDO 

          ju_ghost(pos_rhou,pos_k) = 0.d0
          ju_ghost(pos_rhoE,pos_k) = betaG_ov_betaP - tmp1  

          DO kp = 1,k - 1
             ju_ghost(pos_rhoek + kp - 1,pos_k) = 0.d0
          ENDDO     

          ju_ghost(pos_k,pos_k) = - tmp1

          DO i = k + 1,nb_int_temp
             ju_ghost(pos_rhoek + kp - 1,pos_k) = 0.d0
          ENDDO  

       ENDDO

    ELSE 

       ! Components of species energy depending on translational temperature of heavy-particles
       eiT(pos_em) = 0.d0
       DO i = pos_em + 1,nb_ns 
          tmp1  = phys_data(pos_ei_cell + i - 1)
          DO k = 1,nb_int_temp
             tmp1 = tmp1 - phys_data(pos_eki_cell + (i - 1)*nb_int_temp + k - 1)
          ENDDO
          eiT(i) = tmp1
       ENDDO

       ! Additional common factors 
       rho_emP = u_phys(pos_em)
       rho_emG = u_ghost(pos_em)
       TeP     = phys_data(nb_temp)
       TeG     = ghost_data(nb_temp)
       peP     = rho_emP*Re*TeP

       rho_emP_ov_rho_emG = rho_emP/rho_emG
       rho_emG_ov_rho_emP = rho_emG/rho_emP
       rhoG_pow_ge_m1     = rho**gamma_e_m1
       
       fac1_em = betaG_ov_betaP*peP/rho - rho_emG*ov_rhoP*Re*TeP
       fac2_em = Re/(rho**gamma_e)*gamma_e_m1*rho_emG*(TeG - rho*ov_rhoP*TeP)
       fac3_em = Re*(TeG - TeP*rho_emG_ov_rho_emP)

       ! Column 1 (free electrons)
       j = pos_em 
       DO i = 1,j - 1
          ju_ghost(i,j) = 0.d0
       ENDDO

       ju_ghost(j,j) = - 1.d0

       DO i = j + 1,nb_ns 
          ju_ghost(i,j) = 0.d0
       ENDDO
 
       ju_ghost(pos_rhou,j) = - rho_sum_ov_rhoP_u

       ju_ghost(pos_rhoE,j) = (1.d0 - betaG_ov_betaP)*ek - rho_sum_ov_rhoP_u2 + fac1_em - fac3_em*ov_gamma_e_m1 
       
       DO k = 1,nb_int_temp
          ju_ghost(pos_rhoek + k - 1,j) = 0.d0
       ENDDO

       ju_ghost(pos_rhose,j) = fac2_em - fac3_em/rhoG_pow_ge_m1

       ! Column i,  i = 2,..,nb_ns (heavy particles)
       DO j = pos_em + 1,nb_ns 
 
          DO i = 1,j - 1
             ju_ghost(i,j) = 0.d0
          ENDDO

          ju_ghost(j,j) = - 1.d0

          DO i = j + 1,nb_ns 
             ju_ghost(i,j) = 0.d0
          ENDDO
 
          ju_ghost(pos_rhou,j) = - rho_sum_ov_rhoP_u

          tmp1 = - ei(j) + betaG_ov_betaP*(eiT(j) - ek) + ek - rho_sum_ov_rhoP_u2
          DO k = 1,nb_int_temp
             tmp1 = tmp1 + phys_data(pos_eki_cell + (k - 1)*nb_int_temp + j - 1)*betakG_ov_betakP(k)
          ENDDO
          ju_ghost(pos_rhoE,j) = tmp1

          DO k = 1,nb_int_temp
             ju_ghost(pos_rhoE + k,j) = phys_data(pos_eki_cell + (k - 1)*nb_int_temp + j - 1)*betakG_ov_betakP(k)  &
                                   &  - ei(nb_ns + (k - 1)*nb_int_temp + j) 
          ENDDO

          ju_ghost(pos_rhose,j) = fac2_em

       ENDDO

       ! Column nb_ns + 1
       DO i = 1,nb_ns 
          ju_ghost(i,pos_rhou) = 0.d0
       ENDDO
  
       ju_ghost(pos_rhou,pos_rhou) = rho*ov_rhoP
       ju_ghost(pos_rhoE,pos_rhou) = betaG_ov_betaP*u + u_ghost(pos_rhou)*ov_rhoP

       DO k = 1,nb_int_temp
          ju_ghost(pos_rhoE + k,pos_rhou) = 0.d0 
       ENDDO

       ju_ghost(pos_rhose,pos_rhou) = 0.d0

       ! Column nb_ns + 2
       DO i = 1,nb_ns 
          ju_ghost(i,pos_rhoE) = 0.d0
       ENDDO
  
       ju_ghost(pos_rhou,pos_rhoe) = 0.d0
       ju_ghost(pos_rhoE,pos_rhoE) = - betaG_ov_betaP

       DO k = 1,nb_int_temp
          ju_ghost(pos_rhoE + k,pos_rhoE) = 0.d0 
       ENDDO
 
       ju_ghost(pos_rhose,pos_rhoE) = 0.d0

       ! Column nb_ns + 2 + k,  k = 1,..,nb_int_temp
       DO k = 1,nb_int_temp

          tmp1 = betakG_ov_betakP(k) 

          pos_k = pos_rhoek + k - 1

          DO i = 1,nb_ns 
             ju_ghost(i,pos_k) = 0.d0
          ENDDO 

          ju_ghost(pos_rhou,pos_k) = 0.d0
          ju_ghost(pos_rhoE,pos_k) = betaG_ov_betaP - tmp1  

          DO kp = 1,k - 1
             ju_ghost(pos_rhoek + kp - 1,pos_k) = 0.d0
          ENDDO     

          ju_ghost(pos_k,pos_k) = - tmp1

          DO i = k + 1,nb_int_temp
             ju_ghost(pos_rhoek + kp - 1,pos_k) = 0.d0
          ENDDO  

          ju_ghost(pos_rhose,pos_k) = 0.d0

       ENDDO

       ! Column nb_ns + 2 + nb_int_temp + nb_te  
       DO i = 1,nb_ns 
          ju_ghost(i,pos_rhose) = 0.d0 
       ENDDO

       ju_ghost(pos_rhou,pos_rhose) = 0.d0  
       ju_ghost(pos_rhoE,pos_rhose) = rhoG_pow_ge_m1*ov_gamma_e_m1*(betaG_ov_betaP - rho_emP_ov_rho_emG)

       DO k = 1,nb_int_temp
          ju_ghost(pos_rhoek + k - 1,pos_rhose) = 0.d0
       ENDDO

       ju_ghost(pos_rhose,pos_rhose) = - rho_emG_ov_rho_emP*(rhoP**gamma_e_m1)/rhoG_pow_ge_m1

    ENDIF 
    
  END SUBROUTINE sub_in_neq_rhoiin_T_1DImpl
!------------------------------------------------------------------------------!
!> This subroutine applies a subsonic outlet boundary condition for 1D nonequilibrium flows.
!! The static temperature is imposed.
  SUBROUTINE sub_out_neq_Tout_1DImpl (id, phys_data, u_phys, ghost_data, u_ghost, ju_ghost)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: id                              !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys       !< conservative variables of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data    !< physical properties of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost     !< conservative variables of ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data  !< physical properties of ghost cell
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: ju_ghost  !< ghost state conservative variable Jacobian

    u_ghost    = 0.d0
    ghost_data = 0.d0
    ju_ghost   = 0.d0

    PRINT*,'in "sup_in_neq_1DImpl"'
    PRINT*,'to be implemented yet'
    PAUSE

  END SUBROUTINE sub_out_neq_Tout_1DImpl
!------------------------------------------------------------------------------!
