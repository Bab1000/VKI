!------------------------------------------------------------------------------!
!> This subroutine applies a supersonic inlet boundary conditions for 1D stagnation line nonequilibrium flows.
  SUBROUTINE sup_in_neq_1D_SL_Expl_2nd(id, phys_data1, u_phys1, phys_data2, u_phys2, & 
                                    &  ghost_data1, u_ghost1, ghost_data2, u_ghost2)


    USE mod_general_data,          ONLY: nb_ns, nb_temp, nb_eq, pos_u, pos_v, pos_T, pos_u_cell, & 
                                       & pos_v_cell, pos_T_cell                                         
    USE mod_domain_boundary,       ONLY: boundary, boundary_data, get_boundary_inlet 
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim
 
    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8) :: ov_3
    REAL(KIND=8) :: rho_in, u_in, v_in, p_in
    REAL(KIND=8), DIMENSION(nb_eq) :: inlet_data, prim
    TYPE(boundary) :: bound

    INTEGER, INTENT(IN) :: id                               !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys1       !< conservative variables of the first physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys2       !< conservative variables of the second physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data1    !< physical properties of the first physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data2    !< physical properties of the second physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost1     !< conservative variables of the first ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost2     !< conservative variables of the second ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data1  !< physical properties of the first ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data2  !< physical properties of the second ghost cell    

    ! Inlet values of density, velocity components and pressure
    bound      = boundary_data(id)
    inlet_data = get_boundary_inlet(nb_eq, bound)

    ! Common factor
    ov_3 = 1.d0/3.d0

    ! Ghost state 1
    ! Primitive variables
    ! Species densities
    DO i = 1,nb_ns 
       prim(i) = (4.d0*inlet_data(i) - u_phys2(i))*ov_3
       IF(prim(i)<0.d0) prim(i) = 0.d0
    ENDDO

    ! Velocity components 
    prim(pos_u) = (4.d0*inlet_data(pos_u) - phys_data2(pos_u_cell))*ov_3
    prim(pos_v) = (4.d0*inlet_data(pos_v) - phys_data2(pos_v_cell))*ov_3

    ! Temperature(s)
    DO i = 1,nb_temp
       prim(pos_T + i - 1) = (4.d0*inlet_data(pos_T + i - 1) - phys_data2(pos_T_cell + i - 1))*ov_3
       IF(prim(pos_T + i - 1) < 50.d0) prim(pos_T + i - 1) = 50.d0
    ENDDO

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost1, ghost_data1)

    ! Ghost state 2 
    ! Primitive variables 
    ! Species densities
    DO i = 1,nb_ns 
       prim(i) = 2.d0*inlet_data(i) - u_phys2(i)
       IF(prim(i)<0.d0) prim(i) = 0.d0
    ENDDO

    ! Velocity components 
    prim(pos_u) = 2.d0*inlet_data(pos_u) - phys_data2(pos_u_cell)
    prim(pos_v) = 2.d0*inlet_data(pos_v) - phys_data2(pos_v_cell)

    ! Temperature(s)
    DO i = 1,nb_temp
       prim(pos_T + i - 1) = 2.d0*inlet_data(pos_T + i - 1) - phys_data2(pos_T_cell + i - 1)
       IF(prim(pos_T + i - 1) < 50.d0) prim(pos_T + i - 1) = 50.d0
    ENDDO

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost2, ghost_data2)

  END SUBROUTINE sup_in_neq_1D_SL_Expl_2nd
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
! This subroutine applies a subsonic inlet boundary condition for 1D stagnation line 
! nonequilibrium flows by imposing the radial velocity and the tangential-velocity derivative at the
! inlet (check with Klomfass paper to understand their "physical meaning),
! the inlet static temperature and mixture mass fractions the inlet.
!------------------------------------------------------------------------------!
!> @author
!> Alessandro Turchi
!>
!>@brief
!>Applies subsonic inlet boundary condition with given inlet radial velocity and
! tangential-velocity derivative (check with Klomfass paper to understand their
! "physical meaning), inlet temperature and assuming equilibrium inlet composition.
!> 
!>\b DESCRIPTION: \n 
!> This subroutine applies a subsonic inlet boundary condition for 1D stagnation line 
!> non-equilibrium flows by imposing the radial velocity and tangential-velocityb derivative (check
!> with Klomfass paper to understand their  "physical meaning), the inlet static temperature 
!> and the inlet mass fractions.
!------------------------------------------------------------------------------!
  SUBROUTINE sub_in_vel_neq_1D_SL_Expl_2nd(id, phys_data1, u_phys1, phys_data2, u_phys2, & 
                                    &  ghost_data1, u_ghost1, ghost_data2, u_ghost2)


    USE mod_general_data,          ONLY: nb_ns, nb_temp, nb_eq, pos_u, pos_v, pos_T, pos_u_cell, & 
                                       & pos_v_cell, pos_T_cell, pos_pres_cell, Ri, p_inf, xc, nb_cells, volumes
    USE mod_domain_boundary,       ONLY: boundary, boundary_data, get_boundary_uin, get_boundary_Tin, &
                                       & get_boundary_dv_dyin, get_boundary_yi 
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim
    USE mod_neq_function_pointer,  ONLY: library_compute_eq_composition
 
    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8) :: rho_in, u_in, v_in, p_in
    REAL(KIND=8) :: mdot_in, T_in
    REAL(KIND=8) :: R_mix
    REAL(KIND=8) :: ratio 
    REAL(KIND=8) :: rho_r, rho_l, u_l, v_l, T_l
    REAL(KIND=8) :: frac, frac_inv


    REAL(KIND=8) :: r                                      !< Distance inlet boundary-center of the spherical body
    REAL(KIND=8) :: dv_dy                                  !<Derivative along the direction normal to the stagnation line of 
                                                           ! the v "physical" velocity component (normal to the stagnation line) 
    REAL(KIND=8) :: rho_phys_1, rho_phys_2
    REAL(KIND=8), DIMENSION(nb_eq)   :: prim
    REAL(KIND=8), DIMENSION(nb_ns)   :: xi_in
    REAL(KIND=8), DIMENSION(nb_ns)   :: rhoi_in
    REAL(KIND=8), DIMENSION(nb_ns)   :: yi_in
    REAL(KIND=8), DIMENSION(nb_ns)   :: yi_phys_1, yi_phys_2, yi_l, yi_r, rhoi_l

    TYPE(boundary) :: bound

    INTEGER, INTENT(IN) :: id                               !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys1       !< conservative variables of the first physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys2       !< conservative variables of the second physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data1    !< physical properties of the first physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data2    !< physical properties of the second physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost1     !< conservative variables of the first ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost2     !< conservative variables of the second ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data1  !< physical properties of the first ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data2  !< physical properties of the second ghost cell    

    ! Inlet values of density, velocity components and pressure
    bound      = boundary_data(id)


    ! Inlet velocity, velocity gradient, temperature and composition
    u_in = get_boundary_uin(bound) 
    dv_dy   = get_boundary_dv_dyin(bound)
    T_in    = get_boundary_Tin(bound)
    yi_in = get_boundary_yi(nb_ns,bound)

    ! Inlet pressure is extrapolated from the first two physical states
    ! accounting for the possible cell non-uniformity.
    ratio = volumes(nb_cells+2) / (volumes(nb_cells+2) + volumes(nb_cells+1)) 
    p_in    = phys_data1(pos_pres_cell) - (phys_data2(pos_pres_cell) - phys_data1(pos_pres_cell)) * ratio

    ! Evaluates the specific gas constant of the mixture
    DO i = 1,nb_ns 
       R_mix = R_mix +Ri(i)*yi_in(i) 
    ENDDO

    ! Evaluates the density of the mixtures
    rho_in=p_in/(R_mix*T_in)

    ! Inlet tangential velocity
    r = xc(nb_cells+2)+ volumes(nb_cells+2) * 0.5d0 
    v_in    = dv_dy * r - u_in 
    

    ! Inlet species densities
    DO  i = 1,nb_ns
       rhoi_in(i) = rho_in * yi_in(i)
    ENDDO

    !!!!!!!!!!!!!!!!!
    ! Ghost state 1 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! It is calculated using the definition of the Roe averaged quantities.   !
    ! For this reason the present reconstruction of the ghost cell is really  !
    ! accurate only when the Roe flux splitting with linear reconstruction of !
    ! the fluxes (II order) is used                                           !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Left state (in) of the physical-ghost interface (inlet boundary).
    ! It  is computed using a linear interpolation of the first two physical states. Note that the
    ! only difference wrt the recostruction applied in the flowfield is that no
    ! limiter function is applied here (negligible error considering the small
    ! gradients at the inflow).
    rho_l = 0.d0
    DO i = 1, nb_ns
       rhoi_l(i) = u_phys1(i)-0.5d0*(u_phys2(i)-u_phys1(i))
       rho_l= rho_l + rhoi_l(i)
    ENDDO

    ! Right state (out) of the physical-ghost interface (inlet boundary).
    ! It is calculated usign the definition of the Roe average density
    rho_r =  (rho_in * rho_in)/ rho_l

    ! Calculation of some useful quantities
    frac     = 1.d0/(1.d0+sqrt(rho_r/rho_l))
    frac_inv = 1.d0/frac

    ! Primitive variables

    ! Species densities
    DO i = 1,nb_ns 

       ! Constant extrapolation
!       prim(i) =  rhoi_in(i) 

       ! Linear extrapolation
!       prim(i) = 2.d0* rhoi_in(i) - u_phys1(i)

       ! Roe average state extrapolation
       prim(i)=(yi_in(i) - rhoi_l(i)  * frac / rho_l) * rho_r * (frac_inv / sqrt(rho_r/rho_l))
       yi_r(i)=prim(i)/rho_r

       !Check to avoid negative species densities
       if  (prim(i).lt.1.d-30) prim(i) = 1.d-30

    ENDDO

    ! Velocity components 

    ! Constant extrapolation
!    prim(pos_u) = u_in 
!    prim(pos_v) = v_in 

    ! Linear extrapolation
!    prim(pos_u) = 2.d0*u_in - phys_data1(pos_u_cell)
!    prim(pos_v) = 2.d0*v_in - phys_data1(pos_v_cell)
  
    ! Roe average state extrapolation
    u_l      =  phys_data1(pos_u_cell) - 0.5d0*(phys_data2(pos_u_cell) - phys_data1(pos_u_cell)) 
    v_l      =  phys_data1(pos_v_cell) - 0.5d0*(phys_data2(pos_v_cell) - phys_data1(pos_v_cell)) 
    prim(pos_u) = (u_in - frac * u_l) *  (frac_inv / sqrt(rho_r/rho_l))
    prim(pos_v) = (v_in - frac * v_l) *  (frac_inv / sqrt(rho_r/rho_l))

    ! Temperature(s)
    DO i = 1,nb_temp

       ! Constant extrapolation
!       prim(pos_T + i - 1) = T_in 

       ! Linear extrapolation
!       prim(pos_T + i - 1) = MAX(T_in - (phys_data1(pos_T_cell + i - 1) - T_in), 0.5d0 * T_in) 

       ! Roe average state extrapolation
       T_l = phys_data1(pos_T_cell + i - 1) - 0.5d0*(phys_data2(pos_T_cell + i - 1) - phys_data1(pos_T_cell + i - 1)) 
       prim(pos_T + i - 1) = MIN((T_in - frac * T_l) *  (frac_inv / sqrt(rho_r/rho_l)), 2.d0*T_in)

    ENDDO

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost1, ghost_data1)

    !!!!!!!!!!!!!!!!!
    ! Ghost state 2 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Ghost state 2 primitive variables extrapolation                                      !
    ! Comment all the box below to obtain constant extrapolation form the 1st ghost state  !
    ! In the present form of the boundary condition, there is no need of calculating       !
    ! accurately the second ghost cell since the first ghost cell has been                 !
    ! reconstructed according to the inlet (phisical-ghost interface) values.              !
!----------------------------------------------------------------------------------------------------------------------
    ! Species densities                                                                                               !
!    DO i = 1,nb_ns                                                                                                    !
                                                                                                                      !
       ! Linear extrapolation using 2nd physical state                                                                !
!       prim(i) = 2.d0*rhoi_in(i) - u_phys2(i)                                                                         !
                                                                                                                      !
       ! Linear extrapolation using 1st ghost state                                                                   !
!       prim(i) = 3.d0 * prim(i) - 2.d0 * rhoi_in(i)                                                                  !
                                                                                                                      !
!    ENDDO                                                                                                             !
                                                                                                                      !
                                                                                                                      !
    ! Velocity components                                                                                             !
                                                                                                                      !
    ! Linear extrapolation using 2nd physical state                                                                   !
!    prim(pos_u) = 2.d0*u_in - phys_data2(pos_u_cell)                                                                  !
!    prim(pos_v) = 2.d0*v_in - phys_data2(pos_v_cell)                                                                  !
                                                                                                                      !
    ! Linear extrapolation using 1st ghost state                                                                      !
!    prim(pos_u) = 3.d0 * prim(pos_u) - 2.d0 * u_in                                                                   !
!    prim(pos_v) = 3.d0 * prim(pos_v) - 2.d0 * v_in                                                                   !
                                                                                                                      !
    ! Temperature(s)                                                                                                  !
!    DO i = 1,nb_temp                                                                                                  !
                                                                                                                      !
       ! Linear extrapolation using 2nd physical state                                                                !
!       prim(pos_T + i - 1) = 2.d0*T_in - phys_data2(pos_T_cell + i - 1)                                               !
                                                                                                                      !
       ! Linear extrapolation using 1st ghost state                                                                   !
!       prim(pos_T + i - 1) = MAX(3.d0 * prim(pos_T + i - 1) - 2.d0 * T_in, 0.5d0 * prim(pos_T + i - 1))              !
                                                                                                                      !
!    ENDDO                                                                                                             !
!-------------------------------------------------------------------------------------------------------------------- !
                                                                                                                      
    ! Compute the ghost state conservative variables and physical properties from primitive variables                 
    CALL get_cons_phys_from_prim(prim, u_ghost2, ghost_data2)


  END SUBROUTINE sub_in_vel_neq_1D_SL_Expl_2nd

!------------------------------------------------------------------------------!
!> This subroutine applies a slip wall boundary conditions for 1D stagnation line nonequilibrium flows.
   SUBROUTINE slip_wall_neq_1D_SL_Expl_2nd(id, phys_data1, u_phys1, phys_data2, u_phys2, & 
                                        &  ghost_data1, u_ghost1, ghost_data2, u_ghost2) 

     USE mod_general_data,          ONLY: nb_ns, nb_temp, nb_eq, pos_u, pos_v, pos_T, pos_u_cell, & 
                                        & pos_v_cell, pos_T_cell 
     USE mod_domain_boundary,       ONLY: boundary_data, get_boundary_rhoin, get_boundary_pin
     USE mod_function_pointer,      ONLY: get_cons_phys_from_prim 

     IMPLICIT NONE

     INTEGER :: i
     REAL(KIND=8) :: u, v
     REAL(KIND=8), DIMENSION(nb_temp) :: temp
     REAL(KIND=8), DIMENSION(nb_eq) :: prim

     INTEGER, INTENT(IN) :: id                               !< boundary id
     REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys1       !< conservative variables of the first physical cell
     REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys2       !< conservative variables of the second physical cell
     REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data1    !< physical properties of the first physical cell
     REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data2    !< physical properties of the second physical cell
     REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost1     !< conservative variables of the first ghost cell
     REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost2     !< conservative variables of the second ghost cell
     REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data1  !< physical properties of the first ghost cell
     REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data2  !< physical properties of the second ghost cell

     ! Velocity components 
     u = phys_data1(pos_u_cell)
     v = phys_data1(pos_v_cell)

     ! Ghost state 1
     ! Primitive variables (pressure gradient is set to zero)
     ! Species densities
     DO i = 1,nb_ns 
        prim(i) = u_phys1(i)
     ENDDO

     ! Velocity components 
     prim(pos_u) = -u 
     prim(pos_v) = v 

     ! Temperature(s)
     DO i = 1,nb_temp
        prim(pos_T + i - 1) = phys_data1(pos_T_cell + i - 1)
     ENDDO

     ! Compute the ghost state conservative variables and physical properties from primitive variables 
     CALL get_cons_phys_from_prim(prim, u_ghost1, ghost_data1)

     ! Ghost state 2
     u_ghost2    = u_ghost1
     ghost_data2 = ghost_data1

   END SUBROUTINE slip_wall_neq_1D_SL_Expl_2nd
!------------------------------------------------------------------------------!
!> This subroutine applies an isothermal no slip wall boundary condition for 1D stagnation line 
!! nonequilibrium flows. 
  SUBROUTINE no_slip_iso_Twall_neq_1D_SL_Expl_2nd (id, phys_data1, u_phys1, phys_data2, u_phys2,   & 
                                                &  ghost_data1, u_ghost1, ghost_data2, u_ghost2)

    USE mod_general_data,          ONLY: nb_ns, nb_temp, nb_eq, pos_u, pos_v, pos_T, pos_Te, pos_u_cell, & 
                                      & pos_v_cell, pos_T_cell, pos_pres_cell, Ri, nb_te, pos_em
    USE mod_domain_boundary,       ONLY: boundary, boundary_data, get_boundary_Twall
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim
    USE mod_neq_function_pointer,  ONLY: library_get_molar_fractions

    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8) :: u1, v1, p_ov_T, T1, T2, Twall, Tghost
    REAL(KIND=8) :: u2, v2
    REAL(KIND=8), DIMENSION(nb_ns) :: xi
    REAL(KIND=8), DIMENSION(nb_eq) :: prim
    TYPE(boundary) :: bound

    INTEGER, INTENT(IN) :: id                               !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys1       !< conservative variables of the first physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys2       !< conservative variables of the second physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data1    !< physical properties of the first physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data2    !< physical properties of the second physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost1     !< conservative variables of the first ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost2     !< conservative variables of the second ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data1  !< physical properties of the first ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data2  !< physical properties of the second ghost cell

    ! Wall temperature
    bound = boundary_data(id)
    Twall = get_boundary_Twall(bound)

    ! Velocity components of physical state 1
    u1 = phys_data1(pos_u_cell)
    v1 = phys_data1(pos_v_cell)

    ! Velocity components 
    prim(pos_u) = -u1 
    prim(pos_v) = -v1

    ! Temperature(s)
    DO i = 1,nb_temp
       T1 = phys_data1(pos_T_cell + i - 1)
       T2 = phys_data2(pos_T_cell + i - 1)
       !Tghost = MAX(Twall - 0.5d0*(T2 - T1),0.9d0*Twall)
       Tghost = MAX(Twall - (T1 - Twall),0.8d0*Twall)
       prim(pos_T + i - 1) = Tghost
    ENDDO

    ! Species molar fractions of physical state
    CALL library_get_molar_fractions(u_phys1(1:nb_ns), xi)
   
    ! Species densities of ghost state computed based on the molar fractions 
    ! and pressure of physical state
    IF (nb_te.GE.1) THEN
        p_ov_T = phys_data1(pos_pres_cell)/(prim(pos_T)+xi(pos_em)*(prim(pos_Te)-prim(pos_T)))
    ELSE
    p_ov_T = phys_data1(pos_pres_cell)/prim(pos_T)
    ENDIF
    DO i = 1,nb_ns 
       prim(i) = p_ov_T*xi(i)/Ri(i)
    ENDDO

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost1, ghost_data1)

    ! Ghost state 2 

    ! Velocity components of physical state 1
    u2 = phys_data2(pos_u_cell)
    v2 = phys_data2(pos_v_cell)

    ! Velocity components 
    prim(pos_u) = -u2 
    prim(pos_v) = -v2

    ! Temperature(s)
    DO i = 1,nb_temp
       T1 = phys_data1(pos_T_cell + i - 1)
       T2 = phys_data2(pos_T_cell + i - 1)
       Tghost = MAX(3*Twall - 2*T2,0.9d0*Twall)
      !  Tghost = MAX(Twall - (T2 - Twall),0.8d0*Twall)
       prim(pos_T + i - 1) = Tghost
    ENDDO

    ! Species molar fractions of physical state
    CALL library_get_molar_fractions(u_phys2(1:nb_ns), xi)
   
    ! Species densities of ghost state computed based on the molar fractions 
    ! and pressure of physical state
    p_ov_T = (phys_data2(pos_pres_cell)/prim(pos_T))
    DO i = 1,nb_ns 
       prim(i) = p_ov_T*xi(i)/Ri(i)
    ENDDO

    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost2, ghost_data2)

    ! Comment above and uncomment below to have the second ghost state equal to the first
   !  u_ghost2 = u_ghost1
   !  ghost_data2 = ghost_data1

  END SUBROUTINE no_slip_iso_Twall_neq_1D_SL_Expl_2nd
!------------------------------------------------------------------------------!
SUBROUTINE no_slip_adiabatic_neq_1D_SL_Expl_2nd (id, phys_data1, u_phys1, phys_data2, u_phys2,   & 
                                                &  ghost_data1, u_ghost1, ghost_data2, u_ghost2)

    USE mod_general_data,          ONLY: nb_ns, nb_temp, nb_eq, pos_u, pos_v, pos_T, pos_Te, pos_u_cell, & 
                                      & pos_v_cell, pos_T_cell, pos_pres_cell, Ri, nb_te, pos_em, volumes
    USE mod_domain_boundary,       ONLY: boundary, boundary_data, get_boundary_Twall
    USE mod_function_pointer,      ONLY: get_cons_phys_from_prim 
    USE mod_neq_function_pointer,  ONLY: library_get_molar_fractions, library_get_thermal_cond

    IMPLICIT NONE

    INTEGER :: i
    REAL(KIND=8) :: Twall, p_ov_T, dx, f, df, emis, fac, P, T3
    REAL(KIND=8), DIMENSION(nb_ns) :: xi, rhoi_wall
    REAL(KIND=8), DIMENSION(nb_eq) :: prim
    REAL(KIND=8), DIMENSION(nb_temp) :: Twall_vec, lambda_vec
    REAL(KIND=8), PARAMETER :: sigma = 5.670373D-8
    TYPE(boundary) :: bound

    INTEGER, INTENT(IN) :: id                              !< boundary id
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys1       !< conservative variables of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys2      !< conservative variables of physical cell
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data1    !< physical properties of physical cell 
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: phys_data2   !< physical properties of physical cell 
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost1     !< conservative variables of ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost2    !< conservative variables of ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data1  !< physical properties of ghost cell
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: ghost_data2  !< physical properties of ghost cell

    ! Some constants
    emis  = 0.85d0
    dx = 0.5d0*volumes(3)
    fac = dx*emis*sigma
    
    ! Species molar fractions and pressure of physical state
    CALL library_get_molar_fractions(u_phys1(1:nb_ns), xi)
    P = phys_data1(pos_pres_cell)
    
    ! Initial guess
    Twall = phys_data1(pos_T_cell)
    
    ! Compute the SEB
    Twall_vec = Twall
    rhoi_wall = xi * P / (Ri * Twall)  ! We are enforcing thermal equil. here
    CALL library_get_thermal_cond(rhoi_wall, Twall_vec, lambda_vec)
    T3 = Twall*Twall*Twall
    f = fac*T3*Twall + sum(lambda_vec *(Twall-phys_data1(pos_T_cell:pos_T_cell+nb_temp-1)))
    
    ! Newton loop on Twall
    DO WHILE (abs(f/Twall) > 1.0E-12)
        ! Update Twall
        df = 4.0d0*fac*T3 + sum(lambda_vec) ! Assuming dlam/dTwall*(Twall-Tcell) is small
        Twall = Twall - f/df
        
        ! Compute SEB
        Twall_vec = Twall
        rhoi_wall = xi * P / (Ri * Twall)  ! We are enforcing thermal equil. here
        CALL library_get_thermal_cond(rhoi_wall, Twall_vec, lambda_vec)
        T3 = Twall*Twall*Twall
        f = fac*T3*Twall + sum(lambda_vec *(Twall-phys_data1(pos_T_cell:pos_T_cell+nb_temp-1)))
    ENDDO

    IF (Twall > 3800.0) THEN
      Twall = 3800.0
    ENDIF
    ! Compute the primitive variables at the ghost state
    ! Velocity components 
    prim(pos_u) = -phys_data1(pos_u_cell) 
    prim(pos_v) = -phys_data1(pos_v_cell)
    
    ! Temperatures
    DO i = 1,nb_temp
        prim(pos_T+i-1) = 2.0d0*Twall - phys_data1(pos_T_cell+i-1)
    ENDDO
    
    ! Species densities
    if (nb_te .ge. 1) then
        p_ov_T = P / (prim(pos_T) + xi(pos_em)*(prim(pos_Te)-prim(pos_T)))
    else
        p_ov_T = P / prim(pos_T) 
    endif
    prim(1:nb_ns) = p_ov_T * xi / Ri


    ! Compute the ghost state conservative variables and physical properties from primitive variables 
    CALL get_cons_phys_from_prim(prim, u_ghost1, ghost_data1)

    u_ghost2 = u_ghost1
    ghost_data2 = ghost_data1
    
    !write(*,*) "iso Twall boundary"
    !write(*,*) "Thg = ", prim(pos_T), "  Teg = ", prim(pos_Te)
    !write(*,*) "Thc = ", phys_data(pos_T_cell), "  Tec = ", phys_data(pos_T_cell+1)
    !DO i = 1,nb_ns 
    !   write(*,*) i , u_ghost(i), u_phys(i), xi(i)
    !ENDDO
        
  END SUBROUTINE no_slip_adiabatic_neq_1D_SL_Expl_2nd
  
