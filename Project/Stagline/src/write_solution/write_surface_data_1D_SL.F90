subroutine write_surface_data_1D_SL()
    USE mod_general_data, ONLY: nb_ns, nb_eq, nb_prop, nb_temp, start_u_phys, start_prop_phys, &
        & u, ios_dir, pos_T_cell, pos_pres_cell, pos_u_cell, pos_v_cell, xc, volumes, cell_prop, &
        & pos_ei_cell, simulation_id, Ri, diff_driv, Ji, nb_cells, pos_h0_cell, pos_rho_cell
    USE mod_numerics_data, ONLY: ull, ul, ur, urr, prop_ll, prop_l, prop_r, prop_rr, flag_full_Impl
    USE mod_neq_function_pointer, ONLY: library_get_species_DiffFlux, library_get_molar_fractions, & 
                                         & library_comp_tol, library_get_thermal_cond, &
                                         & library_get_prandtl_number, library_get_frozen_sound_speed
    USE mod_surface_properties
    USE mod_radiation,            ONLY:  wall_Rad
    USE mod_general_ablation


    implicit none
    
    integer, parameter :: out1 = 10
    integer :: i, boundary_edge_index
    real(KIND=8) :: p, velu, velv, du_dr, dv_dr, ov_dr, rc, rc_l, rc_r, vol_l, vol_r
    real(KIND=8) :: pl, pr, velul, velur, velvl, velvr, qconv, qcond, qdiff, qrado, qradi, qrad
    real(KIND=8), dimension(nb_temp) :: Tl, Tr, T, lambda, dT_dr
    real(KIND=8), dimension(nb_ns) :: rhoil, rhoir, rhoi, xil, xir, xi, hi
    real(KIND=8), dimension(nb_cells) :: u_velocity, v_velocity, dv_dy
    real(KIND=8), dimension(nb_cells-1) ::  ddv_ddy
    real(KIND=8), dimension(nb_cells-2) ::  d3v_d3y
    real(KIND=8) :: xc_boundaryEdge, enthalpy_boundaryEdge
    real(KIND=8) :: density_boundaryEdge, velocity_boundaryEdge
    real(KIND=8) :: enthalpy_wall, enthalpy_diff
    real(KIND=8) :: dtau_dx, dp_dxx
    real(KIND=8) :: prandtl, speed_of_sound, heat_transfer_coeff
    real(KIND=8), dimension(nb_ns) :: rhoi_boundaryEdge
    real(KIND=8), dimension(nb_temp) :: temperature_boundaryEdge

    INTERFACE 
      SUBROUTINE apply_bc_1D_SL_Expl(bound_id, phys_data1, u_phys1, phys_data2, u_phys2, &
                                   & ghost_data1, u_ghost1, ghost_data2, u_ghost2)
        INTEGER, INTENT(IN) :: bound_id
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys1, phys_data1, u_phys2, phys_data2 
        REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost1, ghost_data1, u_ghost2, ghost_data2
      END SUBROUTINE apply_bc_1D_SL_Expl
      
      !SUBROUTINE apply_bc_1D_SL_Impl(bound_id, phys_data1, u_phys1, phys_data2, u_phys2, & 
      !                          & ghost_data1, u_ghost1, ghost_data2, u_ghost2, jb)
      !  INTEGER, INTENT(IN) :: bound_id
      !  REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_phys1, phys_data1, u_phys2, phys_data2 
      !  REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u_ghost1, ghost_data1, u_ghost2, ghost_data2
      !  REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: jb
      !END SUBROUTINE apply_bc_1D_SL_Impl
    END INTERFACE 

    ! Get conserved and property variables for right and right-right cells
    ur      = u(start_u_phys+1:start_u_phys+nb_eq)
    urr     = u(start_u_phys+nb_eq+1:start_u_phys+2*nb_eq)
    prop_r  = cell_prop(start_prop_phys+1:start_prop_phys+nb_prop)
    prop_rr = cell_prop(start_prop_phys+nb_prop+1:start_prop_phys+2*nb_prop)

    ! Get the ghost cell values
    !if (flag_full_Impl) then
    !    call apply_bc_1D_SL_Impl(1, prop_r, ur, prop_rr, urr, prop_l, ul, prop_ll, ull)
    !else
        call apply_bc_1D_SL_Expl(1, prop_r, ur, prop_rr, urr, prop_l, ul, prop_ll, ull)
    !endif 

    ! Volumes and positions of left and right state
    vol_l = volumes(2) 
    vol_r = volumes(3)    
    rc_l  = xc(2)
    rc_r  = xc(3)
    rc    = rc_l + 0.5d0*vol_l
    
    ! Left variables
    Tl    = prop_l(pos_T_cell:pos_T_cell+nb_temp-1)
    pl    = prop_l(pos_pres_cell)
    velul    = prop_l(pos_u_cell)
    velvl    = prop_l(pos_v_cell)
    rhoil = ul(1:nb_ns)
    call library_get_molar_fractions(rhoil, xil)
    call library_comp_tol(xil)
    
    ! Right variables
    Tr  = prop_r(pos_T_cell:pos_T_cell+nb_temp-1)
    pr  = prop_r(pos_pres_cell)
    velur  = prop_r(pos_u_cell)
    velvr  = prop_r(pos_v_cell)
    rhoir = ur(1:nb_ns)
    call library_get_molar_fractions(rhoir, xir)
    call library_comp_tol(xir)
    
    ! Surface variables
    T   = (Tl*vol_l+Tr*vol_r)/(vol_l+vol_r)
    DO i= 1,nb_temp
        T_surface(i) = T(1)
    ENDDO
    p   = (pl*vol_l+pr*vol_r)/(vol_l+vol_r)
    velu   = (velul*vol_l+velur*vol_r)/(vol_l+vol_r)
    u_surface = velu
    velv   = (velvl*vol_l+velvr*vol_r)/(vol_l+vol_r)
    v_surface = velv
    rhoi = (rhoil*vol_l+rhoir*vol_r)/(vol_l+vol_r)
    rho_i_surface = rhoi
    call library_get_molar_fractions(rhoi, xi)
    call library_comp_tol(xi)
    
    ! Compute thermal conducivity at the face
    call library_get_thermal_cond(rhoi, T, lambda)
    
    ! Gradients
    ov_dr = 2.d0/(vol_l + vol_r)
    du_dr = ov_dr*(velur - velul)
    dv_dr = ov_dr*(velvr - velvl)
    dT_dr = ov_dr*(Tr - Tl)
    diff_driv = ov_dr*(xir - xil)
    
    ! Compute the species mass diffusion flux
    CALL library_get_species_DiffFlux(p, T(1), T(nb_temp), xi, diff_driv, Ji)
    
    ! Compute the species enthalpies
    hi = (prop_l(pos_ei_cell:pos_ei_cell+nb_ns-1)*vol_l+prop_r(pos_ei_cell:pos_ei_cell+nb_ns-1)*vol_r)/(vol_l+vol_r)
    hi = hi + Ri*T(1)
    hi(1) = hi(1) + Ri(1)*(T(nb_temp)-T(1))
    
    ! Compute convective flux
    qcond = -sum(lambda*dT_dr)
    qdiff = sum(Ji*hi)
    qconv = qcond + qdiff
    
    ! Radiative Flux
    qrado = 0.85d0*5.670373D-8*T(1)**4
    qradi = wall_Rad
    qrad  = qrado - qradi

    CALL library_get_frozen_sound_speed (rhoi, T, speed_of_sound)


    DO i = 1, nb_cells
        u_velocity(i) = cell_prop(start_prop_phys + nb_prop*(i-1) + pos_u_cell)
        v_velocity(i) = cell_prop(start_prop_phys + nb_prop*(i-1) + pos_v_cell)
        dv_dy(i) = (v_velocity(i)+u_velocity(i))/ xc(i+2)
    END DO

    ddv_ddy = 0.0D0
    DO i = 2, nb_cells-1

        ddv_ddy(i) = (dv_dy(i+1) - dv_dy(i-1))/(volumes(2+i+1)+2.0*volumes(2+i)+volumes(2+i-1))
    END DO

    d3v_d3y = 0.0D0
    DO i = 30, nb_cells-2

        d3v_d3y(i) = (ddv_ddy(i+1) - ddv_ddy(i-1))/(volumes(2+i+1)+2.0*volumes(2+i)+volumes(2+i-1))
        if (d3v_d3y(i)>0) then
            boundary_edge_index = i
            xc_boundaryEdge = xc(boundary_edge_index+2)
            EXIT
        end if
    END DO

    enthalpy_wall = cell_prop(start_prop_phys + pos_h0_cell)
    enthalpy_boundaryEdge = cell_prop(start_prop_phys + nb_prop*(boundary_edge_index-1) + pos_h0_cell)
    enthalpy_diff = enthalpy_boundaryEdge - enthalpy_wall
    density_boundaryEdge = cell_prop(start_prop_phys + nb_prop*(boundary_edge_index-1) + pos_rho_cell)
    velocity_boundaryEdge = cell_prop(start_prop_phys + nb_prop*(boundary_edge_index-1) + pos_u_cell)
    temperature_boundaryEdge = cell_prop(start_prop_phys + nb_prop*(boundary_edge_index-1) + pos_T_cell:pos_T_cell+nb_temp-1)

    DO i = 1,nb_ns
        rhoi_boundaryEdge(i) = u(start_u_phys+nb_eq*(boundary_edge_index-1) + i)
    ENDDO

    call library_get_prandtl_number(rhoi_boundaryEdge,temperature_boundaryEdge, prandtl)

    dtau_dx = (abs(qconv)*dv_dy(boundary_edge_index)/enthalpy_diff)*prandtl**(2.D0/3.D0)
    dp_dxx = -2.D0*p/(rc*rc)

    heat_transfer_coeff = qconv/(velocity_boundaryEdge*density_boundaryEdge*enthalpy_diff)

    ! Write surface properties to a file
     IF (ios_dir) THEN
        OPEN(UNIT=out1,FILE='./'//TRIM(simulation_id)//'_surface.dat',STATUS='unknown')
    ELSE
        OPEN(UNIT=out1,FILE='../output/'//TRIM(simulation_id)//'_surface.dat',STATUS='unknown')
    ENDIF
    
    write(out1,*) 'Surface Properties'
    write(out1,100) 'Wall Position [m]', rc
    write(out1,100) 'Temperatures [K]', T
    write(out1,100) 'Pressure [Pa]', P
    write(out1,100) 'Densities [kg/m3]', rhoi
    write(out1,100) 'Mole Fractions', xi
    write(out1,100) 'u [m/s]', velu
    write(out1,100) 'v [m/s]', velv
    write(out1,100) 'q_cond [W/m2]', qcond
    write(out1,100) 'q_diff [W/m2]', qdiff
    write(out1,100) 'q_conv [W/m2]', qconv
    write(out1,100) 'q_rado [W/m2]', qrado
    write(out1,100) 'q_radi [W/m2]', qradi
    write(out1,100) 'q_rad [W/m2]', qrad
    write(out1,100) 'q_tot [W/m2]', qconv + qrad
    write(out1,100) 'h_wall [J/kg]', enthalpy_wall
    write(out1,100) 'mdot [kg/m2 s]',mdot
    write(out1,100) 'mdot_i [kg/m2 s]',mdot_i
    write(out1,100) 'dtau_dx [Pa/m]',dtau_dx
    write(out1,100) 'dp_dxx [Pa/m2]',dp_dxx
    write(out1,100) 'Mach', velu/speed_of_sound

    write(out1,*) 'Edge Properties'
    write(out1,100) 'BL_thickness [m]', xc_boundaryEdge-rc
    write(out1,100) 'h_edge [J/kg]', enthalpy_boundaryEdge
    write(out1,100) 'h_diff [J/kg]', enthalpy_diff
    write(out1,100) 'Temperatures_edge [K]', temperature_boundaryEdge
    write(out1,100) 'Densitiy_edge [kg/m3]', density_boundaryEdge
    write(out1,100) 'u_edge [m/s]', velocity_boundaryEdge
    write(out1,100) 'heat_transfer_coeff [-]',heat_transfer_coeff

    close(out1)

    100 format (A30,': ', 100E15.7)
end subroutine
