!------------------------------------------------------------------------------!
!> This subroutine writes the transport fluxes for 1D stagnation line calorically perfect gas flows. 
  SUBROUTINE write_transp_flux_1D_SL ()

    USE mod_general_data,              ONLY: nb_cells, nb_prop, start_prop_phys, pos_u_cell, pos_v_cell, &
                                          &  pos_T_cell, pos_mu_cell, pos_lambda_cell, volumes, xc,      & 
                                          &  cell_prop, simulation_id, ios_dir
    USE mod_function_pointer,          ONLY: get_stress_tensor_1D_SL

    IMPLICIT NONE

    INTEGER, PARAMETER :: out1 = 20
    INTEGER :: i, j 
    REAL(KIND=8) :: uc, ul, ur, vc, vl, vr, Tl, Tc, Tr
    REAL(KIND=8) :: vol_l, vol_c, vol_r
    REAL(KIND=8) :: ov_dr, du_dr, dv_dr, dT_dr
    REAL(KIND=8) :: tau_rr, tau_rt, tau_tt, q_r
    REAL(KIND=8) :: r, mu, lambda 
    REAL(KIND=8), DIMENSION(nb_prop) :: prop_l, prop_c, prop_r

    IF (ios_dir) THEN
        OPEN(UNIT=out1,FILE='./'//TRIM(simulation_id)//'_transp_fluxes.dat',STATUS='unknown')
    ELSE
        OPEN(UNIT=out1,FILE='../output/'//TRIM(simulation_id)//'_transp_fluxes.dat',STATUS='unknown')
    ENDIF

    ! Flowfield file header
    WRITE(out1,10)'# solver_fvmcc_F90:: 1D stagnation line solution file'
    WRITE(out1,10)'# Calorically perfect gas flow'
    WRITE(out1,10)'# Transport fluxes'
    WRITE(out1,10)'# xc  tau_rr  tau_rt tau_tt  q_r'

    ! Initialization of left and center state
    DO i = 1,nb_prop
       prop_l(i) = cell_prop(nb_prop + i)
       prop_c(i) = cell_prop(2*nb_prop + i)
    ENDDO

    ! Velocity components and temperatures
    ul = prop_l(pos_u_cell)
    vl = prop_l(pos_v_cell)
    Tl = prop_l(pos_T_cell)
    vol_l = volumes(2)   
 
    uc = prop_c(pos_u_cell)
    vc = prop_c(pos_v_cell)
    Tc = prop_c(pos_T_cell)
    vol_c = volumes(3)

    DO i = 1,nb_cells

       ! Right state
       DO j = 1,nb_prop
          prop_r(j) = cell_prop(start_prop_phys + i*nb_prop + j)
       ENDDO 

       ur = prop_r(pos_u_cell)
       vr = prop_r(pos_v_cell)
       Tr = prop_r(pos_T_cell)
       vol_r = volumes(3 + i)

       ! Central state 
       ! Cell centroid location
       r = xc(i + 2)

       ! Dynamic viscosity and thermal conductivity
       mu = prop_c(pos_mu_cell)
       lambda = prop_c(pos_lambda_cell)

       ! Velocity and temperature gradients
       ov_dr = 2.d0/(vol_l + 2.d0*vol_c + vol_r)
       du_dr = ov_dr*(ur - ul)
       dv_dr = ov_dr*(vr - vl)
       dT_dr = ov_dr*(Tr - Tl)

       ! Stress tensor
       CALL get_stress_tensor_1D_SL (mu, r, uc, vc, du_dr, dv_dr, tau_rr, tau_rt, tau_tt)
 
       ! Heat flux
       q_r = - lambda*dT_dr

       WRITE(out1,20)r,tau_rr,tau_rt,tau_tt,q_r

       ! Data exchange
       prop_l = prop_c
       prop_c = prop_r

       ul = uc 
       vl = vc
       Tl = Tc
       vol_l = vol_c

       uc = ur 
       vc = vr
       Tc = Tr
       vol_c = vol_r

    ENDDO   

    CLOSE(out1)

10  FORMAT(A)
20  FORMAT(100(E20.10,1X))

  END SUBROUTINE write_transp_flux_1D_SL
!------------------------------------------------------------------------------!
