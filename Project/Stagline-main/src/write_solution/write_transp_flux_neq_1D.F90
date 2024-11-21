!------------------------------------------------------------------------------!
!> This subroutine writes the transport fluxes for 1D nonequilibrium flows. 
  SUBROUTINE write_transp_flux_neq_1D ()

    USE mod_general_data,              ONLY: nb_cells, nb_prop, nb_eq, nb_ns, nb_int_temp, nb_temp,     & 
                                           & start_u_phys, start_prop_phys, pos_u_cell, pos_T_cell,     & 
                                           & pos_mu_cell, pos_lambda_cell, pos_pres_cell,               & 
                                           & pos_ei_cell, pos_em, Ri, volumes, xc, u, cell_prop,        &
                                           & simulation_id, ios_dir
    USE mod_neq_function_pointer,      ONLY: library_get_molar_fractions, library_get_species_DiffFlux, & 
                                           & library_comp_tol

    IMPLICIT NONE

    INTEGER, PARAMETER :: out1 = 20
    INTEGER :: i, j, k 
    REAL(KIND=8), PARAMETER :: coeff = 4.d0/3.d0
    REAL(KIND=8) :: r, mu, tmp, sum_Ji
    REAL(KIND=8) :: uc, ul, ur, Tc, Tec, pc
    REAL(KIND=8) :: vol_l, vol_c, vol_r
    REAL(KIND=8) :: ov_dx, du_dx, dv_dx, dT_dx
    REAL(KIND=8) :: tau_xx
    REAL(KIND=8) :: q, q_Diff, q_Diff_int, q_Fourier_tr, q_Fourier_int 
    REAL(KIND=8), DIMENSION(nb_prop) :: prop_l, prop_c, prop_r
    REAL(KIND=8), DIMENSION(nb_temp) :: templ, tempc, tempr, grad_T, lambda_vec
    REAL(KIND=8), DIMENSION(nb_ns) :: xil, xic, xir, rhoil, rhoic, rhoir, hic, Vi
    REAL(KIND=8), DIMENSION(nb_ns) :: diff_driv, Ji

    IF (ios_dir) THEN
        OPEN(UNIT=out1,FILE='./'//TRIM(simulation_id)//'_transp_fluxes.dat',STATUS='unknown')
    ELSE
        OPEN(UNIT=out1,FILE='../output/'//TRIM(simulation_id)//'_transp_fluxes.dat',STATUS='unknown')
    ENDIF


    ! Flowfield file header
    WRITE(out1,10)'# solver_fvmcc_F90:: 1D solution file'
    WRITE(out1,10)'# Nonequilibrium gas flow'
    WRITE(out1,10)'# Transport fluxes'
    WRITE(out1,10)'# xc  Vdi  sum(Ji)  tau_xx  q_x  q_x_Ftr  q_x_Fint  q_x_Diff'

    ! Initialization of left and central state
    ! Species densities and molar fractions
    DO i = 1,nb_ns 
       rhoil(i) = u(nb_eq + i)
       rhoic(i) = u(2*nb_eq + i)
    ENDDO

    CALL library_get_molar_fractions (rhoil, xil)
    CALL library_comp_tol(xil) 
   
    CALL library_get_molar_fractions (rhoic, xic)
    CALL library_comp_tol(xic) 

    DO i = 1,nb_prop
       prop_l(i) = cell_prop(nb_prop + i)
       prop_c(i) = cell_prop(2*nb_prop + i)
    ENDDO

    ! Velocity components and temperatures
    ul = prop_l(pos_u_cell)
    DO i = 1,nb_temp
       templ(i) = prop_l(pos_T_cell + i - 1)
    ENDDO
    vol_l = volumes(2)   
 
    uc = prop_c(pos_u_cell)
    DO i = 1,nb_temp
       tempc(i) = prop_c(pos_T_cell + i - 1)
    ENDDO
    vol_c = volumes(3)

    ! Loop over physical cells
    DO i = 1,nb_cells

       ! Right state
       DO j = 1,nb_ns
          rhoir(j) = u(start_u_phys + i*nb_eq + j)
       ENDDO

       CALL library_get_molar_fractions (rhoir, xir)
       CALL library_comp_tol(xir) 

       DO j = 1,nb_prop
          prop_r(j) = cell_prop(start_prop_phys + i*nb_prop + j)
       ENDDO 

       ur = prop_r(pos_u_cell)
       DO j = 1,nb_temp
          tempr(j) = prop_r(pos_T_cell + j - 1)
       ENDDO
       vol_r = volumes(3 + i)

       ! Central state
       ! Cell centroid 
       r = xc(i + 2)

       ! Species specific enthalpies
       Tc  = prop_c(pos_T_cell)
       Tec = prop_c(pos_T_cell + nb_temp - 1)

       DO j = 1,pos_em - 1
          hic(j) = prop_c(pos_ei_cell + j - 1) + Ri(j)*Tc
       ENDDO 

       hic(pos_em) = prop_c(pos_ei_cell + pos_em - 1) + Ri(pos_em)*Tec

       DO j = pos_em + 1,nb_ns
          hic(j) = prop_c(pos_ei_cell + j - 1) + Ri(j)*Tc
       ENDDO  

       ! Dynamic viscosity and thermal conductivity components
       mu = prop_c(pos_mu_cell)
       DO j = 1,nb_temp
          lambda_vec(j) = prop_c(pos_lambda_cell + j - 1)
       ENDDO

       ! Pressure
       pc = prop_c(pos_pres_cell) 
       
       ! Velocity and temperature gradients
       ov_dx = 2.d0/(vol_l + 2.d0*vol_c + vol_r)
       du_dx = ov_dx*(ur - ul)
       DO j = 1,nb_temp
          grad_T(j) = (tempr(j) - templ(j))*ov_dx
       ENDDO

       ! Diffusion driving forces
       diff_driv = 0.d0
       DO j = 1,nb_ns 
          diff_driv(j) = (xir(j) - xil(j))*ov_dx
       ENDDO

       ! Species mass diffusion flux
       CALL library_get_species_DiffFlux(pc, Tc, Tec, xic, diff_driv, Ji)

       ! Species diffusion velocities
       sum_Ji = 0.d0
       DO j = 1,nb_ns
          tmp = Ji(j)
          sum_Ji = sum_Ji + tmp 
          Vi(j)  = tmp/rhoic(j)
       ENDDO

       ! Stress tensor
       tau_xx = coeff*mu*du_dx

       ! Fourier heat flux (component associated to translational energy)  
       q_Fourier_tr = - lambda_vec(1)*grad_T(1)

       ! Fourier heat flux (component associated to internal energy)
       q_Fourier_int = 0.d0
       DO k = 1,nb_int_temp
          q_Fourier_int = q_Fourier_int - lambda_vec(k + 1)*grad_T(k + 1)
       ENDDO

       ! Diffusive heat flux
       q_Diff = 0.d0
       DO j = 1,nb_ns
          q_Diff = q_Diff + Ji(j)*hic(j)
       ENDDO

       ! Total heat flux
       q = q_Fourier_tr + q_Fourier_int + q_Diff

       WRITE(out1,20)r,Vi,sum_Ji,tau_xx,q,q_Fourier_tr,q_Fourier_int,q_Diff

       ! Data exchange
       prop_l = prop_c
       prop_c = prop_r

       ul = uc 
       xil = xic
       rhoil = rhoic
       templ = tempc
       vol_l = vol_c

       uc = ur 
       xic = xir
       rhoic = rhoir
       tempc = tempr
       vol_c = vol_r

    ENDDO   

    CLOSE(out1)

10  FORMAT(A)
20  FORMAT(10000(E20.10,1X))

  END SUBROUTINE write_transp_flux_neq_1D
!------------------------------------------------------------------------------!
