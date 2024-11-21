!------------------------------------------------------------------------------!
!> This subroutine writes the transport fluxes for 1D stagnation line nonequilibrium flows. 
  SUBROUTINE write_transp_flux_neq_1D_SL ()

    USE mod_general_data,              ONLY: nb_cells, nb_dim, nb_prop, nb_eq, nb_ns, nb_int_temp,       & 
                                           & nb_temp, start_u_phys, start_prop_phys, pos_u_cell,         & 
                                           & pos_v_cell, pos_T_cell, pos_mu_cell, pos_lambda_cell,       & 
                                           & pos_pres_cell, pos_ei_cell, pos_em, Ri, volumes, xc,        & 
                                           & u, cell_prop, simulation_id, ios_dir
    USE mod_neq_function_pointer,      ONLY: library_get_molar_fractions, library_get_species_DiffFlux,  & 
                                           & library_comp_tol, library_species_name
    USE mod_function_pointer,          ONLY: get_stress_tensor_1D_SL


    IMPLICIT NONE

    INTEGER, PARAMETER :: out1 = 20
    INTEGER :: i, j, k 
    REAL(KIND=8) :: r, mu, tmp, sum_Ji
    REAL(KIND=8) :: uc, ul, ur, vc, vl, vr, Tc, Tec, pc
    REAL(KIND=8) :: vol_l, vol_c, vol_r
    REAL(KIND=8) :: ov_dr, du_dr, dv_dr, dT_dr
    REAL(KIND=8) :: tau_rr, tau_rt, tau_tt
    REAL(KIND=8) :: q, q_Diff, q_Diff_int, q_Fourier_tr, q_Fourier_int 
    REAL(KIND=8), DIMENSION(nb_prop) :: prop_l, prop_c, prop_r
    REAL(KIND=8), DIMENSION(nb_temp) :: templ, tempc, tempr, grad_T, lambda_vec
    REAL(KIND=8), DIMENSION(nb_ns) :: xil, xic, xir, rhoil, rhoic, rhoir, hic, Vi
    REAL(KIND=8), DIMENSION(nb_ns*nb_dim) :: diff_driv, Ji


    CHARACTER(LEN=20),DIMENSION(nb_ns) :: velo_diff_filename
    CHARACTER(LEN=30):: temp 
    CHARACTER(LEN=1000):: velo_diff_ans
    CHARACTER(LEN=10), DIMENSION(nb_ns):: species_name


    IF (ios_dir) THEN
        OPEN(UNIT=out1,FILE='./'//TRIM(simulation_id)//'_transp_fluxes.dat',STATUS='unknown')
    ELSE
    OPEN(UNIT=out1,FILE='../output/'//TRIM(simulation_id)//'_transp_fluxes.dat',STATUS='unknown')
    ENDIF

    CALL library_species_name(species_name)

    DO i=1,nb_ns
    WRITE (velo_diff_filename(i) ,"('Vdi(',A,') ')")TRIM(species_name(i))
    ENDDO 
    temp=''
    velo_diff_ans=''
   
    DO i=1,int(nb_ns/2) 
    temp = ' "'//TRIM(velo_diff_filename(i*2-1))//'"  "'// TRIM(velo_diff_filename(i*2))//'"'
    velo_diff_ans = TRIM(velo_diff_ans)//TRIM(temp)
    ENDDO
   
   IF (MOD(nb_ns,2).NE.0) THEN   

   velo_diff_ans= TRIM(velo_diff_ans)// ' "' //TRIM(velo_diff_filename(nb_ns))//'"'

   ENDIF


    ! Flowfield file header
    WRITE(out1,10)'# solver_fvmcc_F90:: 1D stagnation line solution file'
    WRITE(out1,10)'# Nonequilibrium gas flow'
    WRITE(out1,10)'# Transport fluxes'
    WRITE(out1,10)'# "xc"  "sum(Ji)"  "tau_rr"  "tau_rt"  "tau_tt"  "q_r"  "q_r_Ftr"  "q_r_Fint"  "q_r_Diff"'&
                                                                                     &//TRIM(velo_diff_ans)//''

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
    vl = prop_l(pos_v_cell)
    DO i = 1,nb_temp
       templ(i) = prop_l(pos_T_cell + i - 1)
    ENDDO
    vol_l = volumes(2)   
 
    uc = prop_c(pos_u_cell)
    vc = prop_c(pos_v_cell)
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
       vr = prop_r(pos_v_cell)
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
       ov_dr = 2.d0/(vol_l + 2.d0*vol_c + vol_r)
       du_dr = ov_dr*(ur - ul)
       dv_dr = ov_dr*(vr - vl)
       DO j = 1,nb_temp
          grad_T(j) = (tempr(j) - templ(j))*ov_dr
       ENDDO

       ! Diffusion driving forces
       diff_driv = 0.d0
       DO j = 1,nb_ns 
          diff_driv(j) = (xir(j) - xil(j))*ov_dr
       ENDDO

       ! Species mass diffusion flux
       CALL library_get_species_DiffFlux(pc, tempc(1), tempc(nb_temp), xic, diff_driv, Ji)

       ! Species diffusion velocities
       sum_Ji = 0.d0
       DO j = 1,nb_ns
          tmp = Ji(j)
          sum_Ji = sum_Ji + tmp 
          Vi(j)  = tmp/rhoic(j)
       ENDDO

       ! Stress tensor
       CALL get_stress_tensor_1D_SL (mu, r, uc, vc, du_dr, dv_dr, tau_rr, tau_rt, tau_tt)

       ! Fourier heat flux (component associated to translational energy)  
       q_Fourier_tr = - lambda_vec(1)*grad_T(1)

       ! Fourier heat flux (component associated to internal energy)
       q_Fourier_int = 0.d0
       DO k = 1,nb_temp-1
          q_Fourier_int = q_Fourier_int - lambda_vec(k + 1)*grad_T(k + 1)
       ENDDO

       ! Diffusive heat flux
       q_Diff = 0.d0
       DO j = 1,nb_ns
          q_Diff = q_Diff + Ji(j)*hic(j)
       ENDDO

       ! Total heat flux
       q = q_Fourier_tr + q_Fourier_int + q_Diff

       WRITE(out1,20)r,sum_Ji,tau_rr,tau_rt,tau_tt,q,q_Fourier_tr,q_Fourier_int,q_Diff,Vi


       ! Data exchange
       prop_l = prop_c
       prop_c = prop_r

       ul = uc 
       vl = vc
       xil = xic
       rhoil = rhoic
       templ = tempc
       vol_l = vol_c

       uc = ur 
       vc = vr
       xic = xir
       rhoic = rhoir
       tempc = tempr
       vol_c = vol_r

    ENDDO   

    CLOSE(out1)

10  FORMAT(A)
20  FORMAT(10000(E20.10,1X))

  END SUBROUTINE write_transp_flux_neq_1D_SL
!------------------------------------------------------------------------------!
