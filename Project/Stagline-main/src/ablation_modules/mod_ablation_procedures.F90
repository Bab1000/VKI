  MODULE mod_ablation_procedures


#include "../config.h"

#ifdef CARBONABLA

    IMPLICIT NONE

    ABSTRACT INTERFACE


      SUBROUTINE proc_compute_recombination(rhoi_wall, T_wall)
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rhoi_wall
        REAL(KIND=8),               INTENT(IN) :: T_wall
      END SUBROUTINE proc_compute_recombination

    END INTERFACE

    PROCEDURE(proc_compute_recombination), POINTER, SAVE :: compute_recombination 

    CONTAINS


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SUBROUTINE ablation_initialize (nb_spec)

        USE mod_physical_model,         
        USE mod_ablation_data,          ONLY: omega_wall, omega_wall_rec, ion_species_matrix,                           &
                                            & C_pos, C3_pos, CN_pos, CO_pos, N_pos, N2_pos, O_pos, O2_pos, e_minus_pos, &
                                            & abla_dir
        USE mod_neq_function_pointer,   ONLY: library_get_species_index, library_get_species_name

        INTEGER, INTENT(IN) :: nb_spec                          !< species number 

        CHARACTER*80 :: mixture, abla_dir_dumm
        CHARACTER*20 :: sp_name, sp_name_2
        INTEGER      :: i, j, k
        INTEGER      :: ion_index, sp_index, tail, ion_num
        INTEGER      :: lenght, var_stat
        CHARACTER(len=5), ALLOCATABLE, DIMENSION(:) :: sp_name_vector

        ! Allocation of wall source term vectors
        ALLOCATE(omega_wall(nb_spec))
        ALLOCATE(omega_wall_rec(nb_spec))
        ! Allocation of species name vector
        ALLOCATE(sp_name_vector(nb_spec))

        omega_wall = 0.d0
        omega_wall_rec = 0.d0

        mixture = get_mixture(physical_model)

        ion_num = 0
        k = 1
        DO i=1,nb_spec
           CALL library_get_species_name(i,sp_name)
           sp_name_vector(i) = sp_name
        END DO

        ! Define the number of ions in the mixture
        DO i=1,nb_spec
           sp_name = sp_name_vector(i)(1:LEN_TRIM(sp_name_vector(i))) 
           IF (IACHAR(sp_name(LEN_TRIM(sp_name):LEN_TRIM(sp_name))).EQ.IACHAR('+')) THEN
              ion_num = ion_num + 1
           ELSE
           ENDIF
        END DO

        ! Allocation of ion-species pair index matrix 
        ALLOCATE(ion_species_matrix(ion_num, 2))

        ! Create the ion-species pair and store it on a line of the matrix 
        DO i=1,nb_spec
           sp_name = sp_name_vector(i)(1:LEN_TRIM(sp_name_vector(i))) 
           IF (IACHAR(sp_name(LEN_TRIM(sp_name):LEN_TRIM(sp_name))).EQ.IACHAR('+')) THEN
              CALL library_get_species_index(sp_name, ion_index) 
              tail = LEN_TRIM(sp_name)-1
              sp_name = sp_name(1:tail)
              DO j=1,nb_spec
                 IF (i.NE.j) THEN
                    sp_name_2=sp_name_vector(j)(1:LEN_TRIM(sp_name_vector(j))) 
                    IF (sp_name_2 == sp_name) THEN
                       CALL library_get_species_index(sp_name_2, sp_index) 
                       ion_species_matrix(k,1) = ion_index
                       ion_species_matrix(k,2) = sp_index
                       k = k+1
                    ELSE
                    ENDIF
                 ELSE
                 ENDIF
              END DO
           ELSE
           ENDIF
        END DO


        ! Set the index of the species involved in the surface ablation and
        ! recombination reactions
        sp_name = 'C'
        CALL library_get_species_index(sp_name, C_pos)
        sp_name = 'C3'
        CALL library_get_species_index(sp_name, C3_pos)
        sp_name = 'CN'
        CALL library_get_species_index(sp_name, CN_pos)
        sp_name = 'CO'
        CALL library_get_species_index(sp_name, CO_pos)
        sp_name = 'N'
        CALL library_get_species_index(sp_name, N_pos)
        sp_name = 'N2'
        CALL library_get_species_index(sp_name, N2_pos)
        sp_name = 'O'
        CALL library_get_species_index(sp_name, O_pos)
        sp_name = 'O2'
        CALL library_get_species_index(sp_name, O2_pos)

        ! Set the index of the electron
        sp_name = 'e-'
        CALL library_get_species_index(sp_name, e_minus_pos)


        IF (ion_num.NE.0) THEN
                compute_recombination => compute_Nitrogen_and_Ions_recombination
        ELSE
                compute_recombination => compute_only_Nitrogen_recombination
        ENDIF

        ! Get the lenght of the path of the directory storing data useful for ablation
        ! using a dummy call to GET_ENVIRONMENT_VARIABLE
        CALL GET_ENVIRONMENT_VARIABLE("ABLA_DATA_DIR",abla_dir_dumm,lenght,var_stat)

        ! Allocate the string of the right dimension for the path 
        ALLOCATE(CHARACTER(LEN=lenght) :: abla_dir)

        ! Get the path of the directory storing data useful for ablation
        CALL GET_ENVIRONMENT_VARIABLE("ABLA_DATA_DIR",abla_dir,lenght,var_stat)

        IF (var_stat.NE.0) THEN 
           print*, 'ATTENTION!'
           print*, 'Environmental variable ABLA_DATA_DIR does not defined!'
           print*, 'Impossible to compute surface chemical heat flux.'
           print*, 'Computation is fine but the blowing heat flux is set to zero in the ablation summary file.'
           print*, 'Set the variable by going to the src/ablation_modules/ablation_data and typing: export ABLA_DATA_DIR=`pwd`.'
           print*, 'Alternatively store the variable in .bashrc file.'
        ENDIF

        END SUBROUTINE ablation_initialize

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SUBROUTINE ablation_finalize ()

        USE mod_ablation_data,          ONLY: omega_wall, omega_wall_rec, ion_species_matrix
        USE mod_ablation_pointers 


        ! Deallocation of wall source term vectors
        IF (ALLOCATED(omega_wall))       DEALLOCATE(omega_wall)
        IF (ALLOCATED(omega_wall_rec))   DEALLOCATE(omega_wall_rec)

        ! Deallocating the ion-species pair index matrix
        IF (ALLOCATED(ion_species_matrix))   DEALLOCATE(ion_species_matrix)

        ! Nullification of ablation procedure pointers
        IF(ASSOCIATED(procedure_ablation_initialize))                    NULLIFY(procedure_ablation_initialize)
        IF(ASSOCIATED(procedure_ablation_write_sol))                     NULLIFY(procedure_ablation_write_sol)
        IF(ASSOCIATED(procedure_ablation_compute_wall_source_terms))     NULLIFY(procedure_ablation_compute_wall_source_terms)
        IF(ASSOCIATED(procedure_ablation_compute_wall_SEB))              NULLIFY(procedure_ablation_wall_diff_flux)
        IF(ASSOCIATED(library_ablation_compute_wall_mass_blowing_rate))  NULLIFY(library_ablation_compute_wall_mass_blowing_rate)
        IF(ASSOCIATED(procedure_ablation_finalize))                      NULLIFY(procedure_ablation_finalize)
        IF(ASSOCIATED(compute_recombination))                            NULLIFY(compute_recombination)
        IF(ASSOCIATED(procedure_ablation_SMB))                           NULLIFY(procedure_ablation_SMB)

        END SUBROUTINE ablation_finalize

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SUBROUTINE ablation_write_sol ()

        USE mod_general_data,           ONLY: simulation_id, ios_dir, nb_temp
        USE mod_ablation_data,          ONLY: mdot_wall, cond_flux_wall, mdot_wall_tot, Twall, seb_stored, N_prob, gamma_rec

        IMPLICIT NONE

        INTEGER, PARAMETER :: out1 = 11
        INTEGER :: j
        REAL(KIND=8) :: q_D,  q_C
        REAL(KIND=8), DIMENSION(nb_temp)  ::  q_F

        CALL ablation_comp_surface_flux (q_D, q_F, q_C)


        IF (ios_dir) THEN
        OPEN(UNIT=out1,FILE='./'//TRIM(simulation_id)//'_ablation_data.dat',STATUS='unknown')
        ELSE
        OPEN(UNIT=out1,FILE='../output/'//TRIM(simulation_id)//'_ablation_data.dat',STATUS='unknown')
        ENDIF

        WRITE(out1,10)'# solver_fvmcc_F90:: 1D stagnation line with ablative boundary wall data'
        WRITE(out1,10)
        WRITE(out1,10)'# SPECIFIC MASS BLOWING RATES (kg/m^2-s)'
        WRITE(out1,10)'# Reaction C_s + O => CO'
        WRITE(out1,20) mdot_wall(1)
        WRITE(out1,10)'# Reaction 2C_s + O_2 => 2CO '
          WRITE(out1,20) mdot_wall(2)
        WRITE(out1,10)'# Reaction C_s + N => CN '
        WRITE(out1,20) mdot_wall(3)
        WRITE(out1,10)'# Reaction 3C_s => C_3'
        WRITE(out1,20) mdot_wall(4)

        WRITE(out1,10)

        WRITE(out1,10)'# TOTAL MASS BLOWING RATE (kg/m^2-s)' 
        WRITE(out1,20) mdot_wall_tot 

        WRITE(out1,10)

        WRITE(out1,10)'# Actual nitridation probability (only for reactivity rebuilding)' 
        WRITE(out1,20) N_prob

        WRITE(out1,10)

        WRITE(out1,10)'# Actual wall recombination probability' 
        WRITE(out1,20) gamma_rec

        WRITE(out1,10)

        WRITE(out1,10)'# WALL TEMPERATURE (K)'
        WRITE(out1,20) Twall

        WRITE(out1,10)

        WRITE(out1,10)'# WALL CONDUCTIVE HEAT FLUX (W/m^2)'
        IF (nb_temp.GT.1) THEN
           DO j = 1,nb_temp
              IF (j.EQ.1) WRITE(out1,10)'   # Translational contribution'
              IF (j.GT.1) WRITE(out1,10)'   # Internal contribution'
              WRITE(out1,20) q_F(j)
           ENDDO
        ELSE
           WRITE(out1,20) q_F(1)
        ENDIF


        WRITE(out1,10)

        WRITE(out1,10)'# WALL DIFFUSIVE HEAT FLUX (W/m^2)'
        WRITE(out1,20) q_D 

        WRITE(out1,10)

        WRITE(out1,10)'# WALL BLOWING (mdot*(hw-hs)) HEAT FLUX (W/m^2)'
        WRITE(out1,20) q_C 

        WRITE(out1,10)

        WRITE(out1,10)'# WALL TOTAL HEAT FLUX (W/m^2)'
        WRITE(out1,20) SUM(q_F) + q_D + q_C

        WRITE(out1,10)

        IF (seb_stored .NE. 0.d0) THEN
                WRITE(out1,10)'# WALL SURFACE ENERGY BALANCE (W/m^2)'
                WRITE(out1,20) seb_stored 
        ELSE
        ENDIF

        CLOSE(out1)

10  FORMAT(A)
20  FORMAT(10000(E20.10,1X))


        END SUBROUTINE ablation_write_sol

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SUBROUTINE ablation_SMB(u_phys1,p_wall,rhoi_wall)

        USE mod_general_data,          ONLY: nb_ns, volumes, nb_temp, nb_dim
        USE mod_algebra,               ONLY: lu_decmp, lu_solver
        USE mod_ablation_data,         ONLY: omega_wall, uwall, Twall, phi_pyro, mdot_wall_tot
        USE mod_ablation_pointers,     ONLY: library_ablation_compute_wall_mass_blowing_rate, &
                                                     & procedure_ablation_compute_wall_source_terms
        USE mod_neq_function_pointer,  ONLY: library_get_molar_fractions, library_get_species_DiffFlux,  & 
                                           & library_compute_eq_composition_pyro, library_comp_tol

        INTEGER ::  i, j, k
        INTEGER ::  smb_it_out
        REAL(KIND=8) :: smb_out
        REAL(KIND=8) :: rhowall, smb_max_old
        REAL(KIND=8) :: u_wall_old
        REAL(KIND=8) :: delta
        REAL(KIND=8) :: vol_p
        REAL(KIND=8) :: ov_dr
        REAL(KIND=8):: d
        INTEGER, DIMENSION(nb_ns) :: indx
        REAL(KIND=8), DIMENSION(nb_ns) :: yi_pyro
        REAL(KIND=8), DIMENSION(nb_ns) :: xi_w, xi_p !< molar fractions of wall and first physical states
        REAL(KIND=8), DIMENSION(nb_ns*nb_dim) :: diff_driv, Ji
        REAL(KIND=8), DIMENSION(nb_ns) :: conv, pyro
        REAL(KIND=8), DIMENSION(nb_temp) :: tempw
        REAL(KIND=8), DIMENSION(nb_ns) :: smb, smb_pert
        REAL(KIND=8), DIMENSION(nb_ns) :: smb_old
        REAL(KIND=8), DIMENSION(nb_ns,nb_ns) :: smb_jac
        REAL(KIND=8), DIMENSION(nb_ns) :: rhoi_wall_old, drhoi_wall

        REAL(KIND=8), DIMENSION(:), INTENT(IN)     :: u_phys1      !< conservative variables of the first physical cell
        REAL(KIND=8), INTENT(IN)                   :: p_wall       !< wall pressure 
        REAL(KIND=8), DIMENSION(:), INTENT(INOUT)  :: rhoi_wall    !< wall species density vector 


        ! Newton's loop exit residual
        smb_out = 1.d-4
        smb_it_out = 10

        ! Some definitions
        vol_p = volumes(3)
        ov_dr = 2.d0/(vol_p)

        ! Inizialization of some variables
        k = 0
        smb_max_old = 1.d10
        smb = 1.d10
        delta =  1.d-2

        ! Evaluation of the pyrolysis gas equilibrium composition at the surface
        ! Evaluated only here since it depends only on P_wall and T_wall
        ! (constant inside this procedure)
        yi_pyro=0.d0
        CALL library_compute_eq_composition_pyro(p_wall, Twall, yi_pyro)
 
newton:  DO

           k = k + 1

           ! Library that compute the total mass flow rate to be used in the evaluation of the wall normal velocity
           CALL library_ablation_compute_wall_mass_blowing_rate(rhoi_wall,Twall,mdot_wall_tot)

           ! Library that compute the source terms
           CALL procedure_ablation_compute_wall_source_terms(rhoi_wall,Twall) 

           ! Wall density to be used in the evaluation of the blowing velocity
           rhowall = 0.d0 
           DO i = 1,nb_ns 
                   rhowall =rhowall +  rhoi_wall(i)
           ENDDO

           ! Evaluation of the blowing velocity
           uwall = mdot_wall_tot * (1.d0+phi_pyro) / rhowall

           ! Get molar fraction of the wall interface
           CALL library_get_molar_fractions(rhoi_wall(1:nb_ns), xi_w)
           CALL library_comp_tol(xi_w)

           ! Get molar fraction of first physical cell
           CALL library_get_molar_fractions(u_phys1(1:nb_ns), xi_p)
           CALL library_comp_tol(xi_p)

           ! Diffusion driving forces
           diff_driv = 0.d0
           DO i = 1,nb_ns 
              diff_driv(i) = (xi_p(i) - xi_w(i)) * ov_dr
           ENDDO

           ! Temperature vector (needed for the call to the diffusive flux calculation)
           tempw = Twall

           ! Species mass diffusion flux
           CALL library_get_species_DiffFlux(p_wall, tempw(1), tempw(nb_temp), xi_w, diff_driv, Ji)

           ! Storage of previous Newton's step residual
           smb_old = smb

           ! Computation of residual of the SMB
           DO i = 1,nb_ns
              conv(i) = rhoi_wall(i) * uwall 
              pyro(i) =  mdot_wall_tot * phi_pyro * yi_pyro(i)
              ! Discern if the species has a source term or not. Removing the species density
              ! from the SMB of the species not involved in the surface processes (ablation
              ! and recombination) solved a lot of problem with the convergence of the
              ! Newton's loop 
              IF (omega_wall(i).EQ.0.d0) THEN
                 smb(i) = (conv(i) + Ji(i) -  pyro(i)) / rhoi_wall(i)
              ELSE
                 smb(i) = conv(i) + Ji(i) - (omega_wall(i) +  pyro(i))
              ENDIF
           ENDDO

           ! Convergence criterion for the Newton's loop (Exit condition 1)
           IF ((MAXVAL(ABS(smb)).LT.smb_out).OR.(k.GT.smb_it_out)) THEN

              ! Before leaving the loop check if species densities are negative. If so, load
              ! previous step values of species densities and blowing velocity.
              IF (MINVAL(rhoi_wall).LT.0.d0) THEN
                 rhoi_wall = rhoi_wall_old
                 uwall = u_wall_old
              ENDIF
           EXIT
           ENDIF

           ! Check if species densities are negative or Newton's loop is diverging. If so
           ! load previous step values of species densities, blowing velocities and residual,
           ! then repeat the Newton's iteration with a reduced perturbation
           IF ((MINVAL(rhoi_wall).LT.0.d0).OR.(MAXVAL(ABS(smb)).GT.smb_max_old)) THEN
              rhoi_wall = rhoi_wall_old
              uwall = u_wall_old
              smb = smb_old
              delta = delta * 0.5d0
           ! If the perturbation is reduced too much without having convergence, exit
           ! (Exit condition 2)
           IF (delta.LT.1.d-10) EXIT
           ENDIF

           ! Store old max value of residual
           smb_max_old = MAXVAL(ABS(smb))

           ! Store old value of species densities
           rhoi_wall_old = rhoi_wall

           ! Store old value of blowing velocity 
           u_wall_old = uwall


jacobian:   DO j = 1,nb_ns

              ! Wall species densities perturbation
              drhoi_wall(j) = delta*MAX(ABS(rhoi_wall(j)),1.d-20)*SIGN(1.d0,rhoi_wall(j))
              rhoi_wall(j) = rhoi_wall(j) + drhoi_wall(j)

              ! Library that compute the total mass flow rate to be used in the evaluation of the wall normal velocity
              CALL library_ablation_compute_wall_mass_blowing_rate(rhoi_wall,Twall,mdot_wall_tot)

              ! Library that compute the source terms
              CALL procedure_ablation_compute_wall_source_terms(rhoi_wall,Twall) 

              ! Wall density to be used in the evaluation of the blowing velocity
              rhowall = 0.d0 
              DO i = 1,nb_ns 
                      rhowall =rhowall +  rhoi_wall(i)
              ENDDO

              ! Evaluation of the blowing velocity
              uwall = mdot_wall_tot * (1.d0+phi_pyro) / rhowall

              ! Get molar fraction of the wall interface
              CALL library_get_molar_fractions(rhoi_wall(1:nb_ns), xi_w)
              CALL library_comp_tol(xi_w)

              ! Diffusion driving forces
              diff_driv = 0.d0
              DO i = 1,nb_ns 
                 diff_driv(i) = (xi_p(i) - xi_w(i))*ov_dr
              ENDDO

              ! Temperature vector (needed for the call to the diffusive flux calculation)
              tempw = Twall

              ! Species mass diffusion flux
              CALL library_get_species_DiffFlux(p_wall, tempw(1), tempw(nb_temp), xi_w, diff_driv, Ji)

              ! Computation of residual of the perturbed SMB
              DO i = 1,nb_ns
                 conv(i) = rhoi_wall(i) * uwall 
                 pyro(i) =  mdot_wall_tot * phi_pyro * yi_pyro(i)
                 ! Discern if the species has a source term or not. Removing the species density
                 ! from the SMB of the species not involved in the surface processes (ablation
                 ! and recombination) solved a lot of problem with the convergence of the
                 ! Newton's loop 
                 IF (omega_wall(i).EQ.0.d0) THEN
                    smb_pert(i) = (conv(i) + Ji(i) -  pyro(i)) / rhoi_wall(i)
                 ELSE
                    smb_pert(i) = conv(i) + Ji(i) - (omega_wall(i) +  pyro(i))
                 ENDIF
                    smb_jac(i,j) = (smb_pert(i) - smb(i)) / drhoi_wall(j)
              ENDDO

              ! Reset ith species density to starting guess value
              rhoi_wall(j) = rhoi_wall_old(j)

           ENDDO jacobian

           ! LU decomposition of the SMB residual Jacobian
           CALL lu_decmp(smb_jac,indx,d)

           ! Solve the system
           CALL lu_solver(smb_jac,indx,smb)

           ! Update the species densities
           DO i = 1, nb_ns
              rhoi_wall(i) = rhoi_wall_old(i) - smb(i)
           ENDDO

        ENDDO newton

        END SUBROUTINE ablation_SMB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SUBROUTINE ablation_compute_source_terms(rhoi_wall,T_wall)

        USE mod_general_data,          ONLY: mi
        USE mod_ablation_data,         ONLY: mdot_wall, omega_wall, omega_wall_rec,             &
                                           & C3_pos, CN_pos, CO_pos, N_pos, O_pos, O2_pos, C_pos


        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: rhoi_wall       !< wall species density vector 
        REAL(KIND=8), INTENT(IN)                :: T_wall          !< wall temperature



        CALL compute_recombination(rhoi_wall,T_wall)

        omega_wall = omega_wall_rec
        omega_wall(C3_pos) = omega_wall(C3_pos) +  mdot_wall(4)                                   !C3
        omega_wall(CN_pos) = omega_wall(CN_pos) +  mdot_wall(3)*mi(CN_pos)/mi(C_pos)              !CN
        omega_wall(CO_pos) = omega_wall(CO_pos) + (mdot_wall(1)+mdot_wall(2))*mi(CO_pos)/mi(C_pos)!CO
        omega_wall(N_pos)  = omega_wall(N_pos)  - (mdot_wall(3)*mi(N_pos)/mi(C_pos))              !N
        omega_wall(O_pos)  = omega_wall(O_pos)  -  mdot_wall(1)*mi(O_pos)/mi(C_pos)               !O
        omega_wall(O2_pos) = omega_wall(O2_pos) -  0.5d0*mdot_wall(2)*mi(O2_pos)/mi(C_pos)        !O2

        END SUBROUTINE ablation_compute_source_terms

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SUBROUTINE compute_only_Nitrogen_recombination(rhoi_w,T_w)

        USE mod_general_data,             ONLY: Ri
        USE mod_ablation_data,            ONLY: gamma_rec, N_pos, N2_pos, omega_wall_rec

        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: rhoi_w                       !< wall species density vector 
        REAL(KIND=8), INTENT(IN)  :: T_w                                        !< wall temperature

        !N + N  => N_2   Surface Nitrogen Recombination
        omega_wall_rec(N_pos)=-rhoi_w(N_pos)*dsqrt(Ri(N_pos)*T_w/(2.d0*acos(-1.d0)))*gamma_rec  
        omega_wall_rec(N2_pos)=-omega_wall_rec(N_pos)

        END SUBROUTINE compute_only_Nitrogen_recombination

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SUBROUTINE compute_Nitrogen_and_Ions_recombination(rhoi_w,T_w)

        USE mod_general_data,             ONLY: Ri, mi
        USE mod_ablation_data,            ONLY: omega_wall_rec, ion_species_matrix, e_minus_pos, &
                                              & gamma_rec, N_pos, N2_pos, omega_wall_rec

        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: rhoi_w                       !< wall species density vector 
        REAL(KIND=8), INTENT(IN)  :: T_w                                        !< wall temperature
        INTEGER :: i, imax

        omega_wall_rec = 0.d0

        !N + N  => N_2   Surface Nitrogen Recombination
        omega_wall_rec(N_pos)=-rhoi_w(N_pos)*dsqrt(Ri(N_pos)*T_w/(2.d0*acos(-1.d0)))*gamma_rec  
        omega_wall_rec(N2_pos)=-omega_wall_rec(N_pos)

        imax = SIZE(ion_species_matrix,1)
        DO i=1,imax 
           ! Compute the ion wall source term (negative since the ions are destroyed)
           omega_wall_rec(ion_species_matrix(i,1))=omega_wall_rec(ion_species_matrix(i,1)) - rhoi_w(ion_species_matrix(i,1))* &
                                                  & dsqrt(Ri(ion_species_matrix(i,1))*T_w/(2.d0*acos(-1.d0)))*1.d0
           ! Compute the species wall source term (positive since the species are created)
           omega_wall_rec(ion_species_matrix(i,2))=omega_wall_rec(ion_species_matrix(i,2)) -                                  &
                                                  & omega_wall_rec(ion_species_matrix(i,1))*                                  &
                                                  & mi(ion_species_matrix(i,2))/mi(ion_species_matrix(i,1))

           ! Compute the electron  wall source term (negative since the ions are destroyed)
           omega_wall_rec(e_minus_pos)=omega_wall_rec(e_minus_pos) + omega_wall_rec(ion_species_matrix(i,1))*     &
                                      &mi(e_minus_pos)/mi(ion_species_matrix(i,1))
        END DO

        END SUBROUTINE compute_Nitrogen_and_Ions_recombination

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SUBROUTINE ablation_compute_SEB(rhoi_wall, T_wall, T_flow, seb)

        USE mod_general_data,          ONLY: nb_ns, nb_temp, volumes
        USE mod_neq_function_pointer,  ONLY: library_get_species_enthalpies, library_get_thermal_cond
        USE mod_ablation_data,         ONLY: mdot_wall_tot, omega_wall, seb_stored, wall_emiss, phi_pyro
        USE mod_ablation_phys_param,   ONLY: SF_const

        INTEGER                 :: i
        REAL(KIND=8)            :: dr_wall                         !< Distance wall-first_cell_center 
        REAL(KIND=8)            :: cond_flux
        REAL(KIND=8)            :: mdot_times_hv
        REAL(KIND=8)            :: hi_times_omega_wall
        REAL(KIND=8)            :: rad_flux_out 
        REAL(KIND=8), DIMENSION(nb_temp)   :: Twall_vec_nb_temp
        REAL(KIND=8), DIMENSION(nb_temp*nb_ns)   :: Twall_vec_nb_temp_ns, h_species
        REAL(KIND=8), DIMENSION(nb_temp) :: lambda_wall_vec

        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: rhoi_wall       !< wall species density vector
        REAL(KIND=8), INTENT(IN)                :: T_wall          !< wall temperature
        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: T_flow          !< temperature(s) of the first physical cell 
        REAL(KIND=8), INTENT(OUT)               :: seb             !< wall temperature


        !1) Conductive heat flux 
        ! Compute thermal conductivity
        Twall_vec_nb_temp = T_wall
        CALL library_get_thermal_cond(rhoi_wall, Twall_vec_nb_temp, lambda_wall_vec) ! << (JBS) I think this is only considering heavy particle thermal cond.

        dr_wall  = volumes(3) * 0.5d0

        cond_flux = 0.d0
        DO i=1, nb_temp
           cond_flux = cond_flux + lambda_wall_vec(i)*(T_flow(i)- T_wall)/dr_wall
        ENDDO




        !2) Enthalpy of the material at the initial state (room temperature) times the total char mass flow rate 
        ! This is the only term of the steady-state conductive flux remaining in 
        ! the SEB after simplification   

        IF (phi_pyro.EQ.0.d0) THEN
           mdot_times_hv = mdot_wall_tot * 0.d0  !Enthalpy of graphite at room temperature is zero
        ELSE
           mdot_times_hv = mdot_wall_tot*(1.d0+phi_pyro) * (-911140.72d0)  !Enthalpy of FM5055 at 297.8K
        ENDIF

        !3) Sum over the species of species enthalpy times species wall source term 
        ! This is the term that substitutes the diffusive heat flux and the blowing
        ! term in the SEB by using the species SMB

        Twall_vec_nb_temp_ns = T_wall
        CALL  library_get_species_enthalpies(Twall_vec_nb_temp_ns, h_species) 

        hi_times_omega_wall = 0.d0
        DO i=1, nb_ns
            hi_times_omega_wall = hi_times_omega_wall+h_species(i)*omega_wall(i)
        ENDDO

        !4) Wall re-radiation
          
        rad_flux_out = SF_const * wall_emiss * T_wall * T_wall * T_wall * T_wall


        !5) Surface Energy Balance (SEB) evaluation 
        seb = cond_flux + mdot_times_hv - rad_flux_out - hi_times_omega_wall

        seb_stored = seb


        END SUBROUTINE ablation_compute_SEB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!! This subroutine computes the stress tensor and the conductive heat
!! flux (for nonequilibrium flows (sphere case - 1 temperature)
!! in case of ablative boundary condition. These values will be used directly as
!! to impose the fluxes at the wall interface, according to the SMB and the SEB. 
!! for nonequilibrium flows (sphere case - 1 temperature).

!NOTE:
! THIS SUBROUTINE IS NO LONGER USED WITH THE NEW VERSION OF THE ABLATIVE
! BOUNDARY CONDITION WORKS SETTING THE STATE OF THE INTERFACE AND THEN
! EXTRAPOLATE TO SET THE STATE OF THE GHOST CELLS
! HOWEVER THE SUBROUTINE IS KEPT HERE FOR POSSIBLE FUTURE USE

        SUBROUTINE ablation_wall_stress_tensor_and_conductive_heat_flux (r_l, r_r, vol_l, vol_r, left_data, right_data, u_left, &
                                                                        & u_right)

        USE mod_general_data,            ONLY: pos_u_cell, pos_v_cell, pos_T_cell,  & 
                                             & pos_mu_cell, pos_lambda_cell
        USE mod_ablation_data,           ONLY: tau_rr_wall, tau_rt_wall, cond_flux_wall 

        IMPLICIT NONE

        INTEGER :: i
        REAL(KIND=8) :: ov_dr
        REAL(KIND=8) :: r, u, v, lambda, mu, T
        REAL(KIND=8) :: dum   
        REAL(KIND=8) :: mul, mur, lambdal, lambdar, Tl, Tr, ul, ur, vl, vr
        REAL(KIND=8) :: du_dr, dv_dr, dT_dr

        REAL(KIND=8), INTENT(IN) :: r_l                       !< radial position of left state
        REAL(KIND=8), INTENT(IN) :: r_r                       !< radial position of rigth state
        REAL(KIND=8), INTENT(IN) :: vol_l                     !< volume of left state
        REAL(KIND=8), INTENT(IN) :: vol_r                     !< volume of right state
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_left      !< conservative variables of left state
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u_right     !< conservative variables of right state
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: left_data   !< physical properties of left state
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: right_data  !< physical properties of right state 

        ! Explicit interface for subroutine stress_tensor_1D_SL_sph
        INTERFACE
            SUBROUTINE stress_tensor_1D_SL_sph (mu, r, u, v, du_dr, dv_dr, tau_rr, tau_rt, dum)
              IMPLICIT NONE
                REAL(KIND=8), INTENT(IN) :: mu, r, u, v, du_dr, dv_dr
                REAL(KIND=8), INTENT(OUT) :: tau_rr, tau_rt, dum   
              END SUBROUTINE stress_tensor_1D_SL_sph
        END INTERFACE

        ! Cell interface position
        r = 0.5d0*(r_l + r_r)

        ! Temperature, velocity components, dynamic viscosity, thermal conductivity of left and right states

        ! Left state
        Tl      = left_data(pos_T_cell)
        ul      = left_data(pos_u_cell)
        vl      = left_data(pos_v_cell)
        mul     = left_data(pos_mu_cell)
        lambdal = left_data(pos_lambda_cell) 


        ! Right state
        Tr      = right_data(pos_T_cell)
        ur      = right_data(pos_u_cell)
        vr      = right_data(pos_v_cell)
        mur     = right_data(pos_mu_cell)
        lambdar = right_data(pos_lambda_cell)


        ! Temperature, velocity, dynamic viscosity and thermal conductivity at cell interface
        T   = 0.5d0*(Tl + Tr)
        u   = 0.5d0*(ul + ur)
        v   = 0.5d0*(vl + vr)
        mu  = 0.5d0*(mul + mur)
        lambda = 0.5d0*(lambdal + lambdar) 


        ! Velocity, temperature gradients
        ov_dr = 2.d0/(vol_l + vol_r)
        du_dr = ov_dr*(ur - ul)
        dv_dr = ov_dr*(vr - vl)
        dT_dr = ov_dr*(Tr - Tl)


        CALL stress_tensor_1D_SL_sph (mu, r, u, v, du_dr, dv_dr, tau_rr_wall, tau_rt_wall, dum)

        ! Wall conductive heat flux 
        cond_flux_wall = - lambda*dT_dr
         
        END SUBROUTINE ablation_wall_stress_tensor_and_conductive_heat_flux 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------------------------------------------------------------------!
!> This subroutine computes the wall diffusive flux to be usde in the ablation
!file summary
  SUBROUTINE ablation_comp_surface_flux (q_Diff,q_Fourier, q_Chem)

    USE mod_general_data,              ONLY: nb_cells, nb_dim, nb_prop, nb_eq, nb_ns, nb_int_temp,       & 
                                           & nb_temp, start_u_phys, start_prop_phys, pos_u_cell,         & 
                                           & pos_v_cell, pos_T_cell, pos_mu_cell, pos_lambda_cell,       & 
                                           & pos_pres_cell, pos_ei_cell, pos_em, Ri, volumes, xc,        & 
                                           & u, cell_prop, simulation_id
    USE mod_neq_function_pointer,      ONLY: library_get_molar_fractions, library_get_species_DiffFlux,  & 
                                           & library_comp_tol, library_get_enthalpy
    USE mod_function_pointer,          ONLY: get_stress_tensor_1D_SL
    USE mod_ablation_data,             ONLY: mdot_wall_tot, abla_dir

    IMPLICIT NONE

    INTEGER, PARAMETER :: in1 = 11
    INTEGER :: i, j, k 
    INTEGER :: n, m, l, ios, mmin, mmax
    REAL(KIND=8) :: mul, mull, tmp, sum_Ji
    REAL(KIND=8) :: ul, ull, vl, vll, Tl, Tll, Tel, Tell, pl, pll
    REAL(KIND=8) :: vol_l, vol_ll
    REAL(KIND=8) :: Tw, muw, pw, rhow
    REAL(KIND=8) :: ov_dr, du_dr, dv_dr, dT_dr
    REAL(kind=8) :: h_graph                               !< graphite enthalpy from file
    REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: T_ar, h_ar
    REAL(KIND=8), DIMENSION(nb_prop) :: prop_l, prop_ll
    REAL(KIND=8), DIMENSION(nb_temp) :: templl, templ, grad_T, lambda_vecl, lambda_vecll
    REAL(KIND=8), DIMENSION(nb_temp) :: lambda_vecw, tempw
    REAL(KIND=8), DIMENSION(nb_temp) :: hw
    REAL(KIND=8), DIMENSION(nb_ns) :: xil, xill, rhoil, rhoill, hil, hill, Vi
    REAL(KIND=8), DIMENSION(nb_ns) :: hiw, rhoiw, xiw 
    REAL(KIND=8), DIMENSION(nb_ns*nb_dim) :: diff_driv, Ji

    REAL(KIND=8), INTENT(OUT)                ::  q_Diff, q_Chem   !< Diffusive and chemical heat fluxes
    REAL(KIND=8), DIMENSION(:), INTENT(OUT)  ::  q_Fourier        !< Conductive heat fluxes


    ! Initialization of left and central state
    ! Species densities and molar fractions
    DO i = 1,nb_ns 
       rhoill(i) = u(2*nb_eq + i)
       rhoil(i) = u(3*nb_eq + i)
    ENDDO

    CALL library_get_molar_fractions (rhoill, xill)
    CALL library_comp_tol(xill) 
   
    CALL library_get_molar_fractions (rhoil, xil)
    CALL library_comp_tol(xil) 

    DO i = 1,nb_prop
       prop_ll(i) = cell_prop(nb_prop + i)
       prop_l(i) = cell_prop(2*nb_prop + i)
    ENDDO

    ! Velocity components and temperatures
    ull = prop_ll(pos_u_cell)
    vll = prop_ll(pos_v_cell)
    DO i = 1,nb_temp
       templl(i) = prop_ll(pos_T_cell + i - 1)
    ENDDO
    vol_ll = volumes(3)   
 
    ul = prop_l(pos_u_cell)
    vl = prop_l(pos_v_cell)
    DO i = 1,nb_temp
       templ(i) = prop_l(pos_T_cell + i - 1)
    ENDDO
    vol_l = volumes(4)


    ! Second phisical cell values

    ! Species specific enthalpies
    Tl  = prop_l(pos_T_cell)
    Tel = prop_l(pos_T_cell + nb_temp - 1)

    DO j = 1,pos_em - 1
       hil(j) = prop_l(pos_ei_cell + j - 1) + Ri(j)*Tl
    ENDDO 

    hil(pos_em) = prop_l(pos_ei_cell + pos_em - 1) + Ri(pos_em)*Tel

    DO j = pos_em + 1,nb_ns
       hil(j) = prop_l(pos_ei_cell + j - 1) + Ri(j)*Tl
    ENDDO 

    ! Dynamic viscosity and thermal conductivity components
    mul = prop_l(pos_mu_cell)
    DO j = 1,nb_temp
       lambda_vecl(j) = prop_l(pos_lambda_cell + j - 1)
    ENDDO

    ! Pressure
    pl = prop_l(pos_pres_cell) 



    ! First phisical cell values

    ! Species specific enthalpies
    Tll  = prop_ll(pos_T_cell)
    Tell = prop_ll(pos_T_cell + nb_temp - 1)

    DO j = 1,pos_em - 1
       hill(j) = prop_ll(pos_ei_cell + j - 1) + Ri(j)*Tll
    ENDDO 

    hill(pos_em) = prop_ll(pos_ei_cell + pos_em - 1) + Ri(pos_em)*Tell

    DO j = pos_em + 1,nb_ns
       hill(j) = prop_ll(pos_ei_cell + j - 1) + Ri(j)*Tll
    ENDDO 

    ! Dynamic viscosity and thermal conductivity components
    mull = prop_ll(pos_mu_cell)
    DO j = 1,nb_temp
       lambda_vecll(j) = prop_ll(pos_lambda_cell + j - 1)
    ENDDO

    ! Pressure
    pll = prop_ll(pos_pres_cell) 


    !Gradients       

    ! Velocity and temperature gradients
    ov_dr = 2.d0/(vol_ll + vol_l)
    du_dr = ov_dr*(ul - ull)
    dv_dr = ov_dr*(vl - vll)
    DO j = 1,nb_temp
       grad_T(j) = (templ(j) - templl(j))*ov_dr
    ENDDO

    ! Diffusion driving forces
    diff_driv = 0.d0
    DO j = 1,nb_ns 
       diff_driv(j) = (xil(j) - xill(j))*ov_dr
    ENDDO


    ! Extrapolating to the surface 
    DO i = 1,nb_temp
       tempw(i) = templl(i) + (templl(i) - templ(i)) * ov_dr * vol_ll / 2.d0
    ENDDO
    DO j = 1,pos_em - 1
       hiw(j) =  hill(j) + (hill(j) - hil(j)) * ov_dr * vol_ll / 2.d0
    ENDDO 

    hiw(pos_em) =  hill(pos_em) + (hill(pos_em) - hil(pos_em)) * ov_dr * vol_ll / 2.d0

    DO j = pos_em + 1,nb_ns
       hiw(j) = hill(j) + (hill(j) - hil(j)) * ov_dr * vol_ll / 2.d0
    ENDDO 

    ! Dynamic viscosity and thermal conductivity components
    muw = mull + (mull - mul) * ov_dr * vol_ll / 2.d0
    DO j = 1,nb_temp
       lambda_vecw(j) = lambda_vecll(j) + (lambda_vecll(j) - lambda_vecl(j)) * ov_dr * vol_ll / 2.d0
    ENDDO

    ! Pressure
    pw = pll + (pll - pl) * ov_dr * vol_ll / 2.d0 

    ! Density and mole fractions
    DO i = 1,nb_ns 
        rhoiw(i) = rhoill(i) + (rhoill(i) - rhoil(i)) * ov_dr * vol_ll / 2.d0
        xiw(i) = xill(i) + (xill(i) - xil(i)) * ov_dr * vol_ll / 2.d0
    ENDDO
 

    ! Species mass diffusion flux
    CALL library_get_species_DiffFlux(pw, tempw(1), tempw(nb_temp), xiw, diff_driv, Ji)

    ! Species diffusion velocities
    sum_Ji = 0.d0
    DO j = 1,nb_ns
       tmp = Ji(j)
       sum_Ji = sum_Ji + tmp 
       Vi(j)  = tmp/rhoiw(j)
    ENDDO

    ! Fourier heat flux (component associated to translational energy)  
    q_Fourier = - lambda_vecw(1)*grad_T(1)

    ! Fourier heat flux (component associated to internal energy)
    IF (nb_temp.GT.1) THEN
       q_Fourier(2) = 0.d0
       DO k = 1,nb_temp-1
          q_Fourier(2) = q_Fourier(2) - lambda_vecw(k + 1)*grad_T(k + 1)
       ENDDO
    ELSE
    ENDIF

    ! Diffusive heat flux
    q_Diff = 0.d0
    DO j = 1,nb_ns
       q_Diff = q_Diff + Ji(j)*hiw(j)
    ENDDO

    ! Open the file with the data of the graphite enthalpy versus temperature
    OPEN(UNIT=in1,FILE=abla_dir//"/graphite_enthalpy.dat",STATUS='old',IOSTAT=ios)


    IF (ios.EQ.0) THEN 

       READ(in1,*) n, m, l 

       k = (m-n)/l+1

       ALLOCATE(h_ar(k))
       ALLOCATE(T_ar(k))

       DO i = 1,k
          READ(in1,*) T_ar(i), h_ar(i)
       ENDDO

       Tw = tempw(1)
       mmin = INT(Tw)-n+1
       mmax = mmin + 1

       h_graph = (h_ar(mmax) - h_ar(mmin)) / (T_ar(mmax) - T_ar(mmin)) * (Tw - T_ar(mmin)) + h_ar(mmin)

       rhow = 0.d0 
       DO j = 1,nb_ns
          rhow = rhow + rhoiw(i) 
       ENDDO

       CALL library_get_enthalpy (rhow, rhoiw, tempw, hw)

       q_Chem = mdot_wall_tot * (hw(1) - h_graph)

    ELSE

      print*, 'ATTENTION!'
      print*, 'Data file graphite_enthalpy.dat not found!' 
      print*, 'Computation is fine but the blowing heat flux is set to zero in the ablation summary file.'

      q_Chem = 0.d0

    ENDIF

    

  END SUBROUTINE ablation_comp_surface_flux
!------------------------------------------------------------------------------!

#endif

  END MODULE mod_ablation_procedures
