!------------------------------------------------------------------------------!
!> This subroutine deallocates common arrays and nullifies pointers to functions/subroutines.
  SUBROUTINE close_solver ()

#include"../config.h"

    USE mod_function_pointer
    USE mod_neq_function_pointer
    USE mod_algebra
    USE mod_numerics_data
    USE mod_domain_boundary
    USE mod_general_data,             ONLY: u, u_old, u_vec, du, ru, delta, geom_source,       & 
                                          & cell_prop, volumes, x, y, xc, yc, Ri,              & 
                                          & rhoi, epsi, sigmai, yi, xi, ei, hi, eiT, eistar,   & 
                                          & temp, rho_eint, beta, betak, ov_betak, eintk, Di,  &
                                          & Dij_mat, chi, lambda_vec, Ji, diff_driv, p_in,     & 
                                          & library_name, iden, s_vec, js_mat, vel,            & 
                                          & log_cfl_ite, log_cfl_number,                       &
                                          & library_name, iden, s_vec, js_mat, vel,            &
                                          & ablation_library_name
#ifdef CARBONABLA
    USE mod_ablation_data,            ONLY: omega_wall
    USE mod_ablation_pointers,        ONLY: procedure_ablation_finalize 
#endif
    USE mod_radiation,                ONLY: rad, rad_old,rad_new
    USE mod_general_ablation,         ONLY: mdot_i
    USE mod_surface_properties,       ONLY: rho_i_surface, T_surface

    IMPLICIT NONE

    WRITE(*,5)'solver_fmvcc_F90:: De-allocating vectors and nullifying pointers'
    PRINT*
    
    ! General data
    IF (ALLOCATED(u))                       DEALLOCATE(u)
    IF (ALLOCATED(u_old))                   DEALLOCATE(u_old)
    IF (ALLOCATED(du))                      DEALLOCATE(du)
    IF (ALLOCATED(ru))                      DEALLOCATE(ru)
    IF (ALLOCATED(rad))                     DEALLOCATE(rad)
    IF (ALLOCATED(rad_old))                 DEALLOCATE(rad_old)
    IF (ALLOCATED(rad_new))                 DEALLOCATE(rad_new)
    IF (ALLOCATED(geom_source))             DEALLOCATE(geom_source)
    IF (ALLOCATED(cell_prop))               DEALLOCATE(cell_prop)
    IF (ALLOCATED(volumes))                 DEALLOCATE(volumes)
    IF (ALLOCATED(x))                       DEALLOCATE(x)
    IF (ALLOCATED(y))                       DEALLOCATE(y)
    IF (ALLOCATED(xc))                      DEALLOCATE(xc)
    IF (ALLOCATED(yc))                      DEALLOCATE(yc)
    IF (ALLOCATED(boundary_data))           DEALLOCATE(boundary_data)   
    IF (ALLOCATED(p_in))                    DEALLOCATE(p_in)
    IF (ALLOCATED(Ri))                      DEALLOCATE(Ri) 
    IF (ALLOCATED(rhoi))                    DEALLOCATE(rhoi)
    IF (ALLOCATED(ei))                      DEALLOCATE(ei)
    IF (ALLOCATED(hi))                      DEALLOCATE(hi)
    IF (ALLOCATED(eiT))                     DEALLOCATE(eiT)
    IF (ALLOCATED(eistar))                  DEALLOCATE(eistar)
    IF (ALLOCATED(yi))                      DEALLOCATE(yi)
    IF (ALLOCATED(xi))                      DEALLOCATE(xi)
    IF (ALLOCATED(Di))                      DEALLOCATE(Di)
    IF (ALLOCATED(Dij_mat))                 DEALLOCATE(Dij_mat)
    IF (ALLOCATED(diff_driv))               DEALLOCATE(diff_driv)
    IF (ALLOCATED(chi))                     DEALLOCATE(chi)
    IF (ALLOCATED(Ji))                      DEALLOCATE(Ji)
    IF (ALLOCATED(lambda_vec))              DEALLOCATE(lambda_vec)
    IF (ALLOCATED(epsi))                    DEALLOCATE(epsi)
    IF (ALLOCATED(sigmai))                  DEALLOCATE(sigmai)
    IF (ALLOCATED(eintk))                   DEALLOCATE(eintk)
    IF (ALLOCATED(betak))                   DEALLOCATE(betak)
    IF (ALLOCATED(ov_betak))                DEALLOCATE(ov_betak)
    IF (ALLOCATED(temp))                    DEALLOCATE(temp)
    IF (ALLOCATED(rho_eint))                DEALLOCATE(rho_eint)
    IF (ALLOCATED(beta))                    DEALLOCATE(beta)
    IF (ALLOCATED(iden))                    DEALLOCATE(iden)
    IF (ALLOCATED(s_vec))                   DEALLOCATE(s_vec)
    IF (ALLOCATED(js_mat))                  DEALLOCATE(js_mat)
    IF (ALLOCATED(vel))                     DEALLOCATE(vel)
    IF (ALLOCATED(u_vec))                   DEALLOCATE(u_vec)
    IF (ALLOCATED(delta))                   DEALLOCATE(delta) 
    IF (ALLOCATED(log_cfl_ite))             DEALLOCATE(log_cfl_ite)
    IF (ALLOCATED(log_cfl_number))          DEALLOCATE(log_cfl_number)
 
 
 
    ! Algebraic data
    IF (ALLOCATED(Ai_mat))                  DEALLOCATE(Ai_mat)
    IF (ALLOCATED(Bi_mat))                  DEALLOCATE(Bi_mat)
    IF (ALLOCATED(Ci_mat))                  DEALLOCATE(Ci_mat)
    IF (ALLOCATED(Yi_vec))                  DEALLOCATE(Yi_vec)
    IF (ALLOCATED(dum_vec1))                DEALLOCATE(dum_vec1)
    IF (ALLOCATED(dum_vec2))                DEALLOCATE(dum_vec2)
    IF (ALLOCATED(dum_mat))                 DEALLOCATE(dum_mat)
    IF (ALLOCATED(alpha_i))                 DEALLOCATE(alpha_i)
    IF (ALLOCATED(invalpha_i))              DEALLOCATE(invalpha_i)

    ! Ablation data
    IF (ALLOCATED(mdot_i))                  DEALLOCATE(mdot_i)

    ! ! Surface properties
    IF (ALLOCATED(rho_i_surface))           DEALLOCATE(rho_i_surface)
    IF (ALLOCATED(T_surface))               DEALLOCATE(T_surface)

    ! Goto BLAS data 
#ifdef GOTO_BLAS 
    IF (ALLOCATED(ipiv_GOTO_BLAS))          DEALLOCATE(ipiv_GOTO_BLAS)
    IF (ALLOCATED(work_GOTO_BLAS))          DEALLOCATE(work_GOTO_BLAS)
#endif

    ! Domain boundary data
    IF (ALLOCATED(yi_bound))                DEALLOCATE(yi_bound)
    IF (ALLOCATED(rhoi_bound))              DEALLOCATE(rhoi_bound)
    IF (ALLOCATED(Tr_bound))                DEALLOCATE(Tr_bound)
    IF (ALLOCATED(Tv_bound))                DEALLOCATE(Tv_bound)
    IF (ALLOCATED(inlet))                   DEALLOCATE(inlet)
    IF (ALLOCATED(inlet_temp))              DEALLOCATE(inlet_temp)
    IF (ALLOCATED(bc))                      DEALLOCATE(bc)
    IF (ALLOCATED(name_bc))                 DEALLOCATE(name_bc)
    IF (ALLOCATED(boundary_data))           DEALLOCATE(boundary_data)
     
    ! Numerics data
    IF (ALLOCATED(stage_coeff))             DEALLOCATE(stage_coeff)
    IF (ALLOCATED(source_terms))            DEALLOCATE(source_terms)
    IF (ALLOCATED(fc))                      DEALLOCATE(fc)
    IF (ALLOCATED(fd))                      DEALLOCATE(fd)
    IF (ALLOCATED(u_left))                  DEALLOCATE(u_left)
    IF (ALLOCATED(u_right))                 DEALLOCATE(u_right)
    IF (ALLOCATED(prim_left))               DEALLOCATE(prim_left)
    IF (ALLOCATED(prim_right))              DEALLOCATE(prim_right)
    IF (ALLOCATED(ull))                     DEALLOCATE(ull)
    IF (ALLOCATED(ul))                      DEALLOCATE(ul)
    IF (ALLOCATED(ur))                      DEALLOCATE(ur)
    IF (ALLOCATED(urr))                     DEALLOCATE(urr)
    IF (ALLOCATED(cons_Roe))                DEALLOCATE(cons_Roe)
    IF (ALLOCATED(phys_data_Roe))           DEALLOCATE(phys_data_Roe)
    IF (ALLOCATED(rhoil))                   DEALLOCATE(rhoil)
    IF (ALLOCATED(rhoir))                   DEALLOCATE(rhoir)
    IF (ALLOCATED(xil))                     DEALLOCATE(xil)
    IF (ALLOCATED(xir))                     DEALLOCATE(xir)
    IF (ALLOCATED(s))                       DEALLOCATE(s)
    IF (ALLOCATED(aux))                     DEALLOCATE(aux)
    IF (ALLOCATED(js))                      DEALLOCATE(js) 
    IF (ALLOCATED(jsl))                     DEALLOCATE(jsl) 
    IF (ALLOCATED(jsr))                     DEALLOCATE(jsr) 
    IF (ALLOCATED(inv_js))                  DEALLOCATE(inv_js)  
    IF (ALLOCATED(left_bound))              DEALLOCATE(left_bound)
    IF (ALLOCATED(right_bound))             DEALLOCATE(right_bound)
    IF (ALLOCATED(jfl))                     DEALLOCATE(jfl)
    IF (ALLOCATED(jfr))                     DEALLOCATE(jfr)
    IF (ALLOCATED(jfdl))                    DEALLOCATE(jfdl)
    IF (ALLOCATED(jfdr))                    DEALLOCATE(jfdr)
    IF (ALLOCATED(Bd))                      DEALLOCATE(Bd)
    IF (ALLOCATED(j_ghost))                 DEALLOCATE(j_ghost)
    IF (ALLOCATED(js_line))                 DEALLOCATE(js_line)
    IF (ALLOCATED(central))                 DEALLOCATE(central)
    IF (ALLOCATED(deltau))                  DEALLOCATE(deltau)
    IF (ALLOCATED(lambda))                  DEALLOCATE(lambda)
    IF (ALLOCATED(diss))                    DEALLOCATE(diss)
    IF (ALLOCATED(right_eig))               DEALLOCATE(right_eig)
    IF (ALLOCATED(left_eig))                DEALLOCATE(left_eig)
    IF (ALLOCATED(absA_eig))                DEALLOCATE(absA_eig)  

    IF (ablation_library_name .NE.'none') THEN
       CALL procedure_ablation_finalize () 
    ELSE
    ENDIF
   
    ! General purpose pointers
    IF (ASSOCIATED(get_prim_from_cons))            NULLIFY(get_prim_from_cons)
    IF (ASSOCIATED(get_cons_from_prim))            NULLIFY(get_cons_from_prim)
    IF (ASSOCIATED(get_cons_phys_from_prim))       NULLIFY(get_cons_phys_from_prim)
    IF (ASSOCIATED(evolve_solution))               NULLIFY(evolve_solution)
    IF (ASSOCIATED(conv_flux_1D))                  NULLIFY(conv_flux_1D)
    IF (ASSOCIATED(conv_flux_1D_SL))               NULLIFY(conv_flux_1D_SL)
    IF (ASSOCIATED(eigsys_neq_1D))                 NULLIFY(eigsys_neq_1D)
    IF (ASSOCIATED(roe_avg_neq_1D))                NULLIFY(roe_avg_neq_1D)
    IF (ASSOCIATED(inv_flux_Jac_neq_1D))           NULLIFY(inv_flux_Jac_neq_1D)
    IF (ASSOCIATED(conv_flux_Jac_1D))              NULLIFY(conv_flux_Jac_1D)
    IF (ASSOCIATED(diff_flux_1D))                  NULLIFY(diff_flux_1D)
    IF (ASSOCIATED(diff_flux_1D_Jac))              NULLIFY(diff_flux_1D_Jac) 
    IF (ASSOCIATED(rec_1D))                        NULLIFY(rec_1D)
    IF (ASSOCIATED(rec_1D_SL))                     NULLIFY(rec_1D_SL)
    IF (ASSOCIATED(limiter))                       NULLIFY(limiter)
    IF (ASSOCIATED(get_phys_from_cons))            NULLIFY(get_phys_from_cons)
    IF (ASSOCIATED(get_transpCoeff))               NULLIFY(get_transpCoeff)
    IF (ASSOCIATED(compute_time_step_1D))          NULLIFY(compute_time_step_1D)
    IF (ASSOCIATED(get_inv_spectral_radius_1D))    NULLIFY(get_inv_spectral_radius_1D)  
    IF (ASSOCIATED(get_visc_spectral_radius_1D))   NULLIFY(get_visc_spectral_radius_1D) 
    IF (ALLOCATED(get_source_term))                DEALLOCATE(get_source_term)
    IF (ALLOCATED(get_source_term_Jac))            DEALLOCATE(get_source_term_Jac)
    IF (ALLOCATED(get_inv_source_term_SL))         DEALLOCATE(get_inv_source_term_SL)
    IF (ALLOCATED(get_inv_source_term_Jac_SL))     DEALLOCATE(get_inv_source_term_Jac_SL)
    IF (ASSOCIATED(get_diff_source_term_SL))       NULLIFY(get_diff_source_term_SL)
    IF (ALLOCATED(get_ghost_state_Expl_1D_1st))    DEALLOCATE(get_ghost_state_Expl_1D_1st)
    IF (ALLOCATED(get_ghost_state_Expl_1D_2nd))    DEALLOCATE(get_ghost_state_Expl_1D_2nd)
    IF (ALLOCATED(get_ghost_state_Impl_1D_1st))    DEALLOCATE(get_ghost_state_Impl_1D_1st)
    IF (ALLOCATED(get_ghost_state_Impl_1D_2nd))    DEALLOCATE(get_ghost_state_Impl_1D_2nd) 
    IF (ALLOCATED(get_ghost_state_Expl_1D_SL_1st)) DEALLOCATE(get_ghost_state_Expl_1D_SL_1st)
    IF (ALLOCATED(get_ghost_state_Expl_1D_SL_2nd)) DEALLOCATE(get_ghost_state_Expl_1D_SL_2nd)
    IF (ALLOCATED(get_ghost_state_Impl_1D_SL_1st)) DEALLOCATE(get_ghost_state_Impl_1D_SL_1st)
    IF (ALLOCATED(get_ghost_state_Impl_1D_SL_2nd)) DEALLOCATE(get_ghost_state_Impl_1D_SL_2nd) 
   
    ! Thermodynamic library pointers
    IF (library_name.NE.'none') THEN
      
       ! De-allocating thermodynamic library data
       CALL library_finalize ()

       ! Nullify thermodynamic library function/subroutine pointers
       IF (ASSOCIATED(library_initialize))                              NULLIFY(library_initialize) 
       IF (ASSOCIATED(library_finalize))                                NULLIFY(library_finalize)
       IF (ASSOCIATED(library_get_Ri))                                  NULLIFY(library_get_Ri)
       IF (ASSOCIATED(library_get_el_data))                             NULLIFY(library_get_el_data) 
       IF (ASSOCIATED(library_get_density))                             NULLIFY(library_get_density)
       IF (ASSOCIATED(library_get_pressure))                            NULLIFY(library_get_pressure)
       IF (ASSOCIATED(library_get_mass_fractions))                      NULLIFY(library_get_mass_fractions)
       IF (ASSOCIATED(library_get_molar_fractions))                     NULLIFY(library_get_molar_fractions) 
       IF (ASSOCIATED(library_get_mass_fractions_from_molar_fractions)) NULLIFY(library_get_mass_fractions_from_molar_fractions)
       IF (ASSOCIATED(library_comp_tol))                                NULLIFY(library_comp_tol)
       IF (ASSOCIATED(library_get_source))                              NULLIFY(library_get_source) 
       IF (ASSOCIATED(library_get_source_Jac))                          NULLIFY(library_get_source_Jac) 
       IF (ASSOCIATED(library_get_temperatures))                        NULLIFY(library_get_temperatures) 
       IF (ASSOCIATED(library_get_energy_densities))                    NULLIFY(library_get_energy_densities) 
       IF (ASSOCIATED(library_get_thermodynamic_data))                  NULLIFY(library_get_thermodynamic_data)
       IF (ASSOCIATED(library_get_data))                                NULLIFY(library_get_data)
       IF (ASSOCIATED(library_get_transpCoeff))                         NULLIFY(library_get_transpCoeff) 
       IF (ASSOCIATED(library_get_species_DiffFlux))                    NULLIFY(library_get_species_DiffFlux) 
       IF (ASSOCIATED(library_compute_eq_composition))                  NULLIFY(library_compute_eq_composition)     
       IF (ASSOCIATED(library_write_fvmcc_solution))                    NULLIFY(library_write_fvmcc_solution)

    ENDIF

5 FORMAT(A)

  END SUBROUTINE close_solver
!------------------------------------------------------------------------------!
