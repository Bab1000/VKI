cmake_minimum_required(VERSION 2.6)

add_sources(stagline
    set_bc_1D.F90
    set_bc_1D_SL.F90
    set_bc_neq_1D.F90
    set_bc_neq_1D_SL.F90
    set_flux_1D.F90
    set_flux_1D_SL.F90
    set_flux_neq_1D.F90
    set_flux_neq_1D_SL.F90
    set_library_pointer.F90
    set_limiter.F90
    set_procedure_pointer.F90
    set_source_term.F90
    set_source_term_neq.F90
    set_source_term_neq_SL.F90
    set_source_term_SL.F90
    set_stress_tensor_1D_SL.F90
    set_time_integration_1D.F90
    set_time_integration_1D_SL.F90
    set_time_step_1D.F90
    set_time_step_1D_SL.F90
    set_transp_coeff.F90
)
