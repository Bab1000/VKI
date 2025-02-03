cmake_minimum_required(VERSION 2.6)

add_sources(stagline
    post_process_solution.F90
    write_sol_1D.F90
    write_sol_1D_SL.F90
    write_sol_2D.F90
    write_sol_neq_1D.F90
    write_sol_neq_1D_SL.F90
    write_sol_neq_2D.F90
    write_solution.F90
    write_surface_data_1D_SL.F90
    write_transp_flux_1D.F90
    write_transp_flux_1D_SL.F90
    write_transp_flux_2D.F90
    write_transp_flux_neq_1D.F90
    write_transp_flux_neq_1D_SL.F90
    write_transp_flux_neq_2D.F90
    write_transport_fluxes.F90
)
