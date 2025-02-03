cmake_minimum_required(VERSION 2.6)

add_sources(stagline
    diff_flux_num_Jac_1D_SL.F90
    ns_flux_1D_SL_sph.F90
    null_ns_flux_1D_SL_Jac.F90
    ns_flux_1D_SL_cyl.F90
    ns_flux_1D_SL_sph_Jac.F90
    stress_tensor_1D_SL_cyl.F90
    ns_flux_1D_SL_cyl_Jac.F90
    null_ns_flux_1D_SL.F90
    stress_tensor_1D_SL_sph.F90
)
