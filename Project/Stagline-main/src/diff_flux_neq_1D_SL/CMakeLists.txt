cmake_minimum_required(VERSION 2.6)

add_sources(stagline
    ns_flux_neq1T_1D_SL_cyl.F90
    ns_flux_neqNT_1D_SL_cyl.F90
    ns_flux_neqNT_Te_1D_SL_cyl.F90
    ns_flux_neq1T_1D_SL_cyl_Jac.F90
    ns_flux_neqNT_1D_SL_cyl_Jac.F90
    ns_flux_neqNT_Te_1D_SL_cyl_Jac.F90
    ns_flux_neq1T_1D_SL_sph.F90
    ns_flux_neqNT_1D_SL_sph.F90
    ns_flux_neqNT_Te_1D_SL_sph.F90
    ns_flux_neq1T_1D_SL_sph_Jac.F90
    ns_flux_neqNT_1D_SL_sph_Jac.F90
    ns_flux_neqNT_Te_1D_SL_sph_Jac.F90
)
