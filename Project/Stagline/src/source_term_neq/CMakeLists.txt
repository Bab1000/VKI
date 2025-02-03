cmake_minimum_required(VERSION 2.6)

add_sources(stagline
    source_term_diff_neq1T_1D_SL_cyl.F90
    source_term_diff_neq1T_1D_SL_cyl_Jac.F90
    source_term_diff_neq1T_1D_SL_sph.F90
    source_term_diff_neq1T_1D_SL_sph_Jac.F90
    source_term_diff_neqNT_1D_SL_cyl.F90
    source_term_diff_neqNT_1D_SL_cyl_Jac.F90
    source_term_diff_neqNT_1D_SL_sph.F90
    source_term_diff_neqNT_1D_SL_sph_Jac.F90
    source_term_diff_neqNT_Te_1D_SL_cyl.F90
    source_term_diff_neqNT_Te_1D_SL_cyl_Jac.F90
    source_term_diff_neqNT_Te_1D_SL_sph.F90
    source_term_diff_neqNT_Te_1D_SL_sph_Jac.F90
    source_term_inv_neq_1D_SL_cyl.F90
    source_term_inv_neq_1D_SL_cyl_Jac.F90
    source_term_inv_neq_1D_SL_sph.F90
    source_term_inv_neq_1D_SL_sph_Jac.F90
    source_term_neq_kinetics.F90
    source_term_neq_kinetics_SL.F90
    source_term_neq_quasi1D.F90
)
