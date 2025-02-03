cmake_minimum_required(VERSION 2.6)

add_sources(stagline
    null_diff_source_term_SL.F90
    null_diff_source_term_SL_Jac.F90
    source_term_diff_1D_SL_cyl.F90
    source_term_diff_1D_SL_cyl_Jac.F90
    source_term_diff_1D_SL_sph.F90
    source_term_diff_1D_SL_sph_Jac.F90
    source_term_diff_num_Jac_SL.F90
    source_term_inv_1D_SL_cyl.F90
    source_term_inv_1D_SL_cyl_Jac.F90
    source_term_inv_1D_SL_sph.F90
    source_term_inv_1D_SL_sph_Jac.F90
    source_term_inv_num_Jac_SL.F90
    source_term_num_Jac.F90
    source_term_quasi1D.F90
)
