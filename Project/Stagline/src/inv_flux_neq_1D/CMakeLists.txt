cmake_minimum_required(VERSION 2.6)

add_sources(stagline
    Pos_Neg_split_Jac_neq1T_1D.F90
    Pos_Neg_split_Jac_neqNT_1D.F90
    Pos_Neg_split_Jac_neqNT_Te_1D.F90
    Steger_Warming_neq1T_1D.F90
    Steger_Warming_neqNT_1D.F90
    Steger_Warming_neqNT_Te_1D.F90
    Yoon_Jameson_Jac_neq_1D.F90
    hlle_neq_1D.F90
    inviscid_flux_Jac_neq1T_1D.F90
    inviscid_flux_Jac_neqNT_1D.F90
    inviscid_flux_Jac_neqNT_Te_1D.F90
    roe_avgState_neq1T_1D.F90
    roe_avgState_neqNT_1D.F90
    roe_avgState_neqNT_Te_1D.F90
    roe_neq_1D.F90
    van_leer_neq_1D.F90
)
