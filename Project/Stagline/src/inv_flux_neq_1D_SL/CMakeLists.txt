cmake_minimum_required(VERSION 2.6)

add_sources(stagline
    inviscid_flux_Jac_neq1T_1D_SL.F90
    Pos_Neg_split_Jac_neqNT_Te_1D_SL.F90
    Steger_Warming_neq1T_1D_SL.F90
    inviscid_flux_Jac_neqNT_1D_SL.F90
    roe_avgState_neq1T_1D_SL.F90
    Steger_Warming_neqNT_1D_SL.F90
    inviscid_flux_Jac_neqNT_Te_1D_SL.F90
    roe_avgState_neqNT_1D_SL.F90
    Steger_Warming_neqNT_Te_1D_SL.F90
    Pos_Neg_split_Jac_neq1T_1D_SL.F90
    roe_avgState_neqNT_Te_1D_SL.F90
    van_leer_neq_1D_SL.F90
    Pos_Neg_split_Jac_neqNT_1D_SL.F90
    roe_neq_1D_SL.F90
    Yoon_Jameson_Jac_neq_1D_SL.F90
    ausm_neq_1D_SL.F90
    pressure_weight_function_1D_SL.F90
    modified_SW_neqNT_Te_1D_SL.F90
    Steger_Warming_1D_SL.F90
    ausmP_up2_neq_1D_SL.F90
    Pos_Neg_split_Jac_1D_SL.F90
)
