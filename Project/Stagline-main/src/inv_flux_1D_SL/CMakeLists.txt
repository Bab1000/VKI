cmake_minimum_required(VERSION 2.6)

add_sources(stagline
    ausmPC_1D_SL.F90
    Pos_Neg_split_Jac_1D_SL.F90
    Steger_Warming_1D_SL.F90
    conv_flux_num_Jac_1D_SL.F90
    roe_1D_SL.F90
    van_leer_1D_SL.F90
    inviscid_flux_Jac_1D_SL.F90
    roe_avgState_1D_SL.F90
    Yoon_Jameson_Jac_1D_SL.F90
)
