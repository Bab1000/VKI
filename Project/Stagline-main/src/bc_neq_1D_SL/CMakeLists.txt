cmake_minimum_required(VERSION 2.6)

add_sources(stagline
    bc_neq_1D_SL_Expl_1st.F90
    bc_neq_1D_SL_Expl_2nd.F90
    bc_neq_1D_SL_Impl_1st.F90
    bc_neq_1D_SL_Impl_2nd.F90
    bc_ablation_neq_1D_SL_Expl_2nd.F90
    bc_ablation_neq_1D_SL_Expl_1st.F90
    bc_ablation_general_neq_1D_SL_Expl_1st.F90
    bc_ablation_general_neq_1D_SL_Expl_2nd.F90
)
