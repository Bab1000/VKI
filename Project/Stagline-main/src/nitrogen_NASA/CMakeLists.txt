cmake_minimum_required(VERSION 2.6)

add_sources(stagline
    mod_function_pointer_NASA.F90
    mod_nitrogen_NASA_CFD_alpha.F90
    mod_nitrogen_NASA_CFD_eq.F90
    mod_nitrogen_NASA_CFD_prop.F90
    mod_nitrogen_NASA_CFD_source_BRVC.F90
    mod_nitrogen_NASA_CFD_source_MT_TTint.F90
    mod_nitrogen_NASA_CFD_source_RVC.F90
    mod_nitrogen_NASA_CFD_source_VC.F90
    mod_nitrogen_NASA_CFD_source_VC_rr.F90
    mod_nitrogen_NASA_CFD_transport.F90
    mod_nitrogen_NASA_zeta.F90
)
