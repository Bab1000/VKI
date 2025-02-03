cmake_minimum_required(VERSION 2.6)

if (${ARGON_CR} MATCHES "ON")
    add_sources(stagline mod_Argon_CR_library.F90)
endif()

if (${MUTATIONPP} MATCHES "ON")
    add_sources(stagline mod_Mutationpp_library.F90)
endif()

if (${NITROGEN_BARI} MATCHES "ON")
    add_sources(stagline mod_Nitrogen_Bari_library.F90)
endif()

if (${NITROGEN_DSMC} MATCHES "ON")
    add_sources(stagline mod_Nitrogen_DSMC_library.F90)
endif()

if (${NITROGEN_FHO} MATCHES "ON")
    add_sources(stagline mod_Nitrogen_FHO_library.F90)
endif()

if (${NITROGEN_NASA} MATCHES "ON")
    add_sources(stagline mod_Nitrogen_NASA_library.F90)
endif()

if (${NITROGEN_PARK} MATCHES "ON")
    add_sources(stagline mod_Nitrogen_Park_library.F90)
endif()

if (${NITROGEN_TTLTH} MATCHES "ON")
    add_sources(stagline mod_Nitrogen_TTLTH_library.F90)
endif()

if (${POLYAT_GAS_WCU} MATCHES "ON")
    add_sources(stagline mod_Polyat_gas_WCU_library.F90)
endif()



