cmake_minimum_required(VERSION 2.6)

# Adds sources to a given target
function (add_sources target)
  # define the <target>_SRCS properties if necessary
  get_property(prop_defined GLOBAL PROPERTY ${target}_SRCS DEFINED)
  if(NOT prop_defined)
    define_property(GLOBAL PROPERTY ${target}_SRCS
      BRIEF_DOCS "Sources for the ${target} target"
      FULL_DOCS "List of source files for the ${target} target")
  endif()
  # create list of sources (absolute paths)
  set(SRCS)
  foreach(src IN LISTS ARGN)
    if(NOT IS_ABSOLUTE "${src}")
      get_filename_component(src "${src}" ABSOLUTE)
    endif()
    list(APPEND SRCS "${src}")
  endforeach()
  # append to global property
  set_property(GLOBAL APPEND PROPERTY "${target}_SRCS" "${SRCS}")
endfunction()

# Source code here
add_sources(stagline
    compute_flowfield.F90
    main.F90
    mesh.F90
    stop_cond.F90
    update_cfl.F90
)

# Generate the config header file
configure_file(
    config.h.in
    ${CMAKE_CURRENT_SOURCE_DIR}/config.h
)

# Subdirectories
add_subdirectory(ablation_modules)
add_subdirectory(algebra)
add_subdirectory(apply_bc_1D)
add_subdirectory(apply_bc_1D_SL)
add_subdirectory(bc_1D)
add_subdirectory(bc_1D_SL)
add_subdirectory(bc_neq_1D)
add_subdirectory(bc_neq_1D_SL)
add_subdirectory(diff_flux_1D)
add_subdirectory(diff_flux_1D_SL)
add_subdirectory(diff_flux_neq_1D)
add_subdirectory(diff_flux_neq_1D_SL)
add_subdirectory(eigensystem_1D)
add_subdirectory(eigensystem_1D_SL)
add_subdirectory(eigensystem_neq_1D)
add_subdirectory(eigensystem_neq_1D_SL)
add_subdirectory(finalize)
add_subdirectory(initialize)
add_subdirectory(inv_flux_1D)
add_subdirectory(inv_flux_1D_SL)
add_subdirectory(inv_flux_neq_1D)
add_subdirectory(inv_flux_neq_1D_SL)
add_subdirectory(library_interface)
add_subdirectory(mesh_adaptor)
add_subdirectory(modules)
add_subdirectory(neq)
add_subdirectory(pg)
add_subdirectory(physical_data)
add_subdirectory(precon_1D)
add_subdirectory(precon_1D_SL)
add_subdirectory(pre_processor)
add_subdirectory(procedure_pointer)
add_subdirectory(radiation)
add_subdirectory(rec_1D)
add_subdirectory(rec_1D_SL)
add_subdirectory(rec_neq_1D)
add_subdirectory(rec_neq_1D_SL)
add_subdirectory(source_term)
add_subdirectory(source_term_neq)
add_subdirectory(time_integration_1D)
add_subdirectory(time_integration_1D_SL)
add_subdirectory(time_step_1D)
add_subdirectory(time_step_1D_SL)
add_subdirectory(transport_coeff)
add_subdirectory(write_solution)

# Chemistry libraries
if (${ARGON_CR} MATCHES "ON")
    add_subdirectory(argon_CR)
endif()

if (${NITROGEN_BARI} MATCHES "ON")
    add_subdirectory(nitrogen_Bari)
endif()

if (${NITROGEN_DSMC} MATCHES "ON")
    add_subdirectory(nitrogen_DSMC)
endif()

if (${NITROGEN_FHO} MATCHES "ON")
    add_subdirectory(nitrogen_FHO)
endif()

if (${NITROGEN_NASA} MATCHES "ON")
    add_subdirectory(nitrogen_NASA)
endif()

if (${NITROGEN_PARK} MATCHES "ON")
    add_subdirectory(nitrogen_Park)
endif()

if (${NITROGEN_TTLTH} MATCHES "ON")
    add_subdirectory(nitrogen_TTLTH)
endif()

if (${POLYAT_GAS_WCU} MATCHES "ON")
    add_subdirectory(Polyat_gas_WCU)
endif()


