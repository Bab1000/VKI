# Basic setup of project and minimum required CMake version
cmake_minimum_required(VERSION 2.6)
project(VKI_Stagnation_Line)
enable_language(Fortran)

# Library options
option(GOTO_BLAS "compile the GOTO BLAS library" OFF)
option(NITROGEN_PARK "include support for the Nitrogen chemistry library of Park" OFF)
option(NITROGEN_NASA "include support for the Nitrogen chemistry library of NASA" OFF)
option(NITROGEN_BARI "include support for the Nitrogen chemistry library of Bari" OFF)
option(NITROGEN_FHO  "include support for the Nitrogen FHO chemistry library" OFF)
option(NITROGEN_TTLTH "include support for the Nitrogen TTLTH chemistry library" OFF)
option(NITROGEN_DSMC "include support for the Nitrogen DSMC chemistry library" OFF)
option(ARGON_CR "include support for the Argon CR library" OFF)
option(POLYAT_GAS_WCU "include support for the POLYAT_GAS_WCU chemistry library" OFF)
option(MUTATIONPP "include support for the Mutation++ library" ON)
option(CARBONABLA "include support for the ablative boundary conditions" ON)

# Descend into the src directory to build all targets and libraries
add_subdirectory(src)
# Mutation++ directories
set(MPP_INSTALL_DIR "" CACHE PATH "Mutation++ installation directory")
if (${MUTATIONPP} MATCHES "ON")
#    if (NOT ${MPP_INSTALL_DIR} MATCHES "")
        include_directories(${MPP_INSTALL_DIR}/include/mutation++)
        include_directories(${MPP_INSTALL_DIR}/include/)
        link_directories(${MPP_INSTALL_DIR}/lib)
#    endif()
endif()
    

# Main executable
get_property(stagline_SRCS GLOBAL PROPERTY stagline_SRCS)
set_source_files_properties(${stagline_SRCS} PROPERTIES COMPILE_FLAGS -cpp)
add_executable(stagline ${stagline_SRCS})

if (${MUTATIONPP} MATCHES "ON")
    target_link_libraries(stagline -lmutation++_fortran)
endif()

# BLAS
set (BLASS_LIB_PATH "" CACHE FILEPATH "The BLAS library file")
if (${GOTO_BLAS} MATCHES "ON")
    target_link_libraries(stagline ${BLASS_LIB_PATH})
endif()

# Installation
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS stagline DESTINATION bin)

# Doxygen generated documentation
option (BUILD_DOCUMENTATION 
    "Use Doxygen to create the HTML based API documentation" OFF)

if (BUILD_DOCUMENTATION)
    FIND_PACKAGE(Doxygen)
    
    if (NOT DOXYGEN_FOUND)
        message(FATAL_ERROR 
            "Doxygen is needed to build the documentation. Please install it 
             correctly")
    endif()
    # Configure the Template Doxyfile for our specific project
    configure_file(Doxyfile.in 
                    ${PROJECT_BINARY_DIR}/Doxyfile  @ONLY IMMEDIATE)
    # Add a custom target to run Doxygen when ever the project is built
    add_custom_target (docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
        SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)
endif()

