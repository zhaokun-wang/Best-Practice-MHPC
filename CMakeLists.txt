cmake_minimum_required(VERSION 3.16)
project(Best_Practice_HPC LANGUAGES Fortran)

# -------------------------------------------------
# Option to select backend: ON = OpenACC, OFF = OpenMP
# -------------------------------------------------
option(USE_OPENACC "Enable OpenACC build (GPU). If OFF, uses OpenMP (CPU)" OFF)
# -------------
# Set compiler
# -------------
set(CMAKE_Fortran_COMPILER nvfortran)

find_program(NF_CONFIG_EXEC nf-config REQUIRED)
find_program(NC_CONFIG_EXEC nc-config REQUIRED)

execute_process(
        COMMAND ${NF_CONFIG_EXEC} --fflags
        OUTPUT_VARIABLE NETCDF_FFLAGS
        OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(
        COMMAND ${NF_CONFIG_EXEC} --flibs
        OUTPUT_VARIABLE NETCDF_FLIBS
        OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(
        COMMAND ${NF_CONFIG_EXEC} --includedir
        OUTPUT_VARIABLE NETCDF_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(
        COMMAND ${NC_CONFIG_EXEC} --libdir
        OUTPUT_VARIABLE NETCDF_LIB_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE)

# -------------------------------------------------
# Logic for OpenMP vs OpenACC vs CUDA
# -------------------------------------------------
if(USE_OPENACC)
    message(STATUS "Building with OpenACC (GPU)")
else()
    find_package(OpenMP REQUIRED)
    message(STATUS "Building with OpenMP (CPU)")
endif()

find_package(MPI REQUIRED)

file(GLOB source_files source/*.f90 source/*.F90)

add_executable(model.x ${source_files})

target_include_directories(model.x PRIVATE ${NETCDF_INCLUDE_DIR})

#-------------------------------------------------
# Compilation flags switch
# -------------------------------------------------
if(USE_OPENACC)
    # Flags for OpenACC
    target_compile_options(model.x PRIVATE
            -Ofast -fPIC
            ${NETCDF_FFLAGS}
            -acc -gpu=cc80 -cuda -Mpreprocess
    )
    target_link_options(model.x PRIVATE -acc -gpu=cc80 -cuda)
else()
    # Flags for OpenMP
    target_compile_options(model.x PRIVATE
            -Ofast -fPIC
            ${NETCDF_FFLAGS}
            ${OpenMP_Fortran_FLAGS}
            -Mpreprocess
    )
    target_link_libraries(model.x PRIVATE ${OpenMP_Fortran_FLAGS})
endif()

# -------------------------------------------------
# Link directories
# -------------------------------------------------
target_link_directories(model.x PRIVATE ${NETCDF_LIB_DIR})

# -------------------------------------------------
# Link libraries switch
# -------------------------------------------------
target_link_libraries(model.x PRIVATE
        ${NETCDF_FLIBS}
        MPI::MPI_Fortran
)

