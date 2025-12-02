cmake_minimum_required(VERSION 3.16)
project(Best_Practice_HPC LANGUAGES Fortran)

# -------------------------------------------------
# Locate nf-config and nc-config
# -------------------------------------------------
find_program(NF_CONFIG_EXEC nf-config REQUIRED)
find_program(NC_CONFIG_EXEC nc-config REQUIRED)

# -------------------------------------------------
# Get compiler flags and library flags
# -------------------------------------------------
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
# Find OpenMP
# -------------------------------------------------
#find_package(OpenMP REQUIRED)

# -------------------------------------------------
# Find MPI (Fortran)
# -------------------------------------------------
find_package(MPI REQUIRED)

# -------------------------------------------------
# Your source files
# -------------------------------------------------
file(GLOB source_files source/*.f90 source/*.F90)

add_executable(model.x ${source_files})

# -------------------------------------------------
# Include *.mod directory for NetCDF-Fortran
# -------------------------------------------------
target_include_directories(model.x PRIVATE ${NETCDF_INCLUDE_DIR})

# -------------------------------------------------
# Compilation flags
# -------------------------------------------------
target_compile_options(model.x PRIVATE
        -Ofast -fPIC
        ${NETCDF_FFLAGS}
)

#${OpenMP_Fortran_FLAGS}

# -------------------------------------------------
# Link directories
# -------------------------------------------------
target_link_directories(model.x PRIVATE ${NETCDF_LIB_DIR})

# -------------------------------------------------
# Link libraries: NetCDF, OpenMP, MPI
# -------------------------------------------------
target_link_libraries(model.x PRIVATE
        ${NETCDF_FLIBS}
        MPI::MPI_Fortran)

#${OpenMP_Fortran_FLAGS}