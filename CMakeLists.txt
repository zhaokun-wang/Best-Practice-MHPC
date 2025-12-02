cmake_minimum_required(VERSION 3.16)
project(Best_Practice_HPC LANGUAGES Fortran)

# -------------------------------------------------
# Set the Fortran compiler to nvfortran
# NOTE: This setting is often better handled outside
# of the CMakeLists.txt file using the -DCMAKE_Fortran_COMPILER=nvfortran
# option when configuring CMake, but this forces it.
# -------------------------------------------------
set(CMAKE_Fortran_COMPILER nvfortran)

# -------------------------------------------------
# Locate nf-config and nc-config
# -------------------------------------------------
find_program(NF_CONFIG_EXEC nf-config REQUIRED)
find_program(NC_CONFIG_EXEC nc-config REQUIRED)

# -------------------------------------------------
# Get compiler flags and library flags for NetCDF
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
# Remove OpenMP finding section
# The original: find_package(OpenMP REQUIRED)
# -------------------------------------------------
find_package(OpenMP REQUIRED)

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
# **Modified** to use OpenACC flags for nvfortran
# -------------------------------------------------
target_compile_options(model.x PRIVATE
    -Ofast -fPIC
    ${NETCDF_FFLAGS}
    # OpenACC flags for nvfortran
    -acc -gpu=managed
)

# -------------------------------------------------
# Link directories
# -------------------------------------------------
target_link_directories(model.x PRIVATE ${NETCDF_LIB_DIR})

# -------------------------------------------------
# Link libraries: NetCDF, MPI
# **Modified** to remove OpenMP linking
# OpenACC linking is handled by the -acc flag and nvfortran
# -------------------------------------------------
target_link_libraries(model.x PRIVATE
        ${NETCDF_FLIBS}
        ${OpenMP_Fortran_FLAGS}
        MPI::MPI_Fortran
)
