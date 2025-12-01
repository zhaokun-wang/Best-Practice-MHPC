cmake_minimum_required(VERSION 3.16)
project(Best_Practice_HPC LANGUAGES Fortran)

execute_process(COMMAND nf-config --fflags
                OUTPUT_VARIABLE NETCDF_FFLAGS
                OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(COMMAND nf-config --flibs
                OUTPUT_VARIABLE NETCDF_FLIBS
                OUTPUT_STRIP_TRAILING_WHITESPACE)

file(GLOB source_files source/*.f90 source/*.F90)
add_executable(module.x ${source_files})

target_compile_options(module.x PUBLIC -Ofast -fPIC ${NETCDF_FFLAGS})
target_link_libraries(module.x PUBLIC ${NETCDF_FLIBS})