cmake_minimum_required(VERSION 3.16)
project(Best_Practice_HPC LANGUAGES Fortran)

# -------------------------------------------------
# Find OpenMP & MPI
# -------------------------------------------------
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

# -------------------------------------------------
# Your source files
# -------------------------------------------------
file(GLOB source_files source/*.f90 source/*.F90)

add_executable(model.x ${source_files})

# -------------------------------------------------
# NetCDF Configuration (FIX MANUALE UBUNTU)
# -------------------------------------------------
# Definiamo il percorso dove risiede la versione MPI della libreria C
set(MPI_NC_ROOT "/usr/lib/x86_64-linux-gnu/netcdf/mpi")

# 1. Include Directories
#    - /usr/include: qui di solito c'è il file netcdf.mod (Fortran)
#    - MPI_NC_ROOT/include: qui ci sono gli header C MPI
target_include_directories(model.x PRIVATE 
    "/usr/include"
    "${MPI_NC_ROOT}/include"
)

# 2. Link Directories (FONDAMENTALE)
#    Diciamo al linker di guardare PRIMA nella cartella MPI
target_link_directories(model.x PRIVATE "${MPI_NC_ROOT}/lib")

# -------------------------------------------------
# Compilation flags
# -------------------------------------------------
target_compile_options(model.x PRIVATE
    -Ofast -fPIC
    ${OpenMP_Fortran_FLAGS}
    # Rimosso NETCDF_FFLAGS perché conteneva path sporchi
)

# -------------------------------------------------
# Link libraries
# -------------------------------------------------
target_link_libraries(model.x PRIVATE
    netcdff               # La lib Fortran (trovata in sistema standard)
    netcdf                # La lib C (trovata in MPI_NC_ROOT grazie a target_link_directories)
    MPI::MPI_Fortran
    OpenMP::OpenMP_Fortran
)
