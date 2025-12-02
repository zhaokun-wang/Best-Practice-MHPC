#!/bin/bash
#SBATCH --job-name=mpi_scaling
#SBATCH --partition=boost_usr_prod
#SBATCH --qos=boost_qos_dbg
#SBATCH --time=00:05:00
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4       
#SBATCH --account=ICT25_MHPC_0
#SBATCH --output=logs/mpi_%j.out
#SBATCH --error=logs/mpi_%j.err

###############################
# Load modules (match what was used for compilation)
###############################
module purge
module load gcc/12.2.0
module load openmpi     
module load netcdf-c
module load netcdf-fortran

###############################
# Build the code
###############################
cd build
make clean
make -j

export OMP_NUM_THREADS=8
export GOTO_NUM_THREADS=1

###############################
# Run with MPI
###############################
srun ./model.x

