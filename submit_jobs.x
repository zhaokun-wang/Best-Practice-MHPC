#!/bin/bash
# submit_scaling_jobs.sh
# This script submits separate SLURM jobs for 1, 2, 4, 8, 16 nodes

###############################
# Define parameters
###############################
NTASKS_PER_NODE=4
NODE_LIST=(1)
JOB_NAME="acc_scaling"
PARTITION="boost_usr_prod"
ACCOUNT="ICT25_MHPC_0"
TIME="00:40:00"        # adjust if needed

###############################
# Load modules (used for compilation)
###############################

module purge
module load gcc
module load netcdf-fortran/4.6.1--hpcx-mpi--2.19--nvhpc--24.5
module load nvhpc

###############################
# Build the code
###############################
cd build
cmake -DUSE_OPENACC=ON ..
make clean
make -j

###############################
# Loop over node counts and submit jobs
###############################
for NNODES in "${NODE_LIST[@]}"; do
    NTASKS=$(( NNODES * NTASKS_PER_NODE ))
    LOGOUT="../logs/accn${NNODES}_%j.out"
    LOGERR="../logs/accn${NNODES}_%j.err"

    echo "Submitting job: $NNODES nodes, $NTASKS tasks"

    sbatch --job-name=${JOB_NAME}_${NNODES} \
           --partition=${PARTITION} \
           --time=${TIME} \
	   --gres=gpu:4 \
           --nodes=${NNODES} \
           --ntasks-per-node=${NTASKS_PER_NODE} \
           --account=${ACCOUNT} \
           --output=${LOGOUT} \
           --error=${LOGERR} \
           --wrap="export GOTO_NUM_THREADS=1; \
		   srun ./model.x"
done

