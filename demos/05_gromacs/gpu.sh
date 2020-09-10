#!/bin/bash -l
#SBATCH --job-name=gpu
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --time=00:15:00
#SBATCH --export=SIFPATH

# so that the srun's below don't inhert --export=
export SLURM_EXPORT_ENV=ALL

# load environment modules
module load Singularity CUDA

# path to SIF
SIF=${SIFPATH}/gromacs_2018.2.sif

# singularity command with required arguments
# "-B /cm/local/apps/cuda" and "-B ${EBROOTCUDA}" are required for the
# container to access the host CUDA drivers and libs
SINGULARITY="$(which singularity) exec --nv -B ${PWD}:/host_pwd \
  -B /cm/local/apps/cuda -B ${EBROOTCUDA} --pwd /host_pwd ${SIF}"

# extend container LD_LIBRARY_PATH so it can find CUDA libs
OLD_PATH=$(${SINGULARITY} printenv | grep LD_LIBRARY_PATH | awk -F= '{print $2}')
export SINGULARITYENV_LD_LIBRARY_PATH="${OLD_PATH}:${LD_LIBRARY_PATH}"

# run Gromacs preliminary step with container
srun ${SINGULARITY} \
    gmx grompp -f pme.mdp

# Run Gromacs MD with container
srun ${SINGULARITY} \
    gmx mdrun -ntmpi 1 -nb gpu -pin on -v -noconfout -nsteps 1000 -s topol.tpr -ntomp 1
