---
title: "Molecular dynamics with GPU containers"
teaching: 5
exercises: 10
questions:
objectives:
- Get started with Nvidia GPU containers for HPC applications
keypoints:
- You can run containerised GPU applications using the flag `--nv`
- Singularity transparently interfaces with HPC schedulers such as Slurm
---


### Run a molecular dynamics simulation on a GPU with containers

For our example we are going to use Gromacs, a quite popular molecular dynamics package, among the ones that have been optimised to run on GPUs through Nvidia containers.

First, let us cd into `demos/05_gromacs`, and ensure that `$SIFPATH` is defined:

```
cd $ERNZ20/demos/05_gromacs
```
{: .bash}

Now, let's verify that the container image `nvcr.io/hpc/gromacs:2018.2` has been correctly downloaded:

```
ls $SIFPATH/gromacs*
```
{: .bash}

```
/nesi/nobackup/nesi99991/nzrse-containers/demos/sif/gromacs_2018.2.sif
```
{: .output}

The current directory has got sample input files picked from the collection of [Gromacs benchmark examples](ftp://ftp.gromacs.org/pub/benchmarks/water_GMX50_bare.tar.gz). In particular, we're going to use the subset `water-cut1.0_GMX50_bare/1536/`. First let's `gunzip` one of the required input files:

```
gunzip conf.gro.gz
```
{: .bash}

In theory, all we should need to do from a Singularity perspective to run a GPU application from a container is to add the runtime flag `--nv`. This will make Singularity look for the Nvidia drivers in the host, and mount them inside the container. 

However, at NeSI (as at many other HPC centres) we have installed the Nvidia drivers and libraries in non-standard locations, so some additional steps are required.

On the host system side, when running GPU applications through Singularity the only requirement consists of the Nvidia driver for the relevant GPU card (the corresponding file is typically called `libcuda.so.<VERSION>`).

Do not execute the next following commands, let us just have a look at them.

* Setup $SINGULARITY variable as convenience
  ```
  SINGULARITY="singularity exec --nv -B ${PWD}:/host_pwd \
    -B /cm/local/apps/cuda -B ${EBROOTCUDA} --pwd /host_pwd ${SIF}"
  ```
  {: .bash}
* Workaround to ensure Singularity can find the host CUDA libraries
  ```
  OLD_PATH=$(${SINGULARITY} printenv | grep LD_LIBRARY_PATH | awk -F= '{print $2}')
  export SINGULARITYENV_LD_LIBRARY_PATH="${OLD_PATH}:${LD_LIBRARY_PATH}"
  ```
  {: .bash}
* Preliminary step
  ```
  ${SINGULARITY} gmx grompp -f pme.mdp
  ```
  {: .bash}
* Production step
  ```
  ${SINGULARITY} gmx mdrun -ntmpi 1 -nb gpu -pin on -v -noconfout -nsteps 5000 -s topol.tpr -ntomp 1
  ```
  {: .bash}

GPU resources are usually made available in HPC systems through schedulers, to which Singularity natively and transparently interfaces. So, for instance let us have a look in the current directory at the Slurm batch script called `gpu.sh`:

```
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
    gmx mdrun -ntmpi 1 -nb gpu -pin on -v -noconfout -nsteps 5000 -s topol.tpr -ntomp 1
```
{: .bash}

Basically, we have just combined the Slurm command `srun` with `singularity exec <..>`, with a few additional steps to enable Singularity to locate the Nvidia drivers and libraries on the host. We can submit the script with:

```
sbatch --account=nesi99991 --reservation=jupyter gpu.sh
```
{: .bash}


> ## Running this example on NeSI
>
> The job will take about 5 minutes to complete and we only have
> a few GPUs available for this workshop, so you may wait in the
> queue for a while. Check the queue with:
> ```
> $ squeue -u $USER
> ```
> {: .bash}
{: .callout}


> ## Nvidia GPU Cloud, a useful resource
>
>The GPU manufacturer *Nvidia* has a dedicated web registry for container images, shipping GPU optimised applications: <https://ngc.nvidia.com>.
>
>You can browse the available containerised packages through the various category boxes. E.g. click on the **High Performance Computing** box, then click on the >**Gromacs** one.
>
> NeSI also has documentation for running NGC containers on our platforms: <https://support.nesi.org.nz/hc/en-gb/articles/360001500156-NVIDIA-GPU-Containers>.
{: .callout}
