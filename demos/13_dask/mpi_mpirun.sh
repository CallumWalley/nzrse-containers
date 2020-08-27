#!/bin/bash

module load Singularity impi
module unload XALT
unset I_MPI_PMI_LIBRARY

mpirun -np 3 singularity run -B $PWD $SIFPATH/dask-mpi_latest.sif dask_example.py
