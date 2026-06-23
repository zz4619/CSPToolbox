#!/bin/bash -l

# stack under SGE with Intel MPI.

# 1. Force bash as the executing shell.

#$ -S /bin/bash

# 2. Wallclock time (format hours:minutes:seconds).

#$ -l h_rt=24:00:00

# 3. RAM per process.

#$ -l mem=4G

# 5. Set the name of the job.

#$ -N Dummy_System

# 6. Select the MPI parallel environment and NUMBER of processes/cores.

#$ -pe mpi 80

# 7. Set the working directory to somewhere in your scratch space.  This is

# a necessary step with the upgraded software stack as compute nodes cannot

# write to $HOME.

# Run at current working dir

#$ -cwd

# 8. Run our MPI job.  GERun is a wrapper that launches MPI jobs on our clusters.

## $ -P AllUsers

#$ -P Gold

#$ -A Imperial_CEng

# RUN VASP
gerun /shared/ucl/apps/vasp/5.4.4-18apr2017/intel-2017/bin/vasp_std > vasp.out
