#!/bin/sh
#PBS -l walltime=0:10:00
#PBS -l select=1:ncpus=4:mem=16gb
#PBS -N SYSTEM_NAME

#script only meant for csofm local minimisation.
# run this script inside the system specific folder!
# LOAD modules

module load tools/prod
module load impi/2021.2.0-intel-compilers-2021.2.0
module load NAGlib/31.1-intel-compilers-2021.4.0
source $EBROOTNAGLIB/scripts/nagvars.sh dynamic int32 vendor
export NAG_KUSARI_FILE=$HOME/NAG_LICENCE/NAG_Keys_25_26.txt
export NAG_HOME=/sw-eb/software/NAGlib/31.1-intel-compilers-2021.4.0/
module load Gaussian/16.C.02-AVX2
eval "$(~/miniforge3/bin/conda shell.bash hook)"


### module load anaconda3/personal

# DEFINE directories and variables
INI=HPC_WORKDIR
EPHEMERAL=/rds/general/user/zz4619/ephemeral/${PBS_JOBNAME}/
OUTPUT=${PBS_O_WORKDIR}/LAM_Integrity
BIN=/rds/general/user/zz4619/home/CSP-Imperial-suite/CSO-FM/CrystalStructureOptimizer-FMV1.3.2/unix_executable/

mkdir -p ${EPHEMERAL}       # DO NOT DELETE. CSO-FM JOBS ARE DONE IN THE SAME PLACE
mkdir -p ${OUTPUT}

# COPY generic input files to EPHEMERAL
cd ${PBS_O_WORKDIR}
cp -r ${PBS_JOBNAME}.res ${EPHEMERAL}/
cp CSO_FM.input *.inter *.input Zmatrix ${EPHEMERAL}/
#cp CSO_FM.input *.inter *.input ${EPHEMERAL}/
cp ${BIN}/LAM_Integrity ${EPHEMERAL}/

# RUN program
cd ${EPHEMERAL}
#./LAM_Integrity CSO_FM.input ${PBS_JOBNAME}.res all > Integrity.output
# Stop after 10 minutes if it is still looping
timeout --signal=TERM --kill-after=20s 5m ./LAM_Integrity CSO_FM.input "${PBS_JOBNAME}.res" DOF > Integrity.output 2>&1
rc=$?

# timeout exit code is 124
if [ "$rc" -eq 124 ]; then
  echo "LAM_Integrity timed out and was terminated" >> Integrity.output
fi

# COPY back outputs
# COPY SPECIFIC FILES ONLY, NOT ENTIRE EPHEMERAL. WILL DOUBLE COUNT CSOFM results
cp ${EPHEMERAL}/Integrity.output ${OUTPUT}/
cp ${EPHEMERAL}/CSO_FM.input ${EPHEMERAL}/mol.input ${EPHEMERAL}/Zmatrix ${OUTPUT}/


