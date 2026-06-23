
#!/bin/sh
#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=16:mem=96gb
#PBS -N SYSTEM_NAME

module load tools/prod
module load impi/2021.2.0-intel-compilers-2021.2.0
module load NAGlib/31.1-intel-compilers-2021.4.0
source $EBROOTNAGLIB/scripts/nagvars.sh dynamic int32 vendor
export NAG_KUSARI_FILE=$HOME/NAG_LICENCE/NAG_Keys_25_26.txt
export NAG_HOME=/sw-eb/software/NAGlib/31.1-intel-compilers-2021.4.0/
module load Gaussian/16.C.02-AVX2

#Defining dirs
EPHEMERAL=/rds/general/user/zz4619/ephemeral/${PBS_JOBNAME}/
OUTPUT=${PBS_O_WORKDIR}/local_minimisation_output
BIN=/rds/general/user/zz4619/home/CSP-Imperial-suite/CSO-RM/CrystalStructureOptimzer-RMV1.4.1/unix_executable

rm -rf ${EPHEMERAL}
mkdir -p ${EPHEMERAL}
rm -rf ${OUTPUT}
mkdir -p ${OUTPUT}

# COPY initialisation structures.res to EPHEMERAL
echo "Conducting CSO-RM for structure ${PBS_JOBNAME}"
cp -r ${PBS_O_WORKDIR}/${PBS_JOBNAME}.res ${PBS_O_WORKDIR}/${PBS_JOBNAME}.info ${EPHEMERAL}/

# COPY generic input files from the corresponding LAM INTEGRITY to CSOFM/EPHEMERAL
cd ${EPHEMERAL}
cp ${PBS_O_WORKDIR}/CSO_RM.input ${PBS_O_WORKDIR}/*.inter ${EPHEMERAL}/
cp ${BIN}/SetupUnitCell ${EPHEMERAL}/
cp ${BIN}/MinimiseLatticeEnergy_RM ${EPHEMERAL}/

cd ${EPHEMERAL}
#./SetupUnitCell CSO_RM.input < 12 > CSO_RM.output
echo 12 | ./SetupUnitCell CSO_RM.input > CSO_RM.output
echo 'SetupUnitCell complete, now running MinimiseLatticeEnergy_RM' >> CSO_RM.output
./MinimiseLatticeEnergy_RM CSO_RM.input ${PBS_JOBNAME} >> CSO_RM.output


# COPY back outputs. Do not copy all files or it will be too large
rm -rf ${EPHEMERAL}/SetupUnitCell
rm -rf ${EPHEMERAL}/MinimiseLatticeEnergy_RM
cp -r ${EPHEMERAL}/* ${OUTPUT}/

