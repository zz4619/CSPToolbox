
#!/bin/sh
#PBS -l walltime=24:00:00
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
PAT=$PBS_O_WORKDIR
INI=HPC_WORKDIR
EPHEMERAL=/rds/general/user/zz4619/ephemeral/${PBS_JOBNAME}/
OUTPUT=${PBS_O_WORKDIR}/local_minimisation_output
BIN=/rds/general/user/zz4619/home/CSP-Imperial-suite/CSO-FM/CrystalStructureOptimizer-FMV1.3.2/unix_executable/

struc_dir=${PAT}/clustered


# Check for reasons to stop from previous runs
#if grep -q "Polished Optimal Lattice energy" ${OUTPUT}/CSO_FM.output
if grep -q "ERROR: Spec types could not be loaded" ${OUTPUT}/CSO_FM.output; then                # Likely erroenous conformation
    echo "This run was erroneous. Skipping"
    exit
fi
if grep -q "ERROR: Error reading g_calc.log. Likely didn't terminate correctly with smaller steps." ${OUTPUT}/CSO_FM.output; then       # Likely erroenous conformation
    echo "This run was erroneous. Skipping"
    exit
fi
if grep -q "Degrees of freedom successfully" ${OUTPUT}/CSO_FM.output; then			            #CSOFM pass flag
    if grep -q "Polished Optimal Lattice energy" ${OUTPUT}/CSO_FM.output; then			        #Polish pass flag
        echo "This run has completed. Skipping"
        exit
    fi
fi

rm -rf ${EPHEMERAL}
mkdir -p ${EPHEMERAL}
rm -rf ${OUTPUT}
mkdir -p ${OUTPUT}

# COPY initialisation structures.res to EPHEMERAL
echo "Conducting CSO-FM for structure ${PBS_JOBNAME}"
cp -r ${INI}/${PBS_JOBNAME}/${PBS_JOBNAME}.res ${EPHEMERAL}/

# COPY generic input files from the corresponding LAM INTEGRITY to CSOFM/EPHEMERAL
cd ${EPHEMERAL}
cp ${PBS_O_WORKDIR}/CSO_FM.input ${PBS_O_WORKDIR}/*.inter ${PBS_O_WORKDIR}/*.input ${PBS_O_WORKDIR}/Zmatrix ${EPHEMERAL}/
cp ${BIN}/MinimiseLatticeEnergy_FM ${EPHEMERAL}/

cd ${EPHEMERAL}
./MinimiseLatticeEnergy_FM CSO_FM.input ${PBS_JOBNAME}.res DOF > CSO_FM.output

# COPY back outputs. Do not copy all files or it will be too large
cp -r ${EPHEMERAL}/CSO_FM.input ${EPHEMERAL}/CSO_FM.output ${EPHEMERAL}/*.out ${EPHEMERAL}/tmp ${OUTPUT}/
cp -r ${EPHEMERAL}/*.inter ${EPHEMERAL}/mol.input ${EPHEMERAL}/Zmatrix ${EPHEMERAL}/results ${OUTPUT}/
cp -r ${EPHEMERAL}/output_LAM.* ${EPHEMERAL}/polished.* ${OUTPUT}/
cp ${PBS_JOBNAME}.res ${OUTPUT}/



