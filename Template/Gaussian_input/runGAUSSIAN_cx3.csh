#!/bin/sh
#PBS -lwalltime=6:00:00
#PBS -l select=1:ncpus=16:mpiprocs=16:mem=40gb
#PBS -N SYSTEM_NAME
#PBS -J 1-2

# Just adding this to bypass slow queues
if [[ $PBS_ARRAY_INDEX -eq 2 ]];
then
    echo "Dummy array job"
    exit
fi

#DEFINE directories
name=SYSTEM_NAME

EPHEMERAL=~/../ephemeral/${PBS_JOBNAME}_${name}/
INPUT=${PBS_O_WORKDIR}/
OUTPUT=${PBS_O_WORKDIR}/

#Load modules
module load tools/prod
module load Gaussian/16.C.02-AVX2

# CREATE directories
rm -rf $EPHEMERAL
mkdir -p $EPHEMERAL

#Execute
cd $EPHEMERAL
pwd
cp ${INPUT}/${name}.com $EPHEMERAL/
g16 ${name}.com ${name}.log 
formchk -3 hess.chk hess.fchk

#Copy back
cp $EPHEMERAL/* ${OUTPUT}/