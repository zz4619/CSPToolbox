#Defining dirs
SYSTEM_NAME=DUMMY_SYSTEM_NAME

OUTPUT=./local_minimisation_output
BIN=/rds/general/user/zz4619/home/CSP-Imperial-suite/CSO-RM/CrystalStructureOptimzer-RMV1.4.1/unix_executable

rm -rf ${OUTPUT}
mkdir -p ${OUTPUT}

echo 12 | ${BIN}/SetupUnitCell CSO_RM.input > CSO_RM.output
echo 'SetupUnitCell complete, now running MinimiseLatticeEnergy_RM' >> CSO_RM.output
${BIN}/MinimiseLatticeEnergy_RM CSO_RM.input ${SYSTEM_NAME} >> CSO_RM.output

