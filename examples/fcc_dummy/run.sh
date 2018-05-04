#Set OpenFOAM folder
OPEN_FOAM_CASE=openfoam/

#Clean CFD
cd ${OPEN_FOAM_CASE}
python clean.py -f
rm -f ../log.openfoam
blockMesh
decomposePar
cd ../

#Run job
CFD_EXE=$(which CPLSediFOAM)
cplexec -c 1 "${CFD_EXE} -case ${OPEN_FOAM_CASE} -parallel > log.openfoam" -m 1 "./MD_dummy_fcc.py"