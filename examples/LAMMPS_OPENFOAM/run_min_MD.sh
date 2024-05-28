
cd openfoam
python clean.py -f
blockMesh
decomposePar
cd ../

#cplexec -c 1 "CPLIcoFoamFE4 -case ./openfoam -parallel" -m 1 "lmp_cpl -i lammps.in" 
cplc++ -o MD minimal_MD.cpp
mpiexec -n 1 CPLIcoFoamFE4 -case ./openfoam -parallel : -n 1 ./MD
