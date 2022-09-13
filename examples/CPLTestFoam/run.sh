blockMesh
decomposePar -force
cplc++ minimal_MD.cpp -o MD
mpiexec -n 1 CPLTestFoam -parallel &
mpiexec -n 1 ./MD
