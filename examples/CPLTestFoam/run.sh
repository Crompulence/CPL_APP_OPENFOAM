decomposePar -force
cplc++ minimal_MD.cpp -o MD
mpiexec -n 2 CPLTestFoam -parallel | grep 'CPLTestFoam' & 
mpiexec -n 2 ./MD
