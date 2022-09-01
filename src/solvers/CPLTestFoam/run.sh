cd Couette/
decomposePar -force
cd ../
wclean
wmake
cplc++ minimal_MD.cpp -o MD
mpiexec -n 2 CPLTestFoam -parallel -case Couette/ | grep 'CPLTestFoam' & 
mpiexec -n 2 ./MD
