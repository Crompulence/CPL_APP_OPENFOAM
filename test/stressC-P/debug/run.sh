#mpirun -n 27 python2 dummyMD_forces.py & PID=$!; mpirun -n 27 CPLIcoFoam -parallel -case test_forces_case ; wait $PID
mpirun -n 27 python2 dummyMD_forces.py : -n 27 CPLIcoFoam -parallel -case test_forces_case
reconstructPar -case test_forces_case
stressComponents -case test_forces_case
writeCellCentres -case test_forces_case
