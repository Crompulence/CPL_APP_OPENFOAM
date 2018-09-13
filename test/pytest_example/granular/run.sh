
cd openfoam
python clean.py -f
blockMesh
decomposePar
cd ../

mpiexec -n 1 CPLSediFOAM -case ./openfoam -parallel : -n 1 py.test -v  ./python_dummy/test_vs_couette_analytical_aspytest.py

