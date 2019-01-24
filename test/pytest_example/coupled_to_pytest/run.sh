if [[ -n "$1" ]]; then
    exe=$1
else
    exe=CPLSediFOAM
fi

cd openfoam
python clean.py -f
blockMesh
decomposePar
cd ../

mpiexec -n 1 $exe -case ./openfoam -parallel : -n 1 pytest -v ./python_dummy/test_vs_couette_analytical_aspytest.py

