
cd openfoam
python clean.py -f
blockMesh
decomposePar
cd ../

cplexec -c 1 "./CPLSediFOAM -case ./openfoam -parallel" -m 1 "./python_dummy/MD_no_particle.py" 

