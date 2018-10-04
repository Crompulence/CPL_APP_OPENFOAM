#!/bin/bash

#Clean MD
cd lammps/
rm -f ../log.lammps vmd_out.dcd thermo_output.txt
cd ../

#Run simulation
cplexec -c 1 "./python_openfoam/CFD_dummy_single.py" -m 1 "lmp_cpl < lammps/single.in"
