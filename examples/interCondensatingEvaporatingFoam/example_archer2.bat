#!/bin/bash

#SBATCH --job-name=my_cpl_demo
#SBATCH --time=0:10:0
#SBATCH --exclusive
#SBATCH --export=none
#SBATCH --account=ecseaf01

#SBATCH --partition=standard
#SBATCH --qos=standard

#SBATCH --nodes=2

# single thread export overriders any declaration in srun
export OMP_NUM_THREADS=1

module load openfoam/com/v2106
module load lammps/13_Jun_2022
module load cray-python
#module load xthi

cd /work/ecseaf01/ecseaf01/gavboi/cpl-library
source SOURCEME.sh
cd /work/ecseaf01/ecseaf01/gavboi/CPL_APP_OPENFOAM
source SOURCEME.sh
cd examples/interCondensatingEvaporatingFoam

blockMesh
decomposePar -force

SHARED_ARGS="--distribution=block:block --hint=nomultithread"

#srun ${SHARED_ARGS} --het-group=0 --nodes=1 --tasks-per-node=2 xthi : --het-group=1 --nodes=1 --tasks-per-node=2 xthi

srun ${SHARED_ARGS} --het-group=0 --nodes=1 --tasks-per-node=2 MD : --het-group=1 --nodes=1 --tasks-per-node=2 CPLinterCondensatingEvaporatingFoam -parallel

