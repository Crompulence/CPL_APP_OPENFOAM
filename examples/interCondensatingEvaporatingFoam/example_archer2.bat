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

module load other-software
module load cpl-openfoam/2106
source $FOAM_CPL_APP/SOURCEME.sh

# using your own installtion: remove the last three lines and use these four lines instead
# remmeber to update the path to the two SOURCEME.sh files
#module load openfoam/com/v2106
#module switch gcc gcc/10.3.0
#module load cray-python
#source /work/y23/shared/cpl-openfoam-lammps/cpl-library/SOURCEME.sh
#source /work/y23/shared/cpl-openfoam-lammps/CPL_APP_OPENFOAM/SOURCEME.sh

# load restor0Dir tool
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
restore0Dir

blockMesh > log.blockMesh 2>&1
decomposePar > log.decomposePar 2>&1

SHARED_ARGS="--distribution=block:block --hint=nomultithread"

srun ${SHARED_ARGS} --het-group=0 --nodes=1 --tasks-per-node=2 MD : --het-group=1 --nodes=1 --tasks-per-node=2 CPLinterCondensatingEvaporatingFoam -parallel

reconstructPar > log.reconstructPar 2>&1


