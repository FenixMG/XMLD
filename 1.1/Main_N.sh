#!/bin/bash
#SBATCH -J Job
#SBATCH -o Job.o%j
#SBATCH -N 7
#SBATCH -p comp
#SBATCH -t 1000:00:00
echo $SLURM_JOB_ID
echo $SLURM_JOB_NAME
echo $SLURM_JOB_NUM_NODES
# Load your modules here
module load compilers/intel/parallel_studio_xe_2015/15.0.1 
module load tools/intel/impi/5.0.2.044
module load applications/anaconda3/4.0.0
# Run your task here
source activate AnacondaEnv

srun --hint=nomultithread -N 1 --ntasks=1 --ntasks-per-node=1 --ntasks-per-socket=1 python -u Job_9/Main_LNS9.py > Job_9/output.out 2>Job_9/error.txt &
srun --hint=nomultithread -N 1 --ntasks=1 --ntasks-per-node=1 --ntasks-per-socket=1 python -u Job_19/Main_LNS19.py > Job_19/output.out 2>Job_19/error.txt &
srun --hint=nomultithread -N 1 --ntasks=1 --ntasks-per-node=1 --ntasks-per-socket=1 python -u Job_26/Main_LNS26.py > Job_26/output.out 2>Job_26/error.txt &
srun --hint=nomultithread -N 1 --ntasks=1 --ntasks-per-node=1 --ntasks-per-socket=1 python -u Job_40/Main_LNS40.py > Job_40/output.out 2>Job_40/error.txt &
srun --hint=nomultithread -N 1 --ntasks=1 --ntasks-per-node=1 --ntasks-per-socket=1 python -u Job_43/Main_LNS43.py > Job_43/output.out 2>Job_43/error.txt &
srun --hint=nomultithread -N 1 --ntasks=1 --ntasks-per-node=1 --ntasks-per-socket=1 python -u Job_55/Main_LNS55.py > Job_55/output.out 2>Job_55/error.txt &
srun --hint=nomultithread -N 1 --ntasks=1 --ntasks-per-node=1 --ntasks-per-socket=1 python -u Job_64/Main_LNS64.py > Job_64/output.out 2>Job_64/error.txt &

wait
source deactivate AnacondaEnv
# matlab -nodisplay -nosplash -nodesktop -r "run('./clusterESRF.m') ;exit;"
                                  
