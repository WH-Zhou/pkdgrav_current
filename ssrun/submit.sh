#!/bin/bash

nodes=$1
runtime=4
checkpoint_file=$HOME/pkdgrav/etc/checkpoint_rub_hist.fdl

bindir=`pwd -P`

echo "submitting:"
echo "#!/bin/bash
#PBS -N $nodes.$2.pkdgrav
#PBS -joe
#PBS -l nodes=$nodes:ppn=16,walltime=$runtime:00:00

cd \$PBS_O_WORKDIR

echo \"job params:

submitted `date`

running on:
\`cat \$PBS_NODEFILE\`

PBS_ENVIRONMENT = \$PBS_ENVIRONMENT
PBS_JOBID = \$PBS_JOBID
PBS_MOMPORT = \$PBS_MOMPORT
PBS_NP = \$PBS_NP
PBS_O_HOME = \$PBS_O_HOME
PBS_O_LOGNAME = \$PBS_O_LOGNAME
PBS_O_QUEUE = \$PBS_O_QUEUE
PBS_O_WORKDIR = \$PBS_O_WORKDIR
PBS_VERSION = \$PBS_VERSION
PBS_GPUFILE = \$PBS_GPUFILE
PBS_JOBNAME = \$PBS_JOBNAME
PBS_NODEFILE = \$PBS_NODEFILE
PBS_NUM_NODES = \$PBS_NUM_NODES
PBS_O_HOST = \$PBS_O_HOST
PBS_O_MAIL = \$PBS_O_MAIL
PBS_O_SERVER = \$PBS_O_SERVER
PBS_QUEUE = \$PBS_QUEUE
PBS_VNODENUM = \$PBS_VNODENUM
PBS_JOBCOOKIE = \$PBS_JOBCOOKIE
PBS_MICFILE = \$PBS_MICFILE
PBS_NODENUM = \$PBS_NODENUM
PBS_NUM_PPN = \$PBS_NUM_PPN
PBS_O_LANG = \$PBS_O_LANG
PBS_O_PATH = \$PBS_O_PATH
PBS_O_SHELL = \$PBS_O_SHELL
PBS_TASKNUM = \$PBS_TASKNUM
\" > queue.out

export PKDGRAV_CHECKPOINT_FDL=$checkpoint_file
export LAUNCH_APP=\`which mpirun\`

cm-launcher \`which mpirun\` -machinefile \$PBS_NODEFILE -np `echo "$OMP_NUM_THREADS*$nodes" | bc -q` ssic               | tee create.out
cm-launcher \`which mpirun\` -machinefile \$PBS_NODEFILE -np `echo "$OMP_NUM_THREADS*$nodes" | bc -q` pkdgrav.mpi ss.par | tee run.out

for i in \`cat \$PBS_NODEFILE\`; do
    ssh $i "killall -9 pkdgrav"
done

" | tee submit.sh

`which qsub` < submit.sh

