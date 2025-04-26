#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --exclusive=user
#SBATCH --switches=1
#SBATCH -p besteffort
#SBATCH --job-name pkdgrav

pkdgrav -sz 20  +overwrite ss.par
