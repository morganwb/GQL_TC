#!/usr/bin/bash
#SBATCH --partition=defq
#SBATCH --time=10:00:00

date
mpiexec -n 8 python3 TaylorCouetteFlow.py cfg_1.cfg
date
mpiexec -n 1 python3 -m dedalus merge_procs snapshots
date