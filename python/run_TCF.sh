#!/usr/bin/bash
#SBATCH -p defq
#SBATCH -t 24:00:00
#SBATCH -N 5

date
#mpirun -np 128 python3 taylor_couette_3d.py --re=1243.81 --eta=0.875 --m=6 --ar=3 --restart=results/TC_3d_re_1243.81_eta_8-7500e-01_Gamma_3/snapshots/snapshots_s1.h5 --mesh_1=16 --mesh_2=8
mpirun -np 128 python3 taylor_couette_3d.py --re=1943.81 --eta=0.875 --m=3 --ar=3 --mesh_1=16 --mesh_2=8
sh ~/GQL_TC/python/analysis_scripts/slurm_analysis.sh
date
