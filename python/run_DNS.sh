#!/usr/bin/bash
#SBATCH -p defq
#SBATCH -t 24:00:00
#SBATCH -N 3

date

mpirun -np 64 python3 taylor_couette_3d.py --re=700 --eta=0.875 --m=6 --ar=3 --mesh_1=8 --mesh_2=8

sh ~/GQL_TC/python/analysis_scripts/slurm_analysis.sh
date
