#!/usr/bin/bash
#SBATCH -p defq
#SBATCH -t 24:00:00
#SBATCH -N 5

date

mpirun -np 128 python3 taylor_couette_3d.py --re=700 --eta=0.875 --m=5 --ar=3 --mesh_1=16 --mesh_2=8 --GQL=True --Lambda_theta=4 --Lambda_z=4

sh ~/GQL_TC/python/analysis_scripts/slurm_analysis.sh

mpirun -np 128 python3 taylor_couette_3d.py --re=700 --eta=0.875 --m=5 --ar=3 --mesh_1=16 --mesh_2=8 --GQL=True --Lambda_z=3 --Lambda_theta=3

sh ~/GQL_TC/python/analysis_scripts/slurm_analysis.sh

date
