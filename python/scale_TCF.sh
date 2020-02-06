#!/usr/bin/bash
#SBATCH -p defq
#SBATCH -t 24:00:00
#SBATCH -N 5

date
mpirun -np 32 python3 taylor_couette_3d.py --re=243.81 --eta=0.875 --m=6 --ar=3 --mesh_1=4 --mesh_2=8
mpirun -np 64 python3 taylor_couette_3d.py --re=243.81 --eta=0.875 --m=6 --ar=3 --mesh_1=8 --mesh_2=8
mpirun -np 128 python3 taylor_couette_3d.py --re=243.81 --eta=0.875 --m=6 --ar=3 --mesh_1=16 --mesh_2=8
date