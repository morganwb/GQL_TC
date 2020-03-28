"""
This script plots the angular momentum of the flow at any given time
and compares the results to the Marcus 2 paper. Could easily be edited
to add more radial variables to plot.

python3 plot_radial_data.py folder_prefix

"""
import h5py
import argparse
import numpy as np
import matplotlib.pyplot as plt

import os
import glob

# Parses filename passed to script
parser = argparse.ArgumentParser(description='Unique folder prefix to identify all runs to plot')
parser.add_argument('folder_prefix', metavar='Rc', type=str, help='.h5 file to plot radial data')
args = parser.parse_args()
prefix = vars(args)['folder_prefix']

prefix_ast = str(prefix) + "*"

eta = 0.875
mu = 0.

# define couette flow
def couette(r):
    A = (1/eta - 1.)*(mu-eta**2)/(1-eta**2)
    B = eta*(1-mu)/((1-eta)*(1-eta**2))
    v0 = A*r + B/r
    return v0
folders = sorted(glob.glob(prefix_ast), key=lambda folder: int(folder.split("_")[-1]))
folders.insert(0,folders[-1])
del folders[-1]

for folder in folders:
    slices = folder + "/slices/slices_s1.h5"
    datafile = h5py.File(slices,'r')
    r_vs_v_plane_avg = datafile['tasks/v_tot'][-1,0,:,:].mean(axis=0)
    r = datafile['scales/r/1.0'][:]
    if "GQL" not in folder:
        plt.plot(r,r_vs_v_plane_avg,label="DNS", linewidth=3)
    else:
        plt.plot(r,r_vs_v_plane_avg,label="GQL, $\Lambda = " + folder.split("_")[-1] + "$")

plt.legend()

plt.xlabel("$r$")
plt.ylabel("Avg. V velocity")
plt.tight_layout()
#plt.style.use('prl')
plt.show()
plot_file_name = 'radial_profile_multiple.png'
plt.savefig(plot_file_name, dpi=300)
