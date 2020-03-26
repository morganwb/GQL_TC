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


for folder in glob.glob(prefix_ast):
    print(folder)
    slices = folder + "/slices/slices_s1.h5"
    print(slices)
    datafile = h5py.File(slices,'r')
    r_vs_v_plane_avg = datafile['tasks/v_tot'][-1,0,:,:].mean(axis=0)
    print(r_vs_v_plane_avg)
    r = datafile['scales/r/1.0'][:]
    v0 = couette(r)
    plt.plot(r,r_vs_v_plane_avg,label=folder[20:])
plt.title('r vs v_total plane average')
plt.legend()
plt.tight_layout()
plt.show()
plot_file_name = 'radial_profile_multiple_prefix.png'
plt.savefig(plot_file_name, dpi=300)
