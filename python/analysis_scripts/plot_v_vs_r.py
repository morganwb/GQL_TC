"""
This script plots the angular momentum of the flow at any given time
and compares the results to the Marcus 2 paper. Could easily be edited
to add more radial variables to plot.

python3 plot_radial_data.py slices_s1.h5

"""
import h5py
import argparse
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (2.6,3.5)

# Parses filename passed to script
parser = argparse.ArgumentParser(description='Passes filename')
parser.add_argument('filename', metavar='Rc', type=str, help='.h5 file to plot radial data')
args = parser.parse_args()
filename = vars(args)['filename']

eta = 0.875
mu = 0.


# Import .h5 file passed to script 
datafile = h5py.File(filename,'r')

r_vs_v_plane_avg = datafile['tasks/v_tot'][-1,0,:,:].mean(axis=0)

r = datafile['scales/r/1.0'][:]


# defined couette flow
A = (1/eta - 1.)*(mu-eta**2)/(1-eta**2)
B = eta*(1-mu)/((1-eta)*(1-eta**2))
v0 = A*r + B/r

# Create plot
plt.title('r vs v plane average')
plt.plot(r,r_vs_v_plane_avg)
#plt.scatter([7.2607308623184235,7.67233758246426],[4.793209944824717,2.6040122651068316],label='Marcus 84 B')
plt.legend()
plt.tight_layout()
plot_file_name = 'radial_profile.png'

plt.savefig(plot_file_name, dpi=300)

