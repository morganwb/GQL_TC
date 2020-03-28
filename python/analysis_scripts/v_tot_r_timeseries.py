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

# Parses filename passed to script
parser = argparse.ArgumentParser(description='Passes filename')
parser.add_argument('filename', metavar='Rc', type=str, help='.h5 file to plot radial data')
args = parser.parse_args()
filename = vars(args)['filename']

# Import .h5 file passed to script 
datafile = h5py.File(filename,'r')
v_tot = datafile['tasks/v_tot'][:,0,0,:].T
t = datafile["scales/sim_time"][:]

# Create plot
plt.imshow(v_tot)
plt.gca().invert_yaxis()
plt.legend()
plt.tight_layout()
plot_file_name = 'v_tot_timeseries.png'

plt.savefig(plot_file_name, dpi=300)
