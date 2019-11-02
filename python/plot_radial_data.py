"""
This script plots both the real and imaginary components of the
eigenvector corresponding to the grid point with the largest
growth rate. To run, pass it a file as follows:

python3 make_eigenvector_plots.py mri_dataset_845.h5

Make sure the .h5 file contains both the eigenvectors and
eigenvalues. The eigenvalues are used to find the grid point
with the largest growth rate.

"""
import h5py
import argparse
import numpy as np
import matplotlib.pyplot as plt

# Parses filename passed to script
parser = argparse.ArgumentParser(description='Passes filename')
parser.add_argument('filename', metavar='Rc', type=str, help='.h5 file to plot scalar data')
args = parser.parse_args()
filename = vars(args)['filename']

# Import .h5 file passed to script 
datafile = h5py.File(filename,'r')

# Create plot
plt.title('Angular Momentum')
plt.plot(range(108),datafile['tasks/Angular Momentum'][:,0,0,0])
plt.tight_layout()
plot_file_name = 'radial_profile.png'
plt.savefig(plot_file_name, dpi=300)

