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
plt.rcParams["figure.figsize"] = (2.6,3.5)

# Parses filename passed to script
parser = argparse.ArgumentParser(description='Passes filename')
parser.add_argument('filename', metavar='Rc', type=str, help='.h5 file to plot scalar data')
args = parser.parse_args()
filename = vars(args)['filename']

eta = 0.875
mu = 0.


# Import .h5 file passed to script 
datafile = h5py.File(filename,'r')

angular_mom_theta_avg = datafile['tasks/Angular Momentum'][-1,0,:,:].mean(axis=0)

r = datafile['scales/r/1.0'][:]


# defined couette flow
A = (1/eta - 1.)*(mu-eta**2)/(1-eta**2)
B = eta*(1-mu)/((1-eta)*(1-eta**2))
v0 = A*r + B/r

# Create plot
plt.title('Angular Momentum')
plot_set = angular_mom_theta_avg
plt.plot(r,plot_set+(v0*r))
plt.plot(r,(v0*r))
plt.scatter([7.2607308623184235,7.67233758246426],[4.793209944824717,2.6040122651068316],label='Marcus 84 B')
plt.legend()
plt.tight_layout()
plot_file_name = 'radial_profile.png'

plt.savefig(plot_file_name, dpi=300)

