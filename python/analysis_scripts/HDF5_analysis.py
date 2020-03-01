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

spn = 1
scalar_vars = ['KE', 'v_tot', 'u_rms', 'v_rms', 'w_rms', 'Re_rms']

# Import .h5 file passed to script 
h = h5py.File(filename,'r')
scales = h['scales']
tasks = h['tasks']
scales_keys = str(scales.keys()).split('\'')
tasks_keys = str(tasks.keys()).split('\'')


for i in range(scales_keys.count(', ')):
	scales_keys.remove(', ')
scales_keys.remove(']>')
scales_keys.remove('r')
scales_keys.remove('theta')
scales_keys.remove('z')
scales_keys.remove('<KeysViewHDF5 [')

for i in range(tasks_keys.count(', ')):
	tasks_keys.remove(', ')
tasks_keys.remove(']>')
tasks_keys.remove('<KeysViewHDF5 [')

print(scales_keys)
print(tasks_keys)

print('-----SCALES-----')

for i in scales_keys:
	print(i)
	dataset = h['scales/'+i]
	print(dataset.shape)

print('-----TASKS-----')

for i in tasks_keys:
	print(i)
	dataset = h['tasks/'+i]
	print(dataset.shape)