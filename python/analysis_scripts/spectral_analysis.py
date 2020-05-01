"""
This script plots the square of the absolute value of the coefficients
from the 3D Taylor-Couette problem.

Usage:

python3 spectral_analysis.py spectra_s1.h5

"""
import h5py
import argparse
import glob
import numpy as np
import matplotlib.pyplot as plt
from dedalus.extras import plot_tools

# Parses filename passed to script
parser = argparse.ArgumentParser(description='Passes filename')
parser.add_argument('filename', metavar='Rc', type=str, help='Find folders with matching prefix')
args = parser.parse_args()
prefix = vars(args)['filename']
print(prefix)

def abs_sqr(xmesh, ymesh, data):
    newdata = np.abs(data)**2
    return xmesh, ymesh, newdata

folders = sorted(glob.glob(prefix), key=lambda folder: int(folder.split("_")[-1]))

for folder in folders:
    filename = folder + "/spectra/spectra_s1.h5"
    datafile = h5py.File(filename,'r')
    coeff = datafile['tasks/vc']
    print(coeff.shape)
    plot_axes = [1, 2]
    slices = [-1, slice(None), slice(None), 0]
    plot_tools.plot_bot(coeff, plot_axes, slices, func=abs_sqr)
    plt.tight_layout()
    plt.savefig("spectral_plot_" + str(folder.split("_")[-1]) + ".png", dpi=300)



"""
# Create plot
plt.imshow(v_tot)
plt.gca().invert_yaxis()
plt.legend()
plt.tight_layout()
plot_file_name = 'v_tot_timeseries.png'

plt.savefig(plot_file_name, dpi=300)
"""
