"""
Plot planes from joint analysis files. For GQL_TC (Morgan's thesis), 

Usage:
    plot_snapshots_3D.py <files>... [--output=<dir>] [--x_axis=<x_axis>] [--y_axis=<y_axis>] [--slice=<slice>]

Options:
    --output=<dir>  Output directory [default: ./img_snapshots]
    --x_axis=<x_axis>  The basis to be plotted on the x axis [default: theta]
    --y_axis=<y_axis> The basis to be plotted on the y axis [default: z]
    --slice=<slice> Where to slice the third basis [default: 7.5 (slicing through r)]
"""

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

from dedalus.extras import plot_tools


def main(filename, start, count, output):
    """Save plot of specified tasks for given range of analysis writes."""

    # Plot settings
    #with h5py.File(filename, mode='r') as file:
    #    tasks = sorted(file['tasks'].keys())
    tasks = ['u','v','w']
    scale = 4
    dpi = 100
    title_func = lambda sim_time: 't = {:.3f}'.format(sim_time)
    savename_func = lambda write: 'write_{:06}.png'.format(write)
    # Layout
    nrows, ncols = 1,3 
    image = plot_tools.Box(2, 1)
    pad = plot_tools.Frame(0.2, 0.2, 0.1, 0.1)
    margin = plot_tools.Frame(0.3, 0.2, 0.1, 0.1)
    plt.style.use('prl')
    # Create multifigure
    mfig = plot_tools.MultiFigure(nrows, ncols, image, pad, margin, scale)
    fig = mfig.figure
    # Plot writes
    with h5py.File(filename, mode='r') as file:
        for index in range(start, start+count):
            for n, task in enumerate(tasks):
                # Build subfigure axes
                i, j = divmod(n, ncols)
                axes = mfig.add_axes(i, j, [0, 0, 1, 1])
                # Call 3D plotting helper, slicing in time
                dset = file['tasks'][task]
                #plot_axes = [3,1]
                new_slices = [index] + run_slices
                plot_tools.plot_bot(dset, plot_axes, new_slices, axes=axes, title=task, even_scale=True)
            # Add time title
            title = title_func(file['scales/sim_time'][index])
            title_height = 1 - 0.5 * mfig.margin.top / mfig.fig.y
            fig.suptitle(title, x=0.48, y=title_height, ha='left')
            # Save figure
            savename = savename_func(file['scales/write_number'][index])
            savepath = output.joinpath(savename)
            fig.savefig(str(savepath), dpi=dpi)
            fig.clear()
    plt.close(fig)


if __name__ == "__main__":

    import pathlib
    from docopt import docopt
    from dedalus.tools import logging
    from dedalus.tools import post
    from dedalus.tools.parallel import Sync

    args = docopt(__doc__)

    output_path = pathlib.Path(args['--output']).absolute()
    # Create output directory if needed
    with Sync() as sync:
        if sync.comm.rank == 0:
            if not output_path.exists():
                output_path.mkdir()
    try:
        slicer = float(args['--slice'])
    except:
        slicer = 7.5
        pass
    x_axis = str(args['--x_axis'])
    y_axis = str(args['--y_axis'])
    plot_axes = []
    for s in [x_axis,y_axis]:
        if s=='z':
            plot_axes.append(1)
        elif s=='theta':
            plot_axes.append(2)
        elif s=='r':
            plot_axes.append(3)
        else:
            print('Bad axis basis specified. Options are z, theta, r. Exiting...')
    run_slices = []
    for i in range(3):
        if i+1 in plot_axes:
            run_slices.append(slice(None))
        else:
            run_slices.append(slicer)

    post.visit_writes(args['<files>'], main, output=output_path)

