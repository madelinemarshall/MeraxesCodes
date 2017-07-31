#!/usr/bin/env python

"""This example script shows how to track the life of a halo starting at a
particular snapshot and going forward in time."""

import sys
import numpy as np
import matplotlib.pyplot as plt
from dragons import nbody
import h5py as h5
from read_one_halo_catalog import read_one_halo_catalog
from tqdm import tqdm

__author__ = "Simon Mutch"
__date__ = "2017-06-21"

def read_life(start_snap, stop_snap, start_index):
    """Read in the life of a single halo.

    Parameters
    ----------
    start_snap : int
        The starting snapshot.
    stop_snap : int
        The stopping redshift (inclusive).
    start_index : int
        The starting index of the halo to trace the evolution of.

    Returns
    -------
    life : ndarray
        The evolution of a selection of the properties for this halo.
    """

    GROUPS_PATH = "/lustre/projects/p070_astro/gpoole/Simulations/Tiamat/catalogs/subfind_{:03d}.catalog_groups_properties"
    HALOS_PATH = "/lustre/projects/p070_astro/gpoole/Simulations/Tiamat/catalogs/subfind_{:03d}.catalog_subgroups_properties"
    TREES_PATH = "/lustre/projects/p070_astro/smutch/input_trees/Tiamat/trees/horizontal_trees_{:03d}.hdf5"
    EXPANSION_FACTOR_PATH = "/lustre/projects/p070_astro/smutch/input_trees/Tiamat/a_list.txt"

    print("Starting from snapshot {:d}, index {:d}...".format(start_snap, start_index))

    n_snaps = stop_snap - start_snap + 1

    # read in the list of expansion factors for each snapshot of the simulation
    # and convert these to redshifts
    a_list = np.loadtxt(EXPANSION_FACTOR_PATH, dtype=float)
    z_list = 1.0/a_list - 1.0

    # define the datatype of what we want to store
    dtype = np.dtype([('id', np.int32),
                      ('redshift', np.float64),
                      ('mvir', np.float32),
                      ('id_MBP', np.float64),
                      ('halo_index',np.float64),
                      ('group_index',np.float64)])

    # create an empty array with one entry for each snapshot
    life = np.zeros(n_snaps, dtype=dtype)

    # loop through each snapshot and store the values
    file_offset = 1
    index = start_index
    for ii, snap in enumerate(range(start_snap, stop_snap+1)):
        # skip this snapshot if the halo is skipped it
        #while (file_offset > 1):# & snap + file_offset < stop_snap:
        if (file_offset > 1):
            file_offset-=1
            continue
        #while (file_offset > 1) & snap + file_offset >= stop_snap:
            #print('break at',snap,file_offset)
            #break

        with h5.File(TREES_PATH.format(snap), "r") as fd:
            entry = fd['trees'][index]

        # Read the halo catalogue files.
        # Note that you could look at the source of this function on the
        # bitbucket website and then adapt it to your needs so that you don't
        # need to read in every catalogue file in a snapshot when you only want
        # one halo that you know the index of. What is done here is incredibly
        # wasteful for the purposes of tracing the history of a single halo.
        # print('Loading snapshot {:d}'.format(snap))
        # halo = nbody.read_halo_catalog(HALOS_PATH.format(snap))[index]
        
        #Loads in data from just one snapshot
        halo=read_one_halo_catalog(HALOS_PATH.format(snap),index)	
        life[ii]['mvir'] = halo['M_vir']
        life[ii]['id'] = entry['id']
        life[ii]['redshift'] = z_list[snap]
        life[ii]['id_MBP'] = halo['id_MBP']
        life[ii]['halo_index'] = index
        life[ii]['group_index'] = entry['group_index'] 
        file_offset = entry['file_offset']
        index = entry['desc_index']
        # if the desc_index is -1 then this halo ceases to exist after this
        # snapshot so break out of the for loop
        if index == -1:
            break

    return life


def plot_mass_evo(life):
    """Plot the evolution of mass for a halo."""

    # create the figure and axis
    fig, ax = plt.subplots(1, 1)

    # plot the log of the mass as a function of redshift
    ax.plot(life[~np.isnan(life['mvir'])]['redshift'], np.log10(life[~np.isnan(life['mvir'])]['mvir']), 'o', lw=4, color='k')

    # reverse the x-axis
    ax.invert_xaxis()

    # set the axis labels
    ax.set_ylabel(r"$\log_{10}(M_{\rm vir}/h)$")
    ax.set_xlabel("redshift")

    # cleanup the whitespace around the figure
    plt.tight_layout()

    # save the plot
    plt.savefig('mass_evo5.pdf')


if __name__ == '__main__':
    # parse the values passed at the command line
    start_snap, stop_snap, start_index = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])

    # read in the life
    life = read_life(start_snap, stop_snap, start_index)

    # save this to a file so that we don't need to calculate it again!
    np.save("life_snap{:d}-{:d}_index{:d}.npy"
            .format(start_snap, stop_snap, start_index),
            life)

    # you could load this file with the e.g. the following
    #life = np.load("life_snap80-94_index0.npy")

    # plot the mass evolution
    plot_mass_evo(life)
