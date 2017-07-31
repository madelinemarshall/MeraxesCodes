import sys
import numpy as np
import matplotlib.pyplot as plt
from dragons import nbody
import h5py as h5
from read_groups import read_group_life

start_snap = 40
stop_snap = 45

GROUPS_PATH = "/lustre/projects/p070_astro/gpoole/Simulations/Tiamat/catalogs/subfind_{:03d}.catalog_groups_properties"
HALOS_PATH = "/lustre/projects/p070_astro/gpoole/Simulations/Tiamat/catalogs/subfind_{:03d}.catalog_subgroups_properties"
TREES_PATH = "/lustre/projects/p070_astro/smutch/input_trees/Tiamat/trees/horizontal_trees_{:03d}.hdf5"
EXPANSION_FACTOR_PATH = "/lustre/projects/p070_astro/smutch/input_trees/Tiamat/a_list.txt"

fname='/home/mmarshal/PhD/id_MBP_MMBH.bin'
#fname='/Users/Maddie/GoogleDrive/PhD/Simulations/history_new_trees/70016445421/id_MBP.bin'
start_id_MBP=np.fromfile(fname,dtype=float,count=-1,sep="")[start_snap]

halo=nbody.read_halo_catalog(HALOS_PATH.format(start_snap))
id_MBPs=halo['id_MBP']

start_halo_index=np.where(id_MBPs==start_id_MBP)
with h5.File(TREES_PATH.format(start_snap), "r") as fd:
             start_index = fd['trees']['group_index'][start_halo_index]

#read_group_life(start_snap,stop_snap,start_index)

########
    HinG_properties=np.load('HinG_properties.npy')
    HinG_file_offsets=np.load('HinG_file_offsets.npy')
    HinG_desc_indices=np.load('HinG_desc_indices.npy')
    HinG_IDs=np.load('HinG_IDs.npy')
    HinG_centrals=np.load('HinG_centrals.npy')
    HinG_flags=np.load('HinG_flags.npy')
    idx_of_halo=np.load('idx_of_halo.npy')
    mass_of_halo=np.load('mass_of_halo.npy')
    list_of_all_halos=np.load('list_of_all_halos.npy')
    n_snaps = stop_snap - start_snap + 1
    snapshots=list(range(start_snap,stop_snap+1))
#    snapshots=np.load('snapshots.npy')
