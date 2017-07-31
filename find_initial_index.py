import numpy as np
from dragons import nbody
import h5py as h5
from read_one_halo_catalog import read_one_halo_catalog


##id_MBP[81&75]=5362546344
#halo=nbody.read_halo_catalog(HALOS_PATH.format(75))
#id_MBPs=halo['id_MBP']
#np.where(id_MBPs==5362546344)


GROUPS_PATH = "/lustre/projects/p070_astro/gpoole/Simulations/Tiamat/catalogs/subfind_{:03d}.catalog_groups_properties"
HALOS_PATH = "/lustre/projects/p070_astro/gpoole/Simulations/Tiamat/catalogs/subfind_{:03d}.catalog_subgroups_properties"
TREES_PATH = "/lustre/projects/p070_astro/smutch/input_trees/Tiamat/trees/horizontal_trees_{:03d}.hdf5"
EXPANSION_FACTOR_PATH = "/lustre/projects/p070_astro/smutch/input_trees/Tiamat/a_list.txt"

start_snap=75
start_index=0
print("Starting from snapshot {:d}, index {:d}...".format(start_snap, start_index))
n_snaps = start_snap

# read in the list of expansion factors for each snapshot of the simulation
# and convert these to redshifts
a_list = np.loadtxt(EXPANSION_FACTOR_PATH, dtype=float)
z_list = 1.0/a_list - 1.0

# define the datatype of what we want to store
dtype = np.dtype([('id', np.int32),
                  ('redshift', np.float64),
                  ('mvir', np.float32)])

# create an empty array with one entry for each snapshot
life = np.zeros(n_snaps, dtype=dtype)

# loop through each snapshot and store the values
index = start_index
skipped=0
for snap in range(start_snap-1,40,-1):
    if skipped!=0:
        skipped-=1
    else:
        found_snap=snap
        with h5.File(TREES_PATH.format(snap), "r") as fd:
            entry = fd['trees']['desc_index']
            offsets = fd['trees']['file_offset']
        progenitor_indices = np.where(entry==index)

        #If progenitor can't be found in the previous snapshot, go back another
        #snapshot and see if it's there
        while (np.size(progenitor_indices[0])==0) & (found_snap>0):
            found_snap=snap-1
            skipped+=1
            with h5.File(TREES_PATH.format(found_snap), "r") as fd:
                entry = fd['trees']['desc_index']
                offsets = fd['trees']['file_offset']
            progenitor_indices = np.where(entry==index)

        masses=np.zeros(len(progenitor_indices[0]))
        for pp in range(0,len(progenitor_indices[0])):
            masses[pp]=read_one_halo_catalog(HALOS_PATH.format(found_snap),progenitor_indices[0][pp])['M_vir']
        index=progenitor_indices[0][np.argmax(masses)]
        offset=offsets[index]
        print(snap,found_snap,index,skipped)
