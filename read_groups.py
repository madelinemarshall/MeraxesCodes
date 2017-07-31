import sys
import numpy as np
import matplotlib.pyplot as plt
from dragons import nbody
import h5py as h5
import parse_tree_flags

def read_group_life(start_snap, stop_snap, start_index):
    """Read in the life of a single halo.
    Parameters
    ----------
    start_snap : int
        The starting snapshot.
    stop_snap : int
        The stopping redshift (inclusive).
    start_index : int
        The starting index of the halo to trace the evolution of.
    """

    GROUPS_PATH = "/lustre/projects/p070_astro/gpoole/Simulations/Tiamat/catalogs/subfind_{:03d}.catalog_groups_properties"
    HALOS_PATH = "/lustre/projects/p070_astro/gpoole/Simulations/Tiamat/catalogs/subfind_{:03d}.catalog_subgroups_properties"
    TREES_PATH = "/lustre/projects/p070_astro/smutch/input_trees/Tiamat/trees/horizontal_trees_{:03d}.hdf5"
    EXPANSION_FACTOR_PATH = "/lustre/projects/p070_astro/smutch/input_trees/Tiamat/a_list.txt"

    n_snaps = stop_snap - start_snap + 1
    snapshots=list(range(start_snap,stop_snap+1))
    # read in the list of expansion factors for each snapshot of the simulation
    # and convert these to redshifts
    a_list = np.loadtxt(EXPANSION_FACTOR_PATH, dtype=float)
    z_list = 1.0/a_list - 1.0

    ################
    HinG_properties =[]
    HinG_file_offsets = []
    HinG_desc_indices = []
    HinG_IDs = []
    HinG_centrals = []
    HinG_flags = []
    desc_group_indices=[]
    grpindices=np.zeros(n_snaps+1)
    grpindices[0]=start_index
    snap=start_snap
    allhalos = nbody.read_halo_catalog(HALOS_PATH.format(snap)) #Read in all
    #halos from the start_snap snapshot
    with h5.File(TREES_PATH.format(snap), "r") as fd:
        group_indices=fd['trees']['group_index'] #Store group_index for all halos
        #in the start_snap snapshot
        central_indices=fd['trees']['central_index'] #Store central_index for all halos
        #in the start_snap snapshot       
    for ii in range(0,n_snaps):
        HinG_indices=[jj for jj in range(0,len(group_indices)) if \
            group_indices[jj]==grpindices[ii]] #Lists the indices of all halos which have the
            #desired (grpidx) group_index
        if ii==0:
            with h5.File(TREES_PATH.format(snap), "r") as fd:
                entry = fd['trees'][HinG_indices] #Store the halo tree properties
                #for all halos in the desired group
        else:
            entry = new_tree[HinG_indices]
            fd.close()
        HinG_properties.append(allhalos[HinG_indices]) #Store the halo catalog properties
        #for all halos in the desired group

        HinG_file_offsets.append(entry['file_offset'])
        HinG_desc_indices.append(entry['desc_index'])
        HinG_IDs.append(entry['id'])
        HinG_centrals.append(entry['central_index'])
        HinG_flags.append(entry['flags'])

        if ii<n_snaps-1:
            snap=snap+1
            allhalos = nbody.read_halo_catalog(HALOS_PATH.format(snap))
            fd = h5.File(TREES_PATH.format(snap), "r")
            new_tree=fd['trees']
            group_indices=new_tree['group_index']

            desc_group_indices.append(group_indices[HinG_desc_indices[ii][HinG_file_offsets[ii]==1]])
            #Only look at the new group index of descendants with file_offset=1
            grpindices[ii+1]=np.argmax(np.bincount(desc_group_indices[ii]))
            #Finds most common desc_group_index and sets this as the grpidx for this snapshot

    #Find halos that are in each snapshot
    set_of_IDs=[]
    for jj in range(0,len(HinG_IDs)):
        set_of_IDs.append(set(HinG_IDs[jj]))
    halos_in_all_snaps=list(set.intersection(*set_of_IDs)) #Halos that are in all snapshots
    list_of_all_halos=list(set.union(*set_of_IDs)) #All halos

    #Store evolving masses of these halos which are found in each snapshot
    idx_of_halo=np.zeros((len(list_of_all_halos),n_snaps))#,dtype='int32')
    mass_of_halo=np.zeros((len(list_of_all_halos),n_snaps))
    for nn, hh in enumerate(list_of_all_halos):
        for ii in range(0,n_snaps):
            idx_of_halo[nn,ii] = next((jj for jj in range(0,len(HinG_IDs[ii])) if HinG_IDs[ii][jj]==hh),'NaN')
            if ~np.isnan(idx_of_halo[nn,ii]):
                mass_of_halo[nn,ii] = HinG_properties[ii][int(idx_of_halo[nn,ii])]['M_vir']
            else:
                mass_of_halo[nn,ii] = 0

    total_mass=np.sum(mass_of_halo,axis=0)

    np.save('HinG_properties',HinG_properties)
    np.save('HinG_file_offsets',HinG_file_offsets)
    np.save('HinG_desc_indices',HinG_desc_indices)
    np.save('HinG_IDs',HinG_IDs)
    np.save('HinG_centrals',HinG_centrals)
    np.save('HinG_flags',HinG_flags)
    np.save('idx_of_halo',idx_of_halo)
    np.save('mass_of_halo',mass_of_halo)
    np.save('list_of_all_halos',list_of_all_halos)
    np.save('snapshots',snapshots)


    ###

    #Plot all of the snapshotsfla
    plt.plot(snapshots,np.log10(total_mass),'--')   
    maxdiff=np.zeros(len(list_of_all_halos))
    flag_alert=np.zeros([n_snaps,len(list_of_all_halos)])
    secondmaxdiff=np.zeros(len(list_of_all_halos))
    differences=np.zeros([n_snaps,len(list_of_all_halos)])
    for nn, hh in enumerate(list_of_all_halos):
    #differences=np.zeros(n_snaps+1)
        for ii in range(1,n_snaps):
            if (mass_of_halo[nn,ii-1]!=0) & (mass_of_halo[nn,ii]!=0):
                differences[ii,nn]=(np.log10(mass_of_halo[nn,ii])- \
                    np.log10(mass_of_halo[nn,ii-1]))
            else:
                differences[ii,nn]=0
        if ~np.isnan(idx_of_halo[nn,ii]):
            flags=parse_tree_flags.TreeFlags()
            halo_flag=flags.flaglist(HinG_flags[ii][int(idx_of_halo[nn,ii])])    
            if any('TREE_CASE_MERGER_PRIMARY' in flg for flg in halo_flag) & \
                ~any('TREE_CASE_MAIN_PROGENITOR' in flg for flg in halo_flag) & \
                ~any('TREE_CASE_MOST_MASSIVE' in flg for flg in halo_flag) & \
                ~any('TREE_CASE_REMNANT' in flg for flg in halo_flag) & \
                ~any('TREE_CASE_BRIDGED' in flg for flg in halo_flag):
                flag_alert[ii,nn]=1

        maxdiff[nn]=np.nanmax(abs(differences[:,nn]))
        differences_list=abs(differences[:,nn]).tolist()
        differences_list.remove(np.nanmax(differences_list))
        secondmaxdiff[nn]=np.nanmax(differences_list)        
        if (maxdiff[nn]>1) & (secondmaxdiff[nn]>0.75):# & (maxdiff[nn]!=1):
            plt.plot(snapshots,differences[:,nn])#np.log10(mass_of_halo[nn,:]))
 #   plt.savefig('mass_evo_MMBHgroup.pdf')
    plt.show()




    #See if the central is the most massive
    #(doesn't seem to be a good indicator of swapping)
    max_mass_loc=idx_of_halo[np.argmax(mass_of_halo,axis=0)]
    central_loc=[0]*(n_snaps+1)
    max_mass_loc=[0]*(n_snaps+1)
    central_match=[0]*(n_snaps+1)
    for ii in range(0,n_snaps+1):
        central_loc[ii]=HinG_centrals[ii][0]
        max_mass_loc[ii]=int(idx_of_halo[np.argmax(mass_of_halo,axis=0)[ii]][ii])
        central_match[ii]=(central_loc[ii]==max_mass_loc[ii])

        
    #See how many large mass changes there are in each snapshot
    problem=[0]*(n_snaps)
    for ii in range(1,n_snaps):
        problem[ii]=len(np.where(differences[ii,:]>1)[0])
    problem_snaps=np.array(snapshots)*(np.array(problem)>1)

    #Print the snapshots in which there is more than 1 big mass change
    for ii in range(1,n_snaps):
        if (np.nanmax(differences[ii,:])>1) & (np.nanmin(differences[ii,:])<-1):
            print(snapshots[ii])
   
    
    #############

if __name__ == '__main__':
    # parse the values passed at the command line
    start_snap, stop_snap, start_index = int(sys.argv[1]), int(sys.argv[2]), \
    int(sys.argv[3])

    # read in the life
    life = read_group_life(start_snap, stop_snap, start_index)

    # save this to a file so that we don't need to calculate it again!
    np.save("life_snap{:d}-{:d}_index{:d}.npy"
            .format(start_snap, stop_snap, start_index),
            life)

    # you could load this file with the e.g. the following
    #  life = np.load("life_snap4-100_index0.npy")
