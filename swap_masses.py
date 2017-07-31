import sys
import numpy as np
import matplotlib.pyplot as plt
from dragons import nbody
import h5py as h5
import parse_tree_flags

##Load data
start_snap=10
stop_snap=60
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


##Find halos with problem

maxdiff=np.zeros(len(list_of_all_halos))
max_change_at_ii=[0]*n_snaps
min_change_at_ii=[0]*n_snaps
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
        if differences[ii,nn]>0.75:
            flag_alert[ii,nn]=1
    #if ~np.isnan(idx_of_halo[nn,ii]):
        #flags=parse_tree_flags.TreeFlags()
        #halo_flag=flags.flaglist(HinG_flags[ii][int(idx_of_halo[nn,ii])])
        #if any('TREE_CASE_MERGER_PRIMARY' in flg for flg in halo_flag) & \
        #    ~any('TREE_CASE_MAIN_PROGENITOR' in flg for flg in halo_flag) & \
        #    ~any('TREE_CASE_MOST_MASSIVE' in flg for flg in halo_flag) & \
        #    ~any('TREE_CASE_REMNANT' in flg for flg in halo_flag) & \
        #    ~any('TREE_CASE_BRIDGED' in flg for flg in halo_flag):
        #    flag_alert[ii,nn]=1

    maxdiff[nn]=np.nanmax(abs(differences[:,nn]))
    differences_list=abs(differences[:,nn]).tolist()
    differences_list.remove(np.nanmax(differences_list))
    secondmaxdiff[nn]=np.nanmax(differences_list)
    #if (maxdiff[nn]>1) & (secondmaxdiff[nn]>0.75):# & (maxdiff[nn]!=1):
    #    plt.plot(snapshots,differences[:,nn])#np.log10(mass_of_halo[nn,:]))

bad_snapshot=[0]*n_snaps
mass_of_halo_new=np.copy(mass_of_halo)
differences_new=np.copy(differences)
for ii in range(1,n_snaps):
     max_change_at_ii[ii]=np.nanargmax(differences_new[ii,:])
     min_change_at_ii[ii]=np.nanargmin(differences_new[ii,:])
     if (differences_new[ii,max_change_at_ii[ii]]>1)&(differences_new[ii,min_change_at_ii[ii]]<-0.75):
         bad_snapshot[ii]=1
         mass_of_halo_new[max_change_at_ii[ii],ii]=mass_of_halo[min_change_at_ii[ii],ii]
         mass_of_halo_new[min_change_at_ii[ii],ii]=mass_of_halo[max_change_at_ii[ii],ii]
         if ii<n_snaps:
             differences_new[ii+1,nn]=(np.log10(mass_of_halo_new[nn,ii+1])- \
                 np.log10(mass_of_halo_new[nn,ii]))

f, (ax1, ax2) = plt.subplots(2, sharex=True,figsize=(8,12))
for nn in range(0,len(list_of_all_halos)):
    if (maxdiff[nn]>1) & (secondmaxdiff[nn]>0.75):# & (maxdiff[nn]!=1):
        ax1.plot(snapshots,np.log10(mass_of_halo[nn,:]),'o')
        ax2.plot(snapshots,np.log10(mass_of_halo_new[nn,:]),'o')
plt.show()
