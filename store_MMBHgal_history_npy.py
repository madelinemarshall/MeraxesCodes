## Determines and stores the history of the host of the most massive BH at z=6
# (ID 70016445421 [adjustable]).
# Input: name of the run folder which holds the Meraxes output to be used
# Output: gal_history_[input_name].npy

import numpy as np
from dragons import meraxes
import os
import sys

cosmo = {'omega_M_0' : 0.308, 
'omega_lambda_0' : 0.692, 
'omega_b_0' : 0.04839912, 
'omega_n_0' : 0.0,
'N_nu' : 0,
'h' : 0.678,
'n' : 0.968,
'sigma_8' : 0.815
}

filename=str(sys.argv[1])

fmeraxes = '/home/mmarshal/data_dragons/'+filename+'/output/meraxes.hdf5'
snapshots, zs, lbts = meraxes.io.read_snaplist(fmeraxes,h=cosmo['h'])
snapshot = 19

gals,sim_props=meraxes.io.read_gals(fmeraxes, sim_props=True,\
                                    snapshot=snapshot,h=cosmo['h'],quiet=True)
props = gals.dtype.names 
ID    = gals['ID']
if len(sys.argv)==2:
  dr = 70016445421
else:
  dr = int(sys.argv[2])
index = np.where(ID==dr)

gal_history = np.zeros(len(snapshots),dtype=gals.dtype)
gal_history[snapshot] = gals[index]

# history
snap = snapshot
tmp = index[0]
while tmp>=0 and snap>0:
    firstprogenitor = meraxes.io.read_firstprogenitor_indices(fmeraxes, snap)[tmp]
    gal_history[snap-1] = meraxes.io.read_gals(fmeraxes,props=props,h =cosmo['h'],\
                                snapshot=snap-1,quiet=True)[firstprogenitor]
    tmp = firstprogenitor
    snap-=1
    
#future
snap = snapshot
tmp = index[0]
while tmp>=0 and snap<len(snapshots)-1:
    descendant = meraxes.io.read_descendant_indices(fmeraxes, snap)[tmp]
    gal_history[snap+1]=meraxes.io.read_gals(fmeraxes,props=props,h =cosmo['h'],\
                                snapshot=snap+1,quiet=True)[descendant]
    tmp = descendant
    snap+=1    

np.save('gal_history'+filename,gal_history)

#save_dir = '/home/mmarshal/PhD/results/bulges'
#if not os.path.exists(save_dir+'/history/%d/'%(dr)):
#    os.makedirs(save_dir+'/history/%d/'%(dr))
#for prop in props:
#    gal_history[prop].astype(np.float).tofile(save_dir+'/history/%d/%s.bin'%(dr,prop))

