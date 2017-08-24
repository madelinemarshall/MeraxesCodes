import numpy as np
from dragons import meraxes
import os
import sys

#####Below need to be the same as the one used in the SAM!!!####
cosmo = {'omega_M_0' : 0.308, 
'omega_lambda_0' : 0.692, 
'omega_b_0' : 0.04839912, 
'omega_n_0' : 0.0,
'N_nu' : 0,
'h' : 0.678,
'n' : 0.968,
'sigma_8' : 0.815
}

filename='bulges_crotonSF'

fmeraxes = '/home/mmarshal/data_dragons/'+filename+'/output/meraxes.hdf5'
snapshots, zs, lbts = meraxes.io.read_snaplist(fmeraxes,h=cosmo['h'])
snapshot = 81

gals,sim_props=meraxes.io.read_gals(fmeraxes, sim_props=True,\
                                    snapshot=snapshot,h=cosmo['h'],quiet=True)
props = gals.dtype.names 
ID    = gals['ID']
dr    = 70016445421
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

plt.plot(gal_history['DiskScaleLength'])
plt.show()        
