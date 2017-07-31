#!/usr/bin/env python

import numpy as np
from dragons import meraxes
import os

cosmo = {'omega_M_0' : 0.308,
'omega_lambda_0' : 0.692,
'omega_b_0' : 0.04839912,
'omega_n_0' : 0.0,
'N_nu' : 0,
'h' : 0.678,
'n' : 0.968,
'sigma_8' : 0.815
}

sim = 'meraxes_on_tiamat/test_1'
fmeraxes = '/home/mmarshal/PhD/results/'+sim+'/meraxes_0.hdf5'
for snapshot in range(78,87):
    gals,sim_props=meraxes.io.read_gals(fmeraxes,\
                                        sim_props=True,\
                                        snapshot=snapshot,\
                                        h=cosmo['h'],quiet=True)
    # You don't want Ghost and You only want massive BHs
    gals = gals[(gals["GhostFlag"]==0) & (gals['BlackHoleMass']>1e-3)]
    for prop in gals.dtype.names:
        gals[prop].astype(np.float).tofile(('/home/mmarshal/PhD/results/massive_bh/Tiamat/%s_z%.2f'%(prop,sim_props['Redshift'])).replace('.','pt')+'.dat')


snapshot = 81
gals,sim_props=meraxes.io.read_gals(fmeraxes,\
                                    sim_props=True,\
                                    snapshot=snapshot,\
                                    h=cosmo['h'],quiet=True)
gals = gals[(gals["GhostFlag"]==0) & (gals['BlackHoleMass']>1e-3)]
biggest_bh_positions = gals['Pos']
biggest_bh_IDs = gals['ID']

gals,sim_props=meraxes.io.read_gals(fmeraxes,\
                                    sim_props=True,\
                                    snapshot=snapshot,\
                                    h=cosmo['h'],quiet=True)
gals = gals[(gals["GhostFlag"]==0)]
Boxsize = sim_props['BoxSize']
width = 0.1 #Mpc
for biggest_bh_ID, biggest_bh_position in zip(biggest_bh_IDs,biggest_bh_positions):
    save_dir = '/home/mmarshal/PhD/environment/Tiamat/%d/'%biggest_bh_ID
    if not os.path.exists(save_dir): os.makedirs(save_dir)
    offsets = gals['Pos'] - biggest_bh_position
    # periodical box
    criteria = ((abs(offsets[:,0])<width) | (abs(offsets[:,0]+Boxsize)<width) | (abs(offsets[:,0]-Boxsize)<width)) &\
               ((abs(offsets[:,1])<width) | (abs(offsets[:,1]+Boxsize)<width) | (abs(offsets[:,1]-Boxsize)<width)) &\
               ((abs(offsets[:,2])<width) | (abs(offsets[:,2]+Boxsize)<width) | (abs(offsets[:,2]-Boxsize)<width))

    # you don't want the massive blackhole again in this catalogue
    nearby_gals = gals[criteria&(gals['ID']!=biggest_bh_ID)]
    for prop in nearby_gals.dtype.names:
        nearby_gals[prop].astype(np.float).tofile('%s/%s.dat'%(save_dir,prop))
