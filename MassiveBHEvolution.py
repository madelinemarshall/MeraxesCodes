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

sim = 'meraxes_on_tiamat_newtrees/Q_100_000_S_006_050'
fmeraxes = '/home/yqin/dragons/results/'+sim+'/meraxes.hdf5'
gals,sim_props=meraxes.io.read_gals(fmeraxes,\
                                    sim_props=True,\
                                    snapshot=range(10,15),\
                                    h=cosmo['h'],quiet=True)
# You don't want Ghost and You only want massive BHs
gals = gals[(gals["GhostFlag"]==0) & (gals['BlackHoleMass']>1e-3)]
#for prop in gals.dtype.names:
#    gals[prop].astype(np.float).tofile(('/home/mmarshal/PhD/results/massive_bh_evolution/Tiamat/%s_z%.2f'%(prop,sim_props['Redshift'])).replace('.','pt')+'.dat')
biggest_bh_IDs = gals['ID']

for biggest_bh_ID in biggest_bh_IDs:
    save_dir = '/home/yqin/shared/environment/Tiamat/%d/'%biggest_bh_ID
    if not os.path.exists(save_dir): os.makedirs(save_dir)
    for prop in gals.dtype.names:
        gals[prop].astype(np.float).tofile(('%s/%s'%(save_dir,prop))+'.dat')
