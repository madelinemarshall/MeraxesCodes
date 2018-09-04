#!/usr/bin/env python
import numpy as np
import sys, os
from dragons import meraxes

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
meraxes.io.set_little_h(cosmo['h'])

name      = 'tuned_reion'
redshift_criteria  = float(sys.argv[1])
redshift_start  = float(sys.argv[2])
redshift_end  = float(sys.argv[3])
cut_property = 'BlackHoleMass'
prop_low     = float(sys.argv[4])
prop_high    = float(sys.argv[5])

fmeraxes = '/home/mmarshal/data_dragons/'+name+'/output/meraxes.hdf5'
props    = ("ID", "Mvir", "GhostFlag", "Type", "StellarMass", "Sfr","BlackHoleMass","BulgeStellarMass")
snapshot = meraxes.io.check_for_redshift(fmeraxes,redshift_criteria)[0]
snapshot_start = meraxes.io.check_for_redshift(fmeraxes,redshift_start)[0]
snapshot_end = meraxes.io.check_for_redshift(fmeraxes,redshift_end)[0]
save_dir = '/home/mmarshal/data_dragons/histories/'+name+'/history_by%s/%03d/%s%02d-%02d/'%(cut_property,snapshot,cut_property,prop_low, prop_high)
if not os.path.exists(save_dir): os.makedirs(save_dir)

gals    = meraxes.io.read_gals(fmeraxes,snapshot=snapshot,quiet=True, props=props)
cut_prop_log    = np.log10(gals[cut_property]*1e10)
indices = np.where((cut_prop_log > prop_low) & (cut_prop_log <= prop_high))[0]

print("Tracked Galaxy IDs: {}".format(gals["ID"][indices]))
dtypes  = gals.dtype
totalN  = len(indices)
print("[Snap%03d] total number galaxy in %s %02d-%02d is %d"%(snapshot,cut_property, prop_low, prop_high, totalN))

mvir_history  = np.zeros([totalN, snapshot_end-snapshot_start+1]) + np.nan
star_history  = np.zeros([totalN, snapshot_end-snapshot_start+1]) + np.nan
bh_history  = np.zeros([totalN, snapshot_end-snapshot_start+1]) + np.nan
bulge_history  = np.zeros([totalN, snapshot_end-snapshot_start+1]) + np.nan
bhbulge_history  = np.zeros([totalN, snapshot_end-snapshot_start+1]) + np.nan
sfr_history   = np.zeros([totalN, snapshot_end-snapshot_start+1]) + np.nan
type_history  = np.zeros([totalN, snapshot_end-snapshot_start+1]) + np.nan
ghost_history = np.zeros([totalN, snapshot_end-snapshot_start+1]) + np.nan
ID_history = np.zeros([totalN, snapshot_end-snapshot_start+1]) + np.nan

for snap in range(snapshot, snapshot_start-1, -1):
    flag_hasprogenitor   = np.where(indices!=-1)[0]
    flag_hasnoprogenitor = np.where(indices==-1)[0]
    print("Past: [Snap%03d] remain number of galaxies is %d"%(snap, len(flag_hasprogenitor)))
    if len(flag_hasprogenitor) == 0:
        break

    gals    = meraxes.io.read_gals(fmeraxes,snapshot=snap,quiet=True, props=props)[indices[flag_hasprogenitor]]
    indices = meraxes.io.read_firstprogenitor_indices(fmeraxes, snapshot=snap)[indices]
    indices[flag_hasnoprogenitor]          = -1 
    mvir_history[flag_hasprogenitor,snap-snapshot_start]  = gals["Mvir"]# *1e10
    star_history[flag_hasprogenitor,snap-snapshot_start]  = gals["StellarMass"]# *1e10
    bh_history[flag_hasprogenitor,snap-snapshot_start]  = gals["BlackHoleMass"]# *1e10
    bulge_history[flag_hasprogenitor,snap-snapshot_start]  = gals["BulgeStellarMass"]# *1e10
    bhbulge_history[flag_hasprogenitor,snap-snapshot_start]  = gals["BlackHoleMass"]/gals["BulgeStellarMass"]
    sfr_history[flag_hasprogenitor,snap-snapshot_start]   = gals["Sfr"]
    type_history[flag_hasprogenitor,snap-snapshot_start]  = gals["Type"]
    ghost_history[flag_hasprogenitor,snap-snapshot_start] = gals["GhostFlag"]
    ID_history[flag_hasprogenitor,snap-snapshot_start] = gals["ID"]


gals    = meraxes.io.read_gals(fmeraxes,snapshot=snapshot,quiet=True, props=props)
cut_prop_log    = np.log10(gals[cut_property]*1e10)
indices = np.where((cut_prop_log > prop_low) & (cut_prop_log <= prop_high))[0]

for snap in range(snapshot, snapshot_end, 1):
    flag_hasprogenitor   = np.where(indices!=-1)[0]
    flag_hasnoprogenitor = np.where(indices==-1)[0]
    print("Future: [Snap%03d] number of galaxies is %d"%(snap, len(flag_hasprogenitor)))
    if len(flag_hasprogenitor) == 0:
        break

    indices = meraxes.io.read_descendant_indices(fmeraxes, snapshot=snap)[indices]
    indices[flag_hasnoprogenitor]          = -1 
    gals    = meraxes.io.read_gals(fmeraxes,snapshot=snap+1,quiet=True, props=props)[indices[flag_hasprogenitor]]

    mvir_history[flag_hasprogenitor,snap+1-snapshot_start]  = gals["Mvir"]# *1e10
    star_history[flag_hasprogenitor,snap+1-snapshot_start]  = gals["StellarMass"]# *1e10
    bh_history[flag_hasprogenitor,snap+1-snapshot_start]  = gals["BlackHoleMass"]# *1e10
    bulge_history[flag_hasprogenitor,snap+1-snapshot_start]  = gals["BulgeStellarMass"]# *1e10
    bhbulge_history[flag_hasprogenitor,snap+1-snapshot_start]  = gals["BlackHoleMass"]/gals["BulgeStellarMass"]
    sfr_history[flag_hasprogenitor,snap+1-snapshot_start]   = gals["Sfr"]
    type_history[flag_hasprogenitor,snap+1-snapshot_start]  = gals["Type"]
    ghost_history[flag_hasprogenitor,snap+1-snapshot_start] = gals["GhostFlag"]
    ID_history[flag_hasprogenitor,snap+1-snapshot_start] = gals["ID"]

mvir_history.tofile(save_dir+'Mvir.bin')
star_history.tofile(save_dir+'StellarMass.bin')
bh_history.tofile(save_dir+'BlackHoleMass.bin')
bulge_history.tofile(save_dir+'BulgeStellarMass.bin')
bhbulge_history.tofile(save_dir+'BHBulge.bin')
sfr_history.tofile(save_dir+'Sfr.bin')
type_history.tofile(save_dir+'Type.bin')
ghost_history.tofile(save_dir+'GhostFlag.bin')
ID_history.tofile(save_dir+'ID.bin')
