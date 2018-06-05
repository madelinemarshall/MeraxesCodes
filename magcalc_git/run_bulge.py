import numpy as np, matplotlib.pyplot as plt
from dragons import meraxes

import sys
import magcalc as mc

# Set up parameters.
direct="tuned_reion_T125"
fname = "/home/mmarshal/data_dragons/"+direct+"/output/meraxes.hdf5"

snapList = []
idxList = []

for z in [0.57,2,3,4,5,6,7,8,9,10]:
    snap = meraxes.io.check_for_redshift(fname, z)[0]
    snapList.append(snap)
    gals = meraxes.io.read_gals(fname, snap, 
                                props = ["GhostFlag", "StellarMass"], 
                                h = .678)
    idxList.append(np.where((gals["GhostFlag"] == 0) & (gals["StellarMass"] > 1e-4))[0])
    #idxList.append([0, 1, 2, 11, 12, 13, 101, 102, 103])

bands = mc.HST_filters(["B435", "V606", "i775", "I814", "z850",
                        "Y098", "Y105", "J125", "H160", "3.6"])
mc.composite_spectra(fname, snapList, idxList, h = 0.678, Om0 = .308, outType = "ph", 
                     sedPath = '/fred/oz013/yqiu/projects/spectra/S99-Salpeter-0.001',
                     restBands = [[1600., 100.], [9000., 100.]],
                     obsBands = bands,
                     prefix = "mags_6",
                     outPath = "/home/mmarshal/results/mags_output/"+direct+"/")

##[1600., 100.] - central, width

###DUST:
#z=2
#MUV=mags.loc[:,['M1600-100','M2000-100']]
#mc.reddening([1600,2000],MUV['M1600-100'],z)
#MUV_dust=MUV+AUV
####To add dust to observed band, convert central wavelength at that z to a rest frame wavelength and use that number in the reddening function.
