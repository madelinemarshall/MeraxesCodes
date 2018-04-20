import numpy as np, matplotlib.pyplot as plt
from dragons import meraxes

import sys
import mag_calc as mc

# Set up parameters.
direct="bulges_correctBHMF_tiamat125"
fname = "/home/mmarshal/data_dragons/"+direct+"/output/meraxes.hdf5"

snapList = []
idxList = []

for z in [0.55, 2, 4, 5, 6, 7, 8, 9, 10]:
    snap = meraxes.io.check_for_redshift(fname, z)[0]
    snapList.append(snap)
    gals = meraxes.io.read_gals(fname, snap, 
                                props = ["GhostFlag", "StellarMass"], 
                                h = .678)
    idxList.append(np.where((gals["GhostFlag"] == 0) & (gals["StellarMass"] > 1e-4))[0])
    #idxList.append([0, 1, 2, 11, 12, 13, 101, 102, 103])

bands = mc.HST_filters(["B435", "V606", "i775", "I814", "z850",
                        "Y098", "Y105", "J125", "H160", "3.6"])
mc.composite_spectra(fname, snapList, idxList, h = 0.678, Om0 = .308, 
                     outType = "ph", sedPath = "/home/yqiu/Projects/Magcalc/input/",
                     restFrame = [[1600., 100.], [9000., 100.]],
                     obsBands = bands,
                     prefix = "mags_6",
                     path = "/home/mmarshal/PhD/results/mags_output/"+direct+"/")

##[1600., 100.] - central, width
#sedPath='/lustre/projects/p113_astro/yqiu/magcalc/input/STARBURST99-Salpeter-0.001' 

##DUST:
#MUV=mags.loc[:,['M1600-100','M2000-100']]
#mc.reddening([1600,2000],MUV['M1600-100'],z)
#MUV_dust=MUV+AUV
##To add dust to observed band, convert central wavelength at that z to a rest frame wavelength and use that number in the reddening function.
