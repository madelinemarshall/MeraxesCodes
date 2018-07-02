import numpy as np, matplotlib.pyplot as plt
from dragons import meraxes

import sys
import magcalc as mc

# Set up parameters.
direct="tuned_reion_T125"
fname = "/home/mmarshal/data_dragons/"+direct+"/output/meraxes.hdf5"

snapList = []
idxList = []

filt1=np.loadtxt('wfc_F625W.dat')
filt2=np.loadtxt('F356W.txt',skiprows=1)
filt3=np.loadtxt('F444W.txt',skiprows=1)
filt4=np.loadtxt('F150W.txt',skiprows=1)
filt5=np.loadtxt('F200W.txt',skiprows=1)
filt6=np.loadtxt('F115W.txt',skiprows=1)
filt7=np.loadtxt('F090W.txt',skiprows=1)
filt8=np.loadtxt('F070W.txt',skiprows=1)
filt9=np.loadtxt('F277W.txt',skiprows=1)
filt10=np.loadtxt('MIRI_F560W.dat')
filt11=np.loadtxt('SLOAN_SDSS_g.dat')
filt12=np.loadtxt('SLOAN_SDSS_i.dat')
filt13=np.loadtxt('SLOAN_SDSS_r.dat')
filt2[:,0]*=1e4
filt3[:,0]*=1e4
filt4[:,0]*=1e4
filt5[:,0]*=1e4
filt6[:,0]*=1e4
filt7[:,0]*=1e4
filt8[:,0]*=1e4
filt9[:,0]*=1e4

for snap in range(213,251):
    snapList.append(snap)
    gals = meraxes.io.read_gals(fname, snap, 
                                props = ["GhostFlag", "StellarMass"], 
                                h = .678)
    idxList.append(np.where((gals["GhostFlag"] == 0) & (gals["StellarMass"] > 1e-4))[0])
    #idxList.append([0, 1, 2, 11, 12, 13, 101, 102, 103])

bands = mc.HST_filters(["B435", "V606", "i775", "I814", "z850",
                        "Y098", "Y105", "J125", "H160", "3.6"])
bands.append(["F625W",filt1.T])
bands.append(["JWST_F356W",filt2.T])
bands.append(["JWST_F444W",filt3.T])
bands.append(["JWST_F150W",filt4.T])
bands.append(["JWST_F200W",filt5.T])
bands.append(["JWST_F115W",filt6.T])
bands.append(["JWST_F090W",filt7.T])
bands.append(["JWST_F070W",filt8.T])
bands.append(["JWST_F277W",filt9.T])
bands.append(["JWST_F560W",filt10.T])
bands.append(['SDSSg',filt11.T])
bands.append(['SDSSi',filt12.T])
bands.append(['SDSSr',filt13.T])

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
