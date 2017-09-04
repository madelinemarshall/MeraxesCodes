import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import sys

filename=[]
#nfiles=len(sys.argv)-1
nfiles=len(sys.argv)-3
prop=str(sys.argv[-2])
redshift=sys.argv[-1]
for ii in range(1,nfiles+1):
  filename.append(sys.argv[ii])

if (redshift==True):
  hist={}
  EXPANSION_FACTOR_PATH = "/lustre/projects/p070_astro/smutch/input_trees/Tiamat/a_list.txt"
  a_list = np.loadtxt(EXPANSION_FACTOR_PATH, dtype=float)
  z_list = 1.0/a_list - 1.0
  for ii in range(0,nfiles):
    hist["{}".format(filename[ii])]=np.load('gal_history_{}.npy'.format(filename[ii]))
    n_snaps=len(hist["{}".format(filename[ii])]['ID'])
    if (prop == 'DiskScaleLength'):
      if ('DiskScaleLength' in hist["{}".format(filename[ii])].dtype.names):
        plt.plot(z_list[0:n_snaps],hist["{}".format(filename[ii])]['DiskScaleLength'],label=filename[ii])
        plt.yscale('log')
      else:
        plt.plot(z_list[0:n_snaps],hist["{}".format(filename[ii])]['StellarDiskScaleLength'],label=filename[ii])
        plt.yscale('log')
    elif (prop == 'BulgeToTotal'):
      if ('BulgeStellarMass' in hist["{}".format(filename[ii])].dtype.names):
        plt.plot(z_list[0:n_snaps],hist["{}".format(filename[ii])]['BulgeStellarMass']/hist["{}".format(filename[ii])]['StellarMass'],label=filename[ii])
    else:
      if (prop in hist["{}".format(filename[ii])].dtype.names):
        plt.plot(z_list[0:n_snaps],hist["{}".format(filename[ii])][prop],label=filename[ii])
        plt.yscale('log')
  plt.gca().invert_xaxis()
  plt.xlim(14)
  plt.xlabel('Redshift')

else:
  hist={}
  for ii in range(0,nfiles):
    hist["{}".format(filename[ii])]=np.load('gal_history_{}.npy'.format(filename[ii]))
    if (prop == 'DiskScaleLength'):
      if ('DiskScaleLength' in hist["{}".format(filename[ii])].dtype.names):
        plt.plot(hist["{}".format(filename[ii])]['DiskScaleLength'],label=filename[ii])
        plt.yscale('log')
      else:
        plt.plot(hist["{}".format(filename[ii])]['StellarDiskScaleLength'],label=filename[ii])
        plt.yscale('log')
    elif (prop == 'BulgeToTotal'):
      if ('BulgeStellarMass' in hist["{}".format(filename[ii])].dtype.names):
        plt.plot(hist["{}".format(filename[ii])]['BulgeStellarMass']/hist["{}".format(filename[ii])]['StellarMass'],label=filename[ii])
    else:
      if (prop in hist["{}".format(filename[ii])].dtype.names):
        plt.plot(hist["{}".format(filename[ii])][prop],label=filename[ii])
        plt.yscale('log')
  plt.xlabel('Snapshot')
plt.ylabel(prop)
plt.legend()
plt.title('Evolution of {} for Host of Most Massive BH Host at z=6'.format(prop))
plt.show()


