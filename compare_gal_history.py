import numpy as np
from dragons import meraxes
import os
import matplotlib.pyplot as plt
import sys

filename=[]
nfiles=len(sys.argv)-1
for ii in range(1,nfiles+1):
  filename.append(sys.argv[ii])

hist={}
for ii in range(0,nfiles):
  hist["{}".format(filename[ii])]=np.load('gal_history_{}.npy'.format(filename[ii]))
  if ('DiskScaleLength' in hist["{}".format(filename[ii])].dtype.names):
    plt.plot(hist["{}".format(filename[ii])]['DiskScaleLength'],label=filename[ii])
  else:
    plt.plot(hist["{}".format(filename[ii])]['StellarDiskScaleLength'],label=filename[ii])

plt.yscale('log')
plt.xlabel('Snapshot')
plt.ylabel('Disk Scale Length')
plt.legend()
plt.title('Evolution of Disk Size for Host of Most Massive BH Host at z=6')
plt.show()


