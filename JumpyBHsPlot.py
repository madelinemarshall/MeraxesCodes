import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['font.size'] = (11)
matplotlib.rcParams['figure.figsize'] = (7.2,5)
#matplotlib.rcParams['font.size'] = (12)
#matplotlib.rcParams['figure.figsize'] = (8.27,6)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

gals_old=np.load('gal_history_default.npy')[30:]
gals_bulges=np.load('gal_history_1102.npy')[30:101]

EXPANSION_FACTOR_PATH = "/lustre/projects/p070_astro/smutch/input_trees/Tiamat/a_list.txt"
a_list = np.loadtxt(EXPANSION_FACTOR_PATH, dtype=float)
z_list = 1.0/a_list - 1.0
z_list=z_list[30:101]

fig, ax=plt.subplots(1,1)
#ax[0].plot(z_list[:101],np.log10(gals_old['BlackHoleMass']/gals_old['StellarMass']))
#ax[0].plot(z_list[:101],np.log10(gals_bulges['BlackHoleMass']/gals_bulges['StellarMass']))
ax.plot(z_list,np.log10(gals_old['DiskScaleLength']*1e3),'orange',linewidth=1.5)
ax.plot(z_list,np.log10(gals_bulges['GasDiskScaleLength']*1e3),'purple',linewidth=1.5)
ax.plot(z_list,np.log10(gals_bulges['Sfr']),'purple',linewidth=1.5)
plt.gca().invert_xaxis()
plt.xlabel('Redshift')
plt.ylabel(r'$\log (R/\mathrm{kpc})$')
plt.legend(['Old Model','New Model'])
plt.show()

