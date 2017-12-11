import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

ID=[60017384810,70016445421,80017150370,80017459020,90016714673,90017453139,100016214064,100016714557,100017474333,
   110015262048,120016389789,130015970208,130016032378,130016214114,130017385486,140016332291,140016389131,
   150015393981,150016663529,160015525754,170015845148,190016663339]#,190017367786,200017221493,230015526011,
   #260015907765,310017213793]

#Load redshift info
fname='/Users/Maddie/GoogleDrive/PhD/Simulations/a_list.txt'
scalefactor=np.fromfile(fname,dtype=float,count=-1,sep="\n")
redshift=1/scalefactor-1
#Set up data structure to hold galaxy info

BHproperties_bulge=np.zeros((13,len(ID),len(redshift)))
addcolumns=['ID','Type','CentralGal','ColdGas','HotGas','StellarMass','StellarDiskScaleLength','Sfr','Rvir','Mvir','Vvir','BulgeStellarMass']
columns=['BlackHoleMass','ID','Type','CentralGal','ColdGas','HotGas','StellarMass','StellarDiskScaleLength','Sfr','Rvir','Mvir','Vvir','BulgeStellarMass']

idx=0
#Loop through the galaxies and import their info
for galID in ID:
    i=0
    for name in columns:
        fname='/Users/Maddie/GoogleDrive/PhD/Simulations/BulgeModelHistory/'+str(int(galID))+'/'+name+'.bin'
        col=np.fromfile(fname,dtype=float,count=-1,sep="")
        BHproperties_bulge[i,idx,:]=col
        i=i+1

    idx=idx+1

for i in range(0,len(ID)):
    y=savgol_filter(np.log10(BHproperties_bulge[0,i,30:]/(BHproperties_bulge[-1,i,30:])), 25, 1)
    #plt.plot(redshift[30:],np.log10(BHproperties_bulge[8,i,30:]/BHproperties_bulge[6,i,30:]))
    if (BHproperties_bulge[-1,i,81]*1e10>10**(10)):
      plt.plot(redshift[30:],y,':r')
    else:
      plt.plot(redshift[30:],y,':b')
y=savgol_filter(np.nanmedian(np.log10(BHproperties_bulge[0,:,30:]/(BHproperties_bulge[-1,:,30:])),axis=0), 25, 1)
plt.plot(redshift[30:],y,'k')
plt.xlabel('Redshift')
plt.gca().invert_xaxis()
plt.xlim([8,2])
plt.ylabel('log(MBH/MBulge)')
plt.tight_layout()
plt.show()
