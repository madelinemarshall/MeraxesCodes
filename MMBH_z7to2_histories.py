import numpy as np
import matplotlib.pyplot as plt

ID=[90017453139,70016445421,280015262506,170012359386,100016714557]

#Load redshift info
fname='/lustre/projects/p070_astro/smutch/input_trees/Tiamat/a_list.txt'
scalefactor=np.fromfile(fname,dtype=float,count=-1,sep="\n")
redshift=1/scalefactor-1
#Set up data structure to hold galaxy info
BHproperties_bulge=np.zeros((13,len(ID),len(redshift)))
addcolumns=['ID','Type','CentralGal','ColdGas','HotGas','StellarMass','StellarDiskScaleLength','Sfr','Rvir','Mvir','Vvir','BulgeStellarMass']
columns=['BHmass','ID','Type','CentralGal','ColdGas','HotGas','StellarMass','StellarDiskScaleLength','Sfr','Rvir','Mvir','Vvir','BulgeStellarMass']

idx=0
#Loop through the galaxies and import their info
for galID in ID:
    #Load data
    fname='/home/mmarshal/PhD/results/MMBHatz7to2/'+str(int(galID))+'/BlackHoleMass.bin'
    BHmass=np.fromfile(fname,dtype=float,count=-1,sep="")
    BHproperties_bulge[0,idx,:]=BHmass
    i=1
    for name in addcolumns:
        fname='/home/mmarshal/PhD/results/MMBHatz7to2/'+str(int(galID))+'/'+name+'.bin'
        col=np.fromfile(fname,dtype=float,count=-1,sep="")
        BHproperties_bulge[i,idx,:]=col
        i=i+1

    idx=idx+1

for i in range(0,len(ID)):
    #y=savgol_filter(np.log10(BHproperties_bulge[0,i,30:]/(BHproperties_bulge[-1,i,30:])), 25, 1)
    y=np.log10(BHproperties_bulge[0,i,25:]*1e10);
    #plt.plot(redshift[30:],np.log10(BHproperties_bulge[8,i,30:]/BHproperties_bulge[6,i,30:]))
    plt.plot(redshift[25:],y)
#y=savgol_filter(np.nanmedian(np.log10(BHproperties_bulge[0,:,30:]/(BHproperties_bulge[-1,:,30:])),axis=0), 25, 1)
#y=np.nanmedian(np.log10(BHproperties_bulge[0,:,30:]/(BHproperties_bulge[-1,:,30:])),axis=0)
#plt.plot(redshift[30:],y,'k')
plt.xlabel('Redshift')
plt.gca().invert_xaxis()
#plt.xlim([8,2])
plt.ylabel('log(MBH)')
plt.legend(['z=7 most massive','z=6 & 5 most massive','z=4 most massive','z=3 most massive','z=2 most massive'],bbox_to_anchor=(1.04,0.5), loc="center left")
plt.tight_layout()
plt.show()

for i in range(0,len(ID)):
    #y=savgol_filter(np.log10(BHproperties_bulge[0,i,30:]/(BHproperties_bulge[-1,i,30:])), 25, 1)
    y=np.log10(BHproperties_bulge[-1,i,25:]*1e10);
    #plt.plot(redshift[30:],np.log10(BHproperties_bulge[8,i,30:]/BHproperties_bulge[6,i,30:]))
    plt.plot(redshift[25:],y)
#y=savgol_filter(np.nanmedian(np.log10(BHproperties_bulge[0,:,30:]/(BHproperties_bulge[-1,:,30:])),axis=0), 25, 1)
#y=np.nanmedian(np.log10(BHproperties_bulge[0,:,30:]/(BHproperties_bulge[-1,:,30:])),axis=0)
#plt.plot(redshift[30:],y,'k')
plt.xlabel('Redshift')
plt.gca().invert_xaxis()
#plt.xlim([8,2])
plt.ylabel('log(MBulge)')
plt.legend(['z=7 most massive','z=6 & 5 most massive','z=4 most massive','z=3 most massive','z=2 most massive'],bbox_to_anchor=(1.04,0.5), loc="center left")
plt.tight_layout()
plt.show()


for i in range(0,len(ID)):
    #y=savgol_filter(np.log10(BHproperties_bulge[0,i,30:]/(BHproperties_bulge[-1,i,30:])), 25, 1)
    y=np.log10(BHproperties_bulge[0,i,25:]/(BHproperties_bulge[-1,i,25:]));
    #plt.plot(redshift[30:],np.log10(BHproperties_bulge[8,i,30:]/BHproperties_bulge[6,i,30:]))
    plt.plot(redshift[25:],y)
#y=savgol_filter(np.nanmedian(np.log10(BHproperties_bulge[0,:,30:]/(BHproperties_bulge[-1,:,30:])),axis=0), 25, 1)
#y=np.nanmedian(np.log10(BHproperties_bulge[0,:,30:]/(BHproperties_bulge[-1,:,30:])),axis=0)
#plt.plot(redshift[30:],y,'k')
plt.xlabel('Redshift')
plt.gca().invert_xaxis()
#plt.xlim([8,2])
plt.ylabel('log(MBH/MBulge)')
plt.legend(['z=7 most massive','z=6 & 5 most massive','z=4 most massive','z=3 most massive','z=2 most massive'],bbox_to_anchor=(1.04,0.5), loc="center left")
plt.tight_layout()
plt.show()
