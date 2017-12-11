import numpy as np
import matplotlib.pyplot as plt

ID=[60017384810,70016445421,80017150370,80017459020,90016714673,90017453139,100016214064,100016714557,100017474333,
   110015262048,120016389789,130015970208,130016032378,130016214114,130017385486,140016332291,140016389131,
   150015393981,150016663529,160015525754,170015845148,190016663339,190017367786,200017221493,230015526011,
   260015907765,310017213793]

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
    fname='/home/mmarshal/PhD/results/bulges/history/'+str(int(galID))+'/BlackHoleMass.bin'
    BHmass=np.fromfile(fname,dtype=float,count=-1,sep="")
    BHproperties_bulge[0,idx,:]=BHmass
    i=1
    for name in addcolumns:
        fname='/home/mmarshal/PhD/results/bulges/history/'+str(int(galID))+'/'+name+'.bin'
        col=np.fromfile(fname,dtype=float,count=-1,sep="")
        BHproperties_bulge[i,idx,:]=col
        i=i+1

    idx=idx+1

for i in range(0,len(ID)):
    f, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(1,5,figsize=(18,3))
    ax1.plot(redshift[30:],np.log10(BHproperties_bulge[4,i,30:]*1e10))
    ax1.set_xlabel('Redshift')
    ax1.invert_xaxis()
    ax1.set_ylabel('Cold Gas Mass')

    ax2.plot(redshift[30:],np.log10(BHproperties_bulge[6,i,30:]*1e10))
    ax2.set_xlabel('Redshift')
    ax2.invert_xaxis()
    ax2.set_ylabel('Stellar Mass')

    ax3.plot(redshift[30:],np.log10(BHproperties_bulge[8,i,30:]))
    ax3.set_xlabel('Redshift')
    ax3.invert_xaxis()
    ax3.set_ylabel('SFR')#SFR')
    ax3.set_title('BH Mass at z=6 = %e' %BHmass[78])

    ax4.plot(redshift[30:][BHproperties_bulge[7,i,30:]*1000<25],np.log10(BHproperties_bulge[7,i,30:][BHproperties_bulge[7,i,30:]*1000<25]*1000))
    ax4.set_xlabel('Redshift')
    ax4.invert_xaxis()
    ax4.set_ylabel('Disk Scale Length (kpc)')

    ax5.plot(redshift[30:],np.true_divide(BHproperties_bulge[-1,i,30:],BHproperties_bulge[6,i,30:]))
    ax5.set_xlabel('Redshift')
    ax5.invert_xaxis()
    ax5.set_ylabel('Bulge To Total Ratio')

    plt.tight_layout()
    plt.show()
