from dragons import meraxes
import numpy as np
import matplotlib.pyplot as plt
cosmo={}; cosmo['h']=0.678
snapshot=63

gals_merge=meraxes.io.read_gals('/home/mmarshal/data_dragons/simon_premerge/output/meraxes.hdf5',\
#gals_merge=meraxes.io.read_gals('/home/mmarshal/data_dragons/merge_counter/output/meraxes.hdf5',\
#gals_merge=meraxes.io.read_gals('/home/mmarshal/data_dragons/simon/output/meraxes.hdf5',\
#gals_merge=meraxes.io.read_gals('/home/mmarshal/data_dragons/maddie_counter/output/meraxes.hdf5',\
snapshot=snapshot,h=cosmo['h'],quiet=True)
gals_merge=gals_merge[gals_merge['Type']==0]
#gals=meraxes.io.read_gals('/home/mmarshal/data_dragons/maddie_counter/output/meraxes.hdf5',\
gals=meraxes.io.read_gals('/home/mmarshal/data_dragons/simon_451/output/meraxes.hdf5',\
snapshot=snapshot,h=cosmo['h'],quiet=True)
gals=gals[gals['Type']==0]
 
IDs=gals['ID']
gals_merge=gals_merge[np.isin(gals_merge['ID'],IDs)]
gals=gals[np.isin(gals['ID'],gals_merge['ID'])]

#props=['BlackHoleMass','StellarMass','BulgeStellarMass','MergerBulgeStellarMass','BlackHoleMass_ID',\
#'BlackHoleMass_MD','BlackHoleMass_Coalescence']#,'InstCount','MergerCount']
#props=['BlackHoleMass','StellarMass']
#props=['StellarDiskScaleLength','GasDiskScaleLength','VGasDisk','VStellarDisk']#'AMstars','AMcold']
#props=['DiskScaleLength','Vmax','StellarMass','BlackHoleMass','ColdGas']
props=['StellarMass']
fig,axes=plt.subplots(3,3)

ii=0
jj=0
for prop in props:
  if prop in ['AMstars','AMcold']:
    for kk in [0,1,2]:
      diff=(((gals_merge[prop][kk]-gals[prop][kk])/gals[prop][kk]))
      axes[jj,ii].hist(diff[np.abs(diff)>0])
      axes[jj,ii].set_xlabel(prop)
      ii+=1
      if ii==3:
        ii=0
        jj+=1
  else:
    diff=(gals_merge[gals[prop]!=0][prop]-gals[gals[prop]!=0][prop])/gals[gals[prop]!=0][prop]
    axes[jj,ii].hist(((diff[(diff!=0)])))
    print(prop+" max: "+str(max(diff))+" min: "+str(min(diff)))
    axes[jj,ii].set_xlabel(prop)
    axes[jj,ii].set_yscale('log', nonposy='clip')
    ii+=1
    if ii==3:
      ii=0
      jj+=1
#mv=(gals['StellarMass']-gals['BulgeStellarMass'])*gals['VStellarDisk']*gals['StellarDiskScaleLength']
#mv_m=(gals_merge['StellarMass']-gals_merge['BulgeStellarMass'])*gals_merge['VStellarDisk']*gals['StellarDiskScaleLength']
#diff=mv_m-mv
#axes[1,2].hist(np.log10(abs(diff[abs(diff)>0])))
#axes[1,2].set_yscale('log', nonposy='clip')
#print(max(diff))
#print(min(diff))


#bf=gals['BulgeStellarMass'][gals['StellarMass']>0]/gals['StellarMass'][gals['StellarMass']>0]
#bf_m=gals_merge['BulgeStellarMass'][gals['StellarMass']>0]/gals_merge['StellarMass'][gals['StellarMass']>0]
#diff=bf_m-bf
#axes[2,0].hist(diff[abs(diff)>0])
#axes[2,0].set_yscale('log', nonposy='clip')


#diff=(gals_merge['InstCount']-gals['InstCount'])/gals['InstCount']
#axes[2,1].set_xlabel('Inst Count')
#axes[2,1].hist(diff)
#axes[2,1].set_yscale('log', nonposy='clip')

diff=(gals_merge['MergerCount']-gals['MergerCount'])
print(gals_merge['MergerCount'])
print(gals['MergerCount'])

axes[2,2].set_xlabel('Merger Count')
axes[2,2].hist(diff)
axes[2,2].set_yscale('log', nonposy='clip')
print(max(diff))
print(min(diff))

plt.show()
