import FitExponential as FE
import numpy as np

mdisk=[1e-6,1e-5,1e-4,1e-3,1e-2,1e-1];
rdisk=[1e-6,1e-5,1e-4,1e-3,1e-2,1e-1];
new_mass=[1e-6,1e-5,1e-4,1e-3,1e-2,1e-1];
new_r=[1e-6,1e-5,1e-4,1e-3,1e-2,1e-1];

final_r=np.zeros((len(mdisk),len(rdisk),len(new_mass),len(new_r)))
R_sq=np.zeros((len(mdisk),len(rdisk),len(new_mass),len(new_r)))
avg_ij=np.zeros((len(new_mass),len(new_r)))
avg_jk=np.zeros((len(new_mass),len(new_r)))
avg_kl=np.zeros((len(new_mass),len(new_r)))

for ii in range(0,len(mdisk)):
  for jj in range(0,len(rdisk)):
    for kk in range(0,len(new_mass)):
      for ll in range(0,len(new_r)):
        final_r[ii,jj,kk,ll],R_sq[ii,jj,kk,ll]=FE.fit_exp(mdisk[ii],rdisk[jj],new_mass[kk],new_r[ll])
    avg_ij[ii,jj]=np.mean(R_sq[ii,jj,:,:])

for jj in range(0,len(rdisk)):
  for kk in range(0,len(new_mass)):
    avg_jk[jj,kk]=np.mean(R_sq[:,jj,kk,:])

for kk in range(0,len(new_mass)):
  for ll in range(0,len(new_r)):
    avg_kl[kk,ll]=np.mean(R_sq[:,:,kk,ll])

