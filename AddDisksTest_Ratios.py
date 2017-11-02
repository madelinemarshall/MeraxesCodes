import FitExponential as FE
import numpy as np
import matplotlib.pyplot as plt

ratr=[1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1,1e1,1e2,1e3,1e4,1e5,1e6,1e7];
ratm=ratr;

rad=1e-3;
m=1e-6;
new_r=(np.ones(len(ratr))*rad)*ratr;
new_mass=(np.ones(len(ratm))*m)*ratm;
ratmonr=np.zeros((len(ratr),len(ratm)))


final_r=np.zeros((len(ratr),len(ratm)))
R_sq=np.zeros((len(ratr),len(ratm)))
avg_ij=np.zeros((len(ratr),len(ratm)))

for ii in range(0,len(ratr)):
  for jj in range(0,len(ratm)):
    final_r[ii,jj],R_sq[ii,jj]=FE.fit_exp(m,rad,new_mass[jj],new_r[ii])
    ratmonr[ii,jj]=(new_mass[jj]/new_r[jj]**2)*(m/rad**2)
    print(final_r[ii,jj],rad,new_r[ii],m,new_mass[ii])
plt.imshow(R_sq, interpolation='nearest');
plt.colorbar()
plt.xlabel('Mass Ratio (1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1,1e1,1e2,1e3,1e4,1e5,1e6,1e7)')
plt.ylabel('Radius Ratio (1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1,1e1,1e2,1e3,1e4,1e5,1e6,1e7)')
plt.show()

