import numpy as np

def MasstoM(Mass):
  LEdd=3.2835e4*Mass
  kb=6.25*(LEdd/1e10)**-0.37+9.00*(LEdd/1e10)**-0.012
  Mbol=4.74-2.5*np.log10(LEdd)
  MB=Mbol+2.5*np.log10(kb)
  MAB=MB-0.09
  M1450=MAB+0.524
  return M1450


def LtoMass(L):
  Mass = np.array(L)/(3.2835e4)
  return Mass

if __name__=='__main__':
  print(MasstoM(5e7))
  print(LtoMass([9.3e13,7e13,2.8e13,0.57e13]))
