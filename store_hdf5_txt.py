import numpy as np  
import pandas as pd
import sys

def save_hdf5_txt(filename):
  dat=pd.read_hdf(filename+'.hdf5') 
  dat=np.array(dat)
  tt=np.array([dat['i775'],dat['M1600']])
  tt=np.transpose(tt)
  np.savetxt(filename+'.txt',tt,fmt='%4.2f',header='i775 M1600')

if __name__=='__main__':
  filename=str(sys.argv[1])
  save_hdf5_txt(filename)
