import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

if __name__=='__main__':
  errs=pd.read_csv('run_errs_2.txt',sep=' ',header=None)
  model_params=pd.read_csv('/home/mmarshal/data_dragons/tune/param_grid.txt',sep=' ',header=None)
  params=np.array(model_params)
  errs_np=errs.as_matrix([0,1,2,3])
  print("Best SMF: {} {} params = {}".format(np.argmin(errs_np[:,0]),errs.get_value(np.argmin(errs_np[:,0]),4),\
    params[np.argmin(errs_np[:,0]),:]))
  print("Best BHMF: {} {} params = {}".format(np.argmin(errs_np[:,1]),errs.get_value(np.argmin(errs_np[:,1]),4),\
    params[np.argmin(errs_np[:,1]),:]))
  print("Best BF: {} {} params = {}".format(np.argmin(errs_np[:,2]),errs.get_value(np.argmin(errs_np[:,2]),4),\
    params[np.argmin(errs_np[:,2]),:]))
  print("Best overall: {} {} params = {}".format(np.argmin(errs_np[:,3]),errs.get_value(np.argmin(errs_np[:,3]),4),\
    params[np.argmin(errs_np[:,3]),:]))
  ## Goodness of fit as function of run number
  plt.plot(errs_np[:,0],'o')
  plt.plot(errs_np[:,1],'o')
  plt.plot(errs_np[:,2],'o')
  plt.plot(errs_np[:,3],'o')
  plt.show()
  ## Goodness of fit of BHMF as function of different parameters
  #plt.plot(params[:,0],errs_np[:,1],'o',label='SFE')
  #plt.plot(params[:,1],errs_np[:,1],'s',label='RME')
  #plt.plot(params[:,2],errs_np[:,1],'^',label='BHGR')
  #plt.plot(params[:,3],errs_np[:,1],'x',label='IDBHGR')
  fig,axes=plt.subplots(2,2)
  plt.title('Goodness of fit of BHMF as function of different parameters')
  axes[0,0].plot(errs_np[:,1],params[:,0],'o',label='SFE')
  axes[1,0].plot(errs_np[:,1],params[:,1],'o',label='RME')
  axes[0,1].plot(errs_np[:,1],params[:,2],'o',label='BHGR')
  axes[1,1].plot(errs_np[:,1],params[:,3],'o',label='IDBHGR')
  plt.legend()
  plt.show()
  ## Goodness of fit of SHMF as function of different parameters
  fig,axes=plt.subplots(2,2)
  plt.title('Goodness of fit of SHMF as function of different parameters')
  axes[0,0].plot(errs_np[:,0],params[:,0],'o',label='SFE')
  axes[1,0].plot(errs_np[:,0],params[:,1],'o',label='RME')
  axes[0,1].plot(errs_np[:,0],params[:,2],'o',label='BHGR')
  axes[1,1].plot(errs_np[:,0],params[:,3],'o',label='IDBHGR')
  plt.legend()
  plt.show()
