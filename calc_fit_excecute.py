import os 
import sys
import numpy as np
import warnings
warnings.filterwarnings("ignore")

for snap in [26,27]:#np.arange(0,56):
  filename='meraxes_{0:0=3d}'.format(snap)
  print("_________ {} __________".format(filename))
  os.system("python /home/mmarshal/PhD/simulation_codes/calc_fit.py {}".format(filename)) 
