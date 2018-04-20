import os 
import sys
import numpy as np

name={26:50,27:51}
for snap in [26,27]:#np.arange(0,56):
  filename='run2/meraxes_{0:0=3d}'.format(snap)
  snap=name[snap]
  os.system("python /home/mmarshal/PhD/simulation_codes/BulgeFraction.py {} {}".format(filename,snap))
  os.system("python /home/mmarshal/PhD/simulation_codes/plotBHMF.py {} {}".format(filename,snap))
  os.system("python /home/mmarshal/PhD/simulation_codes/plotSMF.py {} {}".format(filename,snap))

