import os 
import sys
import numpy as np

for snap in np.arange(0,61):
  filename='run1/meraxes_{0:0=3d}'.format(snap)
  print("____ Analysing file {} ____".format(filename))
  #snap=snap
  #os.system("python /home/mmarshal/PhD/simulation_codes/BulgeFraction.py {} {}".format(filename,snap))
  #os.system("python /home/mmarshal/PhD/simulation_codes/plotBHMF_tuning.py {} {}".format(filename,snap))
  #os.system("python /home/mmarshal/PhD/simulation_codes/plotSMF_tuning.py {} {}".format(filename,snap))
  os.system("python /home/mmarshal/PhD/simulation_codes/MBHMStellarRelation_tuning.py {} {}".format(filename,snap))

for snap in np.arange(0,28):
  filename='run2/meraxes_{0:0=3d}'.format(snap)
  print("____ Analysing file {} ____".format(filename))
  snap=snap+61
  #os.system("python /home/mmarshal/PhD/simulation_codes/BulgeFraction.py {} {}".format(filename,snap))
  #os.system("python /home/mmarshal/PhD/simulation_codes/plotBHMF_tuning.py {} {}".format(filename,snap))
  #os.system("python /home/mmarshal/PhD/simulation_codes/plotSMF_tuning.py {} {}".format(filename,snap))
  os.system("python /home/mmarshal/PhD/simulation_codes/MBHMStellarRelation_tuning.py {} {}".format(filename,snap))

for snap in np.arange(0,19):
  filename='meraxes_{0:0=3d}'.format(snap)
  print("____ Analysing file {} ____".format(filename))
  snap=snap+89
  #os.system("python /home/mmarshal/PhD/simulation_codes/BulgeFraction.py {} {}".format(filename,snap))
  #os.system("python /home/mmarshal/PhD/simulation_codes/plotBHMF_tuning.py {} {}".format(filename,snap))
  #os.system("python /home/mmarshal/PhD/simulation_codes/plotSMF_tuning.py {} {}".format(filename,snap))
  os.system("python /home/mmarshal/PhD/simulation_codes/MBHMStellarRelation_tuning.py {} {}".format(filename,snap))

