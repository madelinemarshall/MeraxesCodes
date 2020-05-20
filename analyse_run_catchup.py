import os 
import sys
import numpy as np

for run in np.arange(0,8):
  filename='meraxes_{0:0=3d}'.format(run)
  print("____ Analysing file {} ____".format(filename))
  #snap=snap
  #os.system("python /home/mmarshal/simulation_codes/BulgeFraction_tuning.py third_attempt/output/{}.hdf5".format(filename))
  #os.system("python /home/mmarshal/simulation_codes/plotBHMF_tuning.py fourth_attempt/output/{}.hdf5".format(filename))
  os.system("python /home/mmarshal/simulation_codes/plotQLF_tuning.py output/{}.hdf5".format(filename))
  os.system("python /home/mmarshal/simulation_codes/plotSMF_tuning.py output/{}.hdf5".format(filename))
  os.system("python /home/mmarshal/simulation_codes/BHGrowthWithRedshift_tuning.py output/{}.hdf5".format(filename))
  #os.system("python /home/mmarshal/simulation_codes/MBHMStellarRelation_tuning.py output/{}.hdf5".format(filename))

#for snap in np.arange(0,28):
#  filename='run2/meraxes_{0:0=3d}'.format(snap)
#  print("____ Analysing file {} ____".format(filename))
#  snap=snap+61
  #os.system("python /home/mmarshal/PhD/simulation_codes/BulgeFraction.py {} {}".format(filename,snap))
  #os.system("python /home/mmarshal/PhD/simulation_codes/plotBHMF_tuning.py {} {}".format(filename,snap))
  #os.system("python /home/mmarshal/PhD/simulation_codes/plotSMF_tuning.py {} {}".format(filename,snap))
#  os.system("python /home/mmarshal/PhD/simulation_codes/MBHMStellarRelation_tuning.py {} {}".format(filename,snap))

#for snap in np.arange(0,19):
#  filename='meraxes_{0:0=3d}'.format(snap)
#  print("____ Analysing file {} ____".format(filename))
#  snap=snap+89
  #os.system("python /home/mmarshal/PhD/simulation_codes/BulgeFraction.py {} {}".format(filename,snap))
  #os.system("python /home/mmarshal/PhD/simulation_codes/plotBHMF_tuning.py {} {}".format(filename,snap))
  #os.system("python /home/mmarshal/PhD/simulation_codes/plotSMF_tuning.py {} {}".format(filename,snap))
#  os.system("python /home/mmarshal/PhD/simulation_codes/MBHMStellarRelation_tuning.py {} {}".format(filename,snap))

