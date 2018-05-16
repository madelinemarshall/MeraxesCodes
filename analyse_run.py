import os 
import sys

snap=int(sys.argv[2])
filename='meraxes_{0:0=3d}'.format(snap)
#os.system("python /home/mmarshal/PhD/simulation_codes/BulgeFraction.py {} {}".format(filename,snap))
os.system("python /home/mmarshal/PhD/simulation_codes/plotBHMF_tuning.py {} {}".format(filename,snap+61))
#os.system("python /home/mmarshal/PhD/simulation_codes/plotSMF.py {} {}".format(filename,snap))

