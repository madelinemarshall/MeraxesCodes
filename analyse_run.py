import os 
import sys

filename=str(sys.argv[1])
snap=int(sys.argv[2])
os.system("python /home/mmarshal/PhD/simulation_codes/BulgeFraction.py {} {}".format(filename,snap))
os.system("python /home/mmarshal/PhD/simulation_codes/plotBHMF.py {} {}".format(filename,snap))
os.system("python /home/mmarshal/PhD/simulation_codes/plotSMF.py {} {}".format(filename,snap))

