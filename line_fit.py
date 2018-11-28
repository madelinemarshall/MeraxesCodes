import numpy as np
import sys

if __name__=="__main__":
  x1=float(sys.argv[1])
  y1=float(sys.argv[2])
  x2=float(sys.argv[3])
  y2=float(sys.argv[4])
  x3=float(sys.argv[5])

  a=(y1-y2)/(x1-x2)
  b=y1-a*x1
  print('y3 = {}'.format(a*x3+b))
