import numpy as np
import matplotlib.pyplot as plt
import sys

def fit_exp(mdisk,rdisk,deltam,rd):
  a = mdisk/rdisk**2 
  b = rdisk
  c = deltam/rd**2
  d = rd

  r=np.logspace(np.log(min([b,d]))-1,np.log(max([b,d]))+1,num=30,base=np.exp(1))
 
  #y1=mdisk*(1-np.exp(-r/rdisk)*(1+r/rdisk))
  #y2=deltam*(1-np.exp(-r/rd)*(1+r/rd)) 
  #plt.plot(np.log(r),y1)
  #plt.plot(np.log(r),y2)
  #plt.plot(np.log(r),y1+y2)
  
  y=a*np.exp(-r/b)+c*np.exp(-r/d)
  plt.plot(np.log(r),np.log(a*np.exp(-r/b)))
  #plt.show()
  plt.plot(np.log(r),np.log(c*np.exp(-r/d)))
  #plt.show()
  plt.plot(np.log(r),np.log(y),'--')
  #plt.show()
 
#  sumr=0
#  sumy=0
#  yy=np.zeros(30);
#  rr=np.zeros(30);
#  maxr=max(b,d)
#  rr[0]=np.log(maxr)-1
#  step=((np.log(maxr)+1)/rr[0])*(1/30)
#  for i in range(0,30):
#    yy[i]=a*np.exp(-np.exp(rr[i])/b)+c*np.exp(-np.exp(rr[i])/d);
#    if ~(np.isinf(np.log(yy[i]))):
#      sumr+=np.exp(rr[i])
#      sumy+=np.log(yy[i]);
#    if i<29:
#      rr[i+1]=rr[i]+step
#  meanr=sumr/30;
#  meany=sumy/30;
#  numer=0
#  denom=0

#  for i in range(0,30):
#    if ~(np.isinf(np.log(yy[i]))):
#      numer+=(np.exp(rr[i])-meanr)*(np.log(yy[i])-meany)
#      denom+=(np.exp(rr[i])-meanr)**2
#  print(-denom/numer)
#  r=np.exp(rr)

  meanr=np.mean(r)
  meany=np.mean(np.log(y))
  numer=sum((r-meanr)*(np.log(y)-meany))
  denom=sum((r-meanr)**2)
  slope=numer/denom
  inter=meany-slope*meanr
  f=np.exp(inter)
  g=-1/slope
  
  #fideal=(mdisk+deltam)/(2*np.pi*g**2)
  #plt.plot(np.log(r),np.log(f*np.exp(-r/g)))
  plt.show()
  ymod=f*np.exp(-r/g)  
 
  ssr=sum((np.log(ymod)-meany)**2)
  sst=sum((np.log(y)-meany)**2)
  rsq=ssr/sst
  return g, rsq

if __name__ == '__main__':
  mdisk=float(sys.argv[1])
  rdisk=float(sys.argv[2])
  deltam=float(sys.argv[3])
  rd=float(sys.argv[4])
  g,rsq=fit_exp(mdisk,rdisk,deltam,rd)
  print(g,rsq)
