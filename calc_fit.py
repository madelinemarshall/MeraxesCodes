import numpy as np
from dragons import meraxes
import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import warnings
warnings.filterwarnings("ignore")

def load_data(filename,meraxes_loc,snapshot,prop,cosmo,bulge=False):
  if bulge:
    gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
        snapshot=snapshot,props=[prop,'GhostFlag','BlackHoleMass','BulgeStellarMass'],\
        h=cosmo['h'],quiet=True)
  else:
    gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
        snapshot=snapshot,props=[prop,'GhostFlag','BlackHoleMass'],\
        h=cosmo['h'],quiet=True)
  gals=gals[(gals["GhostFlag"]==0)]#remove ghosts
  return gals


def calc_SMF(gals,prop,boxwidth):
    maxval=12 
    minval=7.5
    hist, bin_edges = np.histogram(np.log10(gals[prop][gals[prop]>0]*1e10),range=(minval,maxval),bins=30)
    bin_edges=np.array(bin_edges, dtype=np.float128)
    Max=bin_edges[0:-1] + (bin_edges[1]-bin_edges[0])/2.
    SMF=np.log10(hist/(bin_edges[1]-bin_edges[0])/boxwidth**3)
    return SMF


def calc_BHMF(gals,boxwidth):
    maxval=9.7
    minval=4.9
    prop='BlackHoleMass'
    hist, bin_edges = np.histogram(np.log10(gals[prop][gals[prop]>0]*1e10),range=(minval,maxval),bins=24)
    bin_edges=np.array(bin_edges, dtype=np.float128)
    Max=bin_edges[0:-1] + (bin_edges[1]-bin_edges[0])/2.
    #Max=Max[hist>5]
    #hist=hist[hist>5]
    BHMF=np.log10(hist/(bin_edges[1]-bin_edges[0])/boxwidth**3)
    return BHMF


def _h_convertor_y(hubble_h_obs,hubble_h=0.678):
    return (hubble_h/hubble_h_obs)**3.


def calc_disk_frac(gals):
  #Figure 9
  MaxMass=12
  MinMass=10.4
  BinWidth=0.2
  nBins=int(np.ceil((MaxMass-MinMass)/BinWidth))
  TotMinBin=np.zeros(nBins)
  TotNinBin=np.zeros(nBins)
  BulgeMinBin=np.zeros(nBins)
  BinStart=np.zeros(nBins)
  SM=gals['StellarMass']*1e10
  BSM=gals['BulgeStellarMass']*1e10
  for ii in range(0,nBins):
    BinStart[ii]=MinMass+ii*BinWidth
    BinEnd=MinMass+(ii+1)*BinWidth
    TotMinBin[ii]=np.nansum(SM[(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)])
    TotNinBin[ii]=np.size(SM[(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)])
    BulgeMinBin[ii]=np.nansum(BSM[(np.log10(SM)>=BinStart[ii])&(np.log10(SM)<BinEnd)])
    FracInBin=BulgeMinBin/TotMinBin
  BF_model=FracInBin
  return (BF_model)


def calc_disk_frac_SDSS():
#  M=[8.90326,9.102003,9.300743,9.501402,9.703978,9.902666,10.101332,10.301882,10.502413,10.699095,10.901584,11.100241,11.300841,11.499548,11.702157,11.898974]
#  F=[0.8897018,0.88715476,0.88215226,0.8673269,0.8414509,0.7910218,0.7209487,0.6103591,0.48380876,0.3658544,0.26385814,0.18519081,0.120027825,0.085559405,0.08914936,0.09028649]
  M=[10.502413,10.699095,10.901584,11.100241,11.300841,11.499548,11.702157,11.898974]
  F=[0.48380876,0.3658544,0.26385814,0.18519081,0.120027825,0.085559405,0.08914936,0.09028649]
  BF_obs=1-np.array(F)
  return (BF_obs)


if __name__=="__main__":
  #Setup
  cosmo = {'omega_M_0' : 0.308,
  'omega_lambda_0' : 0.692, 'omega_b_0' : 0.04839912,
  'omega_b_0' : 0.04839912,
  'omega_n_0' : 0.0,
  'N_nu' : 0,
  'h' : 0.678,
  'n' : 0.968,
  'sigma_8' : 0.815
  }
  data_folder='/home/mmarshal/data_dragons/'
  meraxes_loc='/output/meraxes.hdf5'
  redshift={63:7,78:6,100:5,116:4,134:3,158:2,194:0.95,213:0.55}
  prop='StellarMass'
  filename='tune'
  meraxes_loc2='/output/'+str(sys.argv[1])+'.hdf5'
  vol=125/cosmo['h']
  filename125=filename
  default='default_fullreion'

  ##SMF
  err=np.zeros(6)
  i=0
  for snapshot in [63,78,100,116,134,158]:
    gals_default=load_data(default,meraxes_loc,snapshot,prop,cosmo)
    gals_bulges=load_data(filename,meraxes_loc2,snapshot,prop,cosmo)
  

    SMF_b=calc_SMF(gals_bulges,prop,vol)
    SMF_d=calc_SMF(gals_default,prop,100)
    diff=np.abs(np.array(SMF_b)-np.array(SMF_d))
    sum_diff=np.sum(diff[np.isfinite(diff)]**2)
    sum_SMF_d=np.sum(SMF_d[np.isfinite(diff)]**2)
    err[i]=np.sqrt(sum_diff)/np.sqrt(sum_SMF_d)
    i+=1  
  err_SMF=(np.max(err))

  ##BHMF
  for snapshot in [213]:
    gals_bulges=load_data(filename,meraxes_loc2,snapshot,prop,cosmo,bulge=True)
    BHMF=calc_BHMF(gals_bulges,vol)
  file = '/home/yqin/dragons/data/bhmf_Shankar2009.txt'
  data0 = np.recfromtxt(file)
  data = data0[data0[:,0]==0.50][:,1:]
  obs_BHMF=np.log10(10**data[:,1]*_h_convertor_y(0.7,cosmo['h']))

  diff=np.abs(np.array(BHMF)-np.array(obs_BHMF))
  sum_diff=np.sum(diff[np.isfinite(diff)]**2)
  sum_obs=np.sum(obs_BHMF[np.isfinite(diff)]**2)
  err_BHMF=np.sqrt(sum_diff)/np.sqrt(sum_obs)

  ##BulgeFraction
  BF_model=calc_disk_frac(gals_bulges)
  BF_obs=calc_disk_frac_SDSS()
  diff=np.abs(np.array(BF_model)-np.array(BF_obs))
  sum_diff=np.sum(diff[np.isfinite(diff)]**2)
  sum_obs=np.sum(BF_obs[np.isfinite(diff)]**2)
  err_BF=np.sqrt(sum_diff)/np.sqrt(sum_obs)
  errs=np.array((err_SMF,err_BHMF,err_BF))


  with open("run_errs.txt", "a+") as myfile:
    myfile.write("\n{} {} {} {} {}".format(err_SMF,err_BHMF,err_BF,np.mean(errs),str(sys.argv[1])))
