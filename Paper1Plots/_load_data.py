from dragons import meraxes


def load_data(filename,snapshot,props):
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
  
  if (props=='All')|(props==None):
    gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
    snapshot=snapshot,\
    h=cosmo['h'],quiet=True)
  else:
    if type(props)==str:
      props=[props]
    if 'GhostFlag' not in props:
      props.append('GhostFlag')
    if 'StellarMass' not in props:
      props.append('StellarMass')
    #props.append('GhostFlag')
    gals=meraxes.io.read_gals(data_folder+filename+meraxes_loc,\
    snapshot=snapshot,props=props,\
    h=cosmo['h'],quiet=True)
  gals=gals[gals['GhostFlag']==0]
  gals=gals[gals['StellarMass']*1e10>1e7]
  if 'BlackHoleMass' in props:
    gals=gals[gals['BlackHoleMass']*1e10>1e5]
  return gals

  '''No cuts
  ##plotSMF

  gals=gals[gals['StellarMass']*1e10>1e7]
  gals=gals[gals['BlackHoleMass']*1e10>1e6]  SAME
  ##MeanBHBulge x2


  gals=gals[gals['BlackHoleMass']*1e10>10**6]
  ##BHGrowth
  ##BHGrowth redshift
  #DiskSize

  gals=gals[gals['BlackHoleMass']*1e10>10**5.5]
  gals=gals[gals['StellarMass']*1e10>1e8]
  return gals
  ##MBHMbulge BT  

  return(gals[(gals["GhostFlag"]==0)&(gals["StellarMass"]>1e-4)])#&(gals["Sfr"]>0)])   OK
  ##MassSfr
  
  gals=gals[gals['StellarMass']*1e10>10**6]
  ##BulgeFraction

  gals=gals[gals['BlackHoleMass']*1e10>10**(6)]
  gals=gals[gals['StellarMass']*1e10>1e8]
  ##MBHM_z0
  ##MBHM_z
  '''
