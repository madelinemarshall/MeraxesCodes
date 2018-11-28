import sys
import numpy as np
from dragons import meraxes
from _function import _function
from _meraxes_util import _quasar_luminosity_boot,\
_Lbol2MB, _Lbol2M1450
import matplotlib.pyplot as plt

eta = 0.06
alpha_q = 1.57 
alpha_q_op = 0.44
cosmo = {'omega_M_0' : 0.308, 
'omega_lambda_0' : 0.692, 
'omega_b_0' : 0.04839912, 
'omega_n_0' : 0.0,
'N_nu' : 0,
'h' : 0.678,
'n' : 0.968,
'sigma_8' : 0.815
}

#HaloTree = sys.argv[1]
#name     = sys.argv[2]
#savename = sys.argv[3]
redshift = float(sys.argv[1])
band = sys.argv[2]

#sim = HaloTree +'/' + name
fmeraxes = '/home/mmarshal/data_dragons/bulges_update1102_full/output/meraxes.hdf5'
snap = meraxes.io.check_for_redshift(fmeraxes, redshift, tol=0.1)[0]

sim_props = meraxes.io.read_input_params(fmeraxes,h=cosmo['h'])
volume = sim_props['Volume']
EddingtonRatio = sim_props['EddingtonRatio']
quasar_open_angle = sim_props['quasar_open_angle']
observed = 1-np.cos(np.deg2rad(quasar_open_angle)/2.)
props = ("GhostFlag","BlackHoleMass","BlackHoleAccretedColdMass",'dt')
bins = np.linspace(-30,-15,16)
Nboot = 100000
fy = np.zeros([Nboot,len(bins)-1])

gals =meraxes.io.read_gals(fmeraxes,props=props,snapshot=snap,h=cosmo['h'],quiet=True)
gals = gals[(gals["GhostFlag"]==0)]
bh = gals["BlackHoleMass"]*1e10
bh_accreted_cold = gals["BlackHoleAccretedColdMass"]*1e10
delta_t = gals["dt"]
Lbol  = _quasar_luminosity_boot(bh[bh_accreted_cold>0],bh_accreted_cold[bh_accreted_cold>0],delta_t[bh_accreted_cold>0],\
                                Nboot=Nboot,eta=eta,EddingtonRatio=EddingtonRatio)
if band == 'B':
    Magn = _Lbol2MB(Lbol)
elif band == 'UV':
    Magn = _Lbol2M1450(Lbol)
else:
    print('What band?')
    sys.exit()

for ii in range(Nboot):
    f,ed = _function(Magn[ii],volume,bins=bins,weights=np.zeros(len(Magn[ii]))+observed,return_edges=True)
    fy[ii] = f

fx = ed[:-1]
fy = np.sort(fy,axis=0)
#fx = ed
#fy = np.sort(fy,axis=0)
alpha=0.95
y_mean = fy[int((1/2.0)*len(fy))]              
y_err =  fy[int((1-alpha)/2.0*len(fy))] 
y_err2 = fy[int((1+alpha)/2.0*len(fy))]     
plt.plot(fx,y_mean,lw=3,color='orange',linestyle='-',zorder=2)
#plt.fill_between(fx, y_err, y_err2, facecolor='blue', alpha=0.5, zorder=1)
#plt.plot(fx[:-1],np.log10(np.mean(fy,axis=0)))
plt.yscale('log')
plt.xlim([-28.4,-15.3])
plt.ylim([10**-10,2*10**-4])
plt.gca().invert_xaxis()

##Compare with original model
#fmeraxes = '/home/mmarshal/data_dragons/default/output/meraxes.hdf5'
#snap = meraxes.io.check_for_redshift(fmeraxes, redshift, tol=0.1)[0]
#
#sim_props = meraxes.io.read_input_params(fmeraxes,h=cosmo['h'])
#volume = sim_props['Volume']
#EddingtonRatio = sim_props['EddingtonRatio']
#quasar_open_angle = sim_props['quasar_open_angle']
#observed = 1-np.cos(np.deg2rad(quasar_open_angle)/2.)
#props = ("GhostFlag","BlackHoleMass","BlackHoleAccretedColdMass",'dt')
#bins = np.linspace(-30,-15,16)
#Nboot = 100
#fy = np.zeros([Nboot,len(bins)-1])
#
#gals =meraxes.io.read_gals(fmeraxes,props=props,snapshot=snap,h=cosmo['h'],quiet=True)
#gals = gals[(gals["GhostFlag"]==0)]
#bh = gals["BlackHoleMass"]*1e10
#bh_accreted_cold = gals["BlackHoleAccretedColdMass"]*1e10
#delta_t = gals["dt"]
#Lbol  = _quasar_luminosity_boot(bh[bh_accreted_cold>0],bh_accreted_cold[bh_accreted_cold>0],delta_t[bh_accreted_cold>0],\
#                                Nboot=Nboot,eta=eta,EddingtonRatio=EddingtonRatio)
#if band == 'B':
#    Magn = _Lbol2MB(Lbol)
#elif band == 'UV':
#    Magn = _Lbol2M1450(Lbol)
#else:
#    print('What band?')
#    sys.exit()
#
#for ii in range(Nboot):
#    f,ed = _function(Magn[ii],volume,bins=bins,weights=np.zeros(len(Magn[ii]))+observed,return_edges=True)
#    fy[ii] = f
#
#fx = ed[:-1]
#fy = np.sort(fy,axis=0)
#
##fx = ed
##fy = np.sort(fy,axis=0)
#alpha=0.95
#y_mean = fy[int((1/2.0)*len(fy))]              
#y_err =  fy[int((1-alpha)/2.0*len(fy))] 
#y_err2 = fy[int((1+alpha)/2.0*len(fy))]     
#plt.plot(fx,y_mean,lw=1,color='orange',linestyle='-',zorder=2)
#plt.fill_between(fx, y_err, y_err2, facecolor='orange', alpha=0.5, zorder=1)
#
##plt.plot(fx[:-1],np.log10(fy[int(Nboot/2)]))
#
plt.show()
#fx.tofile(save_dir+'%s_x.dat'%(savename))
#fy.tofile(save_dir+'%s.dat'%(savename))
