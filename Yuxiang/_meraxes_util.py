#!/usr/bin/env python
import numpy as np
from astropy import constants, units,cosmology
import warnings       
#from cosmolopy import magnitudes

solarM2L = 14729390.502926536 #(1*units.Msun*constants.c**2./units.Myr)

def _AGN_luminosity(BHM, accretedHotBHM, accretedColdBHM,delta_t,eta=0.1, quasarVoL = 180, quasarVoLScaling = 0., EddingtonRatio=1.0,split=False,seed=None,return_fduty=False):
    ''' BHM:              blackhole mass in the end (solar mass)
        accretedHotBHM:   accreted blackhole mass during radio mode in this time step (solar mass)
        accretedColdBHM:  accreted blackhole mass during quasar mode in this time step (solar mass)
        delta_t:          time interval of this snapshot (1e6Myr)
        eta:              accretion efficiency
        quasarVoL:        quasar view of line range (degree)
        quasarVoLScaling: quasar view of line range Scaling (for the purpose of BHM dependent quasarVoL)
        EddingtonRatio:   EddingtonRatio (taken from the simulation!)
        '''
    if len(BHM)!=len(accretedHotBHM) or len(BHM)!=len(accretedColdBHM):
        print("BHM and accretedBHM should have the same length!!")
        return -1
    if seed: np.random.seed(seed=seed)
    glow_time = np.random.random(len(BHM)) * delta_t
    angel = np.random.random(len(BHM))
    solid_angle = 1.-np.cos(np.deg2rad(quasarVoL)/2.) # normalized to 2pi
    m0 = BHM - (1.- eta)*accretedColdBHM
    solid_angle *= (m0*np.exp(EddingtonRatio*glow_time/eta/450.)/1e8)**quasarVoLScaling
    solid_angle[solid_angle>1] = 1.0
    flag_undetected = angel>solid_angle

    # quasar mode
    accretion_timeq = np.log(accretedColdBHM/m0 +1.)*eta*450./EddingtonRatio
    QuasarLuminosity = solarM2L*EddingtonRatio*m0*np.exp(EddingtonRatio*glow_time/eta/450.)/450.
    #QuasarLuminosity/=solid_angle
    if not return_fduty:
        QuasarLuminosity[glow_time>accretion_timeq]=0
        QuasarLuminosity[flag_undetected]=0 

    # radio mode
    m0 -= (1.0 - eta)*accretedHotBHM
    accretion_timer = np.log(accretedHotBHM/m0 +1.)*eta*450./EddingtonRatio
    AGNLuminosity = solarM2L*EddingtonRatio*m0*np.exp(EddingtonRatio*glow_time/eta/450.)/450.
    #AGNLuminosity/=solid_angle
    if not return_fduty:
        AGNLuminosity[glow_time>accretion_timer]=0
        AGNLuminosity[flag_undetected]=0 
    if split:
        if return_fduty:
#!!! INCONSISTENCE of the accretion time, now since AGN luminosity is actually not used, so always use quasar mode accretion time to represent the duty cycle!!!!!
            return [QuasarLuminosity,AGNLuminosity,~flag_undetected, accretion_timeq/delta_t] 
        else:
            return [QuasarLuminosity,AGNLuminosity] 
    else:
        if return_fduty:
            return [(QuasarLuminosity+AGNLuminosity), ~flag_undetected, accretion_timeq/delta_t]
        else:
            return (QuasarLuminosity+AGNLuminosity)

def _quasar_luminosity(BHM, accretedBHM, delta_t, eta=0.1, quasarVoL = 180, quasarVoLScaling = 0., EddingtonRatio=1.0,seed=None,return_fduty=False):
    ''' BHM:              blackhole mass in the end (solar mass)
        accretedBHM:      accreted blackhole mass during quasar mode in this time step (solar mass)
        delta_t:          time interval of this snapshot (1e6Myr)
        eta:              accretion efficiency
        quasarVoL:        quasar view of line range (degree)
        quasarVoLScaling: quasar view of line range Scaling (for the purpose of BHM dependent quasarVoL)
        EddingtonRatio:   EddingtonRatio (taken from the simulation!)
        '''
    result = _AGN_luminosity(BHM, np.zeros_like(BHM), accretedBHM,delta_t,eta=eta,quasarVoL=quasarVoL,quasarVoLScaling=quasarVoLScaling,EddingtonRatio=EddingtonRatio,split=True,seed=seed,return_fduty=return_fduty)
    if return_fduty:
        return [result[0], result[2], result[3]]
    else:
        return result[0]

def _quasar_half_luminosity(BHM, accretedColdBHM, delta_t, eta=0.06, EddingtonRatio=1.0, accretedHotBHM=None):
## at half of the accretion time
    ''' BHM:              blackhole mass in the end (solar mass)
        accretedColdBHM:  accreted blackhole mass during quasar mode in this time step (solar mass)
        delta_t:          time interval of this snapshot (1e6Myr)
        eta:              accretion efficiency
        EddingtonRatio:   EddingtonRatio (taken from the simulation!)
        '''
    m0 = BHM - (1.- eta)*accretedColdBHM
    if accretedHotBHM is not None:
        m0 -= (1.-eta)*accretedHotBHM
    accretion_timeq = np.log(accretedColdBHM/m0 +1.)*eta*450./EddingtonRatio
    QuasarLuminosity = solarM2L*EddingtonRatio*m0*(1.+accretedColdBHM/m0)**0.5 /450.
    return QuasarLuminosity, accretion_timeq/delta_t

def _quasar_luminosity_boot(BHM, accretedColdBHM, delta_t,Nboot=100, eta=0.1, EddingtonRatio=1.0,seed=None):
    if seed: np.random.seed(seed=seed)
    glow_time = np.random.random([Nboot,len(BHM)]) * delta_t
    m0 = BHM - (1.- eta)*accretedColdBHM 
    accretion_timeq = np.log(accretedColdBHM/m0 +1.)*eta*450./EddingtonRatio
    QuasarLuminosity = solarM2L*EddingtonRatio*m0*np.exp(EddingtonRatio*glow_time/eta/450.)/450.
    QuasarLuminosity[glow_time>accretion_timeq]=0 

    return QuasarLuminosity

def _AGN_luminosity_boot(BHM, accretedHotBHM, accretedColdBHM, delta_t,Nboot=100, eta=0.1, EddingtonRatio=1.0,seed=None):
    if len(BHM)!=len(accretedHotBHM) or len(BHM)!=len(accretedColdBHM):
        print("BHM and accretedBHM should have the same length!!")
        return -1
    
    if seed: np.random.seed(seed=seed)
    glow_time = np.random.random([Nboot,len(BHM)]) * delta_t
   
    m0=BHM
    QuasarLuminosity=np.zeros((Nboot,len(BHM)))  
    AGNLuminosity=np.zeros((Nboot,len(BHM)))

    #Quasar mode
    m0 -= (1.- eta)*accretedColdBHM 
    accretion_timeQ = np.log(accretedColdBHM/m0 +1.)*eta*450./EddingtonRatio
    QuasarLuminosity = solarM2L*EddingtonRatio*m0*np.exp(EddingtonRatio*glow_time/eta/450.)/450.
    QuasarLuminosity[glow_time>accretion_timeQ]=0 
    QuasarLuminosity[:,accretedColdBHM<=0]=0
      
    #Radio mode
    m0 -= (1.0 - eta)*accretedHotBHM
    accretion_timeR = np.log(accretedHotBHM/m0 +1.)*eta*450./EddingtonRatio
    AGNLuminosity = solarM2L*EddingtonRatio*m0*np.exp(EddingtonRatio*glow_time/eta/450.)/450.
    AGNLuminosity[glow_time>accretion_timeR]=0
    AGNLuminosity[:,accretedHotBHM==0]=0

    return QuasarLuminosity+AGNLuminosity


from scipy import integrate
def _mass_weighted_electron_optical_depth(z_list, xHII, OmegaM, OmegaB, Hubble_h):

    cosmo = cosmology.FlatLambdaCDM(H0=Hubble_h*100.,
                                    Om0=OmegaM,
                                    Ob0=OmegaB)

    thomson_cross_section = 6.652e-25 * units.cm * units.cm
    density_H = 1.88e-7 * cosmo.Ob0 * cosmo.h**2 / 0.022 * units.cm**-3
    density_He = 0.148e-7 * cosmo.Ob0 * cosmo.h**2 / 0.022 * units.cm**-3 

    cosmo_factor = lambda z: constants.c * (1+z)**2 \
            / cosmo.H(z) * thomson_cross_section

    xHII = xHII[::-1]
    z_list = z_list[::-1]

    def d_te_postsim(z):
        if z <= 4:
            return (cosmo_factor(z) * (density_H + 2.0*density_He)).decompose()
        else:
            return (cosmo_factor(z) * (density_H + density_He)).decompose()
    
    def d_te_sim(z, xHII):
        prefac = cosmo_factor(z)
        return (prefac * (density_H*xHII + density_He*xHII)).decompose()

    post_sim_contrib = integrate.quad(d_te_postsim, 0, z_list[0])[0]

    sim_contrib = np.zeros(z_list.size)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sim_contrib = np.array([integrate.simps(d_te_sim(z_list[:ii+1],
                                                        xHII[:ii+1]),
                                                z_list[:ii+1])
                                for ii in xrange(z_list.size)])

    scattering_depth = sim_contrib + post_sim_contrib

    return z_list, scattering_depth

def _Mbh2Ledd(Mass, EddingtonRatio=1.0):
    '''Black hole mass in Msol to Eddington luminosity in Lsol'''
    return EddingtonRatio*Mass*32731.9788953923

def _Lbol2Mbol(Lbol):
    '''bolometric luminosity to bolometric magnitude'''
    return 4.74 - 2.5*np.log10(Lbol)

def _kBvega(Lbol):
    '''bolometric correction for B band in Vega'''
    return 6.25 * (Lbol/1e10)**-0.37 + 9.00 * (Lbol/1e10)**-0.012

def _dkBdLbol(Lbol):
    return -2.3125e-10 * (Lbol/1e10)**-1.37 - 0.108e-10 * (Lbol/1e10)**-1.012

def _dlogkBdlogLbol(Lbol):
    return Lbol/_kBvega(Lbol)*_dkBdLbol(Lbol)

def _Lbol2MBvega(Lbol):
    '''bolometric luminosity to B magnitude in Vega'''
    return _Lbol2Mbol(Lbol) + 2.5*np.log10(_kBvega(Lbol))

def _k15um(Lbol):
    '''bolometric correction for 15 um'''
    return 7.40 * (Lbol/1e10)**-0.37 + 10.66 * (Lbol/1e10)**-0.014

def _Lbol2M15um(Lbol):
    '''bolometric luminosity to 15um in Vega'''
    return _Lbol2Mbol(Lbol) + 2.5*np.log10(_k15um(Lbol))

def _ksoftX(Lbol):
    '''bolometric correction for 0.5-2keV'''
    return 17.87 * (Lbol/1e10)**0.28 + 10.03 * (Lbol/1e10)**-0.020

def _Lbol2MsoftX(Lbol):
    '''bolometric luminosity to 0.5-2keV'''
    return _Lbol2Mbol(Lbol) + 2.5*np.log10(_ksoftX(Lbol))

def _khardX(Lbol):
    '''bolometric correction for 2-10keV'''
    return 10.83 * (Lbol/1e10)**0.28 + 6.08 * (Lbol/1e10)**-0.020

def _Lbol2MhardX(Lbol):
    '''bolometric luminosity to 2-10keV'''
    return _Lbol2Mbol(Lbol) + 2.5*np.log10(_khardX(Lbol))

def _MBvega2MB(MBvega):
    '''convert B magnitude from Vega to AB system'''
    return MBvega - 0.09

def _MB2M1450(MB):
    '''convert B magnitede to UV 1450 magnitude (all AB system) assuming spetral index is 0.44'''
    return MB + 0.524

def _M14502MB(M1450):
    '''convert UV 1450 magnitede to B magnitude (all AB system) assuming spetral index is 0.44'''
    return M1450 - 0.524

def _MBvega2M1450(MBvega):
    '''convert B magnitede in Vega to UV 1450 magnitude assuming spetral index is 0.44'''
    return _MB2M1450(_MBvega2MB(MBvega))

def _Lbol2MB(Lbol):
    '''bolometric luminosity to B magnitude in AB system'''
    return _MBvega2MB(_Lbol2MBvega(Lbol))

def _Lbol2M1450(Lbol):
    '''bolometric luminosity to UV 1450 magnitude in AB system assuming spetral index is 0.44'''
    return _MB2M1450(_Lbol2MB(Lbol))

#def _sumMag(M1,M2):
#    '''add two absolute AB magnitude together'''
#    return magnitudes.magnitude_AB_from_L_nu(magnitudes.L_nu_from_magAB(M1)+magnitudes.L_nu_from_magAB(M2))

import sys
def red_curve(x):
    ''' k calculation. '''
    Rv = 4.05
    if (x > 2.2):
        slope = (2.659*1.04/2.2 - 2.659*1.04/2.19)/(2.2 - 2.19);
        intcept = 2.659*(-1.857 + 1.04/2.2) + Rv;
        k = slope*(x - 2.2) + intcept
    if (x >= 0.63):
        k = 2.569*(-1.857 + 1.04/x) + Rv
    elif (x >= 0.12):
        k = 2.569*(-2.156 + 1.509/x -0.198/x**2 + 0.011/x**3) + Rv
    elif (x >= 0.05):
        slope = (2.659*(1.509/0.125 - 0.198/(0.125*0.125) + 0.011/0.125**3) -
                 2.659*(1.509/0.120 - 0.198/(0.120*0.120) + 0.011/0.120**3))/(0.125 - 0.120);
        intcept = 2.659*(-2.156 + 1.509/0.12 - 0.198/(0.12*0.12) + 0.011/0.12**3) + Rv;
        k = slope*(x - 0.12) + intcept
    else:
        print(x, "is out of range")
        sys.exit(1)
    return k

from scipy.interpolate import interp1d
dust_correction_slope = interp1d([2.5, 3.8, 5.0, 5.9, 7.0, 8.0],[-0.2, -0.11, -0.14, -0.20, -0.20, -0.15], fill_value='extrapolate')
dust_correction_inter = interp1d([2.5, 3.8, 5.0, 5.9, 7.0, 8.0],[-1.7, -1.85, -1.91, -2.00, -2.05, -2.13], fill_value='extrapolate')

def dust_correction(mags, z):
    amean = np.zeros(len(mags))
    betas = dust_correction_slope(z)*(mags + 19.5) + dust_correction_inter(z)
    for i, beta in enumerate(betas):
        a = 4.43 + 1.99*np.random.normal(beta, 0.35, 1000)
        a[a<0.0] = 0.0
        amean[i] = np.mean(a)
    return amean
        
    #if z == 5 or z == 6:
    #    for m in mag:
    #        intercept = {5: -2.05, 6:-2.22}
    #        slope = {5: -0.17, 6:-0.24}
    #        if m < -18.8:
    #            beta = slope[z]*(m + 18.8) + intercept[z]
    #        else:
    #            beta = -0.08*(m + 18.8) + intercept[z]
    #        a = 4.43 + 1.99*np.random.normal(beta, 0.35, 1000)
    #        a[a < 0.0] = 0.0
    #        amean.append(np.mean(a))
    #else:
    #    intercept = {7: -2.05, 8: -2.13, 9: -2.19, 10: -2.26}
    #    slope = {7: -0.20, 8: -0.15, 9: -0.16, 10: -0.16}
    #    for m in mag:
    #        beta = slope[z]*(m + 19.5) + intercept[z]
    #        a = 4.43 + 1.99*np.random.normal(beta, 0.35, 1000)
    #        a[a < 0.0] = 0.0
    #        amean.append(np.mean(a))
    #
    #return np.array(amean)


def dust_extinction(mags, z, bands):
    brightest = -30 #-26
    faintest = -15 #-15
    mo = np.arange(brightest, faintest, 0.01)
    corr = dust_correction(mo, z)
    mi = ["%.1f" %xx for xx in mo - corr]
    faintest = mo[-1] - corr[-1]
    for ii in xrange(len(bands)):
        extinc = dict(zip(mi, corr*red_curve(bands[ii])/red_curve(0.16)))
        for jj in xrange(len(mags)):
            if mags[jj, -1] >faintest:
                mags[jj, ii] += extinc["%.1f" %faintest]
            else:
                mags[jj, ii] += extinc["%.1f" %(mags[jj, -1])]
    return mags

# WRONG FUNCTION!!!
#def dust_extinction_zf_cosmos_20115(mags, z, bands,distance_modulus=0,lambda0=0.551):
#    brightest = -30 #-26
#    faintest = -10 #-15
#    mo = np.arange(brightest, faintest, 0.01)
#    corr = np.zeros(len(mo))+0.5    
#    mi = ["%.1f" %xx for xx in mo - corr]
#    for ii in xrange(len(bands)):
#        extinc = dict(zip(mi, corr*red_curve(bands[ii])/red_curve(lambda0)))
#        for jj in xrange(len(mags)):
#            if mags[jj, -1]-distance_modulus >faintest:
#                mags[jj, ii] += extinc["%.1f" %faintest]
#            else:
#                mags[jj, ii] += extinc["%.1f" %mags[jj, -1]-distance_modulus]
#    return mags

#
#def _L_B_to_M_B(L_B):
#    '''bolometric luminosity to bolometric magnitude'''
#    M_B = 4.74 - 2.5*np.log10(L_B)
#    return M_B
#
#def _M_B_to_L_B(M_B):
#    '''B band magnitude in Vega to B band luminosity in Vega'''
#    L_B = 10**((4.74-M_B)/2.5)
#    return L_B
#
#def _M_B_to_M_AB_B(M_B):
#    '''B band magnitude in Vega to B band magnitude in AB system'''
#    M_AB_B = M_B-0.09
#    return M_AB_B
#
#def _M_AB_B_to_M_B(M_AB_B):
#    '''B band magnitude in AB system to B band magnitude in Vega'''
#    M_B = M_AB_B+0.09
#    return M_B
#
#def _M_AB_B_to_M1450(M_AB_B):
#    '''B band magnitude in AB system to UVmagnitude M1450'''
#    M1450 = M_AB_B+0.524
#    return M1450
#
#def _M1450_to_M_AB_B(M1450):
#    '''UVmagnitude M1450 to B band magnitude in AB system'''
#    M_AB_B = M1450-0.524 
#    return M_AB_B
#
#def _M_B_to_M1450(M_B):
#    '''B band magnitude in Vega to UVmagnitude M1450'''
#    M_AB_B = _M_B_to_M_AB_B(M_B)
#    M1450 = _M_AB_B_to_M1450(M_AB_B)
#    return M1450
#
#def _M1450_to_M_B(M1450):
#    '''UVmagnitude M1450 to B band magnitude in Vega'''
#    M_AB_B = _M1450_to_M_AB_B(M1450)
#    M_B = _M_AB_B_to_M1450(M_AB_B)
#    return M_B
#
#def _M_AB_B_to_L_B(M_AB_B):
#    ''''B band magnitude in AB system to B band luminosity in Vega'''
#    M_B = _M_AB_B_to_M_B(M_AB_B)
#    L_B = _M_B_to_L_B(M_B)
#    return L_B
#
#def _L_B_to_M_AB_B(L_B):
#    ''''B band luminosity in Vega to B band magnitude in AB system'''
#    M_B = _L_B_to_M_B(L_B)
#    M_AB_B = _M_B_to_M_AB_B(M_B)
#    return M_AB_B
#
#def _M_1450_to_L_B(M1450):
#    '''UVmagnitude M1450 to B band luminosity in Vega'''
#    M_B = _M1450_to_M_B(M1450)
#    L_B = _M_B_to_L_B(M_B)
#    return L_B
#
#def _L_B_to_M1450(L_B):
#    '''B band luminosity in Vega to UVmagnitude M1450'''
#    M_B = _L_B_to_M_B(L_B)
#    M1450 = _M_B_to_M1450(M_B)
#    return M145i0
#
#def _bolometric_to_M_B(L_bol):
#    '''bolometric luminosity to B band magnitude in Vega
#    _bolometric_to_M_B(L_bol) # in units of solar luminosity 
#    '''
#    bolometric_correction_B = 6.25 * (L_bol/1e10)**-0.37 + 9.00 * (L_bol/1e10)**-0.012
#    L_B = L_bol/bolometric_correction_B
#    M_B = 4.74 - 2.5*np.log10(L_B)
#    return M_B
#
#def _bolometric_to_M1450(L_bol):
#    '''bolometric luminosity to UV magnitude M1450
#    _bolometric_to_M1450(L_bol) # in units of solar luminosity 
#    '''
#    M_B = _bolometric_to_M_B(L_bol)
#    M1450 = _M_B_to_M1450(M_B)
#    return M1450
