#!/usr/bin/env python
import numpy as np

def _function_fduty(data, quasarVoL, duty_cycle_acc, volume, bins, poisson_uncert=False):

    duty_cycle_obs  = 1.-np.cos(np.deg2rad(quasarVoL)/2.)
    digitized = np.digitize(data,bins)
    histogram = np.array([np.sum(duty_cycle_acc[digitized == i]) for i in range(1, len(bins))])
    uncert    = np.sqrt(duty_cycle_obs*histogram)

    widths    = bins[1:] - bins[:-1]
    centers   = (bins[1:] + bins[:-1])/2.0
    histogram = duty_cycle_obs * histogram.astype(float) / (volume * widths)

    if not poisson_uncert:
        return  np.dstack((centers, histogram)).squeeze()
    else:
        uncert = duty_cycle_obs * uncert/(volume * widths)
        return  np.dstack((centers, histogram, uncert)).squeeze()

def _function(data, volume, bins, range=None, poisson_uncert=False,
                return_edges=False, weights=None, density=None):
    if weights is not None:
        weights = weights[(~np.isnan(data)) & (~np.isinf(data))] 
    data = data[(~np.isnan(data)) & (~np.isinf(data))]
    if (range is not None and (bins in ['blocks',
                                        'knuth', 'knuths',
                                        'scott', 'scotts',
                                        'freedman', 'freedmans'])):
        if weights is not None:
            weights = weights[(data >=range[0]) & (data <= range[1])] 
        data = data[(data >=range[0]) & (data <= range[1])]

    if isinstance(bins, str):
        log.info("Calculating bin widths using `%s' method..." % bins)
        if bins in ['blocks']:
            bins = bayesian_blocks(data)
        elif bins in ['knuth', 'knuths']:
            dm, bins = knuth_bin_width(data, True)
        elif bins in ['scott', 'scotts']:
            dm, bins = scotts_bin_width(data, True)
        elif bins in ['freedman', 'freedmans']:
            dm, bins = freedman_bin_width(data, True)
        else:
            raise ValueError("unrecognized bin code: '%s'" % bins)
        log.info("...done")

    vals, edges = np.histogram(data, bins, range,  weights=weights, density=density)
    width = edges[1]-edges[0]
    radius = width/2.0
    centers = edges[:-1]+radius
    if poisson_uncert:
        number, edges = np.histogram(data, bins, range,  weights=None, density=density)
        uncert = vals/np.sqrt(number.astype(float))

    vals = vals.astype(float) / (volume * width)

#    if not poisson_uncert:
#        results = np.dstack((centers, vals)).squeeze()
#    else:
#        uncert /= (volume * width)
#        results = np.dstack((centers, vals, uncert)).squeeze()
    results=vals

    if not return_edges:
        return results
    else:
        return results, edges


def _cutinf(x,y,z=None,silence=1):
    if len(x)!=len(y):
        print('x and y have different length!')
        return
    else:
        if silence == 0: print( 'Initially, number is:', len(x))
        if z == None:
            xx=x[(~np.isinf(x)) & (~np.isinf(y))]
            yy=y[(~np.isinf(x)) & (~np.isinf(y))]
            if silence == 0: print('Now, number is:', len(xx))
        else:
            xx=x[(~np.isinf(x)) & (~np.isinf(y)) & (~np.isinf(z))]
            yy=y[(~np.isinf(x)) & (~np.isinf(y)) & (~np.isinf(z))]
            if silence == 0: print('Now, number is:', len(xx))
        return xx,yy

def find_confidence_interval(x, pdf, confidence_level):
    return pdf[pdf > x].sum() - confidence_level

def _density_contour(xdata, ydata, nbins_x=100, nbins_y=100, ax=None, vmin=None, vmax=None, nocontour=False,label=None, **contour_kwargs):
    import matplotlib.pyplot as plt
    import scipy.optimize as so
    from matplotlib.colors import LogNorm
    if ax == None:
        H, xedges, yedges, img = plt.hist2d(xdata, ydata, bins=(nbins_x,nbins_y), vmin=vmin, vmax=vmax, normed=True, norm=LogNorm())
    else:
        H, xedges, yedges, img = ax.hist2d(xdata, ydata, bins=(nbins_x,nbins_y), vmin=vmin, vmax=vmax, normed=True, norm=LogNorm())
    if nocontour: return
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))

    pdf = (H*(x_bin_sizes*y_bin_sizes))

    point_one = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.10))
    point_two = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.20))
    point_three = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.30))
    point_four = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.40))
    point_five = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.50))
    point_six = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.60))
    one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.68))
    two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.95))
    three_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.99))
    levels = [point_one, point_four, one_sigma, two_sigma, three_sigma]
    labels = ['0.1', '0.2', '1_sigma', '2_sigma', '3_sigma']

    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    Z = pdf.T

    if ax == None:
        contour = plt.contour(X, Y, Z, levels=levels, origin="lower",label=labels, **contour_kwargs)
        plt.clabel(contour, iline=1, fontsize=10)
    else:
        contour = ax.contour(X, Y, Z, levels=levels, origin="lower",label=labels, **contour_kwargs)
        ax.clabel(contour, iline=1, fontsize=10)
    return contour

def _sample_select(data, length=100, sample = None):
    import random
    if len(data) == 1:
        x = np.zeros(length)
    else:
        x = np.zeros([len(data),length])
    if sample == None:
        sample = sorted(random.sample(range(len(data[0])), length))
    x[0] = data[0][sample]
    x[1] = data[1][sample]
    return x

def _sample_select_pandas(data, length=100, sample = None):
    import random
    x = np.zeros(length)
    if sample == None:
        sample = sorted(random.sample(data.keys(), length))
    x[j] = data[sample]
    return x

def _sample_select_2d(xdata, ydata, length):
    import random
    x = np.zeros(length)
    y = np.zeros(length)
    try:
        sample = sorted(random.sample(xdata.keys(), length))
    except AttributeError:
        sample = sorted(random.sample(range(len(xdata)), length))
    x = xdata[sample]
    y = ydata[sample]
    return x, y

def _run_hdf5totiamat(sim='REF_EFF_L010N0128', snap_start=0, snap_stop=103,dir='/home/yqin/hydro/tiamat/',move=1):
    from hdf5tobin import subfind
    import os
    
    if move ==1:
        print('transfer subfind HDF5 files of ',sim,'to binary files and save to ',dir)
    else:
        print('transfer subfind HDF5 files of ',sim,'to binary files')
    #prefix = '/lustre/projects/p071_swin/yqin/Smaug/'
    prefix = '/home/yqin/hydro/run_gadget/smaug/'

    for i in range(snap_start,snap_stop+1):
        fname = prefix + sim + '/data/subhalos_%03d/subhalo_%03d'%(i,i)
        result = subfind(fname,subdivide=False,silent=1)

        if(not os.path.exists(dir + sim)): os.makedirs(dir + sim)
        if(not os.path.exists(dir + sim)): os.makedirs(dir + sim)

    if move == 1:
        src = prefix + sim + '/catalogs'
        dst = dir + sim + '/catalogs'
        os.rename(src,dst)

        src = prefix + sim + '/halos'
        dst = dir + sim + '/halos'
        os.rename(src,dst)

    return result


def _run_snapshottobinary(sim='DMONLY_L010N0512', snap_start=0, snap_stop=103,dir='/home/yqin/hydro/tiamat/',move=1):
    from hdf5tobinary import snapshottobinary
    import os

    if move ==1:
        print('transfer snapshots HDF5 files of ',sim,'to binary files and save to ',dir)
    else:
        print('transfer snapshots HDF5 files of ',sim,'to binary files')
    prefix = '/lustre/projects/p071_swin/yqin/Smaug/'
    #prefix = '/home/yqin/hydro/run_gadget/dragons/'

    for i in range(snap_start,snap_stop+1):
        fname = prefix + sim + '/data/snapshot_%03d/snap_%03d.0.hdf5'%(i,i)
        result = snapshottobinary(fname)

        if(not os.path.exists(dir + sim)): os.makedirs(dir + sim)
        if(not os.path.exists(dir + sim)): os.makedirs(dir + sim)

    if move == 1:
        src = prefix + sim + '/snapshots'
        dst = dir + sim + '/snapshots'
        os.rename(src,dst)

    return result

def _sim_match(file1, file2, Nlim=1e10):
    import struct
    f = file(file1,"rb")
    mbr = f.read()
    header1 = struct.unpack('i'*4,mbr[0:16])
    index1 = np.array(struct.unpack('i'*header1[2],mbr[16:16+4*header1[2]]),dtype=np.int)
    #score1 = struct.unpack('f'*header1[2],mbr[16+4*3*header1[2]:16+4*4*header1[2]])
    f.close()

    f = file(file2,"rb")
    mbr = f.read()
    header2 = struct.unpack('i'*4,mbr[0:16])
    index2 = np.array(struct.unpack('i'*header2[2],mbr[16:16+4*header2[2]]),dtype=np.int)
    #score2 = struct.unpack('f'*header2[2],mbr[16+4*3*header2[2]:16+4*4*header2[2]])
    f.close()

    match = np.zeros(min(Nlim,len(index2)),dtype=np.int)

    for i in range(len(match)):
        candidates = np.where(index1 == i)[0]
        if len(candidates) == 0: match[i] = -1
        else:
            if index2[i] in candidates: match[i] = index2[i]
            else: match[i] = -2

    return match

def _diagonal(ax, fit = None, linestyle = '-', diag = False, zoom = False, color ='red', lw=3):
    if diag == True:
        if zoom == True:
            range = (max((ax.get_xlim()[0],ax.get_ylim()[0])),min((ax.get_xlim()[1],ax.get_ylim()[1])))
        else:
            range = (min((ax.get_xlim()[0],ax.get_ylim()[0])),max((ax.get_xlim()[1],ax.get_ylim()[1])))
        ax.set_xlim(range)
        ax.set_ylim(range)
    else:
        ax.set_xlim(ax.get_xlim())
        ax.set_ylim(ax.get_ylim())
    diag = [max((ax.get_xlim()[0],ax.get_ylim()[0])),min((ax.get_xlim()[1],ax.get_ylim()[1]))]
    ax.plot(diag, diag, linestyle = linestyle,color=color,lw=lw)

    if fit:
        xmin = int(ax.get_xlim()[0]) - 1
        xmax = int(ax.get_xlim()[1]) + 1
        x = np.linspace(xmin, xmax, 100)
        ax.plot(x, fit(x), color='blue',linestyle = '--',lw=lw)


def _aligning_axis(ax1,ax2,x=1,y=1):
    if x:
        xlimits = np.array([ax1.get_xlim(),ax2.get_xlim()])
        ax1.set_xlim([xlimits.min(),xlimits.max()])
        ax2.set_xlim([xlimits.min(),xlimits.max()])
    if y:
        ylimits = np.array([ax1.get_ylim(),ax2.get_ylim()])
        ax1.set_ylim([ylimits.min(),ylimits.max()])
        ax2.set_ylim([ylimits.min(),ylimits.max()])

def _bootstrap(data, num_samples = 100000, statistic = np.mean, alpha=0.95, asymmetric=True):
    data = data[~np.isnan(data)]
    n = len(data)
    if n>0:
        samples = np.empty(num_samples)
        for i in range(num_samples):
            idx = np.random.randint(0, n, n)
            samples[i] = statistic(data[idx])
        stat = np.sort(samples)
        y_mean = stat[int((1/2.0)*num_samples)]
        y_err = -stat[int((1-alpha)/2.0*num_samples)]+y_mean
        y_err2 = stat[int((1+alpha)/2.0*num_samples)]-y_mean

        if asymmetric:
            return (y_mean, y_err, y_err2)
        else:
            return (y_mean, (y_err+y_err2)/2.0)
    else:
        if asymmetric:
            return (np.nan, np.nan, np.nan)
        else:
            return (np.nan, np.nan)

def gauss(x, A, mu, sigma):
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def norm_distribution(x, mu, sigma):
    return gauss(x, 0.39894228/mu, mu, sigma)

def _gaussian_distribution_fit_bootstrap(data, num_samples = 10000, statistic = np.mean, alpha=0.95, asymmetric=True):
    from scipy.optimize import curve_fit
    data = data[~np.isnan(data)]
    n = len(data)
    if n>0:
        samples = np.empty(num_samples)
        for i in range(num_samples):
            idx = np.random.randint(0, n, n)
            hist, bin_edges = np.histogram(data[idx], density=True)
            bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
            p0 = [np.max(data[idx]), np.mean(data[idx]), np.var(data[idx])]
            coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)
            samples[i] = coeff[1]
        stat = np.sort(samples)
        y_mean = stat[int((1/2.0)*num_samples)]
        y_err = -stat[int((1-alpha)/2.0*num_samples)]+y_mean
        y_err2 = stat[int((1+alpha)/2.0*num_samples)]-y_mean

        if asymmetric:
            return (y_mean, y_err, y_err2)
        else:
            return (y_mean, (y_err+y_err2)/2.0)
    else:
        if asymmetric:
            return (np.nan, np.nan, np.nan)
        else:
            return (np.nan, np.nan)

def _polyfit_bootstrap(x,y,order,sigma=[], num_samples = 100000, alpha=0.95, asymmetric=True):
    if len(x) != len(y):
        print('x, y have different length!')
        return -1
    elif len(x)<=0:
        print('empty sample!')
        return -2
    else:
        fit = np.empty([num_samples,order+1]) 
        for i in range(num_samples):  
            index = np.random.randint(0, len(x), len(x))
            x_sample = x[index]
            y_sample = y[index]
            if len(sigma)>0:
                sigma_sample = sigma[index]
                fit[i] = np.polyfit(x_sample, y_sample,order, w=sigma_sample)
            else:
                fit[i] = np.polyfit(x_sample, y_sample,order)
        fit = np.transpose(fit)
        stat = np.array([np.sort(fit[i]) for i in range(order+1)])
        fit_mean = stat[:,int(0.5*num_samples)]
        fit_err = -stat[:,int((1-alpha)/2.0*num_samples)]+fit_mean
        fit_err2 = stat[:,int((1+alpha)/2.0*num_samples)]-fit_mean
        if asymmetric:
            return fit_mean, fit_err, fit_err2
        else:
            return fit_mean, (fit_err+fit_err2)/2.0

# Cooling Function
#import h5py
#MIN_TEMP = 4.0                                   
#MAX_TEMP = 8.5                                   
#N_TEMPS = 91                                     
#temp_step = (MAX_TEMP - MIN_TEMP) / (N_TEMPS - 1)
#fcooling = '/home/yqin/bitbucket/meraxes/meraxes/input/cooling_functions/SD93.hdf5'
#group_name = ["mzero","m-30","m-20","m-15","m-10","m-05","m-00","m+05"]            
#cooling_rate = np.zeros((len(group_name),N_TEMPS)) 

#with h5py.File(fcooling,'r') as cooling_function: 
#    for group in range(len(group_name)):
#        for items in range(N_TEMPS): 
#            cooling_rate[group][items] = cooling_function[group_name[group]]['log(lambda_norm)'][items]
#
#def interpolate_temp_dependant_cooling_rate(i_m, logT):
#    global cooling_rate
#    MIN_TEMP = 4.0
#    MAX_TEMP = 8.5
#    N_TEMPS = 91
#    temp_step = (MAX_TEMP - MIN_TEMP) / (N_TEMPS - 1)
#
#    i_t       = ((logT - MIN_TEMP) / temp_step).astype(np.int)
#    i_t[i_t>(N_TEMPS - 1)] = N_TEMPS - 2
#    
#    rate_below = cooling_rate[i_m,i_t]
#    rate_above = cooling_rate[i_m,i_t + 1]
#
#    logT_below = MIN_TEMP + temp_step * i_t
#    rate = rate_below + (rate_above - rate_below) / temp_step * (logT - logT_below)
#    rate[logT<MIN_TEMP] = -27.0
#
#    return rate
#
#def _cooling_function(logT,logZ):
#    import bisect
#
#    metallicities = np.array([-5.0,-3.0,-2.0,-1.5,-1.0,-0.5,+0.0,+0.5])
#    logZ[logZ < metallicities[0]] = metallicities[0]
#    logZ[logZ > metallicities[-1]] = metallicities[-1]
#    
#    i_m = np.array([bisect.bisect_left(metallicities, i) for i in logZ],dtype=np.int)-1
#
#    rate_below = interpolate_temp_dependant_cooling_rate(i_m, logT)
#    rate_above = interpolate_temp_dependant_cooling_rate(i_m + 1, logT)
#
#    rate = rate_below + (rate_above - rate_below) / (metallicities[i_m + 1] - metallicities[i_m]) * (logZ -metallicities[i_m])
#
#    return 10**rate

def _sigmoid(X,w):
    X_w = X**w
    return X_w/(1+X_w)


def _find_redshift(i, desc_indexs, file_offsets, ms, m0, snapshot, max_snapshot=103):
    m_p = None
    while (desc_indexs[snapshot][i] >=0) & (snapshot+file_offsets[snapshot][i] <=max_snapshot) & (ms[snapshot][i]<=m0):
        snapshot_p = snapshot
        m_p = ms[snapshot][i]
        ind = desc_indexs[snapshot][i]
        snapshot +=file_offsets[snapshot][i]
        i = ind
        m = ms[snapshot][i]

    if not m_p:
        if ms[snapshot][i]<=m0:
            return -1, -1, -1, -1 #no desc
        else:
            return -2, -2, -2, -2 #m>m0
    elif ms[snapshot][i]<=m0:
        return -3, -3, -3, -3 #cannot find it until last snapshot
    else:
        return m[0], m_p[0], snapshot, snapshot_p




def _find_redshift_p_2(i, desc_indexs, file_offsets, m0, ms, snapshot):
    m_b, snapshot_b = None, snapshot+1
    m = ms[snapshot][i]
    while (m>=m0) and (snapshot>0) and (snapshot!=snapshot_b):
        m_b, snapshot_b = m, snapshot
        i_prime, m_prime, snapshot_prime, i_tmp, m_tmp, snapshot_tmp = None, 0, None, None, 0, None
        for Nsnap in range(1,min(17,snapshot+1)):
            if len(file_offsets[snapshot-Nsnap]) >0:
                indexs = [a for a,b in enumerate((file_offsets[snapshot-Nsnap] == Nsnap) & (desc_indexs[snapshot-Nsnap] == i)) if b]
                if len(indexs) > 0:
                    samples = ms[snapshot-Nsnap][indexs]
                    indexmax = np.argmax(samples)
                    i_tmp, m_tmp, snapshot_tmp = indexs[indexmax], samples[indexmax], snapshot-Nsnap
                    if m_tmp > m_prime: i_prime, m_prime, snapshot_prime = i_tmp, m_tmp, snapshot_tmp

        if i_prime:
            i, m, snapshot = i_prime, m_prime, snapshot_prime
        elif i_tmp:
            i, m ,snapshot = i_tmp, m_tmp, snapshot_tmp

    if not m_b:
        return -2, -2, -2, -2 #m<m0
    elif snapshot == 0:
        return -3, -3, -3, -3 #cannot find it until snapshot 0
    elif snapshot == snapshot_b:
        return -1, -1, -1, -1 #cannot find it until the earliest
    else:
        return m, m_b, snapshot, snapshot_b


def _find_redshift_p_fast(m0, ms, snapshot):
    m = ms[snapshot]
    m_b = None
    while m > m0:
        m_b = m
        snapshot_b = snapshot
        snapshot -=1
        m = ms[snapshot]

    if m <=0:
        return -1, -1, -1, -1 #cannot find it until the earliest
    elif not m_b:
        return -2, -2, -2, -2 #m<m0
    elif m>m0:
        return -3, -3, -3, -3 #cannot find it until snapshot 0
    else:
        while ms[snapshot-1] == m:
            snapshot-=1
        return m, m_b, snapshot, snapshot_b

def _find_redshift_b_fast(m0, ms, snapshot0,maxsnapshot=103):
    m = ms[0]
    snapshot = 0
    while (m < m0) and (snapshot< maxsnapshot-snapshot0):
        m_p = m
        snapshot_p = snapshot
        snapshot +=1
        m = ms[snapshot]

    if m<m0:                                         
        return -1, -1, -1, -1 #cannot find it until the latest 
    elif snapshot == maxsnapshot+1-snapshot0:
        return -3, -3, -3, -3 #cannot find it until snapshot 103 
    else:
        while (snapshot+1<len(ms)-1) and (ms[snapshot+1] == m):
            snapshot+=1
        if (snapshot+1==len(ms)-1) and (ms[snapshot+1] == m):
            snapshot+=1
        return m, m_p, snapshot+snapshot0, snapshot_p+snapshot0

def _distribution_peak_value_bootstrap(data, num_samples = 100000, statistic = np.mean, alpha=0.95, asymmetric=True):
    data = data[~np.isnan(data)]
    n = len(data)
    if n>0:
        samples = np.empty(num_samples)
        for i in range(num_samples):
            idx = np.random.randint(0, n, n)
            vals, edges = np.histogram(data[idx], bins = max(n/100, 100))
            samples[i] = edges[round(np.median(np.where(vals==vals.max())))]
        stat = np.sort(samples)
        y_mean = stat[int((1/2.0)*num_samples)]
        y_err = -stat[int((1-alpha)/2.0*num_samples)]+y_mean
        y_err2 = stat[int((1+alpha)/2.0*num_samples)]-y_mean

        if asymmetric:
            return (y_mean, y_err, y_err2)
        else:
            return (y_mean, (y_err+y_err2)/2.0)
    else:
        if asymmetric:
            return (np.nan, np.nan, np.nan)
        else:
            return (np.nan, np.nan)

def gauss_up(x, mu, sigma):
    from scipy.special import erfc
    return 0.5*erfc((mu - x)/(1.41421356* sigma))

def gauss_mu_up(x, mu, sigma):  
    from scipy.special import erfc
    return -((np.exp(-((mu - x)**2./(2.*sigma**2.)))* sigma)*0.79788456) + mu

def gauss_sigma2_up(x, mu, sigma):
    from scipy.special import erf, erfc
    return sigma**2. - 2.*sigma*np.exp(-((mu - x)**2./(2.* sigma**2.)))*(x-mu)/2.506628/erfc((mu - x)/(1.41421356*sigma)) - (0.79788456 *sigma*np.exp(-((mu - x)**2./(2.* sigma**2.)))/erfc((mu - x)/(1.41421356*sigma)))**2.

def _malmquist_equations_up(p,x_up,mean,var):
    mu, sigma = p
    return (gauss_mu_up(x_up, mu, sigma) - mean, gauss_sigma2_up(x_up, mu, sigma) - var)

#def gauss_mu(x, mu, sigma):  
#    from scipy.integrate import quad
#    norm = lambda x: np.exp(-0.5*x**2.)*0.39894228
#    return quad(norm, -np.inf, (x-mu)/sigma) *mu - norm((x-mu)/sigma)*sigma

#def gauss_sigma(x, mu, sigma):
#    from scipy.integrate import quad
#    norm = lambda x: np.exp(-0.5*x**2.)*0.39894228
#    return quad(norm, -np.inf, (x-mu)/sigma) *(mu**2. + sigma**2.) - norm((x-mu)/sigma) *sigma*(x+mu)

def gauss_low(x, mu, sigma):
    from scipy.special import erf
    return 0.5*(1+erf((mu - x)/(1.41421356* sigma)))

def gauss_mu_low(x, mu, sigma):  
    from scipy.special import erf
    return 0.5* (mu + np.exp(-((mu - x)**2./(2.* sigma**2.)))* 0.79788456* sigma + mu* erf((mu - x)/(1.41421356*sigma)))/gauss_low(x, mu, sigma) 

def gauss_sigma2_low(x, mu, sigma):
    from scipy.special import erf
    return sigma**2. + 2.*sigma*np.exp(-((mu - x)**2./(2.* sigma**2.)))*(x-mu)/2.506628/(1+erf((mu - x)/(1.41421356*sigma))) - (0.79788456 *sigma*np.exp(-((mu - x)**2./(2.* sigma**2.)))/(1+erf((mu - x)/(1.41421356*sigma))))**2.

def _malmquist_equations_low(p,x_low,mean,var):
    mu, sigma = p
    return (gauss_mu_low(x_low, mu, sigma) - mean, gauss_sigma2_low(x_low, mu, sigma)- var)

def _malmquist_fit_bootstrap_low(data,x_low, num_samples = 100000, statistic = np.mean, alpha=0.95, asymmetric=True):
    from scipy.optimize import fsolve
    data = data[~np.isnan(data)]
    n = len(data)
    if n>0:
        samples = np.empty(num_samples)
        for i in range(num_samples):
            idx = np.random.randint(0, n, n)
            mean = np.mean(data[idx])
            var = np.var(data[idx])
            samples[i] = fsolve(_malmquist_equations_low, (mean, var), args=(x_low,mean,var))[0]

        stat = np.sort(samples)
        y_mean = stat[int((1/2.0)*num_samples)]
        y_err = -stat[int((1-alpha)/2.0*num_samples)]+y_mean
        y_err2 = stat[int((1+alpha)/2.0*num_samples)]-y_mean

        if asymmetric:
            return (y_mean, y_err, y_err2)
        else:
            return (y_mean, (y_err+y_err2)/2.0)
    else:
        if asymmetric:
            return (np.nan, np.nan, np.nan)
        else:
            return (np.nan, np.nan)
