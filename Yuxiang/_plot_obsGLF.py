#!/usr/bin/env python
import numpy as np
#from cosmolopy.luminosityfunction import schechterM as _Schechter

__author__ = ('Yuxiang Qin')                     
__email__ = ("Yuxiang.L.Qin@Gmail.com")          
__all__  = ('plot_obsGLF', 'plot_obsGLF_legend',)

def _h_convertor_x(hubble_h_obs,hubble_h=0.678): 
    return -5. * np.log10(hubble_h/hubble_h_obs) 
                                                 
def _h_convertor_y(hubble_h_obs,hubble_h=0.678): 
    return (hubble_h/hubble_h_obs)**3.           

#def _plot_Schechter(ax,x,M_star,phi_star,alpha,label=None,linestyle='-',lw=2,hubble_h=0.678,markersize=10,legend=True):
#    y = _Schechter(x, phi_star, alpha, M_star)
#    if legend:
#        ax.plot(x+_h_convertor_x(0.7, hubble_h),\
#                y*_h_convertor_y(0.7,hubble_h),\
#                color='gray',linestyle=linestyle,\
#                lw=lw,label = label)
#    else:
#        ax.plot(x+_h_convertor_x(0.7, hubble_h),\
#                y*_h_convertor_y(0.7,hubble_h),\
#                color='gray',linestyle=linestyle,\
#                lw=lw)

def plot_obsGLF_z10pt0(ax,hubble_h=0.678,markersize=10,legend=True):

    #Bouwens et al. 2015
    data = []
    file = '/home/yqin/dragons/programs/galaxy_luminosity/Bouwens15/lfobs_z10.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                    data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='s',color='gray',label='Bouwens et al. 2015')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                    data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='s',color='gray')
    
    data = []
    file = '/home/yqin/dragons/programs/galaxy_luminosity/Bouwens15/lfups_z10.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                yerr=0.9* data[:,1],uplims=True,\
                markersize=markersize,fmt='s',color='gray')
    
    #Bouwens et al. 2015b
    data = []
    file = '/home/yqin/dragons/programs/galaxy_luminosity/Bouwens15/lfobs_z10b.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                    data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='o',color='gray',label='Bouwens et al. 2015b')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                    data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='o',color='gray')
    
    data = []
    file = '/home/yqin/dragons/programs/galaxy_luminosity/Bouwens15/lfups_z10b.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                yerr=0.9* data[:,1],uplims=True,\
                markersize=markersize,fmt='o',color='gray')
    
def plot_obsGLF_z9pt0(ax,hubble_h=0.678,markersize=10,legend=True):

    #Bouwens et al. 2015b
    data = []
    file = '/home/yqin/dragons/programs/galaxy_luminosity/Bouwens15/lfobs_z9b.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                    data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='o',color='gray',label='Bouwens et al. 2015b')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                    data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='o',color='gray')
    
    data = []
    file = '/home/yqin/dragons/programs/galaxy_luminosity/Bouwens15/lfups_z9b.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                yerr=0.9* data[:,1],uplims=True,\
                markersize=markersize,fmt='o',color='gray')
    
def plot_obsGLF_z8pt0(ax,hubble_h=0.678,markersize=10,legend=True):

    #Bouwens et al. 2015
    data = []
    file = '/home/yqin/dragons/programs/galaxy_luminosity/Bouwens15/lfobs_z8.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                    data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='s',color='gray',label='Bouwens et al. 2015')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                    data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='s',color='gray')
    
    data = []
    file = '/home/yqin/dragons/programs/galaxy_luminosity/Bouwens15/lfups_z8.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                yerr=0.9* data[:,1],uplims=True,\
                markersize=markersize,fmt='s',color='gray')

def plot_obsGLF_z7pt0(ax,hubble_h=0.678,markersize=10,legend=True):

    #Bouwens et al. 2015
    data = []
    file = '/home/yqin/dragons/programs/galaxy_luminosity/Bouwens15/lfobs_z7.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                    data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='s',color='gray',label='Bouwens et al. 2015')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                    data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='s',color='gray')
    
    data = []
    file = '/home/yqin/dragons/programs/galaxy_luminosity/Bouwens15/lfups_z7.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                yerr=0.9* data[:,1],uplims=True,\
                markersize=markersize,fmt='s',color='gray')

    #Atek et al. 2015
    data = []
    file = '/home/yqin/dragons/programs/galaxy_luminosity/Bouwens15/Combine_lf_z7.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                    data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='^',color='gray',label='Atek et al. 2015')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                    data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='^',color='gray')
    
def plot_obsGLF_z6pt0(ax,hubble_h=0.678,markersize=10,legend=True):

    #Bouwens et al. 2015
    data = []
    file = '/home/yqin/dragons/programs/galaxy_luminosity/Bouwens15/lfobs_z6.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                    data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='s',color='gray',label='Bouwens et al. 2015')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                    data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='s',color='gray')
    
def plot_obsGLF_z5pt0(ax,hubble_h=0.678,markersize=10,legend=True):

    #Bouwens et al. 2015
    data = []
    file = '/home/yqin/dragons/programs/galaxy_luminosity/Bouwens15/lfobs_z5.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                    data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='s',color='gray',label='Bouwens et al. 2015')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                    data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='s',color='gray')
    
def plot_obsGLF_z4pt0(ax,hubble_h=0.678,markersize=10,legend=True):

    #Bouwens et al. 2015
    data = []
    file = '/home/yqin/dragons/programs/galaxy_luminosity/Bouwens15/lfobs_z4.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    if legend:
        ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                    data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='s',color='gray',label='Bouwens et al. 2015')
    else:
        ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                    data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                    yerr=data[:,2],\
                    markersize=markersize,fmt='s',color='gray')

    x = np.linspace(-23.5,-18.5)
#    _plot_Schechter(ax,x,-20.84,1.36e-3,-1.56,label='van der Burg et al. 2010',linestyle='-',lw=4,hubble_h=hubble_h,markersize=markersize,legend=legend)
#    _plot_Schechter(ax,x,-21.14,1.46e-3,-1.82,label='Yoshida et al. 2006',linestyle='-.',lw=4,hubble_h=hubble_h,markersize=markersize,legend=legend)


#def plot_obsGLF_z3pt0(ax,hubble_h=0.678,markersize=10,legend=True):

#    x = np.linspace(-23.5,-18.5)
#    _plot_Schechter(ax,x,-20.94,1.79e-3,-1.65,label='van der Burg et al. 2010',linestyle='-',lw=4,hubble_h=hubble_h,markersize=markersize,legend=legend)
#    _plot_Schechter(ax,x,-20.97,1.71e-3,-1.73,label='Reddy et al. 2009',linestyle=':',lw=4,hubble_h=hubble_h,markersize=markersize,legend=legend)
#    _plot_Schechter(ax,x,-20.90,1.70e-3,-1.43,label='Sawicki et al. 2006',linestyle='--',lw=4,hubble_h=hubble_h,markersize=markersize,legend=legend)

def plot_obsGLF_legend(ax,markersize=10,color='gray',alpha=1.0):

    #Bouwens et al. 2015
    data = []
    file = '/home/yqin/dragons/programs/galaxy_luminosity/Bouwens15/lfobs_z7.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                yerr=data[:,2],\
                markersize=markersize,fmt='s',color='gray',label='Bouwens et al. 2015')

    #Bouwens et al. 2015b
    data = []
    file = '/home/yqin/dragons/programs/galaxy_luminosity/Bouwens15/lfobs_z10b.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                yerr=data[:,2],\
                markersize=markersize,fmt='o',color='gray',label='Bouwens et al. 2015b')

    #Atek et al. 2015
    data = []
    file = '/home/yqin/dragons/programs/galaxy_luminosity/Bouwens15/Combine_lf_z7.dat'
    for line in open(file,'r'):  data.append(line.split())
    data = np.array(data[1:],dtype=np.float)
    ax.errorbar(data[:,0]+_h_convertor_x(0.69999999,hubble_h),\
                data[:,1]*_h_convertor_y(0.69999999,hubble_h),\
                yerr=data[:,2],\
                markersize=markersize,fmt='^',color='gray',label='Atek et al. 2015')

    #_plot_Schechter(ax,x,-21.14,1.46e-3,-1.82,label='Yoshida et al. 2006',linestyle='-.',lw=4,hubble_h=hubble_h,markersize=markersize,legend=True)
    #_plot_Schechter(ax,x,-20.90,1.70e-3,-1.43,label='Sawicki et al. 2006',linestyle='--',lw=4,hubble_h=hubble_h,markersize=markersize,legend=True)
    #_plot_Schechter(ax,x,-20.97,1.71e-3,-1.73,label='Reddy et al. 2009',linestyle=':',lw=4,hubble_h=hubble_h,markersize=markersize,legend=True)
    #_plot_Schechter(ax,x,-20.94,1.79e-3,-1.65,label='van der Burg et al. 2010',linestyle='-',lw=4,hubble_h=hubble_h,markersize=markersize,legend=True)
    #_plot_Schechter(ax,x,-20.84,2.11e-3,-1.60,label='van der Burg et al. 2010 (Magnification)',linestyle='--',lw=4,hubble_h=hubble_h,markersize=markersize,legend=True)

def plot_obsGLF(ax,z,hubble_h=0.678,markersize=10,legend=True,silent=False,color='gray',alpha=1.0):
    _plot_obsGLF = {10.0: plot_obsGLF_z10pt0,\
                    9.0:  plot_obsGLF_z9pt0,\
                    8.0:  plot_obsGLF_z8pt0,\
                    7.0:  plot_obsGLF_z7pt0,\
                    6.0:  plot_obsGLF_z6pt0,\
                    5.0:  plot_obsGLF_z5pt0,\
                    4.0:  plot_obsGLF_z4pt0}
                   #3.0:  plot_obsGLF_z3pt0}
    try:
        (_plot_obsGLF[z])(ax,hubble_h=hubble_h,markersize=markersize,legend=legend)
    except:
        if not silent:
            print('available GLF at redshift:', _plot_obsGLF.keys())
        return -1
