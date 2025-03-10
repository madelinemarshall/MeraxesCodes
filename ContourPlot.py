from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

def density_estimation(m1, m2):
    index=np.isfinite(m1)&np.invert(np.isnan(m1))&np.isfinite(m2)&np.invert(np.isnan(m2))
    m1=m1[index]
    m2=m2[index]
    xmin = m1.min()
    xmax = m1.max()
    ymin = m2.min()
    ymax = m2.max()
    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]  
    positions = np.vstack([X.ravel(), Y.ravel()])                                                        
    values = np.vstack([m1, m2])                                                                        
    
    kernel = stats.gaussian_kde(values)                                                                 
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z


def contour_plot(x,y,xlab=None,ylab=None,xlims=None,ylims=None,axes=None,colors=None,levels=None,linewidth=2.5,linestyle='-'):
    X,Y,Z=density_estimation(x,y)
    Z=Z/np.sum(Z)
    n = 1000
    t = np.linspace(0, np.amax(Z), n)
    integral = ((Z >= t[:, None, None]) * Z).sum(axis=(1,2))
    f = interpolate.interp1d(integral, t)
    if levels.all():
      levels = f(np.array(np.flip(levels,0)))
    else:  
      levels = f(np.array([0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]))


    if axes==None:
      #if levels==None:
      #  plt.contour(X, Y, Z,colors=colors,linewidths=linewidth)
      #else:
      plt.contour(X, Y, Z,colors=colors,levels=levels,linewidths=linewidth,linestyles=linestyle)
      if xlab is not None:
        plt.xlabel('{}'.format(xlab))
      if ylab is not None:
        plt.ylabel('{}'.format(ylab))
      if xlims==None:
        plt.xlim([min(x),max(x)])
      else:
        plt.xlim(xlims);
      if ylims==None:
        plt.ylim([min(y),max(y)])
      else:
        plt.ylim(ylims);
    else:
      #if levels==None:
      #  CS=axes.contour(X, Y, Z,colors=colors,linewidths=linewidth)
        #axes.clabel(CS,inline=1,fontsize=10)
      #else:
      #plt.scatter(X,Y,c=Z)
      #plt.colorbar()
      #levels=[10,20,30,40,50,60,70,80,90]/max(Z)
      

      CS=axes.contour(X, Y, Z,colors=colors,levels=levels,linewidths=linewidth,linestyles=linestyle)
      #axes.clabel(CS, inline=1, fontsize=10)
      if xlab is not None:
        axes.set_xlabel('{}'.format(xlab))
      if ylab is not None:
        axes.set_ylabel('{}'.format(ylab))
      if xlims==None:
        axes.set_xlim([min(x),max(x)])
      else:
        axes.set_xlim(xlims);
      if ylims==None:
        axes.set_ylim([min(y),max(y)])
      else:
        axes.set_ylim(ylims);
    #    plt.show()
    return



