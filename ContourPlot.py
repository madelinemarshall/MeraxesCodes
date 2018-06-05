from scipy import stats
import numpy as np
import matplotlib.pyplot as plt

def density_estimation(m1, m2):
    index=np.isfinite(m1)&np.invert(np.isnan(m1))&np.isfinite(m2)&np.invert(np.isnan(m2))
    m1=m1[index]
    m2=m2[index]
    xmin = m1.min()
    xmax = m1.max()
    ymin = m2.min()
    ymax = m2.max()
    X, Y = np.mgrid[xmin:xmax:50j, ymin:ymax:50j]                                                     
    positions = np.vstack([X.ravel(), Y.ravel()])                                                        
    values = np.vstack([m1, m2])                                                                        
    
    kernel = stats.gaussian_kde(values)                                                                 
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z

def contour_plot(x,y,xlab=None,ylab=None,xlims=None,ylims=None,axes=None,colors=None,levels=None,linewidth=2.5):
    X,Y,Z=density_estimation(x,y)
    if axes==None:
      #if levels==None:
      #  plt.contour(X, Y, Z,colors=colors,linewidths=linewidth)
      #else:
      plt.contour(X, Y, Z,colors=colors,levels=levels,linewidths=linewidth)
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
      CS=axes.contour(X, Y, Z,colors=colors,levels=levels,linewidths=linewidth)

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
