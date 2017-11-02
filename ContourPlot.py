from scipy import stats
import numpy as np
import matplotlib.pyplot as plt

def density_estimation(m1, m2):
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

def contour_plot(x,y,xlab,ylab,xlims,ylims):
    X,Y,Z=density_estimation(x,y)
    plt.contour(X, Y, Z)
    plt.xlabel('{}'.format(xlab))
    plt.ylabel('{}'.format(ylab))
    if xlims==0:
      plt.xlim([min(x),max(x)])
    else:
      plt.xlim(xlims);
    if ylims==0:
      plt.ylim([min(y),max(y)])
    else:
      plt.ylim(ylims);
    plt.show()
