import numpy as np
import scipy.stats as stats

def chi_sq(sample1hist,sample2hist):

  n_bins=len(sample1hist)
  df = n_bins
  chisq = 0.0
  root2on1=np.sqrt(sum(sample2hist)/sum(sample1hist))
  root1on2=np.sqrt(sum(sample1hist)/sum(sample2hist))

  for ii in range(0,n_bins):
    if (sample1hist[ii] == 0.0) & (sample2hist[ii] == 0.0):
      df -=1   
    else:
      chisq += ((root2on1*sample1hist[ii]-root1on2*sample2hist[ii])**2)/(sample1hist[ii]+sample2hist[ii])

  prob=1-stats.chi2.cdf(chisq,df)
  return chisq, prob
