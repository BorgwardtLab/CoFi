# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 10:03:46 2016

@author: Sir Thomas
"""

import pandas as pd
import os.path as ospath
import matplotlib
import matplotlib.pyplot as plt
import scipy.stats as st
import numpy as np
from scipy.optimize import leastsq


#this diesn't work - WHY???
def freedman_bin_width(data, return_bins=False):
    r"""Return the optimal histogram bin width using the Freedman-Diaconis rule

    The Freedman-Diaconis rule is a normal reference rule like Scott's
    rule, but uses rank-based statistics for results which are more robust
    to deviations from a normal distribution.

    Parameters
    ----------
    data : array-like, ndim=1
        observed (one-dimensional) data
    return_bins : bool (optional)
        if True, then return the bin edges

    Returns
    -------
    width : float
        optimal bin width using the Freedman-Diaconis rule
    bins : ndarray
        bin edges: returned if ``return_bins`` is True

    Notes
    -----
    The optimal bin width is

    .. math::
        \Delta_b = \frac{2(q_{75} - q_{25})}{n^{1/3}}

    where :math:`q_{N}` is the :math:`N` percent quartile of the data, and
    :math:`n` is the number of data points [1]_.

    References
    ----------
    .. [1] D. Freedman & P. Diaconis (1981)
       "On the histogram as a density estimator: L2 theory".
       Probability Theory and Related Fields 57 (4): 453-476

    See Also
    --------
    knuth_bin_width
    scott_bin_width
    bayesian_blocks
    histogram
    """
    data = np.asarray(data)
    if data.ndim != 1:
        raise ValueError("data should be one-dimensional")

    n = data.size
    if n < 4:
        raise ValueError("data should have more than three entries")

    v25, v75 = np.percentile(data, [25, 75])
    dx = 2 * (v75 - v25) / (n ** (1 / 3))

    if return_bins:
        dmin, dmax = data.min(), data.max()
        Nbins = max(1, np.ceil((dmax - dmin) / dx))
        bins = dmin + dx * np.arange(Nbins + 1)
        return dx, bins
    else:
        return dx



def Feature_Distribution(user, data, Feat, bins=0):
    matplotlib.rc('font', family='serif', serif='cm12')
    matplotlib.rc('text', usetex=True)
    matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    #data=pd.read_csv(ospath.join("Features","points.csv"))
    if bins == 0:
        da, bins = freedman_bin_width(data[Feat].values, True)
        bins = len(bins)-1
    if 'loc' in Feat:
        data = data[data[str(Feat)] != 0.]
    histogram = st.histogram(data[Feat].values, numbins=bins)
    x_data = np.linspace(histogram[1],histogram[2]*bins+histogram[1], bins)
    y_data = histogram[0]
    m = np.mean(data[Feat].values)
    md = np.median(data[Feat].values)
    v = np.std(data[Feat].values)
    
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.bar(x_data, y_data, width=histogram[2])
    ax.tick_params(labelsize=20)
    ax.set_ylabel("Frequency", fontsize=22)
    title_string = 'Mean: ' + str(np.round(m,decimals=2)) + r'$\pm$' + str(np.round(v,decimals=3))
    md_string = 'Median: ' + str(np.round(md,decimals=2))
    
    #Put names to the features. The order is important!!!
    if 'd1mom' in Feat:
        Feat = 'detrended mean intensity'
    elif 'd2mom' in Feat:
        Feat = 'detrended variance'
    elif 'd3mom' in Feat:
        Feat= 'detrended skewness'
    elif 'd4mom' in Feat:
        Feat='detrended kurtosis'
    elif 'dereg' in Feat:
        Feat = 'detrended non-linearity'
    elif 'dsens' in Feat:
        Feat = 'detrended chaos'
    elif 'dljungbox' in Feat:
        Feat = 'detrended autocorrelation'
    elif 'length' in Feat:
        Feat= 'length'
    elif 'change' in Feat:
        Feat="change-point"
    elif '2mom' in Feat:
        Feat = 'variance'
    elif '3mom' in Feat:
        Feat= 'skewness'
    elif '4mom' in Feat:
        Feat='kurtosis'
    elif 'inc' in Feat:
        Feat= 'raw increase'
    elif 'end' in Feat:
        Feat = 'end intensity'
    elif 'ereg' in Feat:
        Feat = 'non-linearity'
    elif 'sens' in Feat:
        Feat = 'lyapunov exponent'
    elif 'ljungbox' in Feat:
        Feat = 'autocorrelation'
    elif '1mom' in Feat:
        Feat = 'mean intensity'
    ax.set_title(r'\textbf{ ' + str(Feat) + ' of each sequence}', fontsize=25)
    #plt.title(title_string, fontsize=22)
    #plt.annotate(title_string, xy = (m, 0), xytext = (0, 20), fontsize=22, textcoords = 'offset points', ha = 'right', va = 'bottom',bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.8),arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
    #plt.annotate(md_string, xy = (md, 0), xytext = (80, 60), fontsize=22, textcoords = 'offset points', ha = 'right', va = 'bottom',bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.8),arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
    #plt.suptitle(title_string, fontsize=18)
    #ax.set_xlabel('Values', fontsize=22)
    box = ax.get_position()
    if Feat == 'length':
        ax.set_title('Length',fontsize=25)
        ax.set_xlabel('Time in 30min', fontsize=22)
        ax.set_position([box.x0, box.y0 + box.height * 0.4, box.width, box.height * 0.6])
        r = matplotlib.patches.Rectangle((0,0), 1, 1, fill=False, edgecolor='none',visible=False)
        ax.legend([r,r],[title_string,md_string],prop={'size':22},loc='upper center',bbox_to_anchor=(0.5, -0.2),ncol=1)
    else: 
        ax.set_title(r'\textbf{ ' + str(Feat) + ' of each sequence}', fontsize=25)
        ax.set_position([box.x0, box.y0 + box.height * 0.2, box.width, box.height * 0.8])
        r = matplotlib.patches.Rectangle((0,0), 1, 1, fill=False, edgecolor='none',visible=False)
        ax.legend([r,r],[title_string,md_string],prop={'size':22},loc='upper center',bbox_to_anchor=(0.5, -0.05),ncol=1)
    plt.savefig(ospath.join(user,str(Feat) + ".pdf"), format="pdf")
    plt.clf()
    plt.close()
    
    

def Length_Distribution(init, fitfunc):
    
    errfunc  = lambda p, x, y: (y - fitfunc(p, x))
    
    data=pd.read_csv(ospath.join("points.csv"))
    plt.figure()
    bins = 20
    histogram = st.histogram(data['length'].values, numbins=bins)
    x_data = np.linspace(histogram[1],histogram[2]*bins+histogram[1], bins)
    y_data = histogram[0]
    
    out = leastsq(errfunc, init, args=(x_data,y_data))
    norm = max(fitfunc(out[0],x_data))
    
    return x_data[-1], norm, out[0]

def Draw_Length(maximum, norm,p,fitfunc):
    discard = 0.5
    while discard < 1: #random sampling via Rejection method
        length = int(np.random.uniform(3,maximum))
        discard = (fitfunc(p,length)/norm)*np.random.randn()
    return length
    

def Draw_Histograms(user, bins=0, Features = []):
    matplotlib.rc('font', family='serif', serif='cm12')
    matplotlib.rc('text', usetex=True)
    matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    data=pd.read_csv("points.csv")
    if len(Features) == 0:
        Features = list(data)[1:]
    for F in Features:
        Feature_Distribution(user, data, F, bins)


#Draw_Histograms('test_person',10)
