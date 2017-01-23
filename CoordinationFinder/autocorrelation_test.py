# -*- coding: utf-8 -*-
"""
Created on Mon Aug 08 11:37:32 2016

@author: Sir Thomas
"""

import numpy as np
from statsmodels.stats.diagnostic import acorr_ljungbox
import matplotlib.pyplot as plt
from scipy.stats import chi2

"""
length = 100
acorr =  0.
i = 0
H = int(length/5. - 1)
np.random.seed(42)
runs = 100
while i < runs:
    seq = np.sin(np.linspace(0.,1.,length))
    Abs, p = acorr_ljungbox(seq,lags=H)
    if chi2.cdf(0.95,H)-p[H-1] > 0:
        acorr = acorr + 1
    i = i + 1

print 'TPR on Sine: ' + str(acorr/runs)
"""
"""
runs = 1000
bins = 100
length = np.linspace(10,60,bins)
acorr = []
acorr_ = []
acorr__ = []
sine = []
for l in length:
    H = min([10,int(l/5.)])
    H = 4
    np.random.seed(42)
    i = 0
    acorr.append(0.)
    acorr_.append(0.)
    acorr__.append(0.)
    sine.append(0.)
    while i < runs:
        seq = np.random.randn(l)
        seq_ = np.cumsum(seq)
        seq__ = np.log(seq)
        sq = np.sin(np.linspace(0.,2*np.pi,l))
        Abs, p = acorr_ljungbox(seq,lags=H)
        if p[H-1]-0.05 < 0:
            acorr[-1] = acorr[-1] + 1.
        Abs, p = acorr_ljungbox(seq_,lags=H)
        if p[H-1]-0.05 < 0:
            acorr_[-1] = acorr_[-1] + 1.
        Abs, p = acorr_ljungbox(seq__,lags=H)
        if p[H-1]-0.05 < 0:
            acorr__[-1] = acorr__[-1] + 1.
        Abs, p = acorr_ljungbox(sq,lags=H)
        if p[H-1]-0.05 < 0:
            sine[-1] = sine[-1] + 1.
        i = i + 1
acorr = [x / runs for x in acorr]
acorr_ = [x / runs for x in acorr_]
acorr__ = [x / runs for x in acorr__]
sine = [x / runs for x in sine]

plt.plot(length,acorr, label='Gaussian noise')
plt.plot(length,acorr_, label='Random walk')
plt.plot(length,acorr__, label='log of Random walk')
plt.plot(length, sine,label='smooth sine')
plt.suptitle('Autocorrelation of short random walks', fontsize=20 )
plt.title(str(runs)+' trails', fontsize=18)
plt.xlabel('length', fontsize=20)
plt.ylim(-0.05,1.05)
plt.xticks([10,20,30,40,50,60], fontsize=18)
plt.ylabel('Ljung-Box test rejection rate',fontsize=20)
plt.legend(loc='best', fontsize=18)
plt.show(block=True)

trials = 1000
length = 100
st=[100,200,500,1000]
j = 0
scaling = np.linspace(0,2,length)
ac = [0.]*length
acorr = [0.] * length
acorr_ = [0.] * length
acorr__ = [0.] * length
m = 0
while m < 4:
    np.random.seed(42*m+1)
    steps = np.linspace(0,np.pi,st[m])
    sine = np.diag(np.sin(steps))
    sine = np.dot(sine,np.ones((st[m],length)))
    j = 0
    while j < trials:
        j = j + 1
        noise = np.dot(np.random.randn(st[m],length),np.diag(scaling))
        test = sine + noise #100 columns
        i = 0
        while i < length:
            H = 4
            Abs, p = acorr_ljungbox(test[:,i],lags=H)
            if p[H-1]-0.05 < 0:
                if m == 0:
                    ac[i] = ac[i] + 1.
                if m == 1:
                    acorr[i] = acorr[i] + 1.
                if m == 2:
                    acorr_[i] = acorr_[i] + 1.
                if m == 3:
                    acorr__[i] = acorr__[i] + 1.
            i = i + 1
    m = m +1
acorr = [x / trials for x in acorr]
ac = [x / trials for x in ac]
acorr_ = [x / trials for x in acorr_]
acorr__ = [x / trials for x in acorr__]
plt.ylim(-0.05,1.05)
plt.plot(scaling,ac, label='pi/100')
plt.plot(scaling,acorr, label='pi/200')
plt.plot(scaling,acorr_, label='pi/500')
plt.plot(scaling,acorr__, label='pi/1000')
plt.xticks([0,1,2], fontsize=18)
plt.yticks([0.,0.5,1.], fontsize=18)
plt.suptitle('Critical information density for autocorrelation',fontsize=20)
plt.title(str(trials)+' trails', fontsize=18)
plt.ylabel('Ljung-Box test rejection rate',fontsize=20)
plt.xlabel('noise / amplitude',fontsize=20)
plt.legend(loc='best',fontsize=18)
plt.show(block=True)
"""