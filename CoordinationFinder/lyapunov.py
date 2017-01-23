# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 14:35:33 2016

@author: Sir Thomas
"""

import numpy as np
from attributes import Length_Distribution, Draw_Length, Feature_Distribution
import pandas as pd
from matplotlib import pyplot as plt

def max_lyapunov(test):
    LE = 0
    i = 1
    while i < len(test):
        j = 0
        le = [0]
        while j < len(test):
            if (i != j) and (2*i-j < len(test)) and (i-j > 0):
                le.append(1./(i-j) * np.log(np.abs(test[i]-test[2*i-j])/np.abs(test[i]-test[i-1])))
            j = j + 1
        LE = LE + np.sum(np.asarray(le)) #d np.argmax o we need to take np.abs?
        i = i + 1    
    return float(LE)/len(test)

"""    

fitfunc  = lambda p, x: p[0] + p[1]*x + p[2]*(x**2) + p[3]*(x**3) + p[4]*(x**4)
maximum, norm, p = Length_Distribution([5, 0.5, 0.5, 0.5, 0.2], fitfunc)
LE = []
LECP = []
GS=[]
length = []
i = 0
while i < 1000:
    leng = int(Draw_Length(maximum,norm,p,fitfunc))
    leng = 36
    LE.append(max_lyapunov(np.random.randn((leng)).cumsum().tolist()))
    loc = int(leng/2)
    before = pd.Series(np.random.randn(loc)).cumsum()
    after = pd.Series(np.random.randn(leng-loc)+1., index=np.arange(loc, leng)).cumsum()
    after = after + before[loc-1]
    RW1 = before.append(after)
    RW1.reindex().tolist()
    GS.append(max_lyapunov( np.random.randn((leng)).tolist()))
    LECP.append(max_lyapunov(RW1))
    length.append(leng)
    i = i+1

LE = pd.DataFrame({'Lyapunov exponent random walk':LE, 'Lyapunov exponent RW with CP':LECP, 'Lyapunov exponent Gaussian': GS })

Feature_Distribution('test_person', LE, 'Lyapunov exponent random walk', 10)
Feature_Distribution('test_person', LE, 'Lyapunov exponent RW with CP', 10)
Feature_Distribution('test_person', LE, 'Lyapunov exponent Gaussian', 10)
#LE.hist()
#plt.suptitle('Lyapunov exponent on 1000 test sequences of length 36', fontsize=18)
#plt.show()
"""