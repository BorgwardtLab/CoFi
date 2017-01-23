# -*- coding: utf-8 -*-
"""
Created on Wed Apr 06 11:21:18 2016

@author: Sir Thomas
"""
import pandas as pd
import os.path as ospath
import numpy as np

from changepoint_methods import CP_loc
from statsmodels.stats.diagnostic import acorr_ljungbox
from scipy import stats
from lyapunov import max_lyapunov

import matplotlib.pyplot as plt


from scipy.stats import chi2,combine_pvalues
from scipy.signal import detrend


def Extract_Features(filename, Features=['length', 'fate', 'change', 'mchange', 'osc', 'inc', 'mean', 'sens', 'rate', 'ljungbox'], Columns=['cellNr', 'tree', 'stopReason', 'absoluteTime', ['intNanog','intKlf4']], inclomplete=[0,0], normalization=[0,""]):
    
    data = pd.read_csv(filename)
    it = int(data.describe().ix[0,Columns[0]]) #how many data points
    DF_Columns = ['Start']
    for i in Features:
        if i == 'length' or i == 'fate':
            DF_Columns.append(str(i))
        else:
            for j in Columns[4]:
                DF_Columns.append(str(i)+str(j))
    
    points = pd.DataFrame(columns = DF_Columns) #Data frame for curves as single points
    temp1 = data.loc[0,Columns[0]]
    counter = 0
    points.set_value(counter,'Start',0)
    x = 0
    shift = 0
    if inclomplete[1] == 0:
        shift = -2 #exclude last two data points of every time series
    #features
    length = 0
    
    if normalization[0] == 1:
        print normalization[1]
    
    while x < it: #goes though file column by column
        temp2 = data.loc[x,Columns[0]]
        if ((temp2 != temp1) or (x == it-1)):
            if inclomplete[0] != 0 or int(temp2) != 1: #if the cell number changes, a new curve starts. Make sure that we exclude first cells because of incompleteness of the curves
                x = x + shift
                #location of previous cell starting point
                loc_counter = int(points.loc[counter, 'Start'])
                if normalization[0] == 1:
                    mean_norm = data[loc_counter:x][str(normalization[1])].mean()
                    var_norm = data[loc_counter:x][str(normalization[1])].var()
                    for j in Columns[4]:
                        TF_mean = data[loc_counter:x][str(j)].mean()
                        TF_variance = data[loc_counter:x][str(j)].var()
                        if np.isnan(TF_variance) or np.isnan(var_norm):
                            TF_variance = 1
                            var_norm = 1
                        #data[loc_counter:x][str(j)] = data[loc_counter:x][str(j)].div(data[loc_counter:x][str(normalization[1])])
                        data[loc_counter:x][str(j)] = (  ((data[loc_counter:x][str(j)]-TF_mean)/TF_variance).div((data[loc_counter:x][str(normalization[1])]-mean_norm)/var_norm)  )*TF_variance + TF_mean
                
                if x-loc_counter <= 0:
                    x = loc_counter+1
                length = int(x-loc_counter)
                    
                #get length of curve
                if 'length' in Features:
                    points.set_value(counter,'length', length)
                
                points.set_value(counter, 'fate', data.loc[loc_counter, Columns[2]])
                
                #get mean
                if 'mean' in Features:
                    mean = data[loc_counter:x].mean()
                    for j in Columns[4]:
                        points.set_value(counter,'mean'+str(j),mean.loc[str(j)])
                        
                #detrended Featrues
                if 'd1mom' in Features:
                    for j in Columns[4]:
                        Ser = data.loc[loc_counter:x, str(j)].values.tolist()
                        Ser = np.array(Ser)
                        X = Ser[~np.isnan(Ser)]
                        X = X[~np.isinf(X)]
                        dSer = detrend(X)
                        if len(dSer) == 0:
                            print "Zero length ", Ser
                        points.set_value(counter,'d1mom'+str(j),np.mean(dSer))
                        points.set_value(counter, 'd2mom'+str(j),np.var(dSer))
                        points.set_value(counter, 'd3mom'+str(j),stats.skew(dSer))
                        points.set_value(counter, 'd4mom'+str(j),stats.kurtosis(Ser))
                        y = range(len(dSer))
                        s,i,r,p,e = stats.linregress(dSer,y)
                        points.set_value(counter, 'dereg'+str(j),e)
                    
                        
                        
                #get moments
                if '1mom' in Features:
                    mean = data[loc_counter:x].mean()
                    for j in Columns[4]:
                        points.set_value(counter,'1mom'+str(j),mean.loc[str(j)])
                if '2mom' in Features:
                    osc = data[loc_counter:x].var()
                    for j in Columns[4]:
                        points.set_value(counter, '2mom'+str(j),osc.loc[str(j)])
                if '3mom' in Features:
                    for j in Columns[4]:                        
                        Ser = data.loc[loc_counter:x, str(j)].values.tolist()
                        points.set_value(counter, '3mom'+str(j),stats.skew(Ser))
                if '4mom' in Features:
                    for j in Columns[4]:                        
                        Ser = data.loc[loc_counter:x, str(j)].values.tolist()
                        points.set_value(counter, '4mom'+str(j),stats.kurtosis(Ser))
                #regression and error
                if ('reg' in Features) or ('ereg' in Features):
                    for j in Columns[4]:     
                        Ser = data.loc[loc_counter:x, str(j)].values.tolist()
                        y = range(len(Ser))
                        s,i,r,p,e = stats.linregress(Ser,y)
                        if 'reg' in Features:
                            points.set_value(counter, 'reg'+str(j),s)
                        if 'ereg' in Features:
                            points.set_value(counter, 'ereg'+str(j),e)
                
                
                #endpoint
                if 'end' in Features:
                    for j in Columns[4]:
                        b = data.loc[x-1, str(j)]
                        points.set_value(counter, 'end'+str(j), b)
                
                #Increase
                if 'inc' in Features:
                    for j in Columns[4]:
                        b = data.loc[x-1, str(j)]
                        a = data.loc[loc_counter, str(j)]
                        points.set_value(counter, 'inc'+str(j), b-a)
                
                #rate
                if 'rate' in Features:
                    for j in Columns[4]:
                        a = data.loc[loc_counter, str(j)]
                        points.set_value(counter, 'rate'+str(j), (b-a)/length)
                
                #oscillations
                if 'osc' in Features:
                    osc = data[loc_counter:x].var()
                    for j in Columns[4]:
                        points.set_value(counter, 'osc'+str(j),osc.loc[str(j)])
                
                if (length > 2) and ('mchange' in Features):
                    for j in Columns[4]:
                        diff = []
                        i = 1
                        while i < length:
                            diff.append(data.loc[loc_counter+i:x-1, str(j)].sum()/(length-i+1) - data.loc[loc_counter:loc_counter+i-1, str(j)].sum()/i)    
                            i = i +1
                        diff = pd.Series(diff)
                        steepness = diff.max()
                        points.set_value(counter, 'mchange'+str(j), steepness)
                if (length <= 2) and ('mchange' in Features):
                        points.set_value(counter, 'mchange'+str(j), 0.)
                
                if 'change' in Features:
                    for j in Columns[4]:
                        hit  = 0
                        stop = 0
                        Ser = data.loc[loc_counter:x, str(j)]
                        Serd = [l-k for k, l in zip(Ser.tolist()[:-1], Ser.tolist()[1:])] 
                        #Is a change happening? If T-Test and CUSUM we assume it to happen
                        # Ask T-Test! p=0.05 usual statistical significance
                        if len(Serd) > 10:
                            stop,e,q,hit,b = CP_loc(Serd)
                        #stop = T_Test(0.05,Ser.values)
                        #CUSUM (k,h) = (1.5,1.61), (0.25,8), (0.5,8)
                        #hit = CUSUM([8,0.5],Ser.values)
                        
                        if hit != 0:
                            points.set_value(counter, 'change'+str(j), 1)
                        if hit == 0:
                            points.set_value(counter, 'change'+str(j), 0)
                        
                        if 'loc' in Features:
                            points.set_value(counter, 'loc'+str(j), float(hit)/length)
                            
                    
                #Do the Box-Piecre test for Autocorrelation (Compare to noise)
                #Hyndman http://robjhyndman.com/hyndsight/ljung-box-test/ recommends min(10, T/5) for non-seasonal time series
                #Stata uses min(n/2 − 2, 40) http://www.stata.com/manuals13/tswntestq.pdf
                #It is common to use a Ljung-Box test to check that the residuals from a time series model resemble white noise. if p <0.05 considered noise
                if 'ljungbox' in Features:
                    H = 4 #Gene transcription 10 minutes - 1 data point. Period = 2, Take 2*2 as suggested. 
                    for j in Columns[4]:
                        points.set_value(counter, 'ljungbox'+str(j), 0.)  #maybe make it less arbitrary
                        threshold = 0.05
                        if length > 5:
                            Ser = list(data.loc[loc_counter:x, str(j)].values)
                            H = 7
                            H_ = 4
                            if len(Ser)-2 < H:
                                H = len(Ser)-2
                            if len(Ser)-2 < H_:
                                H_ = len(Ser)-2
                            H = 4
                            threshold = threshold#/H
                            Abs, p = acorr_ljungbox(Ser,lags=H+1)
                            #binary = int(combine_pvalues(p[H:])[1]-threshold < 0)
                            #print 'LJB', j, p[H]
                            binary = int(p[H]-threshold < 0) #make sure null hypothesis is false uncomment for H=4
                            points.set_value(counter, 'ljungbox'+str(j), binary)# p-value based on chi-square distribution - small p value means rejection of random distribution
                            if 'dljungbox' in Features:
                                Abs, p = acorr_ljungbox(dSer,lags=H+1)
                                binary = int(p[H]-threshold < 0)
                                #binary = int(combine_pvalues(p[H_:])[1]-0.0001 < 0)
                                points.set_value(counter, 'dljungbox'+str(j), binary)
                                
                        else:
                            points.set_value(counter, 'ljungbox'+str(j), 0) #assume radnomness if nothing else
                            if 'dljungbox' in Features:
                                points.set_value(counter, 'dljungbox'+str(j), 0)
                
                #Compute the sensitivity to initial condistions (Lyapunov exponent as a measure of chaos)
                #Normalize this to a normally distributed sequence of smae length/mean/var???
                if 'sens' in Features:
                    for j in Columns[4]:
                        val = 1.
                        if length > 10:
                            Ser = data.loc[loc_counter:x, str(j)]
                            val = max_lyapunov(Ser.tolist())
                            if val == -float('Inf'):
                                val = -10
                            #print 'LE', j, val
                        points.set_value(counter, 'sens'+str(j),val)
                
                
                x = x-shift
                points.set_value(counter+1,'Start',x)
                counter = counter +1
        
            if x != it-1:
                temp1 = temp2
        x = x+1

    points = points.dropna(axis=0,how='any')
    points = points[points.fate != 2]
    if inclomplete[0] == 0:
        points = points[points.fate != 0] #exclude last part of tree (incomplete curves)
        points = points[points.fate != 3]
        if 'fate' not in Features:
            points = points[points.fate != 2]
    if 'fate' not in Features:
        points = points.drop('fate', 1)
    points.to_csv(ospath.join('points.csv'))

