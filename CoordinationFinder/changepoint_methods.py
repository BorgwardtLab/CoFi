# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 11:59:19 2016

@author: Sir Thomas
"""
from scipy import stats
import pandas as pd
import numpy as np
import math
import bottleneck as bn

def CUSUM(seq,parameters=[8.,0.5]):
    seq = np.asarray(seq)
    h = parameters[0]
    delta = parameters[1]
    skip = 3
    i = skip
    S = [0]
    hit = 0
    length = seq.size
    while i < length:
        nu1 = seq[0:i].sum()/i
        sig1 = seq[0:i].std()
        s = delta/(sig1**2) *(seq[i] - nu1 - delta/2)
        S.append(max([s  +S[i-skip] , 0 ]))
        if S[-1] > h: #/i to account for multiple hypothesis testing
            hit = i
            h = S[-1]
            #step = np.argmin(S)
            break
        i = i + 1
    return hit
    
def EMWA(seq, parameters, onset=30): #does not work
    l = parameters[0] #recommended 0.05-0.2
    L = parameters[1] # recommended 3 for l<0.1, 2.6-2.8 for l>0.1
    k = onset #onset
    if len(seq) < k:
        print "Sequence too short EWMA"
        return 0
    hit = 0
    mu0 = np.mean(seq[:k])
    z = mu0
    test = [z]
    sig = np.std(seq[:k])
    diff = L * sig * np.sqrt(l/(2-l)) #statistical quality control 9.27 9.28
    for i in xrange(0,k):
        z = l*seq[i] + (1-l)*z #statistical quality control 9.22
        test.append(z)
    for i in xrange(k,len(seq)):
        z = l*seq[i] + (1-l)*z #statistical quality control 9.22
        test.append(z)
        if np.abs(z-mu0)>diff:
            hit = i
            #check = diff-np.mean(seq[:i])-np.abs(z)
            break
    return hit#, test,mu0,diff

def T_Test(seq,p=0.01):
    stop = 0
    A = []
    offset = 5
    i = offset
    length = len(seq)
    while i < length-offset:
        A.append(stats.ttest_ind(seq[0:i], seq[i:length], axis=0, equal_var=True)[1])
        i = i+1
    A = pd.Series(A).abs()
    if A.min() < p:#/length: Bonferroni correction for multiple hypothesis testing most probable cp
        stop = A.idxmin()+offset
        
    return stop

#http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4154476/
#A Statistical Change Point Model Approach for the Detection of DNA Copy Number Variations in Array CGH Data
#returns maximized log-likelihood for an iid normal sample xs
factor = math.log(2*math.pi)
def loglike(xs):
    n = len(xs)
    return -0.5 * n * (factor + math.log(bn.nanstd(xs)+0.0000000001) - 1)
    #n = len(xs)
    #return np.max(np.abs(stats.norm.logpdf(xs, loc=np.mean(xs), scale=np.std(xs))))#-0.5 * n * np.log(2 * np.pi * np.std(xs)) - 0.5 * n
def Binary_Segmentation(xs, maxLike =10., left=None, right=None):
    if left is None:
        left = 0

    if right is None:
        right = len(xs)

    OFFSETPCT = 0.125
    MINNOBS = 5

    ys = xs[left:right]
    offset = min(int(len(ys)*OFFSETPCT), MINNOBS)
    tLeft, tRight = left + offset, right - offset
    if tRight <= tLeft:
        print "Sequence too short"
        print xs
        return 0

    cp = 0
    dataLike = loglike(ys)
    for t in xrange(tLeft, tRight):
        profLike = loglike(xs[left:t]) + loglike(xs[t:right]) #likelihood for both samples added
        lr = 2*(profLike - dataLike)
        if lr > maxLike:
            cp = t
            maxLike = lr
    return cp
    
#QChart, 4of5 test
def QChart(seq,maxlike=5.,parameters_EMWA=(0.2,2.6)):
    if len(seq) < 7:
        print "Sequence too short QChart"
        return 0
    
    mean = float(seq[0]+seq[1])/2
    variance = (seq[0]-mean)**2 + (seq[1]-mean)**2
    stat = []
    chart=[0.,0,0] #sum[0] and score[1] of #elements[2]
    cp = 0
            
    stat_score = [0.]
    for r in xrange(3,len(seq)+1):
        stat.append(stats.norm.ppf(stats.t.cdf(np.sqrt(float(r-1)/r) * (float(seq[r-1]-mean)/np.sqrt(variance)),r-2))) #eq 7, Q Statsistics
        variance = float(r-2)*variance/(r-1) + 1./r*(seq[r-1]-mean)**2 #eq3 in 1991 paper
        mean = float((r-1)*mean+seq[r-1])/r
        chart[0] = chart[0]+ stat[-1]
        chart[1] = chart[1] +int(stat[-1]>0)
        stat_score.append(chart[0])
        if chart[2] < 5:
            chart[2] = chart[2] + 1 #first two: no checks, only increase
        else:
            chart[0] = chart[0]-stat[-5]
            chart[1] = chart[1]-int(stat[-5]>0)
            if (chart[1] < 1 or chart[1] > 4) and np.abs(chart[0]) > maxlike: #resturn most probable location for cp
                cp = r - 6
                maxlike = np.abs(chart[0])
                break #leave this out for making it offline and fpr emwa, enable it for online
    """
    #suggested by paper: EMWA or CUSUM
    cp = EMWA(stat, parameters_EMWA)
    if cp >0:
        cp=cp+2
    """
    return cp#,stat_score,test,mu0,diff
   
def CP_loc(Series): #pass list
    #init = np.array([5.87, 0.57, 0.11, 3.74, 6.5, 0.206, 6.08]) #C1,C2,E1,E2,Q,T,B simulated best p-threshold
    init = np.array([5.87, 0.57, 0.11, 3.74, 6.5, 0.05, 6.08]) #C1,C2,E1,E2,Q,T,B pure T-Test
    return CUSUM(Series,init[0:2]), EMWA(Series, init[2:4], onset=9), QChart(Series,init[4]),T_Test(Series,init[5]),Binary_Segmentation(Series,init[6])
