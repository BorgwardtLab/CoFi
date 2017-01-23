# -*- coding: utf-8 -*-
"""
Created on Wed Jul 06 14:10:50 2016

@author: Sir Thomas
"""
import pandas as pd
import numpy as np

from DTW import KMeans_DTW
from features import Extract_Features


from scipy.cluster.vq import kmeans2, vq, whiten
from sklearn import metrics
from feature_selection import Stability_Analysis
import time as t

from matplotlib import pyplot as plt
import matplotlib

def Simulate_Cell(cp,length,diff=1.,loc=0.): #cp = 1 for diff = 1, else = 0
    if loc == 0:
        loc = length/2
    else:
        loc = int(np.random.uniform(length/2-length/10, length-1))
    before = pd.Series(np.random.randn(loc), index=np.arange(0,loc)).cumsum()
    if cp == 0:
        diff = 0
    after = pd.Series(np.random.randn(length-loc)+diff, index=np.arange(loc,length)).cumsum()
    after = after + before[loc-1]
    testA = before.append(after)
    loc = int(np.random.uniform(length/10+1, length-1))
    before = pd.Series(np.random.randn(loc), index=np.arange(0,loc)).cumsum()
    after = pd.Series(np.random.randn(length-loc)+diff, index=np.arange(loc,length)).cumsum()
    after = after + before[loc-1]
    testB = before.append(after)
    return testA,testB

def tree_(tree, counter,cellNr, cp, maximum, norm, p,gen):
    length1 = int(np.random.uniform(20,40))
    length2 = int(np.random.uniform(30,50))
    test1A,test1B = Simulate_Cell(cp,length1)
    test2A,test2B = Simulate_Cell(cp,length2)
    cell = [cellNr*2]*length1
    cell.extend([cellNr*2 + 1]*length2)
    time = list(np.arange(length1)+counter)
    time.extend(list(np.arange(length2)+counter))
    if cellNr >= 2**(gen-1): #adjust fate
        Simulation = pd.DataFrame({'cellNr': cell, 'tree':[tree]*(length1+length2), 'ProteinA': test1A.append(test2A), 'ProteinB':  test1B.append(test2B), 'absTime': time, 'fate': [0]*(length1+length2) })
    if cellNr < 2**(gen-1):
        Simulation = pd.DataFrame({'cellNr': cell, 'tree':[tree]*(length1+length2), 'ProteinA': test1A.append(test2A), 'ProteinB':  test1B.append(test2B), 'absTime': time, 'fate': [1]*(length1+length2) })
        Simulation = Simulation.append(tree_(tree,counter+length1,cellNr*2,cp,maximum,norm,p,gen),ignore_index=True)
        Simulation = Simulation.append(tree_(tree,counter+length2,cellNr*2+1,cp,maximum,norm,p,gen),ignore_index=True)
    return Simulation



def Create_Simulation(k):
    
    #For drawing from Distribution
    maximum = 90
    norm = 25
    p = [-2.5,-0.003, -0.078, -0.002,0,0.1]
    Simulation = pd.DataFrame()
    tree = 0
    np.random.seed(k*42)
    diff = 1. #Height of cp
    
    while tree < 1:
        tree = tree + 1
        gen = tree + 5
        init = 35
        test1A,test1B = Simulate_Cell(1,init,diff)
        Simulation = pd.DataFrame({'cellNr': 1, 'tree':[str(tree)+'cp']*(init), 'ProteinA': test1A, 'ProteinB':  test1B, 'absTime': list(np.arange(init)), 'fate': [1]*init})
        Simulation = Simulation.append(tree_(str(tree)+'cp',init,1,1,maximum,norm,p,gen))
        test1A,test1B = Simulate_Cell(0,init)
        Simulation = Simulation.append(pd.DataFrame({'cellNr': 1, 'tree':[str(tree)+'nocp']*(init), 'ProteinA': test1A, 'ProteinB':  test1B, 'absTime': list(np.arange(init)), 'fate': [1]*init}))
        Simulation = Simulation.append(tree_(str(tree)+'nocp',init,1,0,maximum,norm,p,gen))
    
    return Simulation.reindex()

def F1_acc(idx):
    Trues = idx[:len(idx)/2]
    Falses = idx[len(idx)/2:]
    number = len(idx)/2
    number = float(number)
    if sum(Trues)/number > 0.5:
        Reference = 1
    else:
        Reference = 0
    i = 0
    TP = 0
    FP = 0
    FN = 0
    a = [0.,0.]
    b = [0.,0.]
    while i < number:
        if Trues[i] == Reference:
            TP = TP + 1
        if Trues[i] != Reference:
            FN = FN + 1
        if Falses[i] == Reference:
            FP = FP + 1
        j = 0
        while j < number: 
            if Trues[i] == Falses[j]:
                a[int(Reference == Trues[i])] = a[int(Reference == Trues[i])] + 1
            elif Trues[i] != Falses[j]:
                b[int(Reference == Trues[i])] = b[int(Reference == Trues[i])] + 1
            if Falses[i] == Trues[j]:
                a[int(Reference == Falses[i])] = a[int(Reference == Falses[i])] + 1
            elif Falses[i] != Trues[j]:
                b[int(Reference == Falses[i])] = b[int(Reference == Falses[i])] + 1
            j = j + 1
        i = i + 1
    s = [(b[0]-a[0])/max(a[0],b[0]), (b[1]-a[1])/max(a[1],b[1])] 
    S = 0.5*sum(s)
    TP = TP/number
    FP = FP/number
    F1 = 2.*TP/(TP+FP+1.)
    
    return F1, S

def Compare_DTW_Char():
    m= 0
    time_DTW = []
    time_CBC = []
    F1_DTW = []
    F1_CBC = []
    S_DTW = []
    S_CBC = []
    while m < 100 :
        m = m + 1
        #Data Ccreation
        Create_Simulation(m).to_csv('Simulation.csv', index=False)
        
        #Data curation
        data=pd.read_csv('Simulation.csv')
        start = t.time()
        Extract_Features("Simulation.csv", Features=['length', 'change', 'mchange', 'osc', 'inc', 'mean', 'sens', 'rate', 'ljungbox'], Columns=['cellNr', 'tree', 'fate', 'absTime', ['ProteinA','ProteinB']])
        points_c = pd.read_csv("points.csv", index_col=1)
        if m == 1:
            F,s = Stability_Analysis('random')
            print F
            print s
        points_c = points_c.loc[:,F[:10]]
        points_c = whiten(points_c.values.astype(float))
        k=2
        centroids_c,labels_c = kmeans2(points_c,k,minit='points')
        idx_c,_ = vq(points_c,centroids_c)
        CBC = t.time()
        
        #Clustering
        idx = KMeans_DTW(data,k)
        DTW = t.time()
        
        
        #Comparison
        print m
        print "Creating a sample dataset out of " + str(len(idx)) + " random walks with length between 20 and 60. Half of them have a changepoint at half of their length."
        F1, S = F1_acc(idx)
        print 'Clustering KMeans according to DTW. '
        print 'F1: ' + str(F1) 
        print 'Silhouette score: ' + str(S) 
        time_DTW.append(float(DTW-CBC))
        S_DTW.append(S)
        F1_DTW.append(F1)
        print 'Time in Minutes: ' + str(float(DTW-CBC)/60)
        F1, S = F1_acc(idx_c)
        print 'Clustering KMeans according to Character.'
        print 'F1: ' + str(F1) 
        print 'Silhouette score: ' + str(S)
        print 'Time in Minutes: ' + str(float(CBC-start)/60) 
        time_CBC.append(float(CBC-start))
        S_CBC.append(S)
        F1_CBC.append(F1)
    
    matplotlib.rc('font', family='serif', serif='cm12')
    matplotlib.rc('text', usetex=True)
    matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
    
    weights = np.ones_like(S_CBC)/len(S_CBC)
    
    #weights = np.ones_like(S_CBC)/float(len(S_CBC)) #bins=10, weights=weights
    plt.hist(S_CBC, bins=[0, 0.1, 0.2, 0.30, 0.40, 0.50, 0.6,0.7,0.8,0.9,1.00], weights=weights, color='green', label=r'global character ' + str(np.round(np.mean(S_CBC),decimals=2)) + '$\pm$' + str(np.round(np.var(S_CBC),decimals=2)))
    #weights = np.ones_like(S_DTW)/float(len(S_DTW))
    plt.hist(S_DTW, bins=[0, 0.1, 0.2, 0.30, 0.40, 0.50, 0.6,0.7,0.8,0.9,1.00], weights=weights, color='blue', label=r'DTW '+ str(np.round(np.mean(S_DTW),decimals=2)) + '$\pm$' + str(np.round(np.var(S_DTW),decimals=2)))
    title_string = 'Comparing performance of feature extraction'
    plt.suptitle(title_string, fontsize=20)  
    plt.xlabel('S', fontsize=20)
    plt.tick_params(labelsize=20)
    plt.legend(loc='best',fontsize=18)
    plt.ylabel('p', fontsize=20)
    plt.show(block=True)
    
    print np.mean(time_CBC), np.std(time_CBC)
    print np.mean(time_DTW), np.std(time_DTW)
    print np.mean(F1_CBC), np.std(F1_CBC)
    print np.mean(F1_DTW), np.std(F1_DTW)
    print np.mean(S_CBC), np.std(S_CBC)
    print np.mean(S_DTW), np.std(S_DTW)
    
    
Compare_DTW_Char()