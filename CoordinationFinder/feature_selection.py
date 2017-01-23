# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 18:42:42 2016

@author: Sir Thomas
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 16:00:07 2016

@author: Sir Thomas
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.vq import kmeans2
import os.path as ospath
from scipy.spatial.distance import pdist

from KMeans import Clean_Data


def Distance_Centroids(Centroids):
    
    i=0
    Distance = 0
    Centroids = np.asarray(Centroids)
    while i < Centroids.shape[1]:
        Distance = Distance + np.linalg.norm(pdist(Centroids[:,i,:],'sqeuclidean')) #Clusters index over 0: Features, 1: clusters, 2: distortions pdist:The points are arranged as m n-dimensional row vectors
        i = i+1
    return Distance.sum()/(i**2) #divide by number of clusters squared to normalize

def Accuracy_Plot(data,Feat,Feat_var,dist):
    score = []
    i = 0
    while i < len(Feat_var):
        rescaled,k = Clean_Data(data,Feat + [Feat_var[i]],0)
        Cent = []
        j = 0
        start = np.random.rand(k,len(Feat)+1)   #when I initialize this for every trial differently I end up having the cluster centers at completely different positions. Now it differs only by a bit
        while j < len(dist):
            Distortion = np.random.randn(rescaled.shape[0],rescaled.shape[1])*dist[j]
            rescaled = rescaled + Distortion

            centroids,labels = kmeans2(rescaled,start, minit='matrix')
            Cent.append(centroids)
            j = j+1
        score.append(Distance_Centroids(np.asarray(Cent)))
        i = i+1
    score = pd.Series(score)
    return_Feat = Feat_var[score.idxmin()]
    Feat_var = Feat_var[:score.idxmin()] + Feat_var[score.idxmin()+1 :]
    return return_Feat, score.min(), Feat_var

def Feature_Search(data,Distortion):
    
    if Distortion == 'linear':
        Distortion = np.arange(0.01,0.10,0.01)
    if Distortion == 'random':
        Distortion = np.abs(np.random.randn(100)*0.05)
    
    Features = list(data)[1:]
        
    F = []
    s = []
    No_Features = len(Features)
    Var_Feat = Features
    i = 0
    while i < No_Features: #Forward search for the Feature which is the most stable under slight distortion of the data
        F_int, s_int, Var_Feat = Accuracy_Plot(data, F, Var_Feat, Distortion)
        F.append(F_int)
        s.append(s_int)
        i = i+1
    return F,s

def Stability_Analysis(Distortion):
    data=pd.read_csv(ospath.join("points.csv"), index_col=1) 
    plt.ion()    
    np.random.seed(42)
    F,s = Feature_Search(data,Distortion)
    #print F,s
    
    """
    #Compare to a test on random data to look for effects of dimensionality:
    s = []
    s_ = []
    i = 0
    m = 100
    while i < m:
        np.random.seed(22*(i+1))
        new = np.random.rand(len(F)+1,len(F)+1)
        new = pd.DataFrame(new)
        F_,s__ = Feature_Search(new,Distortion)
        F,s___ = Feature_Search(data,Distortion)
        print len(F), len(F_), len(s_), len(s__)
        if i == 0:
            s_ = s__
            s = s___
        else:
            s_ = [sum(x) for x in zip(s_,s__)]
            s = [sum(x) for x in zip(s,s___)]
        i = i + 1
    s_ = [float(k)/m for k in s_]
    s = [float(k)/m for k in s]
    
    
    x = np.arange(len(F))
    plt.plot(x,s, label='global features')
    plt.plot(x,s_, label='baseline dimensionality effect')
    #plt.xticks(x, F, rotation=30, fontsize=8)
    plt.semilogy() #x,s)
    plt.xlabel('#features', fontsize=20)
    plt.ylabel("Sensitivity to distortion", fontsize=20)
    plt.tick_params(labelsize=20)
    plt.legend(loc='best')
    plt.suptitle('Greedy forward feature search', fontsize=20)
    plt.title('averaged over ' + str(m) + ' trials.')
    plt.savefig(ospath.join("Stability.png"), format="png")
    plt.show(block=True)
    plt.clf()
    """
    
    return F,s

#print Stability_Analysis('random')