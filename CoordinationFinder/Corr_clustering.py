# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 10:55:48 2016

@author: Sir Thomas
"""

from sklearn.cluster import AgglomerativeClustering
import pandas as pd
import numpy as np
from numpy import linalg as LA
from sklearn.decomposition import PCA

def Feature_Clustering(Corr):
    X = Corr.values
    ward = AgglomerativeClustering(n_clusters=5, linkage='ward').fit(X)
    
    pca = PCA(n_components=5)
    pca.fit(Corr) 
    comp =  pca.components_
    
    label = list(ward.labels_)
    label = pd.Series(label, index=Corr.index.values.tolist(), name="clustering")
    a = [str(i) for i in Corr.keys().tolist()]
    a = a + ['clustering']
    
    Co = pd.concat((Corr,label), axis=1)
    Co = Co.sort(columns='clustering')
    del Co['clustering']
    
    # Compute sum of eigenvalues associated with 
    w, v = LA.eig(Co.values.tolist())
    i = 0
    K = 0.
    K_ = 0.
    N = 0.
    N_ = 0.
    N__ = 0.
    K__ = 0.
    ind = []
    crossN = 0
    crossK = 0
    while i < len(Co.index.values.tolist()):
        if 'Nanog' in Co.index.values[i]:
            for k in Co.loc[:,Co.index.values[i]]:
                k = 0
                while k < len(Co.index.values.tolist()):
                    if 'Nanog' in Co.index.values[k]:
                        N__ = N__ + Co.values[i,k]
                        crossN = crossN + 1
                    k = k + 1
            N = N + abs(w[i])
            N_ = N_ + abs(comp[0,i])
            ind.append(2)
        elif 'Klf4' in Co.index.values[i]:
            for k in Co.loc[:,Co.index.values[i]]:
                k = 0
                while k < len(Co.index.values.tolist()):
                    if 'Klf4' in Co.index.values[k]:
                        K__ = K__ + Co.values[i,k]
                        crossK = crossK + 1
                    k = k + 1
            K = K + abs(w[i])
            K_ = K_ + abs(comp[0,i])
            ind.append(1)
        else :
            ind.append(0)
        i = i + 1
    
    ind = pd.Series(label, index=Corr.index.values.tolist(), name="clustering")
    Co = pd.concat((Corr,label), axis=1)
    Co = Co.sort(columns='clustering')
    del Co['clustering']

    print "Sum of abs eigenvalues of correlation of features associated with Nanog: " + str(N)
    print "Sum of abs eigenvalues of correlation of features associated with Klf4: " + str(K)
    print "Sum of contributions of first PC from Nanog: " + str(N_)
    print "Sum of contributions of first PC from Klf4: " + str(K_)
    print "Total correlation among Nanog: " + str(N__/crossN)
    print "Total correlation among Klf4: " + str(K__/crossK)
    
    
    
    return Co
    

def Test():
    Data = np.random.randn(9,9)
    ind = ['A', 'B','C','D','E','F','G','H','I']
    Data = pd.DataFrame(Data, index=ind)
    print Data
    print Feature_Clustering(Data)

#Test()