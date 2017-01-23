# -*- coding: utf-8 -*-
"""
Created on Tue Jul 05 11:51:39 2016

@author: Sir Thomas
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.cluster.vq import kmeans2, vq, whiten

from dtw import dtw

#This is from http://nbviewer.jupyter.org/github/alexminnaar/time-series-classification-and-clustering/blob/master/Time%20Series%20Classification%20and%20Clustering.ipynb
#Dynamic time wrapping distance without window (slow)
def DTWDistance(s1, s2):
    DTW={}
    
    for i in range(len(s1)):
        DTW[(i, -1)] = float('inf')
    for i in range(len(s2)):
        DTW[(-1, i)] = float('inf')
    DTW[(-1, -1)] = 0

    for i in range(len(s1)):
        for j in range(len(s2)):
            dist= (s1[i]-s2[j])**2
            DTW[(i, j)] = dist + min(DTW[(i-1, j)],DTW[(i, j-1)], DTW[(i-1, j-1)])
		
    return np.sqrt(DTW[len(s1)-1, len(s2)-1])

def curate_matrix(data,cellNr,Proteins,fate):
    it = int(data.describe().ix[0,0])
    temp1 = data.loc[0,str(cellNr)]
    counter = 0
    loc_counter = 0
    lines = pd.DataFrame()
    #if len(Proteins) > 1:
        #lines = lines.astype(object)
    index = 0
    for j in Proteins:
        x = 0
        while x < it: #goes though file column by column
            temp2 = data.loc[x,str(cellNr)]
            if ((temp2 != temp1) or (x == it-1)):
                if int(temp2) != 1 and data.loc[x, str(fate)] != 3 and data.loc[x, str(fate)] != 0:
                    length = int(x-loc_counter)
                    
                    i = 0
                    while i < length:
                        lines.loc[index,i] = data.loc[loc_counter+i,str(j)]
                        i = i + 1
                    index = index + 1
                    
                    counter = counter +1
                    loc_counter = x+1
                if x != it-1:
                    temp1 = temp2
            x = x+1
    return [r.dropna().values.tolist() for a,r in lines.iterrows()]
    
 
def abs_point(data,no):
    lim = len(data)/no
    distance = np.zeros((lim,lim))
    i = 0
    while i < lim:
        j = i
        while j < lim:
            c = 0
            while c < no:
                distance[i,j] = distance[i,j] + DTWDistance(data[i+c*lim],data[j+c*lim])
                c = c + 1
            distance[j,i] = distance[i,j]
            j = j + 1
        i = i + 1
    U,S,V = np.linalg.svd(distance)
    U = U[:,:10]
    S = np.sqrt(S[:10])
    return np.dot(U,np.diag(S))

def KMeans_DTW(data,k,Columns=['cellNr','fate',['ProteinA','ProteinB']]):
    by_column = curate_matrix(data,Columns[0],Columns[2],Columns[1])
    by_column = abs_point(by_column,len(Columns[2]))
    centroids,labels = kmeans2(by_column,k,minit='points')
    idx,_ = vq(by_column,centroids)
    #Show(by_column,idx,k)
    return idx


def Show(points,idx,k):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    colors = ['orange', 'blue', 'azure', 'lightcyan', 'grey', 'lightblue', 'yellow']
    x=0
    while x < k:
        ax.scatter(points[idx==x,0], points[idx==x,1], points[idx==x,2], c=colors[x]) 
        x = x+1
    plt.show(block=True)
    
