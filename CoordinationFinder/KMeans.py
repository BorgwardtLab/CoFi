# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 18:43:25 2016

@author: Sir Thomas
"""

import pandas as pd
import numpy as np
import os
import platform
import matplotlib.pyplot as plt
from scipy.cluster.vq import kmeans2, vq, whiten
from sklearn import metrics
from sklearn.cluster import KMeans
import os.path as ospath
from statsmodels.sandbox.tools.tools_pca import pcasvd
from mpl_toolkits.mplot3d import Axes3D as ax3D
from scipy import stats
import random

from sklearn import manifold,discriminant_analysis

from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

def OptK(rescaled,Feat):
    k = 2
    silhouette = []
    
    while k < 7:
        np.random.seed(k*5)
        centroids = np.random.rand(k,len(Feat))
        centroids,labels = kmeans2(rescaled,centroids, minit='matrix')
        a = labels[0]
        if any(v != a for v in labels): #if there are two or more clusters, i.e. if any value is != 0
            silhouette.append(metrics.silhouette_score(rescaled, labels, metric='sqeuclidean'))
        else:
            silhouette.append(0)
        k = k+1
        
    silhouette = pd.Series(silhouette)
    k = silhouette.idxmax()+2 #choose k with highest silhouette score
    return int(k)

def Clean_Data(data,Feat,Force):
    #Clean up data
    data = data.loc[:,Feat].dropna(axis=1, how='any') #Remove unusable rows (and unused colums?)
    #data = data[(np.abs(stats.zscore(data)) < 3).all(axis=1)] #Exclude outliers if <5 does not work!!!
    for F in Feat:
        data = data[np.abs(data[str(F)]-data[str(F)].mean())<=(3*data[str(F)].std())]   
    p_df = data.values.astype(float) #convert to List
    rescaled = whiten(p_df) #norm to var=1, need to include subtraction of mean!!!
    
    """
    #WHAT FOLLOWS IS TO SUBTRACT THE MEAN - MAYBE WE CAN INCLUDE IT TO whiten??
    #SLOWS DOWN STABILITY ANALYSIS
    """
    analysis = pd.DataFrame(rescaled)
    i = 0
    while i < len(analysis.columns.values):
        meani = analysis.describe().iloc[1,i]
        for j,rows in analysis.iterrows():
            rescaled[j,i] = rescaled[j,i] - meani
        i = i + 1
    #print pd.DataFrame(rescaled).describe()
    #Choose k
    if Force == 0:
        k = OptK(rescaled,Feat)
    if Force != 0:
        k = Force #User Forces amount of clusters
    
    return rescaled, k


def Show_KMeans(Space='PCA', Force=0, filename=0, Feat=0, Bias =0, Tree = 1, Columns=['cellNr', 'tree', 'stopReason', 'absoluteTime', ['intNanog','intKlf4']]):
    data=pd.read_csv(ospath.join("points.csv"), index_col=1)
    if Feat == 0:
        Feat = list(data)[1:]
    if Bias == 'nanog':
        i = 0
        while i < len(Feat):
            if 'Klf4' in Feat[i]:
                Feat.pop(i)
            else:
                i = i+1
    plt.ion()  
    rescaled, k = Clean_Data(data, Feat, Force)
    #RUN KMEANS
    np.random.seed(k*5)
    centroids = np.random.rand(k,len(Feat))
    centroids,labels = kmeans2(rescaled,centroids)
    idx,_ = vq(rescaled,centroids)
    
    colors = ['orange', 'blue', 'azure', 'lightcyan', 'grey', 'lightblue', 'yellow']
    

    
    if Space == 'PCA' and len(Feat) > 2:
        ax = ax3D(plt.gcf())
        #RUN PCA
        pca = pcasvd(rescaled, keepdim=0, demean=False)
        
        xpc, ypc, zpc = (0, 1, 2)
        xreduced, factors, evals, evecs = pca
        singvals = np.sqrt(evals)
        scale = 1
        
        # data
        xs = factors[:, xpc] * singvals[xpc]**(1. - scale)
        ys = factors[:, ypc] * singvals[ypc]**(1. - scale)
        zs = factors[:, zpc] * singvals[zpc]**(1. - scale)
        
        cluster=0
        while cluster < k:
            ax.scatter(xs[idx==cluster], ys[idx==cluster], zs[idx==cluster], c=colors[cluster]) 
            cluster = cluster+1
        
        tvars = np.dot(np.eye(factors.shape[0], factors.shape[1]),evecs) * singvals**scale
    
        i = 0
        while i < len(Feat): # enumerate(xreduced.columns.values):
            x, y, z = tvars[i][xpc], tvars[i][ypc], tvars[i][zpc]
            a = Arrow3D([0,x*4], [0,y*4], [0,z*4], mutation_scale=20, lw=1, arrowstyle="-|>", color="k")
            ax.add_artist(a)
            ax.text(x* 4.4, y * 4.4, z* 4.4, Feat[i], color='k', fontsize=14)
            i = i+1
    
        ax.set_xlabel('PC1', fontsize=16)
        ax.set_ylabel('PC2', fontsize=16)
        ax.set_ylabel('PC3', fontsize=16)
                
        if Force == 0:
            ax.text2D(0.05, 0.95, 'Maximal mean silhouette coefficient ' + str(round(metrics.silhouette_score(rescaled, labels, metric='sqeuclidean'),3))+' using k=' + str(int(k)) + ' clusters.', transform=ax.transAxes, fontsize=14)
        if Force != 0:
            ax.text2D(0.05, 0.95, 'Mean silhouette coefficient ' + str(round(metrics.silhouette_score(rescaled, labels, metric='sqeuclidean'),3))+' using k=' + str(int(k)) + ' clusters.', transform=ax.transAxes)
    
    elif Space == '2D':
        X = rescaled
        n_samples, n_features = X.shape
        #TSNE
        #print "Is it"
        #tsne = manifold.TSNE(n_components=3, init='pca', random_state=5)
        #print "coming?"
        #X = tsne.fit_transform(X)
        #print "Now."
        
        #Locally linear embedding
        n_neighbors=30
        clf = manifold.LocallyLinearEmbedding(n_neighbors, n_components=2, method='standard')
        X = clf.fit_transform(X)

        
        x_min, x_max = np.min(X, 0), np.max(X, 0)
        X = (X - x_min) / (x_max - x_min)
        
        plt.figure()
        ax = plt.subplot(111)
        
        xs = X[:, 0]
        ys = X[:, 1] 
        cluster=0
        lb = ['Stationary','Changing']
        while cluster < k:
            plt.scatter(xs[idx==cluster], ys[idx==cluster], c=colors[cluster], label=lb[cluster]) 
            cluster = cluster+1
        plt.xlim((-0.05,1.05))
        plt.ylim((-0.05,1.05))
        plt.legend(loc='best',fontsize=18)
        plt.xticks([]) 
        plt.yticks([])
        plt.suptitle("Locally Linear Embedding for k-means clustering visualization",fontsize=18)
        #plt.title('Mean silhouette score ' + str(round(metrics.silhouette_score(rescaled, labels, metric='sqeuclidean'),3))+' using ' + str(int(k)) + ' clusters.', fontsize=16)
        plt.show()
        
    
    data = data.loc[:,Feat].dropna(axis=0, how='any')
    for F in Feat:
        data = data[np.abs(data[str(F)]-data[str(F)].mean())<=(3*data[str(F)].std())] 
    p_df = data.values.astype(float)
    if Space == 'raw' and len(Feat) > 3:
        ax = ax3D(plt.gcf())
        #Determine the three Features with most variation and use them as the axes
        var = pd.DataFrame(p_df).var(axis=0).order(ascending=False).index.tolist()
        x=0
        while x < k:
            ax.scatter(p_df[idx==x,var[0]], p_df[idx==x,var[1]], p_df[idx==x,var[2]], c=colors[x]) 
            x = x+1

        ax.set_xlabel(Feat[var[0]])
        ax.set_ylabel(Feat[var[1]])
        ax.set_zlabel(Feat[var[2]])
      
    if platform.system() == 'Windows':    
        plt.show(block=True)
    if platform.system() == 'Linux':
        plt.show(block=False) #not showing anything
    rep_cellNr = []
    rep_tree = []
    if filename != 0:
        data['Cluster'] = labels
        #from sklearn.metrics import pairwise_distances_argmin_min,pairwise_distances
        
        kmeans = KMeans(n_clusters=2, random_state=0).fit(rescaled)
        ind_closest = []
        for i in range(k):
            d = kmeans.transform(rescaled)[:, i]
            ind_closest.extend(list(np.argsort(d)[::-1][:6]))
        
        #ind_closest, no = pairwise_distances_argmin_min(centroids, rescaled)
        # The array closest contains the index of the point in rescaled that is closest to each centroid.
        #Feat = 'Cluster'
        #from coordination import Coordination_Plots
        #Coordination_Plots(data,[Feat],filename,"Coordination of clusters.png")
    
        #save total result to .csv
        if Tree == 1:
            curves=pd.read_csv(filename)
            from save_tree import save_tree
            rep_cellNr, rep_tree = save_tree(data,curves,'tree.csv',Columns,ind_closest)
    print "Closest traces: ", rep_cellNr, rep_tree
    return k, rep_cellNr, rep_tree

#Show_KMeans(Space='2D', Force=2, filename='NanogKlf4nucmemcomma.csv', Feat=0, Bias =0, Tree = 1, Columns=['cellNr', 'tree', 'stopReason', 'timepoint', ['intNanog']])