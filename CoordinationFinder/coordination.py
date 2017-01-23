# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 08:37:14 2016

@author: Sir Thomas
"""

import pandas as pd
import os.path as ospath
import numpy as np
from matplotlib import pyplot as plt

def Coordination(filename, data, Feat):
    curves=pd.read_csv(ospath.join(filename))
    index = data.index
    md = data.median(axis=0)[Feat] #We split the different features into two categories according to the median
    if md == 1:
        md = 0 #correct for binary features
    daughters = 0 #how many
    uncoordinated = 0
    coordinated = 0
    partial = [] #start ofmother cells
    coord = []
    uncoord = []
    a = curves.loc[index, ['cellNr','tree']]
    for row in zip(a['cellNr'], a['tree']):
        daughter1 = a.loc[(a['cellNr'] == 2*row[0]) & (a['tree'] == row[1])]
        daughter2 = a.loc[(a['cellNr'] == 2*row[0]+1) & (a['tree'] == row[1])]
        mother = a.loc[(a['cellNr'] == row[0]) & (a['tree'] == row[1])]
        if len(daughter1.index) != 0 & len(daughter2.index != 0):
            daughter1_ = (data.loc[daughter1.index.values, Feat].values > md)
            daughter2_ = (data.loc[daughter2.index.values, Feat].values > md)
            mother_ = (data.loc[mother.index.values, Feat].values > md)
            if daughter1_ == daughter2_:
                if daughter1_ == mother_:
                    coordinated = coordinated + 1 #all fall into same category
                    coord.append(mother.index.values)
                if daughter1_ != mother_:
                    daughters = daughters + 1 #only daughters behave similar
                    partial.append(mother.index.values)
            if daughter1_ != daughter2_:
                uncoordinated = uncoordinated + 1 #no correlation
                if mother_ == daughter1_:
                    uncoord.append(daughter1.index.values)
                if mother_ == daughter2_:
                    uncoord.append(daughter2.index.values)
        del daughter1, daughter2, mother
    return daughters, uncoordinated, coordinated, uncoord, partial, coord
    

def Coordination_Plots(data, user, filename, Feat=[], title="Coordination of features.png"):
    if len(Feat)==0:
        Feat = list(data)[1:]
    plt.ioff()
    Analysis = pd.DataFrame(np.zeros((len(Feat),3)), index=Feat, columns = ['Coordination','Partial','Divergence'])
    for F in Feat:
        p,n,c, u_, p_, c_ =  Coordination(filename, data, F)
        Analysis.loc[F,'Partial'] = p
        Analysis.loc[F,'Divergence'] = n
        Analysis.loc[F,'Coordination'] = c
    u_ = pd.Series(u_)
    p_ = pd.Series(p_)
    c_ = pd.Series(c_)
    #print u_
    #print p_
    #print c_
    ax = Analysis.plot(kind='bar', stacked='true', fontsize=6)
    fig = ax.get_figure()    
    #fig.title("Coordination in smallest trees", fontsize=20) #CHECK FOR CORRECT SYNTAX!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    fig.savefig(ospath.join(user,title), format="png", dpi=500) 
    fig.clf()


def Coordination_Table(filename, user, Feat=[]):
    data = pd.read_csv(ospath.join("points.csv"), index_col=1)
    title="Coordination of Features.png"
    Coordination_Plots(data,filename,user,Feat,title)
