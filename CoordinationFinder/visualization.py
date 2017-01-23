# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 16:00:07 2016

@author: Sir Thomas
"""

import pandas as pd
import matplotlib.pyplot as plt
import os.path as ospath
import numpy as np
import matplotlib

from attributes import Draw_Histograms


def Plots():
    data=pd.read_csv(ospath.join("points.csv"), index_col=0)
    Feat = list(data)[1:]
    plt.ioff()
    for a in Feat:
        for b in Feat:
            if a != b:
                plt.scatter(data.loc[:, a], data.loc[:,b])
                plt.xlabel(a)
                plt.ylabel(b)
                plt.savefig(ospath.join("Plots",a + " vs " + b + ".png"), format="png") 
                plt.clf()
    
from KMeans import Clean_Data
            
def BoxPlots(user, Feat = []):
    data=pd.read_csv(ospath.join("points.csv"), index_col=0)
    if len(Feat) == 0:
        Feat = list(data)[1:]
    plt.ioff()
    matplotlib.rc('font', size=3)
    data,k = Clean_Data(data,Feat,0)
    x = np.arange(len(Feat)+1)
    plt.boxplot(data)
    plt.xticks(x, [0]+Feat, rotation=30, fontsize=7)
    plt.xlabel("Features", fontsize=12)
    plt.ylabel("whitened range of Data", fontsize=12)
    plt.title("Boxplots of Features", fontsize=20)
    plt.savefig(ospath.join(user, "Boxplots of Features.png"), format="png", dpi=500) 
    plt.clf()
    matplotlib.rc('font', size=12)

def BarPlot(user, Heatc, Heatd, Heatp, Feat=[]):
    C = list(np.diag(Heatc.loc[Feat,Feat]))
    P = list(np.diag(Heatp.loc[Feat,Feat]))
    D = list(np.diag(Heatd.loc[Feat,Feat]))
    print [float(i) for i in C]
    C = [float(i) for i in C]
    P = [float(i) for i in P]
    D = [float(i) for i in D]
    x = np.arange(len(Feat))
    
    plt.bar(x,C, color='g', label='coordinated')
    plt.bar(x,P, bottom=C,color='b', label='partial')
    plt.bar(x,D, bottom=[i+j for i,j in zip(C,P)], color='r', label='divergent')
    plt.xticks(x+0.1, Feat, rotation=30, fontsize=7)
    plt.legend(loc='best')
    plt.savefig(ospath.join(user, "BarplotFeatures.png"), format="png", dpi=500) 
    plt.clf()
    
def Save_Histograms(user, bins=0, Feat = []):
    plt.ioff()
    Draw_Histograms(user, bins, Feat)

