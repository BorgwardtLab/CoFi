# -*- coding: utf-8 -*-
"""
Created on Mon May 02 18:43:56 2016

@author: Sir Thomas
"""

import sys
from PyQt4.QtGui import *
from PyQt4.QtCore import *
from PyQt4.uic import *

import matplotlib as mpl
import matplotlib.patches as mpatches
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)


import pandas as pd
import numpy as np
import os as os
import platform
import time
import math
import csv

#from dtw import dtw

from features import Extract_Features #Extract_Features(filename)
from KMeans import Show_KMeans
from DTW import DTWDistance #Show_KMeans('PCA', 0, filename,Feat, Bias,Tree=1), DTWDistance(s1, s2)
from feature_selection import Stability_Analysis #F,s = Stability_Analysis('random') 
from alltrees import *
from changepoint_methods import CP_loc #returns  CUSUM,T_Test
from visualization import Save_Histograms, BoxPlots, BarPlot
from ARL_CoFi import ARL #return ARL_C, ARL_T according to users input
from Corr_clustering import Feature_Clustering

#Global variables for whole data
file_path = "NanogKlf4nucmemcomma.csv"
cells = [] #Feature space, columns allowed, D1, D2, ntree, ptree label indices to the resprective trees
data = [] #Original data
access =[] #name of columns of absulute time and protein intensities in data
Columns = [] #all columns of data which are taken into account by user
votes = [] #coordination votes by users, contains list of trees
Feat = [] #List of names of Features
Summary = [] #Summary of coordination of trees wrt all single features
Heat = pd.DataFrame() #heat plot of pairwise coordination/partiallity/divergence
Heatc = pd.DataFrame() #only coordination
Heatp = pd.DataFrame()
Heatd = pd.DataFrame()
heatdisp = 0 #which heat map is currently on display
user = "test_person"

#global varibale for location of mother cell in the data
counter=0 #the index for cells an d votes are the start indices from the origrinal data
incomplete = [0,0] #1 for taking incomplete curves into account (surface cells in the tree), 1 for excluding last two data points

#representative tracnes for clusters
rep_cellNr = []
rep_tree = []
        
#For plotting
patch = []
legend = ['Nanog','Klf4'] #blue before red
label = ['CUSUM', 'T-Test', 'EWMA', 'QChart', 'BinSeg'] 
lc = ['b-', 'r-', 'g-']
l__c = ['bo', 'ro', 'go']
l_c= ['blue','red','green']

#For CP voting
TF = "" #current TF
subcurve = 1 #1:M, 2:D1, 3:D2
CP = 0
normalization = 0
np.random.seed(42)
np.random.RandomState(0)

def Wipe_all():
    #Global variables for whole data
    global file_path
    global cells
    global data
    global access
    global Columns
    global votes
    global Feat
    global Summary
    global Heat
    global Heatc
    global Heatp
    global Heatd
    global heatdisp
    global user
    global normalization
    global rep_cellNr
    global rep_tree
    rep_cellNr = []
    rep_tree = []
    normalization = 0
    file_path = "NanogKlf4nucmemcomma.csv"
    cells = [] #Feature space, columns allowed, D1, D2, ntree, ptree label indices to the resprective trees
    data = [] #Original data
    access =[] #name of columns of absulute time and protein intensities in data
    Columns = [] #all columns of data which are taken into account by user
    votes = [] #coordination votes by users, contains list of trees
    Feat = [] #List of names of Features
    Summary = [] #Summary of coordination of trees wrt all single features
    Heat = pd.DataFrame() #heat plot of pairwise coordination/partiallity/divergence
    Heatc = pd.DataFrame() #only coordination
    Heatp = pd.DataFrame()
    Heatd = pd.DataFrame()
    heatdisp = 0 #which heat map is currently on display
    user = "test_person"
    global counter
    global incomplete
    #global varibale for location of mother cell in the data
    counter=0 #the index for cells an d votes are the start indices from the origrinal data
    incomplete = [0,1] #1 for taking incomplete curves into account (surface cells in the tree)
    #For plotting
    global patch
    global legend
    global label
    global lc
    global l_c
    global l__c
    patch = []
    legend = ['Nanog','Klf4'] #blue before red
    label = ['CUSUM', 'T-Test', 'EWMA', 'QChart', 'BinSeg']
    lc = ['b-', 'r-', 'g-']
    l__c = ['bo', 'ro', 'go']
    l_c= ['blue','red','green']
    global TF
    global subcurve
    global CP
    TF = "" #current TF
    subcurve = 1 #1:M, 2:D1, 3:D2
    CP = 0
    np.random.seed(42)
    np.random.RandomState(0)
    
#From here on there are no system relevant functions
def talk(string): #talk to user
    t.label.setText(string)
    t.activateWindow()
    t.show()
    QApplication.processEvents()
def DisplaySpecs():
    Text = "<table border=\"0\" style=\"width:50%\"><tr>"
    for i in Feat:
        coord= Summary.loc[counter,str(i)]
        Text = Text + "<tr><td>" + synonym_to_word(i) + "</td><td><font color=\"" + color(coord) + "\">" + coord + "</font></td></tr>"
    Text = Text + "</table></br>"
    #Text = Text + str(Summary.loc[counter,'cellNr']) + " " + str(Summary.loc[counter,'tree'])
    MD1,MD2,D1D2 = Inner_prod()
    #Text = Text + "Dynamic time wrapping inner distance:<br/>Mother-Daughter2: " + str(np.round(DTWDistance(np.array(m)[:,1:].reshape(-1,1), np.array(d1)[:,1:].reshape(-1,1)), decimals=2)) + "<br/> Mother-Daughter1: " + str(np.round(dtw(np.array(m)[:,1:].reshape(-1,1), np.array(d2)[:,1:].reshape(-1,1), dist=lambda x, y: np.linalg.norm(x - y, ord=1))[0], decimals=2))+ "<br/> Daughter1-Daughter2: " + str(np.round(dtw(np.array(d1)[:,1:].reshape(-1,1), np.array(d2)[:,1:].reshape(-1,1), dist=lambda x, y: np.linalg.norm(x - y, ord=1))[0], decimals=2))    
    #Text = Text + "Dynamic time wrapping inner distance:<br/>Mother-Daughter2: " + str(MD1) + "<br/> Mother-Daughter1: " + str(MD2)+ "<br/> Daughter1-Daughter2: " + str(D1D2)
    talk(Text)
def synonym_to_word(Name,rep=True):
    if rep == True:
        if 'dereg' in Name:
            return Name.replace('dereg','Detrended non-linearity ')
        elif 'ereg' in Name:
            return Name.replace('ereg','Non-linearity ')
        elif 'dljungbox' in Name:
            return Name.replace('dljungbox','Detrended autocorrelation ')
        elif 'ljungbox' in Name:
            return Name.replace('ljungbox','Autocorrelation ')
        elif 'd1mom' in Name:
            return Name.replace('d1mom','Detrended mean ')
        elif '1mom' in Name:
            return Name.replace('1mom','Mean ')
        elif 'd2mom' in Name:
            return Name.replace('d2mom','Detrended variance ')
        elif '2mom' in Name:
            return Name.replace('2mom','Variance ')
        elif 'd3mom' in Name:
            return Name.replace('d3mom','Detrended skewness ')
        elif '3mom' in Name:
            return Name.replace('3mom','Skewness ')
        elif 'd4mom' in Name:
            return Name.replace('d4mom','Detrended kurtosis ')
        elif '4mom' in Name:
            return Name.replace('4mom','Kurtosis ')
        elif 'sens' in Name:
            return Name.replace('sens','Sensitivity ')
        elif 'end' in Name:
            return Name.replace('end','Endpoint ')
        elif 'inc' in Name:
            return Name.replace('inc','Increase ')
        elif 'length' in Name:
            return Name.replace('length','Length ')
        elif 'change' in Name:
            return Name.replace('change','Change-point ')
        else:
            return Name
    if rep != True:
        if 'dereg' in Name:
            return 'Detrended non-linearity '
        elif 'ereg' in Name:
            return 'Non-linearity '
        elif 'dljungbox' in Name:
            return 'Detrended autocorrelation '
        elif 'ljungbox' in Name:
            return 'Autocorrelation '
        elif 'd1mom' in Name:
            return 'Detrended mean '
        elif '1mom' in Name:
            return 'Mean '
        elif 'd2mom' in Name:
            return 'Detrended variance '
        elif '2mom' in Name:
            return 'Variance '
        elif 'd3mom' in Name:
            return 'Detrended skewness '
        elif '3mom' in Name:
            return 'Skewness '
        elif 'd4mom' in Name:
            return 'Detrended kurtosis '
        elif '4mom' in Name:
            return 'Kurtosis '
        elif 'sens' in Name:
            return 'Sensitivity '
        elif 'end' in Name:
            return 'Endpoint '
        elif 'inc' in Name:
            return 'Increase '
        elif 'length' in Name:
            return 'Length '
        elif 'change' in Name:
            return 'Change-point '
        else:
            return Name
def get_checked_Feat(): #generate list of features the user has selected
    global normalization
    F = []
    if d.length.isChecked() == True:
        F.append('length')
    #if d.fate.isChecked() == True:
    #    F.append('fate')
    if d.change.isChecked() == True:
        F.append('change')
    #if d.mchange.isChecked() == True:
    #    F.append('mchange')
    if d.osc.isChecked() == True:
        F.append('2mom')
    if d.skew.isChecked() == True:
        F.append('3mom')
    if d.kurt.isChecked() == True:
        F.append('4mom')
    if d.end.isChecked() == True:
        F.append('end')
    if d.inc.isChecked() == True:
        F.append('inc')
    #if d.trend.isChecked() == True:
    #    F.append('reg')
    if d.nolin.isChecked() == True:
        F.append('ereg')
    if d.norm.isChecked() == True:
        normalization = 1
    #if d.rate.isChecked() == True:
    #    F.append('rate')
    if d.mean.isChecked() == True:
        F.append('1mom')
    if d.ljungbox.isChecked() == True:
        F.append('ljungbox')
    if d.sens.isChecked() == True:
        F.append('sens')
    #if d.rel_loc.isChecked() == True:
    #    if 'change' in F:
    #        F.append('loc')
    if d.detrend.isChecked() == True:
        F.append('d1mom')
        F.append('d2mom')
        F.append('d3mom')
        F.append('d4mom')
        F.append('dljungbox')
        F.append('dereg')
    return F, [int(d.include_surface.isChecked()), int(d.LP.isChecked())]
def Out_Analysis():
    text = "Previous users checked this: "
    var = get_votes(votes,counter)
    text = text+str(var[0]) +" coordination | " + str(var[1]) +" partial | " + str(var[2]) + " divergent <br/>"
    return text
def color(coord):
    if coord == "coordinated":
        return 'green'
    if coord == "partial":
        return 'blue'
    if coord == "divergent":
        return "red"
def Display(): #Use HTML
    coord = Summary.loc[counter,'Cluster']
    cellNr = int(cells.loc[counter, 'cellNr'])
    tree = cells.loc[counter,'tree']
    Text = "<font color=\"" + color(coord) + "\" size =6> Cell Nr " + str(cellNr) + " in line " + str(tree) + " is " + coord + ". </font><br/>"
    w.Description.setText(Text)
def Analysis():
    #extract all users and run precision and recall on it
    #meta analysis on all "true" true positives and all "possible" true_positives and the average
    if len(data) != 0:
        splash.show()
        app.processEvents()
        Text = "<p>Precision and Recall of the algorithm with respect to divergent subtrees <table border=\"0\" ><tr>"
        Text = Text + "<tr><td>" + "user" + "</td><td>" + "precision" + "</td><td>" +  "Recall" + "</td></tr>"
        all_users = get_voters(votes, access[1:])
        for x in all_users:
            pre, rec = Precision_Recall([x])
            Text = Text + "<tr><td>" + str(x) + "</td><td>" +  str(np.round(pre,decimals=2)) + "</td><td>" +  str(np.round(rec,decimals=2)) + "</td></tr>"
        pre, rec = Precision_Recall(all_users)
        Text = Text + "<tr><td>" + "<b>average</b>" + "</td><td>" +  str(np.round(pre,decimals=2)) + "</td><td>" +  str(np.round(rec,decimals=2)) + "</td></tr>"
        all_u, true_u = get_uncoord(votes, access[1:])
        pre, rec =Precision_Recall(all_users,all_u)
        Text = Text + "<tr><td>" + "<b>OR</b>" + "</td><td>" +  str(np.round(pre,decimals=2)) + "</td><td>" +  str(np.round(rec,decimals=2)) + "</td></tr>"
        pre, rec =Precision_Recall(all_users,true_u)
        Text = Text + "<tr><td>" + "<b>AND</b>" + "</td><td>" +  str(np.round(pre,decimals=2)) + "</td><td>" +  str(np.round(rec,decimals=2)) + "</td></tr></table><\p>"
        
        #cp output
        ARL_C,ARL_T,CF_C,CF_T = ARL(votes,access[1:],data,Columns,access[0],all_users)
        #insert h-line and more
        Text = Text + "<p><br/> ARL1 and FP&FN for CUSUM and T-test according to different users Inputs<table border=\"0\" ><tr>"
        Text = Text + "<tr><td>" + "user" + "</td><td>" + "ARL1 T-Test" + "</td><td>"+ "FPFN T-Test" + "</td><td>" + "ARL1 T-Test" + "</td><td>"+ "FPFN T-Test" + "</td></tr>"
        i = 0
        while i < len(all_users):
            Text = Text + "<tr><td>" + str(all_users[i]) + "</td><td>" +  str(np.round(ARL_T[i],decimals=2)) + "</td><td>" +  str(np.round(CF_T[i],decimals=2)) + "</td><td>"  +  str(np.round(ARL_C[i],decimals=2)) + "</td><td>" +  str(np.round(CF_C[i],decimals=2)) + "</td></tr>"
            i = i + 1
        app.processEvents()
        splash.finish(w)
        talk(Text+"</table><\p>") 
def About():
    html = "<p><b>Coordination finder 1.2</b></p>"
    html = html+"<p>by Thomas Gumbsch</p>"
    html = html+ "<p>This program was part of a master thesis in the <a href=\"www.bsse.ethz.ch/mlcb \">MLCB lab</a> of ETHZ. <br> Supervisor: Dr. Dean Bodenham, Prof. Dr. Karsten Borgwardt <br> The full thesis can be found on the <a href=\"www.tgumbsch.ethz.ch\">Homepage</a></p>"
    a.label.setText(html)
    a.label.setTextFormat(Qt.RichText)
    a.label.setTextInteractionFlags(Qt.TextBrowserInteraction)
    a.label.setOpenExternalLinks(True)
    a.exec_()
def Documentation():
    if platform.system() == 'Windows':  
        os.system("Documentation.pdf")
    if platform.system() != 'Windows':  
        os.system("open Documentation.pdf")
def Inner_prod(): #DTW using sum of all proteins
    m, d1, d2 = Curves()
    k = 0
    MD1 = 0
    MD2 = 0
    D1D2 = 0
    while k < len(access)-1:
        MD1 = MD1 + np.round(DTWDistance(list(m.values[:,k+1]),list(d1.values[:,k+1])), decimals=2)
        MD2 = MD2 + np.round(DTWDistance(list(m.values[:,k+1]), list(d2.values[:,k+1])), decimals=2)
        D1D2 = D1D2 + np.round(DTWDistance(list(d2.values[:,k+1]), list(d1.values[:,k+1])), decimals=2)
        k = k + 1
    return MD1, MD2,D1D2
def Curves(): #return curves for plotting according to the conter
    cellNr = int(cells.loc[counter, 'cellNr'])
    tree = cells.loc[counter,'tree']
    cellm = data.loc[(data[Columns[0]]==cellNr) & (data[Columns[1]]==str(tree)), access]
    celld1 = data.loc[(data[Columns[0]]==cellNr*2) & (data[Columns[1]]==str(tree)), access]
    celld2 = data.loc[(data[Columns[0]]==cellNr*2+1) & (data[Columns[1]]==str(tree)), access]
    return cellm, celld1, celld2  
def Plotting(): #plot in widget, mother has one widget, daughters share one
    m, d1, d2 = Curves()
    TFs = list(m.columns[1:])
    try:
        w.daughterl.removeWidget(w.canvasd) #detach so a new canvas can appear
        w.motherl.removeWidget(w.canvasm)
    except:
        test = 1
    c,c1,c2 = Feat_value('Cluster')
    bgcolors = ['papayawhip', 'lavender', 'azure', 'lightcyan', 'grey', 'lightblue', 'yellow']
    fig1 = Figure() 
    fig1.set_facecolor('white')
    ax1f1 = fig1.add_subplot(111, axisbg = bgcolors[c])
    ax1f1.set_xlabel('Time')
    ax1f1.set_ylabel('Intensity')
    ax1f1.set_title('Mother')
    fig2 = Figure()
    fig2.set_facecolor('white')
    ax1f2 = fig2.add_subplot(212,axisbg = bgcolors[c1])
    ax1f2.set_xlabel('Time')
    ax1f2.set_ylabel('Intensity')
    ax1f2.set_title('Daughter2')
    ax2f2 = fig2.add_subplot(211, axisbg = bgcolors[c2])
    i = 0
    while i < len(TFs): #plot all TFs in the axis
        ax1f1.plot(m[access[0]].values,m[TFs[i]].values, lc[i])
        ax1f2.plot(d1[access[0]].values,d1[TFs[i]].values, lc[i])
        ax2f2.plot(d2[access[0]].values,d2[TFs[i]].values, lc[i])
        i = i +1
    ax2f2.set_xlabel('Time')
    ax2f2.set_ylabel('Intensity')
    ax2f2.set_title('Daughter1')
    ax1f1.legend(patch, TFs, loc='best')
    ax1f2.legend(patch, TFs, loc='best')
    ax2f2.legend(patch, TFs, loc='best')
    i = 0
    while i < len(TFs):
        if 'change'+str(TFs[i]) in Feat: #plot changepoints
            
            Diff = [l-k for k, l in zip(np.array(m[str(TFs[i])], dtype=pd.Series).tolist()[:-1], np.array(m[str(TFs[i])], dtype=pd.Series).tolist()[1:])] 
            mcn, e,q, mtn,b = CP_loc(Diff)
            ax1f1.plot((m.iloc[mtn,0]), (m.iloc[mtn, i+1]), l__c[i])
            if mtn != 0:
                ax1f1.axes.annotate(label[1], xy = (m.iloc[mtn,0], m.iloc[mtn, i+1]), xytext = (-20, 20), textcoords = 'offset points', ha = 'right', va = 'bottom',bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))        
            
            Diff = [l-k for k, l in zip(np.array(d1[str(TFs[i])], dtype=pd.Series).tolist()[:-1], np.array(d1[str(TFs[i])], dtype=pd.Series).tolist()[1:])] 
            mcn, e,q, mtn,b = CP_loc(Diff)
            ax1f2.plot((d1.iloc[mtn,0]), (d1.iloc[mtn, i+1]), l__c[i])
            if mtn != 0:
                ax1f2.axes.annotate(label[1], xy = (d1.iloc[mtn,0], d1.iloc[mtn, i+1]), xytext = (-20, 20), textcoords = 'offset points', ha = 'right', va = 'bottom',bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))        

            Diff = [l-k for k, l in zip(np.array(d2[str(TFs[i])], dtype=pd.Series).tolist()[:-1],np.array(d2[str(TFs[i])], dtype=pd.Series).tolist()[1:])]             
            mcn, e,q, mtn,b = CP_loc(Diff)
            ax2f2.plot((d2.iloc[mtn,0]), (d2.iloc[mtn, i+1]), l__c[i])
            if mtn != 0:
                ax2f2.axes.annotate(label[1], xy = (d2.iloc[mtn,0], d2.iloc[mtn, i+1]), xytext = (-20, 20), textcoords = 'offset points', ha = 'right', va = 'bottom',bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))        
        i = i + 1
    fig2.tight_layout()
    w.canvasd = FigureCanvas(fig2)
    w.daughterl.addWidget(w.canvasd)
    w.canvasd.setParent(w.daughter)
    w.canvasd.setFocusPolicy( Qt.ClickFocus )
    w.canvasd.setFocus()
    cid = w.canvasd.mpl_connect('button_press_event', daughter_click) #move to one of the daughter cells
    w.canvasm = FigureCanvas(fig1)
    w.motherl.addWidget(w.canvasm)  
    w.canvasm.setParent(w.mother)
    w.canvasm.setFocusPolicy( Qt.ClickFocus )
    w.canvasm.setFocus()
    cim = w.canvasm.mpl_connect('button_press_event', mother_click)#move to top of tree
    w.canvasm.draw()
    w.canvasd.draw()
    Draw_Tree()
def Plot_curve(also_show=[],other_user=[]):
    m, d1, d2 = Curves()
    try:
        v.Changepointl.removeWidget(v.canvasc) #detach so a new canvas can appear
    except:
        test = 1
    fig = Figure()
    fig.set_facecolor('white')
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Intensity')
    ax1.set_title('Click on location of Changepoint and press Next')
    if subcurve == 1:
        ax1.plot(m[str(access[0])].values,m[str(TF)].values) 
        a = 0
        for i in also_show:
            if i != m.index.values[0]:
                ax1.annotate(label[a], xy = (m.loc[i,str(access[0])], m.loc[i,str(TF)]), xytext = (-20, 20), textcoords = 'offset points', ha = 'right', va = 'bottom',bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))    
            a = a + 1
    elif subcurve == 2:
        ax1.plot(d1[str(access[0])].values,d1[str(TF)].values)
        a = 0
        for i in also_show:
            if i != d1.index.values[0]:
                ax1.annotate(label[a], xy = (d1.loc[i,str(access[0])], d1.loc[i,str(TF)]), xytext = (-20, 20), textcoords = 'offset points', ha = 'right', va = 'bottom',bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))    
            a = a + 1
    elif subcurve == 3:
        ax1.plot(d2[str(access[0])].values,d2[str(TF)].values)
        a = 0
        for i in also_show:
            if i != d2.index.values[0]:
                ax1.annotate(label[a], xy = (d2.loc[i,str(access[0])], d2.loc[i,str(TF)]), xytext = (-20, 20), textcoords = 'offset points', ha = 'right', va = 'bottom',bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))    
            a = a + 1
    for i in other_user:
        ax1.axvline(i, color='b', linestyle='-')
    if CP != 0:
        ax1.axvline(CP, color='k', linestyle='--')
    v.canvasc = FigureCanvas(fig)
    v.Changepointl.addWidget(v.canvasc)
    v.canvasc.setParent(v.Changepoint)
    v.canvasc.setFocusPolicy( Qt.ClickFocus )
    v.canvasc.setFocus()
    ci1 = v.canvasc.mpl_connect('button_press_event', cp_click) #move to one of the daughter cells
    v.canvasc.draw()
def cp_click(event):
    global CP
    if CP != 0:
        CP = 0
    elif CP == 0:
        CP = event.xdata 
    Plot_curve()
def mother_click(event):
    UP()
def daughter_click(event):
    if event.y > float(w.canvasd.height())/2: #half of the widget size
        D2()
    elif event.y <=float(w.canvasd.height())/2:
        D1()
def vPlotting(): #plot cuves on the voting dialog, no changepoints and bg color accoding to cluster here!
    m, d1,d2 = Curves()
    TFs = list(m.columns[1:])
    try:
        v.daughterl.removeWidget(v.canvasdv)
        v.motherl.removeWidget(v.canvasmv)
    except:
        m, d1,d2 = Curves()
        
    fig1v = Figure()
    fig1v.set_facecolor('white')
    ax1f1v = fig1v.add_subplot(111)
    ax1f1v.set_xlabel('Time')
    ax1f1v.set_ylabel('Intensity')
    ax1f1v.set_title('Mother')
    fig2v = Figure()
    fig2v.set_facecolor('white')
    ax1f2v = fig2v.add_subplot(212)
    ax1f2v.set_xlabel('Time')
    ax1f2v.set_ylabel('Intensity')
    ax1f2v.set_title('Daughter1')
    ax2f2v = fig2v.add_subplot(211)
    ax2f2v.set_xlabel('Time')
    ax2f2v.set_ylabel('Intensity')
    ax2f2v.set_title('Daughter2')
    i = 0
    while i < len(TFs):
        ax1f1v.plot(m[access[0]].values,m[TFs[i]].values, lc[i])
        ax1f2v.plot(d1[access[0]].values,d1[TFs[i]].values, lc[i])
        ax2f2v.plot(d2[access[0]].values,d2[TFs[i]].values, lc[i])
        i = i +1
    ax1f1v.legend(patch, TFs, loc='best')
    ax1f2v.legend(patch, TFs, loc='best')
    ax2f2v.legend(patch, TFs, loc='best')
    v.canvasdv = FigureCanvas(fig2v)
    v.daughterl.addWidget(v.canvasdv)
    v.canvasmv = FigureCanvas(fig1v)
    v.motherl.addWidget(v.canvasmv)
    v.canvasmv.draw()
    v.canvasdv.draw()
def Next_Plot(event): #cycle through the heat plots
    global heatdisp
    if heatdisp != 0:
        try:
            w.plotl.removeWidget(w.canvash)
        except:
            print("This should not be possible")
        if heatdisp == "divergent":
            heatdisp = "Correlation"
        elif heatdisp == "Correlation":
            heatdisp = "overview"
        elif heatdisp == "partial":
            heatdisp = "divergent" 
        elif heatdisp == "coordinated":
            heatdisp = "partial"
        elif heatdisp == "overview":
            heatdisp = "coordinated"
        Plot_Heat(heatdisp)  
        
#Navigation
def next_tree(): #neglecting current position, in the tree, find the next tree
    global counter
    if counter >= votes.index.values[-1]:
        counter = votes.index.values[0]
    else:
        counter = int(cells.loc[counter, 'ntree'])
def prev_tree(): #neglecting current position in the tree, find the previous tree
    global counter
    if counter <= votes.index.values[0]:
        counter = votes.index.values[-1]
    else:
        counter = int(cells.loc[counter, 'ptree'])
def Forward(): #Display statring node of next tree
    next_tree()
    Display()
    Plotting()
    DisplaySpecs()
def Back(): #Display starting node of previous tree
    prev_tree() 
    Display()
    Plotting()
    DisplaySpecs()
def D1(): #down through daughter1
    global counter
    counter = int(cells.loc[counter,'D1'])
    Plotting()
    Display()
    DisplaySpecs()
def D2(): #down through daughter2
    global counter
    counter = int(cells.loc[counter,'D2'])
    Plotting()
    Display()
    DisplaySpecs()
def UP(): #go up 
    test = cells.loc[counter, 'tree']
    while cells.loc[counter,'tree'] == test:
        next_tree()
    prev_tree()
    Plotting()
    Display()
    DisplaySpecs()

#Analysis
def Check_Coordination(Feat): #return coordination and colors according to median. TO DO: Account for clustering in dirrerent algorithm
    md = cells.median(axis=0)[Feat] #We split the different features into two categories according to the median
    if md == 1:
        md = 0 #correct for binary features
    mother_, daughter1_, daughter2_ = Feat_value(Feat)
    daughter1_ = (float(daughter1_) >= md)
    daughter2_ = (float(daughter2_) >= md)
    mother_ = (float(mother_) >= md)
    if daughter1_ == daughter2_:
        if daughter1_ == mother_:
            return 'coordinated'
        if daughter1_ != mother_:
            return 'partial'
        if daughter2_ == mother_:
            return 'coordinated'
        if daughter2_ != mother_:
            return 'partial'
    if daughter1_ != daughter2_:
        return 'divergent'
def Check_ExactCoordination(Feat): #Not only looking at median
    mother_, daughter1_, daughter2_ = Feat_value(Feat)
    if daughter1_ == daughter2_:
        if daughter1_ == mother_:
            return 'coordinated'
        if daughter1_ != mother_:
            return 'partial'
        if daughter2_ == mother_:
            return 'coordinated'
        if daughter2_ != mother_:
            return 'partial'
    if daughter1_ != daughter2_:
        return 'divergent'
def Feat_value(Feat):
    cellNr = int(cells.loc[counter,'cellNr'])
    tree = cells.loc[counter,'tree']
    return cells.loc[counter,Feat], cells.loc[(cells[Columns[0]]==cellNr*2) & (cells[Columns[1]]==str(tree)),Feat].values[0], cells.loc[(data[Columns[0]]==cellNr*2+1) & (data[Columns[1]]==str(tree)),Feat].values[0]
def ExportOpen(): #show export Dialog and save files
    if len(data) != 0:
        e.label.setText("The files will be saved in the subfolder <br/>" + str(os.path.abspath('./'+str(user))))
        r = e.exec_()
        QApplication.processEvents()
        splash.show()
        splash.raise_()
        QApplication.processEvents()
        if r == QDialog.Accepted:
            if not os.path.exists(user): #make new folder for new user
                os.makedirs(user)
            if e.Votes.isChecked() == True:
                data.to_csv(os.path.join(user,"Data.csv")) #the original Data
                cells.to_csv(os.path.join(user,"Featrues.csv")) #Feature space
                votes.to_csv(os.path.join(user,"Votes.csv")) #all votes
                Summary.to_csv(os.path.join(user,"Coordination.csv")) #summary of coordination   
            if e.Graphics.isChecked() == True:
                plt.savefig(os.path.join(user,"Clustering.pdf"))  #plt has currently the clustering
                plt.clf()
                """
                #Save plots and corresponding cluster index
                for index,rows in cells.iterrows():
                    single_trace = data.loc[(data[Columns[0]]==int(cells.loc[index, 'cellNr'])) & (data[Columns[1]]==cells.loc[index, 'tree']), access]
                    i = 0                
                    TFs = list(single_trace.columns[1:])
                    while i < len(TFs):
                        x,y =zip(*sorted(zip(single_trace[access[0]].values, single_trace[TFs[i]].values))) 
                        plt.plot(x,y, label=str(TFs[i]))
                        i = i + 1
                    plt.legend(loc='best',fontsize=16)
                    plt.title('Cell number ' + str(rows['cellNr'])+' of tree ' +str(rows['tree']), fontsize=18)
                    plt.xlabel('Time', fontsize=16)
                    plt.ylabel('Intensity', fontsize=16)
                    plt.xticks(fontsize=14)
                    plt.yticks(fontsize=14)
                    if not os.path.exists(os.path.join(user, 'Cluster_' + str(rows['Cluster']))): #make new folder for new cluster
                        os.makedirs(os.path.join(user, 'Cluster_' + str(rows['Cluster'])))
                    if len(single_trace[access[0]].values) > 0:
                        plt.savefig(os.path.join(user, 'Cluster_' + str(rows['Cluster']), str(rows['tree'])+str(rows['cellNr'])+'.pdf'))
                    plt.clf()
                """
                
                #Save representative plots
                check_rep = [[i,j] for i,j in zip(rep_cellNr,rep_tree)]
                print check_rep
                for index,rows in cells.iterrows():
                    if [int(cells.loc[index, 'cellNr']),cells.loc[index, 'tree']] in check_rep:
                        single_trace = data.loc[(data[Columns[0]]==int(cells.loc[index, 'cellNr'])) & (data[Columns[1]]==cells.loc[index, 'tree']), access]
                        i = 0                
                        TFs = list(single_trace.columns[1:])
                        for i in range(len(TFs)):
                            x,y =zip(*sorted(zip(single_trace[access[0]].values, single_trace[TFs[i]].values))) 
                            plt.plot(x,y, label=str(TFs[i]))
                            i = i + 1
                        plt.legend(loc='best',fontsize=16)
                        plt.suptitle("Representative of cluster" + str(rows['Cluster']),fontsize=18)
                        plt.title('Cell number ' + str(rows['cellNr'])+' of tree ' +str(rows['tree']), fontsize=16)
                        plt.xlabel('Time', fontsize=16)
                        plt.ylabel('Intensity', fontsize=16)
                        plt.xticks(fontsize=14)
                        plt.yticks(fontsize=14)
                        if not os.path.exists(os.path.join(user, 'Cluster_' + str(rows['Cluster']))): #make new folder for new cluster
                            os.makedirs(os.path.join(user, 'Cluster_' + str(rows['Cluster'])))
                        if len(single_trace[access[0]].values) > 0:
                            plt.savefig(os.path.join(user, 'Cluster_' + str(rows['Cluster']), str(rows['tree'])+str(rows['cellNr'])+'.pdf'))
                        plt.clf()
                
                Association_Plot()
                Plot_Heat("divergent",1)
                Plot_Heat("partial",1)
                Plot_Heat("coordinated",1)
                Plot_Heat("overview",1)
                #Coordination_Table(file_path,user) #heat maps are better
                Save_Histograms(user,10)
                BoxPlots(user)
                BarPlot(user, Heatc, Heatd, Heatp, Feat)
                
                
                #Save Correlation
                Correlation = cells.loc[:,Feat].corr()
                fig1 = Figure()
                fig1.set_facecolor('white')
                ax1f1 = fig1.add_subplot(111)
                im = ax1f1.imshow(Correlation.values, interpolation='none', cmap=cm.seismic) 
                ax1f1.tick_params(axis = 'both', which = 'all', labelsize = 8)
                x = np.arange(len(Feat))
                ax1f1.set_xticklabels(Feat, minor=False,rotation=30)
                ax1f1.xaxis.tick_bottom()
                ax1f1.xaxis.set_label_position('bottom') 
                ax1f1.xaxis.set_visible(False)
                ax1f1.set_yticks(x, minor=False)
                ax1f1.set_yticklabels(Feat, minor = False)
                for tick in ax1f1.yaxis.get_major_ticks(): #display all ticks on the y-axis (DO NOT CHANGE!!!)
                    tick.label1On = False
                    tick.label2On = True
                ax1f1.set_title('Correlation of Features') #display no ticks on the x-axis (no space)
                cbar_ax = fig1.add_axes([0.05, 0.15, 0.05, 0.7]) #colorbar is on the left side
                fig1.colorbar(im, cax=cbar_ax)
                canvas = FigureCanvas(fig1)
                canvas.print_figure(os.path.join(user, "Correlation_features.pdf"))
                #fig1.savefig(os.path.join(user, "Correlation_features.png"), format="png", dpi=500) 
                
        splash.finish(w)

def Association_Plot():
    for F in Feat:
        if 'change' in F:
            cluster1 = cells.loc[cells['Cluster']==0]
            cluster2 = cells.loc[cells['Cluster']==1]
            a = [int(i) for i in cluster1.loc[:,str(F)].values.tolist()]
            b = [int(i) for i in cluster2.loc[:,str(F)].values.tolist()]
            weightsa = np.ones_like(a)/float(len(a))
            weightsb = np.ones_like(b)/float(len(b))
            plt.hist([a,b], label=['stationary','changing'], color=['orange','blue'], bins=[0,0.5,1], weights=[weightsa,weightsb])
            plt.legend(loc='best',fontsize=18)
            plt.suptitle('Change-point frequency by cluster',fontsize=20)
            plt.xticks([0.3,0.7],['no change-point', 'change-point'], fontsize=18)
            plt.yticks(fontsize=18)
            plt.ylabel('Frequency of change-points', fontsize=18)
            plt.savefig(os.path.join(user,'Cluster_association_'+str(F)+'.pdf'))
            plt.clf()
            plt.close()
        else:
            bp =  cells.boxplot(column= str(F), by='Cluster', patch_artist=True, return_type='dict', whis=[16,84])
            colors = ['Orange', 'Blue', 'azure', 'lightcyan', 'grey', 'lightblue', 'yellow']
            ticks = ['Stationary', 'Changing']
            k = 0
            for i in bp.keys():
                for box in bp[i]['boxes']:
                    box.set(color='#7570b3',linewidth=2)
                    box.set(facecolor='#1b9e77')
                for whisker in bp[i]['whiskers']:
                    whisker.set(color='#7570b3',linewidth=2)
                for cap in bp[i]['caps']:
                    cap.set(color='#7570b3',linewidth=2)
                for median in bp[i]['medians']:
                    median.set(color='#b2df8a',linewidth=2)
                k = k + 1
                
            plaintext = synonym_to_word(F,False)
            plt.suptitle(str(plaintext) + "association with Cluster", fontsize=18)
            plt.title('', fontsize=18)
            plt.ylabel(str(plaintext), fontsize=18)
            plt.xlabel('Cluster', fontsize=18)
            plt.yticks(fontsize=18)
            k = 0
            plt.xticks([1,2],['Stationary', 'Changing'] , fontsize=18)
            for i in  plt.gca().get_xticklabels():
                i.set_color(colors[k]) 
                k = k + 1
            plt.savefig(os.path.join(user,'Cluster_association_'+str(F)+'.pdf'))
            plt.clf()
            plt.close()
    
    

def Make_Coord_DF(): #Create summary. This is a major time issue
    global Summary
    ##################################################
    ##BUG: IF VOTES LOADED WITH DIFFERENT FEATURES####
    ##INCOMPLETE TREES MIGHT MAKE IT TO SUMMARY#######
    ##################################################
    Summary = votes.loc[:,['cellNr','tree']]
    global counter
    save = counter
    
    for index, rows in Summary.iterrows():
        counter = index
        Summary.set_value(index, 'Cluster', Check_ExactCoordination('Cluster'))
        for F in Feat:    
            Summary.set_value(index, str(F), Check_Coordination(F))
        
    counter = save

def prepare_LoadList(columns):
     o.CellNr.clear()
     o.tree.clear()
     o.stopReason.clear()
     o.absoluteTime.clear()
     o.TFs.clear()
     o.norml.clear()
     if normalization == 0:
         o.widget.hide()
     for i in columns:
         itemC = QListWidgetItem(str(i))
         o.CellNr.addItem(itemC)
         itemt = QListWidgetItem(str(i))
         o.tree.addItem(itemt)
         items = QListWidgetItem(str(i))
         o.stopReason.addItem(items)
         itemT = QListWidgetItem(str(i))
         o.absoluteTime.addItem(itemT)
         itemP = QListWidgetItem(str(i))
         o.TFs.addItem(itemP)
         itemN = QListWidgetItem(str(i))
         o.norml.addItem(itemN)
         if (i == 'intNanog' ) and len(access) == 0:
             itemP.setSelected(True)
         elif (i == 'ProteinA' or i == 'ProteinB') and len(access) == 0:
             itemP.setSelected(True)
         elif i == 'cellNr'and len(access) == 0:
             itemC.setSelected(True)
         elif i == 'tree'and len(access) == 0:
             itemt.setSelected(True)
         elif i == 'stopReason'and len(access) == 0:
             items.setSelected(True)
         elif i == 'fate' and len(access) == 0:
             items.setSelected(True)
         elif i == 'absoluteTime'and len(access) == 0:
             itemT.setSelected(True)
         elif i == 'absTime'and len(access) == 0:
             itemT.setSelected(True)
         elif i == 'timepoint'and len(access) == 0:
             itemT.setSelected(True)
         elif i == 'intNucmem'and len(access) == 0:
             itemN.setSelected(True)
def Load(): #Set up new Data
    r = d.exec_() #open windows and let the user specify the way of performing the computation
    if len(str(d.lineEdit.text())) != 0 and '.csv' in str(d.lineEdit.text()) and r == QDialog.Accepted: #only start when we have data and the user choses ok
        Wipe_all()
        global file_path
        global cells
        global data
        global Feat
        global user
        global Summary
        global votes
        global Columns
        global access
        global patch
        global incomplete
        global counter
        global normalization
        global rep_cellNr
        global rep_tree
        Features, incomplete = get_checked_Feat()
        file_path = str(d.lineEdit.text())
        data = pd.read_csv(os.path.join(file_path))
        prepare_LoadList(list(data.columns.values))
        r_ = o.exec_()
        print("-------------1. Reading in Data-------------")
        Columns = []
        Columns = [str(o.CellNr.selectedItems()[0].text()), str(o.tree.selectedItems()[0].text()), str(o.stopReason.selectedItems()[0].text()), str(o.absoluteTime.selectedItems()[0].text()), [str(x.text()) for x in o.TFs.selectedItems()]]
        if r_ == QDialog.Accepted and len(Columns) > 4:
            access = [Columns[3]]
            access.extend(Columns[4])
            i = 1
            while i < len(access): 
                patch.append(mpatches.Patch(color=l_c[i-1], label=access[i])) 
                i = i +1
            QApplication.processEvents()
            splash.show() #create loading splash screen
            splash.raise_()
            splash.activateWindow()
            QApplication.processEvents()
            n_channel = ""
            if normalization == 1:
                n_channel = str(o.norml.selectedItems()[0].text())
            if len(str(d.user.text())) != 0:
                user = str(d.user.text())            
            if d.StabilityAnalysis.isChecked() == True and len(d.no_Feat.text()) != 0:
                No = int(d.no_Feat.text())
                print("-------------2. Extracting Features-------------")
                Extract_Features(file_path,Features,Columns,incomplete,[normalization, n_channel])
                print("-------------3. Perfoming Stability Analysis-------------")
                F,s = Stability_Analysis('random')
                Features = F[:No] #take the first 10 features
                print("-------------4. Preparing For Show-------------")
                k, rep_cellNr, rep_tree = Show_KMeans('PCA', 0, file_path, Features)
            if d.StabilityAnalysis.isChecked() == True and len(d.no_Feat.text()) == 0:
                print("-------------2. Extracting Features-------------")
                Extract_Features(file_path,Features,Columns,incomplete,[normalization, n_channel])
                F = []
                for i in Features:
                    if i == 'length' or i== 'fate':
                        F.append(i)
                    else:
                        for j in Columns[4]:
                            F.append(str(i)+str(j))
                Features = F
                print("-------------3. Preparing For Show-------------")
                k, rep_cellNr, rep_tree = Show_KMeans('PCA', 0, file_path,Columns)
            if d.StabilityAnalysis.isChecked() != True:
                print("-------------2. Extracting Features-------------")
                Extract_Features(file_path,Features,Columns,incomplete,[normalization, n_channel])
                F = []
                for i in Features:
                    if i == 'length' or i== 'fate':
                        F.append(i)
                    else:
                        for j in Columns[4]:
                            F.append(str(i)+str(j))
                Features = F
                print("-------------3. Preparing For Show-------------")
                k, rep_cellNr, rep_tree = Show_KMeans('PCA', 0, file_path, Features, Columns)
            cells = pd.read_csv("tree.csv", index_col=0)
            Feat = list(cells)[:]
            i = 0
            while 'Cluster' != Feat[i] or i > 1000: #Take all Features but the first one and not CLuster and what is located after it (cellNr, tree, )
                i = i+1
            while i < len(Feat):
                Feat.pop(i) 
            TFs = ""
            for i in Columns[4]: #name of the voting file depends on how many TFs the user selected. This has to be done because coordination depends on the number of selected proteins.
                TFs = TFs+str(i)
            if os.path.isfile('votes'+str(TFs)+str(incomplete[0])+str(incomplete[1])+'.csv') == False: #name it also according to features and original filename
                votes = return_min_trees(cells) 
            if os.path.isfile('votes'+str(TFs)+str(incomplete[0])+str(incomplete[1])+'.csv') == True: #read in the votes.csv
                votes = pd.read_csv('votes'+str(TFs)+str(incomplete[0])+str(incomplete[1])+'.csv', index_col=0)
                votes = votes[np.isfinite(votes['cellNr'])]
            print("-------------5. Creating Lookup Dataframe For Coordination-------------")
            Make_Coord_DF()
            Create_Heat(k) #Heat map of feature coordination     
            #counter = int(cells.loc[counter, 'ntree']) ##Delete it!!!!!!!!!!!!!
            counter = votes.index.values[0]
            Forward()
            Draw_Tree() #tree structure in Lines
            splash.finish(w)
def Open_CVS():
    new_file_path = QFileDialog.getOpenFileName()
    if ".csv" in new_file_path:
        d.lineEdit.setText(new_file_path)
#Plot tree on w.plot
def generation(cellnumber):
    return int(math.log(cellnumber,2))
def plotlocation(cellnumber):
    offset = 0
    g = generation(cellnumber)
    if g != 0:
        i = 1
        while i <= g:
            i = i + 1
            offset = offset + 1/float(2**(i))
    return (cellnumber-(2**g))/float(2**g) - offset #evenly spaces cells of one generation between 0 and 1. The offset depends on the generation and makes it look "treelike"
def Line_Subtree(): #draws line for eache subtree starting at lifetime of mothercell and then moving along each branch to the death of the daughter cells
    m, d1,d2 = Curves()
    M, D1,D2 = Feat_value('cellNr')
    y = [plotlocation(M), plotlocation(M),plotlocation(D1),plotlocation(D1),plotlocation(D1),plotlocation(M),plotlocation(D2) ,plotlocation(D2)]
    x = [m[access[0]].values[0], m[access[0]].values[-1], d1[access[0]].values[0],d1[access[0]].values[-1],d1[access[0]].values[0],m[access[0]].values[-1], d2[access[0]].values[0], d2[access[0]].values[-1]]
    coord = Summary.loc[counter,'Cluster']
    c= color(coord)
    return x,y,c
def Draw_Tree(): #plot the current tree in w.tree. has to be redrawn for every change in subtree because of the location of the annotation (current subtree)
    x = [] #collect all points to get min/max for plot
    y = []
    try:
        w.Treel.removeWidget(w.canva)
    except:
        test = 1
    global counter
    temp = counter
    a,b,c = Line_Subtree()
    currenty = b[1]
    currentx = a[1]
    look = cells.loc[counter,'tree']
    f = Figure()
    f.set_facecolor('white')
    ax = f.add_subplot(111)
    for index,rows in cells.iterrows():
        counter = index
        if cells.loc[counter,'allowed'] == True and cells.loc[counter,'tree']== look:
            a,b,c = Line_Subtree()
            ax.plot(a,b,str(c[0])+'-')
            x.extend(a)
            y.extend(b)
    maxx = max(x)
    minx = min(x)
    ax.set_xlim(minx-(maxx-minx)/10, maxx+(maxx-minx)/10)
    ax.set_ylim(min(y)-0.25, max(y)+0.25)
    ax.annotate('you are here', xy = (currentx, currenty), xytext = (-20, 20), textcoords = 'offset points', ha = 'right', va = 'bottom',bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
    ax.get_yaxis().set_visible(False)
    ax.set_xlabel('absolute time')
    if platform.system() != 'Windows':    
        ax.axis('off')
    ax.set_title('Current tree')
    colors = ['green', 'blue', 'red']
    Curvename = ['coordinated','partial','divergent']
    i = 0
    patchc = []
    while i < len(colors): #create legend
        patchc.append(mpatches.Patch(color=colors[i], label=colors[i])) 
        i = i +1
    ax.legend(patchc, Curvename, loc='best')
    w.canva = FigureCanvas(f)
    w.Treel.addWidget(w.canva)
    w.canva.draw()  
    counter = temp
def Create_Heat(k): #create all four heat plots
    global Heat
    global Heatc
    global Heatp
    global Heatd
    global heatdisp
    heatdisp = "overview"
    Heatc = pd.DataFrame(np.zeros((len(Feat),len(Feat))), columns=Feat, index=Feat)
    Heatd = pd.DataFrame(np.zeros((len(Feat),len(Feat))), columns=Feat, index=Feat)
    Heatp = pd.DataFrame(np.zeros((len(Feat),len(Feat))), columns=Feat, index=Feat)
    Heat = pd.DataFrame(np.zeros((len(Feat),len(Feat))), columns=Feat, index=Feat)
    for index,rows in Summary.iterrows():
        for i in Feat:
            for j in Feat:
                if rows[i] == rows[j]:
                    Heat.set_value(i,j,Heat.loc[i,j]+1) 
                    if rows[i] == 'coordinated':
                        Heatc.set_value(i,j,Heatc.loc[i,j]+1) 
                    if rows[i] == 'partial':
                        Heatp.set_value(i,j,Heatp.loc[i,j]+1) 
                    if rows[i] == 'divergent':
                        Heatd.set_value(i,j,Heatd.loc[i,j]+1) 
    Plot_Heat(heatdisp)
    #create figure with legend according to culsters only
    try:
        w.Legend_Clustersl.removeWidget(w.canvasl)
    except:
        test = 1
    bgcolors = ['papayawhip', 'lavender', 'azure', 'lightcyan', 'grey', 'lightblue', 'yellow']
    Clustername = ['Cluster 1','Cluster 2','Cluster 3','Cluster 4', 'Cluster 5','Cluster 6', 'Cluster 7']
    i = 0
    patchc = []
    print k
    while i < k: #k is return value from KMeans, number of clusters
        patchc.append(mpatches.Patch(color=bgcolors[i], label=bgcolors[i])) 
        i = i +1
    fig0 = Figure()#figsize=(3,0.5))
    fig0.legend(patchc, Clustername, loc='center',ncol=7)
    fig0.set_facecolor('white')
    w.canvasl = FigureCanvas(fig0)
    w.Legend_Clustersl.addWidget(w.canvasl)
    w.canvasl.draw()
def Plot_Heat(what, save=0): #plot it in w.plot what=pointer to one of the four
    try:
        w.plotl.removeWidget(w.canvash)
    except:
        test = 1
    
    fig1 = Figure()
    fig1.set_facecolor('white')
    ax1f1 = fig1.add_subplot(111)
    if what == "divergent":
            im = ax1f1.imshow(Heatd.values, interpolation='none', cmap=cm.seismic) 
    elif what == "partial":
            im = ax1f1.imshow(Heatp.values, interpolation='none', cmap=cm.seismic) 
    elif what == "coordinated":
            im = ax1f1.imshow(Heatc.values, interpolation='none', cmap=cm.seismic) 
    elif what == "overview":
            im = ax1f1.imshow(Heat.values, interpolation='none', cmap=cm.seismic) 
    elif what == "Correlation":
            Co = cells.loc[:,Feat] #Feature_Clustering(cells.loc[:,Feat].corr()) #uncomment for more analysis on correlation on real dataset
            im = ax1f1.imshow(Co.corr().values, interpolation='none', cmap=cm.seismic)
    ax1f1.tick_params(axis = 'both', which = 'all', labelsize = 8)
    copy_Feat = Feat
    ticklabels = [synonym_to_word(i) for i in copy_Feat]
    x = np.arange(len(ticklabels))
    """
    if what == "Correlation":
        ax1f1.set_xticklabels(Co.index.tolist(), minor=False,rotation=30)
    else:
    """
    ax1f1.set_xticklabels(ticklabels, minor=False,rotation=30)
    ax1f1.xaxis.tick_bottom()
    ax1f1.xaxis.set_label_position('bottom') 
    ax1f1.xaxis.set_visible(False)
    ax1f1.set_yticks(x, minor=False)
    """
    if what == "Correlation":
        ax1f1.set_yticklabels(Co.index.tolist(), minor = False)
    else:
    """
    ax1f1.set_yticklabels(ticklabels, minor = False)
    for tick in ax1f1.yaxis.get_major_ticks(): #display all ticks on the y-axis (DO NOT CHANGE!!!)
        tick.label1On = False
        tick.label2On = True
        
    if what != "Correlation":
        ax1f1.set_title('Heat score ' + str(what)) #display no ticks on the x-axis (no space)
    elif what == "Correlation":
        ax1f1.set_title('Peason correlation for Features')
    cbar_ax = fig1.add_axes([0.05, 0.15, 0.05, 0.7]) #colorbar is on the left side
    fig1.colorbar(im, cax=cbar_ax)
    w.canvash = FigureCanvas(fig1)
    w.plotl.addWidget(w.canvash)
    w.canvash.setParent(w.plot)
    w.canvash.setFocusPolicy( Qt.ClickFocus )
    w.canvash.setFocus()
    ci = w.canvash.mpl_connect('button_press_event', Next_Plot)
    if save != 0:
        fig1.savefig(os.path.join(user, "Correlation"+str(what)+".png"), format="png", dpi=1200) 
    w.canvash.draw()  
#display voting screen and and collect users votes
def Vote():
    if len(data) != 0:
        splash.show()
        splash.raise_()
        app.processEvents()
        vPlotting()
        global votes
        TFs = ""
        for i in Columns[4]:
            TFs = TFs+str(i)
        if os.path.isfile('votes'+str(TFs)+str(incomplete[0])+str(incomplete[1])+'.csv') == True:
            votes = pd.read_csv('votes'+str(TFs)+str(incomplete[0])+str(incomplete[1])+'.csv', index_col=0)
            print('open existing file')
        if os.path.isfile('votes'+str(TFs)+str(incomplete[0])+str(incomplete[1])+'.csv') == False:
            votes = return_min_trees(cells)
        global TF
        global subcurve
        TF = access[1]
        subcurve = 1
        Plot_curve()
        app.processEvents()
        splash.finish(w)
        app.processEvents()
        r = v.exec_()
        if r == QDialog.Accepted:
            votes.to_csv(os.path.join('votes'+str(TFs)+str(incomplete[0])+str(incomplete[1])+'.csv'))
            talk("Dear "+str(user)+",<br> your votes have been saved to " + 'votes'+str(TFs)+str(incomplete[0])+str(incomplete[1])+'.csv'+"'.") #only save dataframe when suer exits with ok
        if r != QDialog.Accepted:
            talk("Dear "+str(user)+", your votes have been discarded.")
def Vote_Coord():
    global votes
    votes, stuff1,stuff2 = save_tree(votes, counter, user, 1)
    c.vote.setText("coordinated")
    vProgress()
def Vote_Partial():
    global votes
    votes, stuff1, stuff2 = save_tree(votes, counter, user, 2)
    c.vote.setText("partial")
    vProgress()
def Vote_Uncoord():
    global votes
    votes, stuff1, stuff2 = save_tree(votes, counter, user, 3)
    c.vote.setText("divergent")
    vProgress()
def Output_Labels():
    var = get_votes(votes,counter)
    c.ana.setText(Check_ExactCoordination('Cluster'))
    c.coord.setText(str(var[0])+" coordinated")
    c.partial.setText(str(var[1])+" partial")
    c.uncoord.setText(str(var[2])+" divergent")
def vProgress(): #next subtree in voting screen
    global counter
    Output_Labels()
    c.exec_()
    counter = get_tree(votes, user)
    vPlotting()
#Metaanalysis of all votes
def Precision_Recall(users, index =[]): #pass a list!!!
    global counter
    temp = counter
    precision = 0
    recall = 0
    true_pos = 0
    for u in users:
        for i,row in cells.iterrows():
            counter = i
            if row['allowed'] == 1:
                if math.isnan(votes.loc[counter, str(u)]) == False:
                    if Check_ExactCoordination('Cluster')[0] == 'divergent':
                        precision = precision + 1
                    if int(votes.loc[counter, str(u)]) == 3:
                        recall = recall + 1
                        if Check_ExactCoordination('Cluster')[0] == 'divergent':
                            true_pos = true_pos + 1
    if len(index) > 0:
        recall = 0
        true_pos = 0
        for i in index:
            counter = i
            if i in cells.index.values:
                recall = recall + 1
                if Check_ExactCoordination('Cluster')[0] == 'divergent':
                    true_pos = true_pos + 1
    if precision == 0:
        precision = 1
    if recall == 0:
        recall = 1
    precision = float(true_pos)/precision
    recall = float(true_pos)/recall
    counter = temp
    return precision, recall

def Vote_CP():
    global votes
    global subcurve
    global TF
    global counter
    global CP
    m, d1, d2 = Curves()
    #other_users = get_user_votes(votes, TF, access[1:], subcurve,counter) #does not work if he dont have any previous user
    other_users=[]
    if subcurve == 1:
        Diff = [l-k for k, l in zip(np.array(m[str(TF)], dtype=pd.Series).tolist()[:-1],np.array(m[str(TF)], dtype=pd.Series).tolist()[1:])] 
        cusum, e, q, ttest, b = CP_loc(Diff)
        cusum = cusum + m.index.values[0] #CP detection yield only relative coordinates
        e = e + m.index.values[0]
        ttest = ttest + m.index.values[0] #CP detection yield only relative coordinates
        q = q + m.index.values[0]
        b = b + m.index.values[0]
    elif subcurve == 2:
        Diff = [l-k for k, l in zip(np.array(d1[str(TF)], dtype=pd.Series).tolist()[:-1],np.array(d1[str(TF)], dtype=pd.Series).tolist()[1:])] 
        cusum, e, q, ttest, b = CP_loc(Diff)
        cusum = cusum + d1.index.values[0]
        ttest = ttest + d1.index.values[0]
        e = e + d1.index.values[0]
        q = q + d1.index.values[0]
        b = b + d1.index.values[0]
    elif subcurve == 3:
        Diff = [l-k for k, l in zip(np.array(d2[str(TF)], dtype=pd.Series).tolist()[:-1],np.array(d2[str(TF)], dtype=pd.Series).tolist()[1:])] 
        cusum, e, q, ttest, b = CP_loc(Diff)
        cusum = cusum + d2.index.values[0]
        ttest = ttest + d2.index.values[0]
        e = e + d2.index.values[0]
        q = q + d2.index.values[0]
        b = b + d2.index.values[0]
    Plot_curve([int(cusum),int(ttest), int(e),int(q),int(b)], other_users) #plot the reference votes, then pause for 5 s
    votes, stuff1, stuff2 = save_tree(votes,counter,str(user)+str(TF)+str(subcurve),CP)
    v.Next.setEnabled(False)
    app.processEvents()
    time.sleep(5)
    v.Next.setEnabled(True)
    counter, TF, subcurve = get_single_curve(votes,counter,access[1:], user)
    CP = 0
    Plot_curve()
    
def Exit_():
    app.quit()

app = QApplication(sys.argv)
splash_pix = QPixmap('conti.png')
splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)  
splash.setMask(splash_pix.mask())
splash.show()
app.processEvents()
w = loadUi("show.ui")
d = loadUi("load.ui")
v = loadUi("vote.ui")
c = loadUi("compare_analysis.ui")
t = loadUi("talk.ui")
a = loadUi("about.ui")
e = loadUi("export.ui")
o = loadUi("open.ui")

w.connect(w.prev_tree, SIGNAL("clicked()"), Back)
w.connect(w.next_tree, SIGNAL("clicked()"), Forward)
w.connect(w.actionOpen, SIGNAL("triggered()"), Load)
w.connect(w.actionVote, SIGNAL("triggered()"), Vote)
w.connect(w.actionAbout, SIGNAL("triggered()"), About)
w.connect(w.actionDocumentation, SIGNAL("triggered()"), Documentation)
w.connect(w.actionAnalysis, SIGNAL("triggered()"), Analysis)
w.connect(w.actionExport, SIGNAL("triggered()"), ExportOpen)
w.connect(w.actionClose, SIGNAL("triggered()"), Exit_)
#w.showFullScreen()

d.connect(d.Open_CVS, SIGNAL("clicked()"), Open_CVS)
d.lineEdit.setText(file_path)
o.TFs.setSelectionMode(QAbstractItemView.ExtendedSelection)

v.connect(v.coord, SIGNAL("clicked()"), Vote_Coord)
v.connect(v.parcoord, SIGNAL("clicked()"), Vote_Partial)
v.connect(v.uncoord, SIGNAL("clicked()"), Vote_Uncoord)
v.connect(v.Next, SIGNAL("clicked()"), Vote_CP)

app.processEvents()
w.showMaximized()
w.raise_()
splash.finish(w)
sys.exit(app.exec_())