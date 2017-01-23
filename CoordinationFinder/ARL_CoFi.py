# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 16:35:37 2016

@author: Sir Thomas
"""

from changepoint_methods import CP_loc
import pandas as pd
import numpy as np
import math

def Check_Difference(cellNr,tree,CP,data, absTime,TF,Columns):
    m = data.loc[(data[Columns[0]]==cellNr) & (data[Columns[1]]==str(tree)), [absTime,str(TF)]]
    cusum, a, b, ttest, c = CP_loc(np.array(m[str(TF)], dtype=pd.Series).tolist())
    returnc = 0
    returnt = 0
    if cusum > 0:
        if CP-m.iloc[0,0]-cusum >= 0:
            returnc = CP-m.iloc[0,0]-cusum #true positive, delay
        if CP-m.iloc[0,0]-cusum < 0:
            returnc= -1 #false positive
    if cusum == 0 and CP > 0:
        returnc = -2 #false negative
    if ttest > 0: #we tretat T-Test here the same as CUSUM. I.E. when TTest detects a change before the user indicates a CP, it will be treated as a FP
        if CP-m.iloc[0,0]-ttest >= 0:
            returnt = CP-m.iloc[0,0]-ttest
        if CP-m.iloc[0,0]-ttest < 0:
            returnt = -1
    if ttest == 0 and CP > 0:
        returnt = -2
    return returnc,returnt #gives distance, if positive. -1 for false positive, -2 for false negative
        
            

def ARL(list_of_trees,TFs,cells,Columns,absTime,voters):
    ARL1_C = [] #ARL1
    ARL1_T = []
    CF_C =[] #combined false positive and false negative
    CF_T = []
    for v in voters:
        c = 0
        t = 0
        normalization_C = 1
        normalization_T = 1 #to not devide by zero
        FPRC = 0
        FPRT = 0
        NT = 1
        NC = 1
        for j in TFs: 
            test1 = list_of_trees.loc[:,['cellNr','tree',str(v)+str(j)+str(1)]]
            test2 = list_of_trees.loc[:,['cellNr','tree',str(v)+str(j)+str(2)]]
            test3 = list_of_trees.loc[:,['cellNr','tree',str(v)+str(j)+str(3)]]
            for i,row in test1.iterrows():
                if math.isnan(np.array(row)[2]) != True:
                    tempc, tempt = Check_Difference(int(np.array(row)[0]),np.array(row)[1],np.array(row)[2],cells,absTime,j,Columns)
                    if tempc >= 0:
                        c = c + tempc
                        normalization_C = normalization_C + 1
                    if tempt >= 0.:
                        t = t + tempt
                        normalization_T = normalization_T + 1
                    if tempc < 0:
                        FPRC = FPRC + 1
                        NC = NC + 1
                    if tempt < 0:
                        FPRT = FPRT + 1
                        NT = NT + 1
            for i,row in test2.iterrows():
                if math.isnan(np.array(row)[2]) != True:
                    tempc, tempt = Check_Difference(int(np.array(row)[0]*2),np.array(row)[1],np.array(row)[2],cells,absTime,j,Columns)
                    if tempc >= 0:
                        c = c + tempc
                        normalization_C = normalization_C + 1
                    if tempt >= 0.:
                        t = t + tempt
                        normalization_T = normalization_T + 1
                    if tempc < 0:
                        FPRC = FPRC + 1
                        NC = NC + 1
                    if tempt < 0:
                        FPRT = FPRT + 1
                        NT = NT + 1
            for i,row in test3.iterrows():
                if math.isnan(np.array(row)[2]) != True:
                    tempc, tempt = Check_Difference(int(np.array(row)[0]*2+1),np.array(row)[1],np.array(row)[2],cells,absTime,j,Columns)
                    if tempc >= 0:
                        c = c + tempc
                        normalization_C = normalization_C + 1
                    if tempt >= 0.:
                        t = t + tempt
                        normalization_T = normalization_T + 1
                    if tempc < 0:
                        FPRC = FPRC + 1
                        NC = NC + 1
                    if tempt < 0:
                        FPRT = FPRT + 1
                        NT = NT + 1
        CF_C.append(float(FPRC)/NC)
        CF_T.append(float(FPRT)/NT)
        ARL1_C.append(float(c)/normalization_C)
        ARL1_T.append(float(t)/normalization_T)
    return ARL1_C,ARL1_T,CF_C,CF_T
            
    