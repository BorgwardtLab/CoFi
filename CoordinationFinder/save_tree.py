# -*- coding: utf-8 -*-
"""
Created on Mon May 02 13:53:00 2016

@author: Sir Thomas
"""
import pandas as pd
import numpy as np
import os as os

def save_tree(data,curves, filename,Columns=['cellNr', 'tree', 'stopReason', 'absoluteTime', ['intNanog','intKlf4']],closest=[]):
    data['cellNr'] = curves.loc[data.index, Columns[0]]
    data['tree'] = curves.loc[data.index, Columns[1]]
    rep_cellNr = data.iloc[closest].cellNr
    rep_tree = data.iloc[closest].tree
    counter = 0
    trees = data
    sLength = len(trees['cellNr'])
    trees['D1'] = pd.Series(np.zeros(sLength), index=trees.index)
    trees['D2'] = pd.Series(np.zeros(sLength), index=trees.index)
    trees['allowed'] = pd.Series(np.zeros(sLength), index=trees.index)
    trees['ntree'] = pd.Series(np.zeros(sLength), index=trees.index).apply(str) #index of the enxt tree
    trees['ptree'] = pd.Series(np.zeros(sLength), index=trees.index).apply(str) #index of the previous tree
    
    #mark trees as allwoed
    for index,row in trees.iterrows():
        cellNr = row['cellNr']
        tree = row['tree']
        trees.set_value(index,'allowed',allowed_tree(data,cellNr,tree))
    #indexing for navigation
    catch = pd.Series(0,index = trees.tree.unique())
    for index,row in trees.iterrows(): 
        cellNr = row['cellNr']
        tree = row['tree']
        if len(data.loc[(data['cellNr']==cellNr*2) & (data['tree']==str(tree))]) != 0:
            trees.set_value(index, 'D1', data.loc[(data['cellNr']==cellNr*2) & (data['tree']==str(tree))].index.values[0]) #what happens if you click on D1
        if len(data.loc[(data['cellNr']==cellNr*2+1) & (data['tree']==str(tree))]) != 0:
            trees.set_value(index, 'D2', data.loc[(data['cellNr']==cellNr*2+1) & (data['tree']==str(tree))].index.values[0]) #what happens if you click on D2
        if trees.loc[index,'allowed']  == True or index == 0: #only do the computation for allowed trees
            counter = index
            #look for the next tree
            while counter < trees.index.values[-1]:
                counter = counter + 1
                if counter in trees.index.values:
                    if not (trees.loc[counter,'cellNr'] == cellNr and trees.loc[counter,'tree'] == tree) and not (trees.loc[counter,'cellNr'] >= cellNr*2 and trees.loc[counter,'tree'] == tree): #try to catch paralell but disconnected trees with the same tree index
                        if catch.loc[str(trees.loc[counter,'tree'])] < 2:
                            if tree == trees.loc[counter,'tree']:  
                                catch.set_value(str(trees.loc[counter,'tree']), catch.loc[str(trees.loc[counter,'tree'])] +1) #catch the top two disconnected subtrees 
                            break
            #keep looking for the next full subtree
            if counter < trees.index.values[-1]:
                while trees.loc[counter,'allowed'] != True:
                    if counter >= trees.index.values[-1]:
                        break
                    counter = counter + 1
                    while counter not in trees.index.values:
                        if counter > trees.index.values[-1]:
                            break
                        counter = counter + 1
            if counter >= trees.index.values[-1]:
                counter = trees.index.values[0]
            #update entry
            trees.set_value(index, 'ntree',counter)
            counter = index
            #look for a change in tree while browsing backwards
            while counter > 0:
                if counter in trees.index.values:
                    if trees.loc[counter,'tree'] != tree:
                        break
                counter = counter - 1
            #keep looking until you find a complete subtree
            if counter > 0:
                while trees.loc[counter,'allowed'] != True:
                    if counter <= 1:
                        counter = trees.index.values[-1]
                    counter = counter - 1
                    while counter not in trees.index.values:
                        if counter <= 1:
                            break
                        counter = counter - 1
                    if counter <= 1: #
                        counter = trees.index.values[-1]
            #now go back until the tree changes again
            temp = trees.loc[counter,'tree']
            while counter > 0:
                if counter in trees.index.values:
                    if trees.loc[counter,'tree'] != temp:
                        break
                counter = counter - 1
            #go now forward for the next subtree
            if counter > 0:
                while trees.loc[counter,'allowed'] != True:
                    if counter <= 1:
                        break
                    counter = counter + 1
                    while counter not in trees.index.values:
                        if counter <= 1:
                            break
                        counter = counter + 1    
                    if counter <= 1: #
                        break #
            if counter <= 0:
                counter = trees.index.values[0]
            #update entry
            trees.set_value(index, 'ptree',counter)
    
    #udpate D1 and D2 for jumping incomplete trees
    for index,row in trees.iterrows():
        if trees.loc[index,'allowed'] == True:
            if int(trees.loc[index,'D1']) == index:
                if trees.loc[int(trees.loc[int(trees.loc[index,'D1']), 'D1']), 'allowed'] == True:
                    trees.set_value(index, 'D1', trees.loc[int(trees.loc[index,'D1']), 'D1'])
                if trees.loc[int(trees.loc[int(trees.loc[index,'D1']), 'D2']), 'allowed'] == True:
                    trees.set_value(index, 'D1', trees.loc[int(trees.loc[index,'D1']), 'D2'])
            if trees.loc[int(trees.loc[index,'D1']), 'allowed'] != True:
                trees.set_value(index, 'D1', index)
                
            if int(trees.loc[index,'D2']) == index:
                if trees.loc[int(trees.loc[int(trees.loc[index,'D2']), 'D1']), 'allowed'] == True:
                    trees.set_value(index, 'D2', trees.loc[int(trees.loc[index,'D2']), 'D1'])
                if trees.loc[int(trees.loc[int(trees.loc[index,'D2']), 'D2']), 'allowed'] == True:
                    trees.set_value(index, 'D2', trees.loc[int(trees.loc[index,'D2']), 'D2'])
            if trees.loc[int(trees.loc[index,'D2']), 'allowed'] != True:
                trees.set_value(index, 'D2', index) 
    #TO DO: delete all rows that are not part of the tree structure  
    trees.to_csv(os.path.join(filename))
    return rep_cellNr, rep_tree

def allowed_tree(cells, cellNr, tree): #VERY IMPORTANT check for validity of tree
    if len(cells.loc[(cells['cellNr']==cellNr*2) & (cells['tree']==str(tree))]) != 0:
        if len(cells.loc[(cells['cellNr']==cellNr*2+1) & (cells['tree']==str(tree))]) != 0:
            return True
        if len(cells.loc[(cells['cellNr']==cellNr*2+1) & (cells['tree']==str(tree))]) == 0:
            return False
    if len(cells.loc[(cells['cellNr']==cellNr*2) & (cells['tree']==str(tree))]) == 0:
        return False
