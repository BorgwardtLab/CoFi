# -*- coding: utf-8 -*-
"""
Created on Mon May 09 14:31:47 2016

@author: Sir Thomas
"""
import pandas as pd
import numpy as np
import math


def return_min_trees(cells): #create list of subtrees
    list_of_trees = pd.DataFrame(index=cells.index)
    for counter,rows in cells.iterrows():
        if rows['allowed'] == 1:
            list_of_trees.set_value(counter,'cellNr', rows['cellNr'])
            list_of_trees.set_value(counter,'tree', rows['tree'])
        if rows['allowed'] != 1:
            list_of_trees.drop(counter, inplace=True) #what does inplace mean?
    #print list_of_trees
    return list_of_trees
def save_tree(list_of_trees, counter, user, vote): #save vote of tree
    list_of_trees.set_value(counter, str(user), int(vote))
    return list_of_trees
def get_tree(list_of_trees, user): #move to the next tree whcih the user has not voted on yet
    if str(user) in str(list_of_trees.columns.values):
        for i,r in list_of_trees.iterrows():
            if math.isnan(r[str(user)]) == True:
                counter = i
                break
        return counter
def get_votes(list_of_trees, counter): #collect the results (#of votes c,p,d) of the subtree at location=counter
    votes = list_of_trees.loc[counter, :]
    votes = np.asarray(votes)[2:]
    votes = votes[2:]
    coord = 0
    partial=0
    uncoord=0
    for i in votes:
        try:
            int(i)
            if math.isnan(i) != True:
                if int(i) == 1:
                    coord = coord + 1
                if int(i) == 2:
                    partial = partial + 1
                if int(i) == 3:
                    uncoord = uncoord +1
        except:
            print i
    return coord,partial,uncoord
def get_uncoord(list_of_trees, TFs): #how many trees are AND or OR divergent
    users = get_voters(list_of_trees, TFs)
    all_un =[]
    true_un=[]
    for i, r in list_of_trees.iterrows():
        k = 0
        for x in users:
            if list_of_trees.loc[i,x] == 3:
                k = k + 1
        if k > 0: #at least one user saays that this tree is uncoordinated
            all_un.append(i)
        if k >= len(users): #all users agree, that this is uncoordinated
            true_un.append(i)
    return all_un, true_un
def get_voters(list_of_trees, TFs):
    voters = list(list_of_trees.columns.values[2:])
    k = 0
    while k < len(voters):
        print voters[k]
        for i in TFs:
            print i
            if str(i) in str(voters[k]):
                print voters[k]
                del voters[k]
                print voters
                k = k - 1
        k = k+ 1
    return voters
def get_user_votes(list_of_trees,TF,TFs,subcurve,counter):
    voters = get_voters(list_of_trees,TFs)
    user_votes = []
    for v in voters:
        if str(v)+str(TF)+str(subcurve) in str(list_of_trees.columns.values):
            if math.isnan(list_of_trees.loc[counter,str(v)+str(TF)+str(subcurve)]) != True:
                user_votes.append(list_of_trees.loc[counter,str(v)+str(TF)+str(subcurve)])
                print "Here it is: " + str(list_of_trees.loc[counter,str(v)+str(TF)+str(subcurve)].values)
    return user_votes
def get_single_curve(list_of_trees, counter, TFs,user): #always jump to the next curve which has not been votes on before. If there is none left in the current subtree, jump to the next subtree
    for i in TFs:
        for j in [1,2,3]:
            if str(user)+str(i)+str(j) not in str(list_of_trees.columns.values):
                return counter,i,j
                break
            if math.isnan(list_of_trees.loc[counter,str(user)+str(i)+str(j)]) == True:
                return counter,i,j
                break
    counter = counter + 1
    while counter not in list_of_trees.index.values:
        counter = counter + 1
    return counter, TFs[0], 1
                
