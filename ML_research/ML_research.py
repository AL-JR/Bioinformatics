# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 02:13:40 2018

@author: aalco
"""

from sklearn import svm
from sklearn.model_selection import KFold
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from argparse import ArgumentParser



#implement argparse



def otucolarray(idname,otu_table):
    otucol = otu_table.loc[:,idname].tolist()
    array = np.array(otucol)
    return array
   
def XYarraymaker(otu_table,type_list, totalsamplecount,featurecount):
    orderedlist = []
    type_array = np.zeros(totalsamplecount)
    type_index = 1
    index = 0
    for typelist in type_list:
        for idname in typelist:
            array = np.array(otu_table.loc[:,idname])
            orderedlist.append(array)
            type_array[index] = type_index
            index += 1
        type_index += 1
    orderedlist = np.asarray(orderedlist)
    return orderedlist , type_array
#creating type_list and typecount_list and featurecount and samplecount
def dataparser(otu_table,map_table):
    type_table = map_table.loc[: , ['#SampleID', 'Group']]
    ID_list = type_table['#SampleID'].tolist()
    group_list = type_table['Group'].tolist()
    type_dict = dict(zip(ID_list,group_list))
    dicot_list,monocot_list,soil_list = [],[],[]
    for sample in type_table.loc[:,'#SampleID']:
        if type_dict[sample] == 'Dicot':
            dicot_list.append(sample)
        elif type_dict[sample] == 'Monocot':
            monocot_list.append(sample)
        else:
            soil_list.append(sample)
    type_list = [dicot_list,monocot_list,soil_list]
    typecount_list = [len(dicot_list),len(monocot_list),len(soil_list)]
    totalsamplecount = len(dicot_list) + len(monocot_list) + len(soil_list)
    
    soilsample = otu_table.loc[:, type_list[0][0]].tolist()
    featurecount =len(soilsample)
    
    return type_list, typecount_list,totalsamplecount, featurecount , type_dict
def accuracy(testY,prediction_list):
    numcorrect = 0
    scount = 0
    while scount < len(testY):
        if testY[scount] == prediction_list[scount]:
            numcorrect+=1
        scount +=1
    accuracy = numcorrect/(scount +1)
    return accuracy
    

#Main Code
#extracting data
otu_table = pd.read_csv('all_bact_otu_table_e.txt', sep = '\t')
map_table = pd.read_csv('allBact_mappings.txt', sep = '\t')
type_list, typecount_list,totalsamplecount, featurecount , type_dict = dataparser(otu_table,map_table)
X,Y = XYarraymaker(otu_table, type_list, totalsamplecount,featurecount)
trainX_list , trainY_list = [],[]
testX_list, testY_list = [],[]
k = 20
cumulative_accuracy = True
clf = svm.SVC(gamma= 'auto')
kaccuracy_list = []
kfolds = 2
if cumulative_accuracy == True:
    while kfolds <= k:
    #kfold partitioning of data
        kf = KFold(n_splits = kfolds,shuffle = True)
        kf.split(X)
        data_split = kf.split(X)
        for train_index,test_index in kf.split(X):
            X_train = X[train_index]
            X_test = X[test_index]
            Y_train,Y_test = Y[train_index],Y[test_index]
            trainX_list.append(X_train)
            trainY_list.append(Y_train)
            testX_list.append( X_test)
            testY_list.append(Y_test)
    
    #SVM prediction and error calculation
    
        k_count = 0
        accuracylist = []
        predictionarrays = []
        while k_count < kfolds:
            samples,features = trainX_list[k_count],trainY_list[k_count]
            clf.fit(samples,features)
            predictionarray = clf.predict(testX_list[k_count])
            acc = accuracy(testY_list[k_count], predictionarray)
            accuracylist.append(acc)
            predictionarrays.append(predictionarray)
            k_count += 1
        avg_err = sum(accuracylist)/len(accuracylist)
        kaccuracy_list.append(avg_err)
        kfolds += 1
#plotting error rates for each fold
    folds =range(2,k+1)
    y_pos = np.arange(len(kaccuracy_list))
    plt.bar(y_pos, kaccuracy_list, align='center', alpha=0.5)
    plt.xticks(y_pos, folds)
    plt.ylabel('Accuracy %')
    plt.title('Accuracy for folds 2 to k')
     
    plt.show()
            
else:
#kfold partitioning of data
    kf = KFold(n_splits = kfolds,shuffle = True)
    kf.split(X)
    data_split = kf.split(X)
    for train_index,test_index in kf.split(X):
        X_train = X[train_index]
        X_test = X[test_index]
        Y_train,Y_test = Y[train_index],Y[test_index]
        trainX_list.append(X_train)
        trainY_list.append(Y_train)
        testX_list.append( X_test)
        testY_list.append(Y_test)

#SVM prediction and error calculation

    k_count = 0
    accuracylist = []
    predictionarrays = []
    while k_count < kfolds:
        samples,features = trainX_list[k_count],trainY_list[k_count]
        clf.fit(samples,features)
        predictionarray = clf.predict(testX_list[k_count])
        acc = accuracy(testY_list[k_count], predictionarray)
        accuracylist.append(acc)
        predictionarrays.append(predictionarray)
        k_count += 1
    avg_acc = sum(accuracylist)/len(accuracylist)
    kaccuracy_list.append(avg_err)
    kfolds += 1
    
    print('Average Error: ', avg_err)
    print('Number of folds: ', k)


