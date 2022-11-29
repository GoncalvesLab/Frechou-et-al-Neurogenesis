# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 23:11:36 2022

@author: Agus
"""

'''
Imports
'''

#matplotlib.use('Qt5Agg')
import os
import numpy as np
import glob
import sklearn.decomposition 
from sklearn import metrics
from sklearn.linear_model import LogisticRegression

'''
FUNCTIONS
'''
def neurons_per_file(filelist, filename):
    
    #get number of neurons in each file
    neurons = []
    for file in filelist:
        traces = np.load(os.path.join(file, filename))
        neurons.append([file, traces.shape[0]])
    neurons = np.array(neurons)
    return neurons

def subsample(filelist, filename, size):
    #subsample traces to the appropiate size
    ftraces_good_matched_all = []
    for file1 in filelist:
        ft = []
        traces = np.load(os.path.join(file1, filename))
        for index in range(0, 10):
            idx = np.random.randint(traces.shape[0], size=size)
            ftraces_good_matched = traces[idx,:]
            ft.append(ftraces_good_matched)
            index +=1
        ftraces_good_matched_all.append(ft)
    return ftraces_good_matched_all

def deltaf_f(traces):
    #calculate delta F/F for every trace
    #traces is array containing traces for all neurons in one mouse
    delta_all = []
    for index in range(len(traces)):
        ftrace = traces[index]
        delta_rolling = []
        for index in range(0, len(ftrace), 100):
            x = ftrace[index:index+100]
            percentile = np.percentile(x, 20)
            delta = (ftrace[index:index+100] - percentile)/percentile
            delta_rolling.extend(delta)
        delta_all.append(delta_rolling)
    delta_all = np.array(delta_all)
    return delta_all

def bin_data(spatial_dir, nbins=21):    
    #bin positions
    bins = np.linspace(0, 100, nbins, dtype='int')
    spatial_binned = np.zeros(len(spatial_dir))
    count = 0
    for index in range(len(bins)-1):
        spatial_binned[np.logical_and(spatial_dir>bins[index], spatial_dir<=bins[index+1])] = count
        count += 1
    return spatial_binned

#%%

'''
Linear Decoder: Population Spatial Information Content

'''

'''
Load Files
'''

'''
preprocessing
'''

#load groups to subsample
filelist_group1 = glob.glob(r"filepath")
filelist_group2 = glob.glob(r"filepath")


#get number of neurons per mouse    

neurons_group1 = neurons_per_file(filelist_group1, 'fluorescence.npy')[:,1]
neurons_group2 = neurons_per_file(filelist_group2, 'fluorescence.npy')[:,1]

#remove the mice that dont have enough cells   

# Determine Information Content for first group
w = filelist_group1

#subsample
ftraces_good_matched_all = subsample(w, 'fluorescence.npy', 48) #size = amount of cells you're sampling from your data

positions_all = []
for filename in w:    
    
    #make list with all positions in filename
    positions = np.load(os.path.join(filename, 'positions.npy'))
    positions_all.append(positions)


filelist = [ftraces_good_matched_all, positions_all]

#run for all pseudo animals that have each 10 animals 

average_all = []
for iiii in range(len(filelist[0])):
    
    # bin positions 
    spatial_binned = bin_data(filelist[1][iiii])
        
    accuracy_all_finn = []
    accuracytrain_all_finn = []
        
    for index in range(len(filelist[0][iiii])):
        
        ftraces_good = filelist[0][iiii][index]
        
        base_all = deltaf_f(ftraces_good)        
        
        a = list(zip(spatial_binned, base_all.T))
        
        x=[]
        for i_a in range(len(a)):
            for i_x in range(len(x)):
                if a[i_a][0]==x[i_x][0][0]:
                    x[i_x].append(a[i_a])
                    break
            else:
                x.append([a[i_a]])
        
        '''
        logistic regression
        '''

        LR_scores = []
        y_pred = []
        acc_score = []

        intercept_all = []
        accuracy_all = []
        accuracytrain_all = []
        dprime_all = []
        
        for i in range(len(x)):
            
            if i < 19:
            
                X1 = []
                y1 = []
                for ind in range(len(x[i])):  
                    X1.append(x[i][ind][1])   
                    y1.append(x[i][ind][0])   
                        
                X2 = []
                y2 = []
                for ii in range(len(x[i+1])):  
                    X2.append(x[i+1][ii][1])   
                    y2.append(x[i+1][ii][0])
            
            elif i == 19:
                X1 = []
                y1 = []
                for indm in range(len(x[i])):  
                    X1.append(x[i][indm][1])   
                    y1.append(x[i][indm][0])   
                        
                X2 = []
                y2 = []
                for iim in range(len(x[0])):  
                    X2.append(x[0][iim][1])   
                    y2.append(x[0][iim][0])
            
            
            
            rng = np.random.RandomState(0)
            train_size = 0.75
            
            X1_train, X1_test, y1_train, y1_test = sklearn.model_selection.train_test_split(X1, y1, train_size=train_size, random_state=rng)
            X2_train, X2_test, y2_train, y2_test = sklearn.model_selection.train_test_split(X2, y2, train_size=train_size, random_state=rng)
            
            X_train = X1_train + X2_train
            y_train = y1_train + y2_train
            X_test = X1_test + X2_test
            y_test = y1_test + y2_test
            
            #define model
            log_reg= LogisticRegression()
            
            
            #fit to data
            log_reg.fit(X_train, y_train)
            
            #score = log_reg.score(X_test, y_test)
            y_pred=log_reg.predict(X_test)
            #acc_score = metrics.accuracy_score(y_test, y_pred)
            #LR_scores.append(acc_score)
            weights = log_reg.coef_
            
            intercept = log_reg.intercept_
            
            X1aa = np.array(X1_test).T
            X2aa = np.array(X2_test).T
            
            n_1_t = np.dot(weights,X1aa) + intercept
            n_2_t = np.dot(weights,X2aa) + intercept
            var = np.mean([np.var(n_1_t),np.var(n_2_t)])
            
            dprime = (np.square(np.mean(n_1_t)-np.mean(n_2_t)))/var
            dprime = dprime/len(weights[0])
            dprime_all.append(dprime)
        
            train_accuracy = log_reg.score(X_train, y_train)
            #print(f"Accuracy on the training data: {train_accuracy:.2%}")    
            
            test_accuracy = metrics.accuracy_score(y_test, y_pred)
            #print("Accuracy on the testing data:", test_accuracy)
            
            accuracy = [train_accuracy, test_accuracy]
            accuracy_all.append(test_accuracy)
            accuracytrain_all.append(train_accuracy)

    print(np.mean(accuracy_all))  #print test
    print(np.mean(accuracytrain_all)) #print train
    print('file processed: ' + w[iiii])
    decoder_mouse = np.mean(dprime_all)
    np.savetxt(os.path.join(w[iiii], 'mouse_decoder_subsample'), [np.mean(accuracy_all)])
    np.savetxt(os.path.join(w[iiii], 'mouse_decoder_dprime'), [decoder_mouse])