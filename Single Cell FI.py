# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 23:11:35 2022

@author: Agus
"""

'''
Imports
'''
#matplotlib.use('Qt5Agg')
import os
import numpy as np
import glob


'''
FUNCTIONS
'''

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
Single-Cell Spatial Information Content using Fisher Information
'''

'''
Load Files
'''

filelist = glob.glob(r"filepath")

fi_all_mice = []
for filename in filelist:
    
    # Load Files 
    spatial_dir = np.load(os.path.join(filename, 'positions.npy'))
    ftraces_good = np.load(os.path.join(filename, 'fluorescence.npy'))
    
    '''
    Data Preprocessing
    '''
    
    # bin positions
    spatial_binned = bin_data(spatial_dir)
    spatial_binned = [l.tolist() for l in spatial_binned]
    
    # Determine DeltaF/F of calcium data
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
    

    fi_all = []
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
        
        '''
        FI per cell
        '''
        
        m1 = np.array(X1).T
        m2 = np.array(X2).T
        
        mean1 = np.mean(m1, axis=1) #mean response per cell in position 1
        mean2 = np.mean(m2, axis=1) #mean response per cell in position 2
        
        m1var = np.var(m1,axis=1)
        m2var = np.var(m2,axis=1)
        
        fi = np.square(mean1-mean2)/((m1var+m2var)/2)
        
        fi_all.append(fi)

    fe = np.array(fi_all).T
    fi_per_cell = np.mean(fe, axis=1)
    
    np.savetxt(os.path.join(filename, 'FI per cell'), fi_per_cell)
        
    fi_all_mice.append(fi_per_cell)
    
    print('file processed: ' + filename)
    
    
    
    
    
    
    
    