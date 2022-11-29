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
from scipy import signal


'''
FUNCTIONS
'''

def remove_noise(traces):
    #apply filter to remove noise
    ftraces_filtered = []
    for index in range(len(traces)):
        b, a = signal.butter(3, 0.05) #Create an order 3 lowpass butterworth filter and return the filter coefficients.
        y = signal.filtfilt(b, a, traces[index]) #This function applies a linear digital filter twice, once forward and once backwards. The combined filter has zero phase and a filter order twice that of the original.
        ftraces_filtered.append(y)
        
    [l.tolist() for l in ftraces_filtered]
    ftraces_filtered = np.array(ftraces_filtered)
    
    return ftraces_filtered

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

def a_filter(trace, N=100):
        np.ndarray.tolist(trace)
        through = False
        for index in range(len(trace)):
            if index > N:
                if trace[index] > np.std(trace[index-1-N:index-1]) * 2:
                    through = True
                elif trace[index] < np.std(trace[index-1-N:index-1]) * 0.5:
                    through = False
                yield trace[index] if through else 0
                if trace[index] <= np.std(trace[index-1-N:index-1]) * 0.5:
                    through = False
            else:
                yield 0
                
#function using the traces without peaks to calculate the rolling mean
def a_filter2(trace, peaks, N=100):
    np.ndarray.tolist(trace)
    through = False
    for index in range(len(trace)):
        if index > N:
            if trace[index] > np.std(peaks[index-1-N:index-1]) * 2:
                through = True
            elif trace[index] < np.std(peaks[index-1-N:index-1]) * 0.5:
                through = False
            yield trace[index] if through else 0
            if trace[index] <= np.std(peaks[index-1-N:index-1]) * 0.5:
                through = False
        else:
            yield 0

#%%

'''
Area under the curve normalized
'''

'''
Load Files
'''

filelist = glob.glob(r"filepath")


activity_all_mice = []
for filename in filelist:
    
    # Load files
    spatial_dir = np.load(os.path.join(filename, 'positions.npy'))
    ftraces_good = np.load(os.path.join(filename, 'fluorescence.npy'))
    
    # Bin positions
    spatial_binned = bin_data(spatial_dir)

    # Determine DeltaF/F of calcium data
    delta = deltaf_f(ftraces_good)
    
    '''
    removing noise
    '''   

    ftraces_filtered = remove_noise(delta) 
    
    activity_all = []
    for i in range(len(ftraces_filtered)):
        
        #no iteration
        ca_traces = list(a_filter(ftraces_filtered[i]))
        no_peaks = ftraces_filtered[i] - ca_traces
        activity1 = np.sum(ca_traces)
        
        #iteration 1
        ca_traces2 = list(a_filter2(ftraces_filtered[i], no_peaks))
        no_peaks2 = ftraces_filtered[i] - ca_traces2
        activity2 = np.sum(ca_traces2)
    
        #iteration 2
        ca_traces3 = list(a_filter2(ftraces_filtered[i], no_peaks2))
        no_peaks3 = ftraces_filtered[i] - ca_traces3
        activity3 = np.sum(ca_traces3)
        
        #iteration 3
        ca_traces4 = list(a_filter2(ftraces_filtered[i], no_peaks3))
        activity4 = np.sum(ca_traces4)
        
        #activity normalized to distance traveled in cm:
        #1) Extract all the end of laps
        l = []
        for index in range(len(spatial_dir)-1):
            if spatial_dir[index] > spatial_dir[index+1]:
                l.append(spatial_dir[index])
            else:
                l.append(0)
        l.append(0)
                
        #2) Sum end of laps to get final position
        p = []
        adds = 0
        for index in range(len(spatial_dir)-1):
            if spatial_dir[index] > spatial_dir[index+1]:
                p.append(spatial_dir[index] + adds)
                adds += l[index]
            else:
                p.append(spatial_dir[index] + adds)
                adds += l[index]
        p.append(spatial_dir[index] +1 + adds)
    
        #get the last value of the cummulative position equal to the distance traveled
        distance_traveled = p[-1]
        #convert distance_traveled to cm (1 lap (100) = 180cm)
        dist_trav_cm = distance_traveled*180/100
        
        #normalize activity to distance traveled in cm      
        activity_all.append(activity4/dist_trav_cm)
    print('dist traveled total' + str(dist_trav_cm))    
    np.savetxt(os.path.join(filename, 'activity_area_under_curve_normed'), activity_all)
    activity_all_mice.append(activity_all)
    print('file processed: ' + filename)
    
            
            
