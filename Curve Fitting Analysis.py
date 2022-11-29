# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 23:11:36 2022

@author: Agus
"""

import scipy
#matplotlib.use('Qt5Agg')
import os
import numpy as np
import glob
import sklearn.decomposition
from scipy.optimize import curve_fit

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
Curve Fitting Analysis

'''

# unimodal vonmises function M(θ) = A exp{k[cos 2(θ − φ) − 1]}
def vonmises_unimodal(x, k, phi, A, B): #stds
    return B + A*np.exp(k*(np.cos(x - phi) - 1)) 

rsq_group = []
rng = np.random.RandomState(0)
train_size = 0.75

'''
Load Files
'''

filelist = glob.glob(r"filepath")

pars_alls = []
rsq_alls = []
for filename in filelist:
    
    positions = np.load(os.path.join(filename, 'positions.npy'))
    ftraces_good = np.load(os.path.join(filename, 'fluorescence.npy'))
    
    base_all = deltaf_f(ftraces_good)
    
    [l.tolist() for l in base_all] 
    base_all = np.array(base_all) 
    
    positions_degrees = []    
    for value in positions:
        x = (value*360)/100
        positions_degrees.append(x)
    
    positions_radians = np.radians(positions_degrees) 
    
    positions_radians_2 = []
    for value in positions_radians:
        x = value - np.pi
        positions_radians_2.append(x)
    
    rsq_all = []
    rsq1_all = []
    rsq_all_bin= []
    rsq1_all_bin = []
    pars_data = []
    pars_data_bin = []
    bins_per_cell_all = []
    circ_all = []
    
    for index in range(len(base_all)):
         
        x = positions_radians_2
        y = np.ndarray.tolist(base_all[index])
        
        init_vals = [1, 0, np.max(y), np.percentile(y, 5)] # for [k, phi, A, B]
        
        #split data into train and test sets    
        x_train, x_test, y_train, y_test = sklearn.model_selection.train_test_split(x, y, train_size=train_size, random_state=rng)
        #calculate optimal parameters using the train data
        pars1, cov1 = curve_fit(vonmises_unimodal, x_train, y_train, p0=init_vals, maxfev=100000, bounds=([0,-3, 0, np.min(y)], [np.inf, 3, np.inf, np.max(y)]))
        #print('parameters: {}'.format(pars1))
        
        #calculate R2 of curve on train data    
        residuals = y_train - vonmises_unimodal(x_train, *pars1)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y_train - np.mean(y_train))**2)
        rsq = 1 - (ss_res / ss_tot)
        
        #calculate R2 on test data using vonmises with optimal parameters calculated on the train data    
        residuals = y_test - vonmises_unimodal(x_test, *pars1)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y_test - np.mean(y_train))**2)
        rsq1 = 1 - (ss_res / ss_tot)
        
        rsq_all.append(rsq)
        rsq1_all.append(rsq1)
        pars_data.append(pars1)
        
        #bin and average data
        y_binned_mean = scipy.stats.binned_statistic(x, y, statistic= 'mean', bins=20)
        y_binned_var = scipy.stats.binned_statistic(x, y, statistic=lambda l: np.var(l), bins=20)
        
        # get list of values of x and y for each bin
        
        bin_lists_x = []
        bin_lists_y = []
        for ind in range(len(y_binned_mean[1])-1):
            lx = []
            ly = []
            for indexw in range(len(x)):
                if(x[indexw] >= y_binned_mean[1][ind] and x[indexw] < y_binned_mean[1][ind+1]):
                    lx.append(x[indexw])
                    ly.append(y[indexw])
            bin_lists_x.append(lx)
            bin_lists_y.append(ly)   
        
        #separate train and test in each bin
        x_train_bi = []
        x_test_bi = []
        y_train_bi = []
        y_test_bi = []
        for indexx in range(len(bin_lists_y)):
            
            x_train, x_test, y_train, y_test = sklearn.model_selection.train_test_split(bin_lists_x[indexx], bin_lists_y[indexx], train_size=train_size, random_state=rng)
            x_train_bi.append(x_train)
            y_train_bi.append(y_train)
            x_test_bi.append(x_test)
            y_test_bi.append(y_test)
        
        #determine mean/std per bin per group    
        groups = [x_train_bi, y_train_bi, x_test_bi, y_test_bi]
        groups_bin = []
        for g in groups:
            bins_group = []
            for l in g:
                xl = np.mean(l) #/np.std(l)
                bins_group.append(xl)
            groups_bin.append(bins_group)
            
        x_train_bin = groups_bin[0]
        y_train_bin = groups_bin[1]
        x_test_bin = groups_bin[2]
        y_test_bin = groups_bin[3]
        
        #determine the mean/std for x and y data (for plotting see below)
        groupss = [bin_lists_x, bin_lists_y]
        groups_binn = []
        for gr in groupss:
            bins_groups = []
            for li in gr:
                xls = np.mean(li)#/np.std(li)
                bins_groups.append(xls)
            groups_binn.append(bins_groups)
            
        x_bin = groups_binn[0]
        y_bin = groups_binn[1]
        
        #calculate optimal parameters using the train bin data
        pars2, cov2 = curve_fit(vonmises_unimodal, x_train_bin, y_train_bin, p0=init_vals, maxfev=1000000, bounds=([0,-3, 0, np.min(y)], [np.inf, 3, np.inf, np.max(y)]))
        #print('parameters: {}'.format(pars2))
          
        #cross validate
        
        #calculate R2 of curve on train data    
        residuals = np.array(y_train_bin - vonmises_unimodal(x_train_bin, *pars2))
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y_train_bin - np.mean(y_train_bin))**2)
        rsq_bin = 1 - (ss_res / ss_tot)
        #calculate R2 on test data using vonmises with optimal parameters calculated on the test data    
        residuals = np.array(y_test_bin - vonmises_unimodal(x_test_bin, *pars2))
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y_test_bin - np.mean(y_test_bin))**2)
        rsq1_bin = 1 - (ss_res / ss_tot) 
        

        #circular variance
        besseli1 = scipy.special.i1(pars2[0])
        besseli0 = scipy.special.i0(pars2[0])
        
        circ_var = 1 - (besseli1/besseli0)
        circ_var_deg = 360*np.sqrt(circ_var)
        pars3 = np.append(pars2, circ_var_deg)
        
        #x_bin in degrees           
        x_bin_degrees = []
        for x in x_bin:
            x_bin_degrees.append(np.degrees(x))
            
        x_test_bin_deg = []
        for xu in x_test_bin:
            x_test_bin_deg.append(np.degrees(xu))
        
        rsq_all_bin.append(rsq_bin) #train r2 per lap
        rsq1_all_bin.append(rsq1_bin) #test r2 per lap 
        pars_data_bin.append(pars3) 
        #bins_per_cell_all.append(bins_per_cell)
        
        pars2_to_print = ['{:f}'.format(item) for item in pars3]
        
#            plt.figure()
#            plt.plot(x_test_bin_deg, y_test_bin)
#            plt.plot(x_bin_degrees, vonmises_unimodal(x_bin, *pars2), 'r')
#            plt.text(-3,np.min(y_bin), 'R2test = {0:.3f}  R2train = {1:.3f}\n[k, phi, A, B, std] = {2}'.format(rsq1_bin, rsq_bin, pars2_to_print))
#            #plt.savefig(os.path.join(filename, 'curve3_fit_cell_test{:03d}.png'.format(index)), dpi=150)
        #plt.close('all')
        
        
        pars2_to_print = ['{:f}'.format(item) for item in pars3]
        
#        plt.figure()
#        plt.plot(x_test_bin_deg, y_test_bin)
#        plt.scatter(x_bin_degrees, vonmises_unimodal(x_bin, *pars2))
#        plt.text(-3,np.min(y_bin), 'R2test = {0:.3f}  R2train = {1:.3f}\n[k, phi, A, B, std] = {2}'.format(rsq1_bin, rsq_bin, pars2_to_print))
#        plt.savefig(os.path.join(filename, 'curve3_fit_cell_test{:03d}.png'.format(index)), dpi=150)
#        plt.close('all')
        
    pars_comp = list(zip(pars_data,pars_data_bin))   
    
    pars_data_bin_list = []
    for v in pars_data_bin:
        n = np.ndarray.tolist(v)  
        pars_data_bin_list.append(n)
        
    #transform data parameters (not binned) to list to save    
    pars_data_list = []
    for v in pars_data:
        n = np.ndarray.tolist(v)  
        pars_data_list.append(n)    
        
    np.savetxt(os.path.join(filename, 'pars_data_list2'), pars_data_list)
    np.savetxt(os.path.join(filename, 'rsq_train2'), rsq_all)
    np.savetxt(os.path.join(filename, 'rsq_test2'), rsq1_all)
    
    np.savetxt(os.path.join(filename, 'pars_data_bin2'), pars_data_bin_list) 
    np.savetxt(os.path.join(filename, 'rsq_train_bin2'), rsq_all_bin)
    np.savetxt(os.path.join(filename, 'rsq_test_bin2'), rsq1_all_bin)
    
    pars_alls.append(pars_data_bin_list)
    rsq_alls.append(rsq1_all_bin)
    
    rsq_group.append(rsq1_all_bin)
    print(filename)








