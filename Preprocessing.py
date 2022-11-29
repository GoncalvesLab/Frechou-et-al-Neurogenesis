# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 10:47:26 2022

@author: M. Agustina Frechou
"""

'''
Imports
'''

#matplotlib.use('Qt5Agg')
import os
import numpy as np
import glob
import h5py 
from skimage.transform import resize

'''
Functions
'''

def extract_zones(textures_trace, t_voltages = [0.6, 1.2, 1.8, 2.4]):
        
        delay = 15000
        transitions = [0]
        textures = []
        mean_values = []
        thre_values = []
        t_voltages = np.array(t_voltages)
        
        while True:
        
            last_transition = transitions[-1]
            
            ref = textures_trace[last_transition+delay:last_transition+delay+1000, 0]
            if ref.size == 0:
                #print("Last transition at the end, assuming next texture")
                textures.append((textures[-1] + 1) % len(t_voltages))
                #print(f"There are {len(transitions)} zones.")
                break
            current_mean = np.mean(ref)
            textures.append(np.argmin(np.abs(t_voltages - current_mean)))
            mean_values.append(current_mean)
            current_threshold = 0.8 * (np.max(ref) - np.min(ref))
            thre_values.append(current_threshold)
            outside_this_texture = np.where(np.abs(textures_trace[last_transition+delay:, 0] - current_mean) > current_threshold)[0]
            if outside_this_texture.size == 0:
                #print(f"There are {len(transitions)} zones.")
                break
            transitions.append(outside_this_texture[0] + last_transition+delay)
        # Return transitions and textures as a tuple of two lists:
        return np.array(transitions), np.array(textures)

#%%

'''
Load Files
'''

filelist = glob.glob(r'filepath') 
brain_region = "DG3*" #add brain region name with * to determine the tif lengths for every session recorded of that brain region

activity_all = []

for filename in filelist:

    '''
    LOAD FILES
    '''
    mouse_dir = filename
    cells = np.load(os.path.join(mouse_dir, 'suite2p\plane0\iscell.npy'))
    ftraces = np.load(os.path.join(mouse_dir, 'suite2p\plane0\F.npy'))
    deconvtraces = np.load(os.path.join(mouse_dir, 'suite2p\plane0\spks.npy'))
    
    # Frames per Second:
    fps = 15.253
    
    #Total frames of video
    tf = len(deconvtraces[1])
    
    '''
    SELECT USEFUL ROIs
    '''
    #select good ROI from iscell and get only good traces
    
    ftraces_good = ftraces[cells[:,0]==1]
    dtraces_good = deconvtraces[cells[:,0]==1]
       
#%%
        
    '''
    spatial information
    '''
    
    #using the length of the video folder
    a = os.path.split(mouse_dir)[0]
    #aprime = os.path.split(a)[0] #when you have 2 folders to go back
    b = glob.glob(os.path.join(a, brain_region)) #if above uncommented add aprime here
    lenexp = []
    for file in b:
        c = glob.glob(os.path.join(file, "*.tif"))
        lens = len(c)
        lenexp.append(lens)
        
    print(lenexp)
    
    spatial_dir = glob.glob(os.path.join(a, brain_region, 'SyncData*'))   #aprime and add the second folder here too in the syncdata
    
    '''
    clean files
    '''
    normedlist_all = []
    trans1_all = []
    indexes_final = []
    count = 0
    
    mset_all = []
    tset_all = []
    for index in range(len(spatial_dir)):
        
        f = h5py.File(os.path.join(spatial_dir[index], 'Episode001.h5'), 'r+')
        list(f.keys())
        a_group_key = list(f.keys())[0]
    
        mset = np.array(f['AI/movement']) #converts the h5 file to numpy array
        tset = np.array(f['AI/texture'])
    
        t = np.arange(mset.size) / 30000 # 
    
        #plt.figure()
        #plt.plot(t, mset)
        #plt.plot(t, tset)
        #plt.savefig(os.path.join(filename, 'mset_tset_file{}.png'.format(index + 1)), dpi=150)
        #plt.close('all')
    
        t2 = np.arange(len(mset)) / fps
        keep_spatial = (t <= t2[-1]) #-1 starts from the last value of t2. This variable keeps the values of t2 that are equal or smaller to t (so we cut that extra part that is saved with the data that we dont need)
        mset = mset[keep_spatial] #cuts mset and leaves only the part that we need
        tset = tset[keep_spatial] #cuts tset and leaves only the part that we need 
        t = t[keep_spatial] #cleans t
    
        # Unpack tuple into two lists: get transitions and zones
        transitions, zones = extract_zones(tset)
        
        mset_all.extend(mset)
        tset_all.extend(tset)
        
        cum = []
        # resample mset 
        res1 = resize(mset, (lenexp[index],1), order=1, preserve_range=True)
        trans = (transitions/len(mset))*lenexp[index]
        #convert floats in trans to integers so you can use them to slice
        trans1 = [int(i) for i in trans]
        #do the cumsum for every zone
        for k in range(1, len(trans1)-2): # this selects all values from the transitions except the first and last
            cum.append(np.cumsum(res1[trans1[k]:trans1[k+1]]))
        
        [l.tolist() for l in cum]
        
        #normalize
        normed = []
        for i in range(len(cum)):
            if zones[i] == 3:
                normed.append((cum[i]/np.max(cum[i]))*25 + 75)
            if zones[i] == 2:
                normed.append((cum[i]/np.max(cum[i]))*25 + 50)
            if zones[i] == 0:
                normed.append((cum[i]/np.max(cum[i]))*25 + 25)
            if zones[i] == 1:
                normed.append((cum[i]/np.max(cum[i]))*25)
    
        normed = [l.tolist() for l in normed]
        normed = np.array(normed)
    
        flatten = lambda l: [item for sublist in l for item in sublist]
        normedlist = flatten(normed)
        normedlist = np.array(normedlist)
    
        normedlist_all.extend(normedlist)
        trans1_all.append(trans1)
    
        lens = list(range(lenexp[index] + 1))
        vid = lens[trans1[1]:trans1[-2]]
        video = [x + count for x in vid]
        count += lenexp[index] + 1
        indexes_final.extend(video)
    
    
    
    final_f = ftraces_good[:,indexes_final]
    
    final_spks = dtraces_good[:,indexes_final]
    
    np.save(os.path.join(mouse_dir, 'positions.npy'), normedlist_all)
    
    np.save(os.path.join(mouse_dir, 'spks_final.npy'), final_spks)
    
    np.save(os.path.join(mouse_dir, 'fluorescence.npy'), final_f)
    
    np.save(os.path.join(mouse_dir, 'mset'), mset_all)
    np.save(os.path.join(mouse_dir, 'tset'), tset_all)
    
    print('Processed file:' + filename)

    # check everything makes sense
    #plt.plot(ftraces_good[0])
    #plt.plot(mset_all)
    
    '''
    preprocessing for tuning code
    '''
    
    # Subtract mean:
    dataf_mean = np.mean(ftraces_good, axis=1)[:, None]
    dataf2 = ftraces_good - dataf_mean
    # Threshold of "dataf2"
    threshold2 = (np.std(ftraces_good, axis=1) * 2)[:, None]
    threshold = threshold2 + dataf_mean
    # Fluorescence signal above mean+2*sigma:
    above_threshold = (dataf2 > threshold2) & (dtraces_good > 0) #binary data == (data2 > 0). True are cells with activity above the threshold
    
    data2 = dtraces_good.copy()
    # Use calculated filter on "Suite2P" cleaned signal:
    # set all values to 0 where NOT above_theshold:
    data2[~above_threshold] = 0
    
    # Plot each cleaned&re-thresholded trace in one figure:
    
    dataf_mean = np.mean(ftraces_good, axis=1)[:, None]
    threshold2 = (np.std(ftraces_good, axis=1) * 2)[:, None]
    threshold = threshold2 + dataf_mean
    
    np.save(os.path.join(mouse_dir, 'cleaned_traces_binary.npy'), above_threshold) #number of active cells (yes/no type of data)
    np.save(os.path.join(mouse_dir, 'cleaned_traces_data.npy'), data2) #actual values of active cleaned cells 
    
    
    
    
    
    
    
    
    
    
    
    
    
