#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 16:11:24 2021

@author: np3
"""

import pandas as pd
import numpy as np
import time


#csv = "/home/np3/Documentos/LNBio/Trem_tes/a.csv"
csv = "/home/np3/Documentos/LNBio/Dev/np3_ms_workflow/src/spectra_viewer/static/a_test.csv"

a = pd.read_csv(csv)

#csv = "/home/np3/Documentos/LNBio/Trem_tes/b.csv"

csv = "/home/np3/Documentos/LNBio/Dev/np3_ms_workflow/src/spectra_viewer/static/b_test.csv"
b = pd.read_csv(csv)


def merge_and_share(dfA, dfB, 
                    mz_col_name="mz",
                    int_col_name="int",
                    bin_size=0.025):
    '''
    Merge mz rows in a bin size window, sliding it on a bin size step
    Also inform all the coincident peaks

    Parameters
    ----------
    dfA : Pandas dataframe
        Peak list in pandas dataframe.
    dfB : Pandas dataframe
        Peak list in pandas dataframe.
    mz_col_name : str, optional
        Column header for m/z peak list. The default is "mz".
    int_col_name : TYPE, optional
        Column header for intensity peak list. The default is "int".
    bin_size : float, optional
        Bin size window for merge. The default is 0.025.

    Returns
    -------
    dfA_merge : Pandas dataframe
        dfA with peaks merged.
    dfB_merge : Pandas dataframe
        dfB with peaks merged.
    dfA_inter : Pandas dataframe
        Shared peaks from dfA.
    dfB_inter : Pandas dataframe
        Shared peaks from dfA..

    '''
    # From pandas to numpy
    a=dfA.to_numpy()
    b=dfB.to_numpy()

    
    # initial value for slider is the min of mzs
    slider = int(min(dfA[mz_col_name].min(), dfB[mz_col_name].min()))
    max_v = max(dfA[mz_col_name].max(), dfB[mz_col_name].max())
    
    
    # Store values for outputs
    final_mz_a = []
    final_mz_b = []
    
    final_int_a = []
    final_int_b = []
    
    inter_mz_a = []
    inter_mz_b = []
    
    inter_int_a = []
    inter_int_b = []
    
    # Agregating    
    while slider < max_v:
        
        aux_a = np.where((a[:,0] >= slider) &
                         (a[:,0] <= slider+bin_size))
        
        if len(aux_a[0]) > 0:
            new_mz_a = np.average(a[aux_a][:,0], weights=a[aux_a][:,1])
            new_int_a = np.sum(a[aux_a][:,1])
            
            final_mz_a.append(new_mz_a)
            final_int_a.append(new_int_a)
    
        aux_b = np.where((b[:,0] >= slider) &
                         (b[:,0] <= slider+bin_size))
        
        if len(aux_b[0]) > 0:
            new_mz_b = np.average(b[aux_b][:,0], weights=b[aux_b][:,1])
            new_int_b = np.sum(b[aux_b][:,1])
            
            final_mz_b.append(new_mz_b)
            final_int_b.append(new_int_b)
            
        # Colect shared peaks
        if len(aux_a[0]) > 0 and len(aux_b[0]) > 0:
            inter_mz_a.append(new_mz_a)
            inter_mz_b.append(new_mz_b)
    
            inter_int_a.append(new_int_a)
            inter_int_b.append(new_int_b)
        
        slider += bin_size
    
    # Create output dataframes
    dfA_merge = pd.DataFrame({mz_col_name: final_mz_a,
                              int_col_name: final_int_a})
    
    dfB_merge = pd.DataFrame({mz_col_name: final_mz_b,
                              int_col_name: final_int_b})
    
    dfA_inter = pd.DataFrame({mz_col_name: inter_mz_a,
                              int_col_name: inter_int_a})
    
    dfB_inter = pd.DataFrame({mz_col_name: inter_mz_b,
                              int_col_name: inter_int_b})
    
    return dfA_merge, dfB_merge, dfA_inter, dfB_inter


def aggregate_mz(nparray, bin_size=0.025, step=0.01):
    
    work_arr = nparray.copy()
    slider = int(nparray[:,0].min())
    max_v = nparray[:,0].max()
    
    while slider < max_v:
        aux = np.where((work_arr[:,0] >= slider) &
                       (work_arr[:,0] <= slider+bin_size))
        
        aux_neg = np.where(~((work_arr[:,0] >= slider) &
                             (work_arr[:,0] <= slider+bin_size)))
        
        if len(aux[0]) > 0:
            new_mz = np.average(work_arr[aux][:,0],
                                weights=work_arr[aux][:,1])
            new_int = np.sum(work_arr[aux][:,1])
            
            work_arr = work_arr[aux_neg].copy()
            
            work_arr = np.append(work_arr, [[new_mz, new_int]], axis=0)
            
        slider += step
    
    return work_arr

def shared_mz(npa, npb, bin_size=0.025):
    
    slider = min(npa[:,0].min(), npb[:,0].min())
    max_v = max(npa[:,0].max(), npb[:,0].max())
    
    out_a = []
    out_b = []
    
    while slider < max_v:
        aux_a = np.where((npa[:,0] >= slider) &
                         (npa[:,0] <= slider+bin_size))
        aux_b = np.where((npb[:,0] >= slider) &
                         (npb[:,0] <= slider+bin_size))
        
        if len(aux_a[0]) > 0 and len(aux_a[0]) > 0:
            if len(out_a) == 0:
                out_a = npa[aux_a].copy()
            else:
                np.append(out_a, npa[aux_a])
            
            if len(out_b) == 0:
                out_b = npb[aux_b].copy()
            else:
                np.append(out_b, npb[aux_b])
        
        slider += bin_size
    
    return out_a, out_b

def merge_and_share2(dfA, dfB, 
                    mz_col_name="mz",
                    int_col_name="int",
                    bin_size=0.025,
                    step=0.01,
                    inter=True):
    '''
    Merge mz rows in a bin size window, sliding it on a bin size step
    Also inform all the coincident peaks

    Parameters
    ----------
    dfA : Pandas dataframe
        Peak list in pandas dataframe.
    dfB : Pandas dataframe
        Peak list in pandas dataframe.
    mz_col_name : str, optional
        Column header for m/z peak list. The default is "mz".
    int_col_name : TYPE, optional
        Column header for intensity peak list. The default is "int".
    bin_size : float, optional
        Bin size window for merge. The default is 0.025.

    Returns
    -------
    dfA_merge : Pandas dataframe
        dfA with peaks merged.
    dfB_merge : Pandas dataframe
        dfB with peaks merged.
    dfA_inter : Pandas dataframe
        Shared peaks from dfA.
    dfB_inter : Pandas dataframe
        Shared peaks from dfA..

    '''
    # From pandas to numpy
    a=dfA[[mz_col_name, int_col_name]].to_numpy()
    b=dfB[[mz_col_name, int_col_name]].to_numpy()

    # Agregating    
    a = aggregate_mz(a, bin_size=bin_size, step=step)
    b = aggregate_mz(b, bin_size=bin_size, step=step)
    
    if inter:
        shr_a, shr_b = shared_mz(a,b,bin_size=bin_size)
        
   
    # Create output dataframes
    dfA_merge = pd.DataFrame({mz_col_name: a[:,0],
                              int_col_name: a[:,1]})
    
    dfB_merge = pd.DataFrame({mz_col_name: b[:,0],
                              int_col_name: b[:,1]})
    
    if inter:
        dfA_inter = pd.DataFrame({mz_col_name: shr_a[:,0],
                                  int_col_name: shr_a[:,1]})
    
        dfB_inter = pd.DataFrame({mz_col_name: shr_b[:,0],
                                  int_col_name: shr_b[:,1]})
    
        return dfA_merge, dfB_merge, dfA_inter, dfB_inter
    
    else:
        return dfA_merge, dfB_merge


print("Sem steps")
t0 = time.time()
A_merge, B_merge, A_int, B_int = merge_and_share(a,b)
t1 = time.time()
print(round(t1-t0,2), end=" secs\n")



print("Com steps, sem inter")
t0 = time.time()
A_mergenp, B_mergenp = merge_and_share2(a,b, inter=False)
t1 = time.time()
print(round(t1-t0,2), end=" secs\n")

print("Com steps, com inter")
t0 = time.time()
A_mergenp, B_mergenp, A_intnp, B_intnp = merge_and_share2(a,b)
t1 = time.time()
print(round(t1-t0,2), end=" secs\n")