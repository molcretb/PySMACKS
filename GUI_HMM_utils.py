# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 15:05:18 2026

@author: molcre0000
"""

import numpy as np
from hmmlearn import hmm, vhmm

from tkinter.filedialog import askopenfilenames
from tkinter import Tk
import matplotlib.pyplot as plt
import pandas as pd
from scipy import optimize

def TbT_init(DD_channel, DA_channel, nb_hidden_states = 2, dim = 2, HMM_iter = 10, return_dyn_param = 0):
    
    pi_init_0 = (1/nb_hidden_states) * np.ones(nb_hidden_states)
    a_0 = (1/nb_hidden_states) * np.ones((nb_hidden_states, nb_hidden_states))
        
    max_DD = np.max(DD_channel)
    max_DA = np.max(DA_channel)
    min_DD = np.min(DD_channel)
    min_DA = np.min(DA_channel)
    
    DD_init_means = np.linspace(min_DD, max_DD, nb_hidden_states)
    DA_init_means = np.linspace(min_DA, max_DA, nb_hidden_states)
    
    init_means = []
    
    for i in range(nb_hidden_states):
        # states are going from low FRET (state 0) to high FRET (state = nb_hidden_states -1)
        
        init_means.append([DD_init_means[nb_hidden_states-i-1], DA_init_means[i]])
        
    X = np.stack([DD_channel, DA_channel], axis = 1)
        
    best_score = best_model = None
        
    for idx in range(HMM_iter):
        
        model = hmm.GaussianHMM(n_components= nb_hidden_states, covariance_type="full", init_params='sct', params='smct', random_state=idx)  # params='stmc'
        model.n_features = dim
        model.means_ = init_means

        model.fit(X)
            
        score = model.score(X)
            
        if best_score is None or score > best_score:
            best_model = model
            best_score = score

    result = best_model.predict(X)
        
    HMM_TbT_param = {'start_prob': best_model.startprob_, 'trans_mat': best_model.transmat_,'means': best_model.means_,'covar_mat': best_model.covars_, 'predict': result}
    if return_dyn_param == 0:
        return result
    else:
        return HMM_TbT_param

def TbT_HMM_pipeline(traces_dict, nb_hidden_states = 2, dim = 2, HMM_iter = 10):
    
    # we assume traces_dict to be an OpenFRET dataset
    
    #traces_dict_ID = list(traces_dict.keys())
    
    nb_traces = len(traces_dict.traces)
    
    compt = 0
    
    stack_rate_mat = np.zeros((nb_traces, nb_hidden_states, nb_hidden_states))
    
    stack_prob_init = np.zeros((nb_traces, nb_hidden_states))
    
    for i in range(nb_traces):
        print(compt+1)
        # DD_channel = traces_dict[i]['Intensity_DD']
        # DA_channel = traces_dict[i]['Intensity_DA']
        DD_channel = traces_dict.traces[i].channels[0].data
        DA_channel = traces_dict.traces[i].channels[1].data
        
        # AA_channel = traces_dict[i]['Intensity_AA']
        
        HMM_TbT_param_i = TbT_init(DD_channel, DA_channel, nb_hidden_states = nb_hidden_states, dim = dim, HMM_iter = HMM_iter, return_dyn_param = 1)
    
        #traces_dict[i]['HMM_param'] = HMM_TbT_param_i
        
        traces_dict.traces[i].metadata['HMM_param'] = HMM_TbT_param_i
        
        stack_rate_mat[compt,:,:] = HMM_TbT_param_i['trans_mat']
        
        stack_prob_init[compt,:] = HMM_TbT_param_i['start_prob']
        
        compt += 1
    
    
    
    return traces_dict, stack_rate_mat, stack_prob_init
    

def plot_TbT_kin_rates(stack_rate_mat, framerate = 10): #framerate corresponds to the number of frame per seconds, 10 fps by default
    nb_hidden_states = stack_rate_mat.shape[1]
    
    # low FRET is state 0, high FRET is state 1
    
    compt = 0
    rate_dict = {}
    
    fig, ax = plt.subplots()
    
    for i in range(nb_hidden_states):
        for j in range(nb_hidden_states):
            if i != j:
                label_i = r'$k_{' + str(i) + str(j) + r'}$'
                rate_dict[label_i] = stack_rate_mat[:,i,j]*framerate
                #ax[compt].boxplot(stack_rate_mat[:,i,j], label=label_i)
                #ax[compt].legend()
                compt += 1
    ax.boxplot(rate_dict.values(), tick_labels=rate_dict.keys())
    plt.ylabel(r'Kinetic rate ($s^{-1}$)')
    plt.show()

def plot_trace_with_HMM(traces_dict, ID = 0, scale = 1):
    
    traces_dict_ID = list(traces_dict.keys())
    
    i = traces_dict_ID[ID]
    
    DD_channel = traces_dict[i]['Intensity_DD']
    DA_channel = traces_dict[i]['Intensity_DA']
    
    fig, ax = plt.subplots()
    
    ax.plot(DD_channel, 'orange', label='DD')
    ax.plot(DA_channel, 'red', label='DA')
    ax.plot(traces_dict[i]['HMM_param']['predict']*np.max([DD_channel, DA_channel])*scale, 'k--', label='HMM')
    ax.legend()
    
def calc_sum_neg_scores(a, pi, traces_dict, nb_hidden_states = 2, dim = 2):
    # nb_traces = len(dataset_traces)
    #traces_dict_ID = list(traces_dict.keys())
    nb_traces = len(traces_dict.traces)
    sum_neg_scores = 0
    #for k in traces_dict_ID:
    for k in range(nb_traces):
        try:
            model = hmm.GaussianHMM(n_components= nb_hidden_states, covariance_type="full", init_params='mc', params='mc')
            model.startprob_ = pi
            model.transmat_ = a
            model.n_features = dim
            # model.means_ = traces_dict[k]['HMM_param']['means']
            # model.covars_ = traces_dict[k]['HMM_param']['covar_mat']
            model.means_ = traces_dict.traces[k].metadata['HMM_param']['means']
            model.covars_ = traces_dict.traces[k].metadata['HMM_param']['covar_mat']
            
            # DD_channel = np.array(traces_dict[k]['Intensity_DD'])
            # DA_channel = np.array(traces_dict[k]['Intensity_DA'])
            DD_channel = np.array(traces_dict.traces[k].channels[0].data)
            DA_channel = np.array(traces_dict.traces[k].channels[1].data)
            # AA_channel = np.array(traces_dict[k]['Intensity_AA'])
            
            # X = np.stack([DD_channel, DA_channel, AA_channel], axis = 1)
            
            X = np.stack([DD_channel, DA_channel], axis = 1)
            
            sum_neg_scores = sum_neg_scores - model.score(X)
        except ValueError:
            sum_neg_scores = sum_neg_scores
    return sum_neg_scores

def Ens_HMM(traces_dict, trans_mat_0, prob_init_0):
    """
    Minimize target_func with respect to a matrix variable where each row sums to 1.

    Parameters:
    - target_func: Function that takes a matrix (2D numpy array) as input and returns a scalar.
    - initial_matrix: Initial guess for the matrix, shape (m, n).

    Returns:
    - result_matrix: The optimized matrix satisfying the constraints.
    """

    m, n = trans_mat_0.shape
    
    # Flatten the initial matrix to a vector for the optimizer
    x0 = np.concat([trans_mat_0.flatten(), prob_init_0])

    # Define the constraint: each row sums to 1
    # For each row i, sum of variables x[i*n : (i+1)*n] == 1
    constraints = []
    for i in range(m):
        def row_sum_constraint(x, row=i):
            return np.sum(x[row*n:(row+1)*n]) - 1
        constraints.append({'type': 'eq', 'fun': row_sum_constraint})
    
    def prob_init_sum_constraint(x):
        return np.sum(x[-n:]) - 1
    constraints.append({'type': 'eq', 'fun': prob_init_sum_constraint})

    # Optional: bounds to keep variables between 0 and 1 (if needed)
    #bounds = [(0.0,1.0) for k in range(m * (n+1))]

    # Define the objective function for the optimizer
    def objective(x):
        trans_mat = x[0:m*n].reshape((m, n))
        pi = x[-n:]
        return calc_sum_neg_scores(trans_mat, pi, traces_dict, nb_hidden_states = n, dim = 2) # dim = 2 or 3

    # Run the optimizer
    #result = optimize.minimize(objective, x0, constraints=constraints, bounds=bounds, method='SLSQP')
    
    result = optimize.minimize(objective, x0, constraints=constraints, method='SLSQP')

    # Check if optimization was successful
    if not result.success:
        raise RuntimeError(f"Optimization failed: {result.message}")

    # Reshape the result back into a matrix
    Ens_trans_mat = result.x[0:m*n].reshape((m, n))
    Ens_start_prob = result.x[-n:]
    return Ens_trans_mat, Ens_start_prob