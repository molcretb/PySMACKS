# -*- coding: utf-8 -*-
"""
Created on Fri Jan  2 10:05:32 2026

@author: molcre0000
"""

import numpy as np

def gauss_em_prob_dens(mu, V, x):
    """
    Function that generates the Gaussian emission probability densities over time

    Parameters
    ----------
    mu : Numpy array, dimension d
        Mean values for state i
    V : Numy matrix, dimension d*d
        Covariance matrix for state i
    x : Numpy array, dimension d*T
        Observable vector of length T

    Returns
    -------
    gauss_dens_time : Numpy array, dimension T
        Gaussian emission probability densities for state i at every time

    """
    gauss_dens_time = np.exp(-0.5*np.diag(np.matmul(np.transpose(x-mu),np.matmul(np.linalg.inv(V),x-mu)))) \
    / ((2*np.pi)**(len(mu)/2)*np.linalg.det(V)**0.5)
    
    return gauss_dens_time

def forward_backward(pi_init, a, b, x):
    """
    Forward-Backward algorithm to generate the forward and backward probabilities

    Parameters
    ----------
    pi_init : Numpy array, dimension d*d
        Start probabilities
    a : Numpy array, dimension d*d
        Transition probabilities
    b : list of Numpy arrays
        List containing d 2-elements sublists with mean vectors and covariance matrices for eahc states
    x : Numpy array, dimension d*T
        Observable vector of length T

    Returns
    -------
    alpha : Numpy array
        Forward probabilities
    beta : Numpy array
        Backward probabilities
    prod_prob_alpha : float
        Production probabilities
    gauss_prob_list : Numpy array, dimension T*d
        Time series of Gauss emission prob. dens. for each state i

    """
    
    d = len(pi_init) # number of states
    
    T = x.shape[1] # length of the time series
    
    gauss_prob_list = np.zeros((T, d)) # list of all time series of Gauss emission prob. dens. for each state i
    
    for i in range(d):
        gauss_prob_list[:,i] = gauss_em_prob_dens(b[i][0], b[i][1], x) # b[i][0] and b[i][1] are respectively the mean values vector and covariance matrix for state i 
    
    # forward variable calculation
    
    alpha = np.zeros((T, d))
    
    # Initiation
    
    alpha[0,:] = np.matmul(pi_init, gauss_prob_list[0,:])
    
    # Recursion steps
    
    for k in range(1,T):
        alpha[k,:] = np.matmul(alpha[k-1,:], a) * gauss_prob_list[k,:]
        
    # Calculation of the production proability
        
    prod_prob_alpha = np.sum(alpha[T-1,:])
    
    # backward variable calculation
    
    beta = np.zeros((T, d))
    
    # Initiation
    
    beta[T-1,:] = 1
    
    # Recursion steps
    
    for k in range(T-2, -1, -1):
        beta[k,:] = np.matmul(a, beta[k+1,:]) * gauss_prob_list[k+1,:]
        
    # Calculation of the production proability, (should be the same as the alpha one)
    
    prod_prob_beta = np.sum(np.matmul(pi_init,gauss_prob_list[0,:]) * beta[0,:])
    
    return alpha, beta, prod_prob_alpha, gauss_prob_list

def Baum_Welch_algo(a, b, alpha, beta, prod_prob, gauss_prob_list, x, FRET_const = True):
    
    d = len(a) # number of states
    T = len(alpha) # length of the time series
    
    gamma = alpha * beta / prod_prob
    
    alpha_stack = np.moveaxis(np.stack(([alpha for i in range(d)]), axis = 0), 1, 0).T
    
    beta_stack = np.moveaxis(np.stack(([beta for i in range(d)]), axis = 0), 1, -1)
    
    beta_stack[:,:,0:-1] = beta_stack[:,:,1:]
    
    gauss_stack = np.moveaxis(np.stack(([gauss_prob_list for i in range(d)]), axis = 0), 1, -1)
    
    gauss_stack[:,:,0:-1] = gauss_stack[:,:,1:]
    gauss_stack[:,:,-1] = np.ones(gauss_stack[:,:,-1].shape)
    
    a_stack = np.moveaxis(np.stack(([a for i in range(T)]), axis = 0), 0, -1)
    
    gamma_trans = alpha_stack * a_stack * gauss_stack * beta_stack / prod_prob
    
    # Update parameters
    
    sum_gamma_i = np.sum(gamma, axis = 0)
    
    pi_update = gamma[0,:]
    
    a_update = np.sum(gamma_trans, axis = 2) / np.stack(([sum_gamma_i for i in range(d)]), axis = 0).T
    
    if FRET_const == False:
        gamma_stack = np.stack(([gamma for i in range(d)]), axis = 0).T
        mu_update = np.sum(gamma_stack*x.T, axis = 1) / sum_gamma_i
    else:
        # we assume column 0 and 1 of x are donor and acceptor channels
        I_tot = np.mean(np.sum(x.T[:,0:2], axis = 1).reshape((T,1))) # Average total donor + acceptor intensity
        gamma_stack = np.stack(([gamma for i in range(d)]), axis = 0).T
        mu_update = I_tot * np.sum(gamma_stack*x.T, axis = 1) / np.sum(np.sum(x.T[:,0:2], axis = 1).reshape((T,1)) * gamma, axis = 0)
    
    b_update = []
    
    for k in range(d):
        x_stack = np.zeros((d, d, T))
        for j in range(T):
            x_i = x.T[j,:].reshape((d,1))
            x_stack[:,:,j] = np.matmul(x_i, x_i.T)
        V_i = np.sum(gamma[:,k] * x_stack, axis = 2) / sum_gamma_i[k] - np.matmul(b[k][0], b[k][0].T)
        
        b_update.append([mu_update[k], V_i])