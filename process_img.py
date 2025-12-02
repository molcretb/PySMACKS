# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 16:07:52 2025

@author: molcre0000
"""

import numpy as np
from tkinter.filedialog import askopenfilenames
from tkinter import Tk
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image
from bigfish import detection, plot
from scipy.optimize import minimize
import ruptures as rpt
import cv2
#%matplotlib





def load_data():
    root = Tk(className='Open trajectories', )
    file_path = askopenfilenames()
    root.destroy()
    try:
    # Assuming data is in a structured format like CSV or similar
        img = Image.open(file_path[0])
    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
    #X = np.concatenate((df['<d_d>'].tolist(), df['<a_d>'].tolist()))
    #X = np.array(df['<d_d>'].tolist(), df['<a_d>'].tolist())
    nb_frame = img.n_frames
    img_stack = np.zeros((img.height, img.width, nb_frame))
    for i in range(nb_frame):
        img.seek(i)
        img_stack[:,:,i] = np.array(img)
    return img_stack, file_path

def load_submovies(filename):
    img = Image.open(filename)
    nb_frame = img.n_frames
    img_stack = np.zeros((img.height, img.width, nb_frame))
    for i in range(nb_frame):
        img.seek(i)
        img_stack[:,:,i] = np.array(img)
    return img_stack

def detect_spot(img_stack, plot_detect = 0, kernel_size=2, min_distance=2): #### !!!!!!!! img_stack is now (img.height, img.width, nb_frame), need to adapt the code
    # Z-projection of the stack with sum
    img_sum = np.sum(img_stack,axis=2)
    # rough background removal with average median of the Z-projection
    img_sum_filter = img_sum - np.median(img_sum)
    # spots detection using a LoG filter, need to specify the expected standard deviation of the spots as kernel size \
        # and the minimal distance between two spots
    coord_spots = detection.detect_spots(img_sum_filter, log_kernel_size=kernel_size, minimum_distance=min_distance)
    # plot of detected spots
    if plot_detect == 1 :
        plot.plot_detection(img_sum_filter*(img_sum_filter>=0),coord_spots, contrast=True)
    
    return coord_spots

def filter_close_prox_spots(coord_spots, min_dist = 7):
    xv, xh = np.meshgrid(coord_spots[:,0], coord_spots[:,0])
    yv, yh = np.meshgrid(coord_spots[:,1], coord_spots[:,1])
    dist = np.sqrt((xv - xh)**2+(yv - yh)**2)
    dist_filter = dist >= min_dist
    coords_filter = np.where(np.sum(dist_filter,axis = 1)==len(coord_spots)-1)
    # return the filtered spot coordinates
    return coord_spots[coords_filter]

def spot_trace_extractor(img_stack, coord_spots, filenames):  # to be completed
    traces_data = {}
    for k in range(len(coord_spots)):
        xl = np.max([coord_spots[k,0]-6, 0])
        xr = np.min([coord_spots[k,0]+ 6, img_stack.shape[0]])
        yu = np.max([coord_spots[k,1]-6, 0])
        yd = np.min([coord_spots[k,1]+ 6, img_stack.shape[1]])
        trace = np.max(np.max(img_stack[xl:xr,yu:yd,:], axis = 0),axis=0)  # example
        traces_data[str(k)] = {'x_coord': coord_spots[k,0], 'y_coord': coord_spots[k,1]}
        traces_data[str(k)]['Intensity'] = trace
    
    for j in range(len(filenames)-1):
        print('Processing movies ' + str(j))
        img_stack = load_submovies(filenames[j+1])
        for k in range(len(coord_spots)):
            xl = np.max([coord_spots[k,0]-6, 0])
            xr = np.min([coord_spots[k,0]+ 6, img_stack.shape[0]])
            yu = np.max([coord_spots[k,1]-6, 0])
            yd = np.min([coord_spots[k,1]+ 6, img_stack.shape[1]])
            trace = np.max(np.max(img_stack[xl:xr,yu:yd,:], axis = 0),axis=0)
            traces_data[str(k)]['Intensity'] = np.concat((traces_data[str(k)]['Intensity'], trace))
    return traces_data

def detect_bleaching_event(traces_data, percentage_value = 0.05):
    for i in range(len(traces_data)):
        print(i)
        one_perc = int(len(traces_data[str(i)]['Intensity'])*percentage_value)   # X% of the total movie duration, in frames
        trace_start_value = np.median(traces_data[str(i)]['Intensity'][0:one_perc])
        trace_final_value = np.median(traces_data[str(i)]['Intensity'][-one_perc:])
        std_trace = np.std(traces_data[str(i)]['Intensity'])
        if (trace_start_value - trace_final_value) > std_trace:
        # trace_start_value = np.std(traces_data[str(i)]['Intensity'][0:one_perc])
        # trace_final_value = np.std(traces_data[str(i)]['Intensity'][-one_perc:])
        # if trace_start_value <= 1.2 * trace_final_value:
            algo_c = rpt.KernelCPD(kernel="linear", min_size=2).fit(traces_data[str(i)]['Intensity'])
            result = algo_c.predict(n_bkps=1)
            traces_data[str(i)]['bleaching_event'] = result[0]
        else:
            traces_data[str(i)]['bleaching_event'] = 0
        #plt.plot(traces_data[str(i)]['Intensity']);
        #plt.axvline(x = result[0], color = 'r', label = 'predicted breakpoint');
        #plt.show()
    return traces_data

def plot_traces_by_10(traces_data, traces_ID, nb_plot = 10):
    """Function used to plot traces with detected bleach point by set of 10 traces (can be customized) """
    fig, axs = plt.subplots(nb_plot)
    fig.suptitle('Traces ' + str(traces_ID) + ' to ' + str(traces_ID+nb_plot-1))
    for k in range(nb_plot):
        axs[k].plot(traces_data[str(traces_ID + k)]['Intensity']);
        axs[k].axvline(x = traces_data[str(traces_ID + k)]['bleaching_event'], color = 'r', label = 'predicted breakpoint');
    plt.show()
        

def get_SNR(traces_data):
    for i in range(len(traces_data)):
        if traces_data[str(i)]['bleaching_event'] == 0:
            traces_data[str(i)]['SNR'] = 0
        else:
            bleach_event = traces_data[str(i)]['bleaching_event']
            trace_i = traces_data[str(i)]['Intensity']
            traces_data[str(i)]['SNR'] = (np.mean(trace_i[0:bleach_event]) - np.mean(trace_i[bleach_event:])) / np.std(trace_i[bleach_event:])
    return traces_data
    

def filter_traces_on_SNR(traces_data, SNR_thresh = 8):
    traces_high_SNR = {}
    for i in range(len(traces_data)):
        if traces_data[str(i)]['SNR'] >= SNR_thresh:
            traces_high_SNR[str(i)] = traces_data[str(i)]
    return traces_high_SNR, list(traces_high_SNR.keys())

def plot_traces_high_SNR(traces_high_SNR, high_SNR_IDs, pick_ID):
    fig, axs = plt.subplots(5)
    for k in range(5):
        back_noise = np.median(traces_high_SNR[str(high_SNR_IDs[pick_ID + k])]['Intensity'])
        # plot high SNR traces with background removed (taken as median of bleached section)
        axs[k].plot(traces_high_SNR[str(high_SNR_IDs[pick_ID + k])]['Intensity']-back_noise)
        axs[k].axvline(x = traces_high_SNR[str(high_SNR_IDs[pick_ID + k])]['bleaching_event'], color = 'r', label = 'predicted breakpoint')
    plt.show()

def logP_chrom_corr(x, img1, img2):
    H, W = img1.shape[0:2]
    # xmesh, ymesh = np.meshgrid(range(0,img1.shape[1]), range(0,img1.shape[0]))
    # zmesh = np.ones((img1.shape[0], img1.shape[1]))
    # #ori_coord = np.stack((xmesh, ymesh, zmesh))   might be needed later, keep it in case
    # ori_coord_flat = np.reshape(np.stack((xmesh, ymesh, zmesh)), shape = (3, img1.shape[0]*img1.shape[1]))
    
    # list_px_value_img1 = np.reshape(img1, shape = (1, img1.shape[0]*img1.shape[1]))
    # list_px_value_img2 = np.reshape(img2, shape = (1, img2.shape[0]*img2.shape[1]))
    
    matA = np.array([[x[0], 0, (1-x[0])*W/2],  \
                     [0, x[1], (1-x[1])*H/2], \
                     [0, 0, 1]])
        
    matB = np.array([[1, 0, x[2] + x[4]],  \
                     [0, 1, x[3] + x[5]], \
                     [0, 0, 1]])
        
    matC = np.array([[np.cos(x[6]), np.sin(x[6]), 0],  \
                     [-np.sin(x[6]), np.cos(x[6]), 0], \
                     [0, 0, 1]])
        
    matD = np.array([[1, 0, -x[4]],  \
                     [0, 1, -x[5]], \
                     [0, 0, 1]])
        
    omega = np.matmul(matA, np.matmul(matB, np.matmul(matC, matD)))
    
    
    trans_img1 = cv2.warpPerspective(img1,omega, (W, H))
    
    # trans_coord_flat = np.matmul(omega, ori_coord_flat).astype(int)
    
    # filter_x = [0 <= x < W for x in trans_coord_flat[0,:]]
    # filter_y = [0 <= x < H for x in trans_coord_flat[1,:]]
    
    # comb_filter = [a and b for a, b in zip(filter_x, filter_y)]
    # # trans_coord_flat_filter contains the coordinates of the pixels we need to extract from img2
    # trans_coord_flat_filter = trans_coord_flat[:,comb_filter]  #list of all transformed pixel that fit into the image (coordinates comrpised in a W * H image)
    
    # intensity_trans_img1 = list_px_value_img1[0,comb_filter]
    # intensity_corres_img2 = img2[trans_coord_flat_filter[0,:],trans_coord_flat_filter[1,:]]
    
    # neglogP = 0.5 * W * H * np.log(np.sum((intensity_trans_img1 - intensity_corres_img2)**2))   # negqtive logP, so we can use the minimize function (instead of maximize) 
    
    neglogP = 0.5 * W * H * np.log(np.sum((trans_img1-img2)**2))
    
    return neglogP

def minimize_logP(img1, img2):
    H, W = img1.shape[0:2]
    # we want to minimize -logP with regard to sx, sy, dx, dy, cx, cy and theta
    res = minimize(logP_chrom_corr, [1, 1, 1, 1, int(W/2), int(H/2), 0.01], method='Nelder-Mead', args=(img1, img2), options={'maxfev': 500})
    return res

def generate_chrom_ab_corr_movie(img1, x, plot_img = 0):
    H, W = img1.shape[0:2]
    
    matA = np.array([[x[0], 0, (1-x[0])*W/2],  \
                     [0, x[1], (1-x[1])*H/2], \
                     [0, 0, 1]])
        
    matB = np.array([[1, 0, x[2] + x[4]],  \
                     [0, 1, x[3] + x[5]], \
                     [0, 0, 1]])
        
    matC = np.array([[np.cos(x[6]), np.sin(x[6]), 0],  \
                     [-np.sin(x[6]), np.cos(x[6]), 0], \
                     [0, 0, 1]])
        
    matD = np.array([[1, 0, -x[4]],  \
                     [0, 1, -x[5]], \
                     [0, 0, 1]])
        
    omega = np.matmul(matA, np.matmul(matB, np.matmul(matC, matD)))
    
    trans_img1 = cv2.warpPerspective(img1,omega, (W, H))
    
    if plot_img == 1:
        plt.figure();plt.imshow(trans_img1)
    
    return trans_img1 
    
    
    