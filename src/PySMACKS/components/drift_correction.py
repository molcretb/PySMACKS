# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 14:12:05 2026

@author: molcre0000
"""
import numpy as np
from bigfish import detection
from PIL import Image
from PySMACKS.components.chromatic_aberrations_correction import *
from PySMACKS.components.utils import *
import trackpy as tp
import pandas as pd
from tkinter.filedialog import askopenfilenames
from tkinter import Tk
from scipy.signal import savgol_filter
from PySMACKS.components.traces_extractor import filter_close_prox_spots

def correct_drift(file_path, kernel_size=2, min_distance=2, ref_AA = 0, perc_frames = 0.05, get_coord = 0, Z_project = 1):
    
    time_frames = 0
    
    comp_frame_drift = 0
    
    for k in range(len(file_path)):
        print('Processing movie ' + str(k))
        img = Image.open(file_path[k])
        nb_frame = img.n_frames
        if k < len(file_path) - 1:
            nb_frame_old = img.n_frames
        if Z_project == 0:
            list_frames = [int(i) for i in np.floor(np.linspace(0,nb_frame-1,int(nb_frame*perc_frames)))]
            for j in list_frames:
                img.seek(j)
                img_raw = np.array(img)
                coord_spots = detection.detect_spots(img_raw, log_kernel_size=kernel_size, minimum_distance=min_distance)
                if (k+j) == 0:
                    df_coords = pd.DataFrame(data={'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots)), 'frame':np.zeros(len(coord_spots))})
                    time_frames = nb_frame_old*k + j
                    comp_frame_drift = comp_frame_drift + 1
                else:                    
                    time_frames = nb_frame_old*k + j
                    new_coords = {'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots))+time_frames, 'frame':np.zeros(len(coord_spots)) + comp_frame_drift}
                    df_coords = pd.concat([df_coords, pd.DataFrame(data=new_coords)]).reset_index(drop=True)
                    comp_frame_drift = comp_frame_drift + 1
        else:
            img.seek(0)
            img0 = np.array(img)
            dim_img = img0.shape
            img_raw = np.zeros(((dim_img[0], dim_img[1], nb_frame)))
            img_raw[:,:,0] = img0
            for j in range(1,nb_frame):
                img.seek(j)
                img_raw[:,:,j] = np.array(img)
            img_raw_projZ = np.sum(img_raw, axis = 2)
            coord_spots = detection.detect_spots(img_raw_projZ, log_kernel_size=kernel_size, minimum_distance=min_distance)
            
            if k ==0:
                df_coords = pd.DataFrame(data={'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots)), 'frame':np.zeros(len(coord_spots))})
                time_frames = time_frames + nb_frame
            else:
                new_coords = {'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots))+time_frames, 'frame':np.zeros(len(coord_spots)) + k}
                df_coords = pd.concat([df_coords, pd.DataFrame(data=new_coords)]).reset_index(drop=True)
                time_frames = time_frames + nb_frame
        
    df_linked = tp.link(df_coords, 2, memory=time_frames)  # 5
    
    
    list_frames_drift = df_linked['frame_real'].unique()
    drift = tp.compute_drift(df_linked)
    
    spot_IDs = df_linked['particle'].unique()
    
    coord_spots_track = np.zeros((len(spot_IDs), 2))
    
    if get_coord == 1:
    
        for j in spot_IDs:
            df_i = df_linked[df_linked['particle']==j]
            if len(df_i) >  2:
                min_frame_i = np.min(df_i['frame'])
                if min_frame_i == 0:
                    drift_corr_i = [0, 0]
                else:
                    drift_corr_i = [drift['y'][min_frame_i], drift['x'][min_frame_i]]
                    
                y_spot = int(df_i[df_i['frame']==min_frame_i]['y'].tolist()[0] - drift_corr_i[0])
                x_spot = int(df_i[df_i['frame']==min_frame_i]['x'].tolist()[0] - drift_corr_i[1])
                
                coord_spots_track[j, 0] = y_spot
                coord_spots_track[j, 1] = x_spot
            
        mask = ~np.any(coord_spots_track <= 0, axis=1)
        coord_spots_track = coord_spots_track[mask]
    
    return drift, coord_spots_track, list_frames_drift


def correct_drift_ALEX(file_path, kernel_size=2, min_distance=2, ref_AA = 0, perc_frames = 0.05, get_coord = 0, Z_project = 1):
    
    time_frames = 0
    
    comp_frame_drift = 0
    
    for k in range(len(file_path)):
        print('Processing movie ' + str(k))
        img = Image.open(file_path[k])
        nb_frame = img.n_frames
        if k < len(file_path) - 1:
            nb_frame_old = img.n_frames
        if Z_project == 0:
            list_frames = [int(i)*2+ref_AA for i in np.floor(np.linspace(0,int(nb_frame/2)-1,int(nb_frame/2*perc_frames)))]
            for j in list_frames:
                img.seek(j)
                img_raw = np.array(img)
                coord_spots = detection.detect_spots(img_raw, log_kernel_size=kernel_size, minimum_distance=min_distance)
                if (k+j) == 0:
                    df_coords = pd.DataFrame(data={'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots)), 'frame':np.zeros(len(coord_spots))})
                    time_frames = nb_frame_old*k + j
                    comp_frame_drift = comp_frame_drift + 1
                else:                    
                    time_frames = nb_frame_old*k + j
                    new_coords = {'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots))+time_frames, 'frame':np.zeros(len(coord_spots)) + comp_frame_drift}
                    df_coords = pd.concat([df_coords, pd.DataFrame(data=new_coords)]).reset_index(drop=True)
                    comp_frame_drift = comp_frame_drift + 1
        else:
            img.seek(0)
            img0 = np.array(img)
            dim_img = img0.shape
            img_raw = np.zeros(((dim_img[0], dim_img[1], int(nb_frame/2))))
            img_raw[:,:,0] = img0
            for j in range(1,int(nb_frame/2)):
                img.seek(2*j)
                img_raw[:,:,j] = np.array(img)
            img_raw_projZ = np.sum(img_raw, axis = 2)
            coord_spots = detection.detect_spots(img_raw_projZ, log_kernel_size=kernel_size, minimum_distance=min_distance)
            
            if k ==0:
                df_coords = pd.DataFrame(data={'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots)), 'frame':np.zeros(len(coord_spots))})
                time_frames = time_frames + int(nb_frame/2)
            else:
                new_coords = {'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots))+time_frames, 'frame':np.zeros(len(coord_spots)) + k}
                df_coords = pd.concat([df_coords, pd.DataFrame(data=new_coords)]).reset_index(drop=True)
                time_frames = time_frames + int(nb_frame/2)
        
    df_linked = tp.link(df_coords, 2, memory=time_frames)  # 5
    
    
    list_frames_drift = df_linked['frame_real'].unique()
    drift = tp.compute_drift(df_linked)
    
    spot_IDs = df_linked['particle'].unique()
    
    coord_spots_track = np.zeros((len(spot_IDs), 2))
    
    
    if get_coord == 1:
        for j in spot_IDs:
            df_i = df_linked[df_linked['particle']==j]
            if len(df_i) > 2:
                min_frame_i = np.min(df_i['frame'])
                if min_frame_i == 0:
                    drift_corr_i = [0, 0]
                else:
                    drift_corr_i = [drift['y'][min_frame_i], drift['x'][min_frame_i]]
                    
                y_spot = int(df_i[df_i['frame']==min_frame_i]['y'].tolist()[0] - drift_corr_i[0])
                x_spot = int(df_i[df_i['frame']==min_frame_i]['x'].tolist()[0] - drift_corr_i[1])
                
                coord_spots_track[j, 0] = y_spot
                coord_spots_track[j, 1] = x_spot
            
        mask = ~np.any(coord_spots_track <= 0, axis=1)
        coord_spots_track = coord_spots_track[mask]
    
    return drift, coord_spots_track, (list_frames_drift-ref_AA)/2

def extract_DA_spots_drift(file_path, drift, kernel_size=3, min_distance=2, ref_AA = 1, perc_frames = 0.05, Z_project = 1):
    
    time_frames = 0
    
    comp_frame_drift = 0
    
    df_coords = []
    
    for k in range(len(file_path)):
        print('Processing movie ' + str(k))
        img = Image.open(file_path[k])
        nb_frame = img.n_frames
        if k < len(file_path) - 1:
            nb_frame_old = img.n_frames
        if Z_project == 0:
            list_frames = [int(i)*2+ref_AA for i in np.floor(np.linspace(0,int(nb_frame/2)-1,int(nb_frame/2*perc_frames)))]
            for j in list_frames:
                img.seek(j)
                img_raw = np.array(img)
                coord_spots = detection.detect_spots(img_raw, log_kernel_size=kernel_size, minimum_distance=min_distance)
                if (k+j-ref_AA) == 0:
                    df_coords = pd.DataFrame(data={'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots)), 'frame':np.zeros(len(coord_spots))})
                    time_frames = nb_frame_old*k + j
                    comp_frame_drift = comp_frame_drift + 1
                else:
                    
                    time_frames = nb_frame_old*k + j
                    new_coords = {'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots))+time_frames, 'frame':np.zeros(len(coord_spots)) + comp_frame_drift}
                    df_coords = pd.concat([df_coords, pd.DataFrame(data=new_coords)]).reset_index(drop=True)
                    comp_frame_drift = comp_frame_drift + 1
        else:
            img.seek(1)
            img0 = np.array(img)
            dim_img = img0.shape
            img_raw = np.zeros(((dim_img[0], dim_img[1], int(nb_frame/2))))
            img_raw[:,:,0] = img0
            for j in range(1,int(nb_frame/2)):
                img.seek(2*j+1)
                img_raw[:,:,j] = np.array(img)
            img_raw_projZ = np.sum(img_raw, axis = 2)
            coord_spots = detection.detect_spots(img_raw_projZ, log_kernel_size=kernel_size, minimum_distance=min_distance)
            
            if k ==0:
                df_coords = pd.DataFrame(data={'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots)), 'frame':np.zeros(len(coord_spots))})
                time_frames = time_frames + int(nb_frame/2)
            else:
                new_coords = {'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots))+time_frames, 'frame':np.zeros(len(coord_spots)) + k}
                df_coords = pd.concat([df_coords, pd.DataFrame(data=new_coords)]).reset_index(drop=True)
                time_frames = time_frames + int(nb_frame/2)
                        
    df_linked = tp.link(df_coords, 3, memory=time_frames, adaptive_stop = 0.5, adaptive_step = 0.95)  # 5
    
    
    list_frames_drift = df_linked['frame_real'].unique()
    
    spot_IDs = df_linked['particle'].unique()
    
    coord_spots_track = np.zeros((len(spot_IDs), 2))
    
    for j in spot_IDs:
        df_i = df_linked[df_linked['particle']==j]
        if len(df_i) > 2:
            min_frame_i = np.min(df_i['frame'])
            if min_frame_i == 0:
                drift_corr_i = [0, 0]
            else:
                drift_corr_i = [drift['y'][min_frame_i-1], drift['x'][min_frame_i-1]]
                
            y_spot = int(df_i[df_i['frame']==min_frame_i]['y'].tolist()[0] - drift_corr_i[0])
            x_spot = int(df_i[df_i['frame']==min_frame_i]['x'].tolist()[0] - drift_corr_i[1])
            
            coord_spots_track[j, 0] = y_spot
            coord_spots_track[j, 1] = x_spot
        
    mask = ~np.any(coord_spots_track <= 0, axis=1)
    coord_spots_track = coord_spots_track[mask]
    
    return coord_spots_track

def extract_DA_spots_drift_No_ALEX(file_path, drift, kernel_size=3, min_distance=2, perc_frames = 0.05, Z_project = 1):
    
    time_frames = 0
    
    comp_frame_drift = 0
    
    df_coords = []
    
    for k in range(len(file_path)):
        print('Processing movie ' + str(k))
        img = Image.open(file_path[k])
        nb_frame = img.n_frames
        if k < len(file_path) - 1:
            nb_frame_old = img.n_frames
        if Z_project == 0:
            list_frames = [int(i) for i in np.floor(np.linspace(0,nb_frame-1,int(nb_frame*perc_frames)))]
            for j in list_frames:
                img.seek(j)
                img_raw = np.array(img)
                coord_spots = detection.detect_spots(img_raw, log_kernel_size=kernel_size, minimum_distance=min_distance)
                if (k+j) == 0:
                    df_coords = pd.DataFrame(data={'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots)), 'frame':np.zeros(len(coord_spots))})
                    time_frames = nb_frame_old*k + j
                    comp_frame_drift = comp_frame_drift + 1
                else:                    
                    time_frames = nb_frame_old*k + j
                    new_coords = {'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots))+time_frames, 'frame':np.zeros(len(coord_spots)) + comp_frame_drift}
                    df_coords = pd.concat([df_coords, pd.DataFrame(data=new_coords)]).reset_index(drop=True)
                    comp_frame_drift = comp_frame_drift + 1
        else:
            img.seek(0)
            img0 = np.array(img)
            dim_img = img0.shape
            img_raw = np.zeros((dim_img[0], dim_img[1], nb_frame))
            img_raw[:,:,0] = img0
            for j in range(1,nb_frame):
                img.seek(j)
                img_raw[:,:,j] = np.array(img)
            img_raw_projZ = np.sum(img_raw, axis = 2)
            coord_spots = detection.detect_spots(img_raw_projZ, log_kernel_size=kernel_size, minimum_distance=min_distance)
            
            if k ==0:
                df_coords = pd.DataFrame(data={'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots)), 'frame':np.zeros(len(coord_spots))})
                time_frames = time_frames + nb_frame
            else:
                new_coords = {'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots))+time_frames, 'frame':np.zeros(len(coord_spots)) + k}
                df_coords = pd.concat([df_coords, pd.DataFrame(data=new_coords)]).reset_index(drop=True)
                time_frames = time_frames + nb_frame
                
        
    df_linked = tp.link(df_coords, 3, memory=time_frames, adaptive_stop = 0.5, adaptive_step = 0.95)  # 5
    
    
    list_frames_drift = df_linked['frame_real'].unique()
    
    spot_IDs = df_linked['particle'].unique()
    
    coord_spots_track = np.zeros((len(spot_IDs), 2))
    
    for j in spot_IDs:
        df_i = df_linked[df_linked['particle']==j]
        if len(df_i) > 2:
            min_frame_i = np.min(df_i['frame'])
            if min_frame_i == 0:
                drift_corr_i = [0, 0]
            else:
                drift_corr_i = [drift['y'][min_frame_i-1], drift['x'][min_frame_i-1]]
                
            y_spot = int(df_i[df_i['frame']==min_frame_i]['y'].tolist()[0] - drift_corr_i[0])
            x_spot = int(df_i[df_i['frame']==min_frame_i]['x'].tolist()[0] - drift_corr_i[1])
            
            coord_spots_track[j, 0] = y_spot
            coord_spots_track[j, 1] = x_spot
        
    mask = ~np.any(coord_spots_track <= 0, axis=1)
    coord_spots_track = coord_spots_track[mask]
    
    return coord_spots_track


def add_drift_to_traces_dict(traces_dict, drift):
    """
    Function that add the XY microscope stage drift to the traces dictionary

    Parameters
    ----------
    traces_dict : dict
        Dictionnary containing all the traces results.
    drift : Numpy array
        YX drift of the microscope stage

    Returns
    -------
    traces_dict : dict
        Dictionnary containing all the traces results with updated drift.

    """
    
    traces_dict['drift'] = drift
    
    return traces_dict

def correct_drift_single_sub_movie(img_stack, kernel_size=1.5, min_distance=2, ref_AA = 0):
    
    nb_frame = img_stack.shape[2]
    
    for k in range(nb_frame):
        print('Processing frame ' + str(k))
        coord_spots = detection.detect_spots(img_stack[:,:,k], log_kernel_size=kernel_size, minimum_distance=min_distance)
        if k == 0:
            df_coords = pd.DataFrame(data={'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame':np.zeros(len(coord_spots))})
        else:
            new_coords = {'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame':np.zeros(len(coord_spots))+k}
            df_coords = pd.concat([df_coords, pd.DataFrame(data=new_coords)]).reset_index(drop=True)
        
    df_linked = tp.link(df_coords, 5, memory=nb_frame)  # 5
    drift = tp.compute_drift(df_linked)
    
    spot_IDs = df_linked['particle'].unique()
    
    coord_spots_track = np.zeros((len(spot_IDs), 2))
    
    for j in spot_IDs:
        df_i = df_linked[df_linked['particle']==j]
        min_frame_i = np.min(df_i['frame'])
        if min_frame_i == 0:
            drift_corr_i = [0, 0]
        else:
            drift_corr_i = [drift['y'][min_frame_i], drift['x'][min_frame_i]]
            
        y_spot = int(df_i[df_i['frame']==min_frame_i]['y'].tolist()[0] - drift_corr_i[0])
        x_spot = int(df_i[df_i['frame']==min_frame_i]['x'].tolist()[0] - drift_corr_i[1])
        
        coord_spots_track[j, 0] = y_spot
        coord_spots_track[j, 1] = x_spot
        
    mask = ~np.any(coord_spots_track < 0, axis=1)
    coord_spots_track = coord_spots_track[mask]
    
    return drift, coord_spots_track


def correct_drift_syn_movie(movie_syn, kernel_size=1.5, min_distance=2, ref_AA = 0, perc_frames = 0.05, get_coord = 1, Z_project = 0):
    
    time_frames = 0
    
    comp_frame_drift = 0
    
    for k in range(movie_syn.shape[2]):
        print('Processing frame ' + str(k))
        
        if Z_project == 0:
            if k == 0:
                img_raw = np.sum(movie_syn[:,:,0:2], axis = 2).astype(float)
            elif k == movie_syn.shape[2]-1:
                img_raw = np.sum(movie_syn[:,:,k-2:k], axis = 2).astype(float)
            else:
                img_raw = np.sum(movie_syn[:,:,k-1:k+1], axis = 2).astype(float)
            coord_spots = detection.detect_spots(img_raw, log_kernel_size=kernel_size, minimum_distance=min_distance)
            if k == 0:
                df_coords = pd.DataFrame(data={'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots)), 'frame':np.zeros(len(coord_spots))})
            else:
                new_coords = {'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots))+time_frames, 'frame':np.zeros(len(coord_spots)) + k}
                df_coords = pd.concat([df_coords, pd.DataFrame(data=new_coords)]).reset_index(drop=True)
                comp_frame_drift = comp_frame_drift + 1
        else:
            img.seek(0)
            img0 = np.array(img)
            dim_img = img0.shape
            img_raw = np.zeros(((dim_img[0], dim_img[1], nb_frame)))
            img_raw[:,:,0] = img0
            for j in range(1,nb_frame):
                img.seek(j)
                img_raw[:,:,j] = np.array(img)
            img_raw_projZ = np.sum(img_raw, axis = 2)
            coord_spots = detection.detect_spots(img_raw_projZ, log_kernel_size=kernel_size, minimum_distance=min_distance)
            
            if k ==0:
                df_coords = pd.DataFrame(data={'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots)), 'frame':np.zeros(len(coord_spots))})
                time_frames = time_frames + nb_frame
            else:
                new_coords = {'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots))+time_frames, 'frame':np.zeros(len(coord_spots)) + k}
                df_coords = pd.concat([df_coords, pd.DataFrame(data=new_coords)]).reset_index(drop=True)
                time_frames = time_frames + nb_frame
        
    df_linked = tp.link(df_coords, 5, memory=time_frames)  # 5
    
    
    list_frames_drift = df_linked['frame_real'].unique()
    drift_estimated = tp.compute_drift(df_linked)
    
    spot_IDs = df_linked['particle'].unique()
    
    coord_spots_track = np.zeros((len(spot_IDs), 2))
    
    if get_coord == 1:
    
        for j in spot_IDs:
            df_i = df_linked[df_linked['particle']==j]
            if len(df_i) >  2:
                min_frame_i = np.min(df_i['frame'])
                if min_frame_i == 0:
                    drift_corr_i = [0, 0]
                else:
                    drift_corr_i = [drift_estimated['y'][min_frame_i], drift_estimated['x'][min_frame_i]]
                    
                y_spot = int(df_i[df_i['frame']==min_frame_i]['y'].tolist()[0] - drift_corr_i[0])
                x_spot = int(df_i[df_i['frame']==min_frame_i]['x'].tolist()[0] - drift_corr_i[1])
                
                coord_spots_track[j, 0] = y_spot
                coord_spots_track[j, 1] = x_spot
            
        mask = ~np.any(coord_spots_track <= 0, axis=1)
        coord_spots_track = coord_spots_track[mask]
    
    return drift_estimated, coord_spots_track, list_frames_drift