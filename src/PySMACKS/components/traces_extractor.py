# -*- coding: utf-8 -*-
"""
Created on Thu Dec  4 11:01:07 2025

@author: Bastien Molcrette, Schmid lab, Chemistry department, University of Basel
"""
import numpy as np
from bigfish import detection, plot
import ruptures as rpt
from scipy.ndimage import median_filter
from PySMACKS.components.chromatic_aberrations_correction import *
from PySMACKS.components.utils import *
import cv2
from TVDCondat2013 import tvd_2013
from openfret import Dataset, Trace, Channel, Metadata, write_data
from datetime import date
from uuid_extensions import uuid7


def detect_spot(img_DD, img_DA, kernel_size=2, min_distance=2, proj_oper = 'sum', size_median_filt = 10):
    """
    Function used to detect the spots xy localizations from both the DD and DA channels

    Parameters
    ----------
    img_DD : Numpy array
        DD channel.
    img_DA : Numpy array
        DA channel.
    kernel_size : int, optional
        standard deviation of the LoG kernel to be used as spot detector, in pixels. The default is 2.
    min_distance : int, optional
        minimal distance that is considered to distinguish two adjacent spots, in pixels. The default is 2.
    proj_oper : str, optional
        detection method that is used for the Z-projection, default is sum projection, anything else would be considered as max projection
    size_median_filt : int, optional
        size in pixels of the median filter radius for the rough background correction. Default is 10.

    Returns
    -------
    coord_spots : Numpy array
        List of coordinates of the detected spots.

    """
    print('Start of the spots detection step...')
    
    # get the min and max pixel intensity values from the DD and DA channels respectively, to normalize each channel during the Z-projection
    min_DD = np.min(np.min(img_DD, axis = 1), axis = 0)
    min_DA = np.min(np.min(img_DA, axis = 1), axis = 0)
    max_DD = np.max(np.max(img_DD, axis = 1), axis = 0)
    max_DA = np.max(np.max(img_DA, axis = 1), axis = 0)
    
    # Spot detection is performed either by sum or max Z-projection\
        # sum is likely to be more robust toward high noisy pixels but it might lose the short traces (spots that bleach quickly), especially with long movies
    if proj_oper == 'sum':
        # Z-projection of the stack with sum
        print('Spot detection performed with sum Z-projection')
        img_sum = np.sum((img_DD-min_DD)/(max_DD - min_DD)+(img_DA-min_DA)/(max_DA - min_DA),axis=2)
    else:
        # Z-projection of the stack with max
        print('Spot detection performed with max Z-projection')
        img_sum = np.max((img_DD-min_DD)/(max_DD - min_DD)+(img_DA-min_DA)/(max_DA - min_DA),axis=2)
    
    # rough background removal with local median of the Z-projection
    img_back = median_filter(img_sum, size = size_median_filt)
    img_sum_filter = img_sum - img_back
    
    # spots detection using a LoG filter, need to specify the expected standard deviation of the spots as kernel size \
        # and the minimal distance between two spots
    coord_spots = detection.detect_spots(img_sum_filter*(img_sum_filter>=0), log_kernel_size=kernel_size, minimum_distance=min_distance)
    
    print(str(len(coord_spots)) + ' spots detected')
    
    # plot of detected spots
    plot_detect_spots(img_sum_filter*(img_sum_filter>=0),coord_spots)
    
    #return the set of spots xy-coordinates
    return coord_spots

def plot_detect_spots(img, coord_spots, factor_con_max = 0.5, factor_con_min = 1):
    """
    Function used to plot the results of the spots detector, from the bigfish library

    Parameters
    ----------
    img : 2D Numpy array
        Original frame.
    coord_spots : Numpy array
        List of coordinates of the detected spots.

    Returns
    -------
    None.

    """
    
    plot_detection_modif_BM(img, coord_spots, factor_con_max = factor_con_max, factor_con_min = factor_con_min)
    
    return

def filter_close_prox_spots(coord_spots, min_dist = 7):
    """
    Function used to filter the spots that are too close from each other, so it might mess with the traces extraction

    Parameters
    ----------
    coord_spots : Numpy array
        List of coordinates of the detected spots.
    min_dist : int, optional
        Minimal distance in pixels between spots; every pair of spots whose inter-spots distance is less than this value will be removed. The default is 7.

    Returns
    -------
    filter_coord_spots  : Numpy array
        List of coordinates of the distance-filtered detected spots.

    """
    
    print('Removing the clusters of spots (interspots distance less than ' + str(min_dist) + ' pixels)')
    
    # matrices used to build an inter-spots distances matrix
    xv, xh = np.meshgrid(coord_spots[:,0], coord_spots[:,0])
    yv, yh = np.meshgrid(coord_spots[:,1], coord_spots[:,1])
    
    # matrix with distances between spots i and j (i = row, j = column)
    dist = np.sqrt((xv - xh)**2+(yv - yh)**2)
    
    # filter the spots with interspots distance more than minimal distance
    dist_filter = dist >= min_dist
    
    # get the spots IDs that are isolated
    coords_filter = np.where(np.sum(dist_filter,axis = 1)==len(coord_spots)-1)
    
    filter_coord_spots = coord_spots[coords_filter]
    
    print(str(len(filter_coord_spots)) + ' remaining spots (' + str(len(coord_spots)-len(filter_coord_spots)) + ' spots removed)')
    
    # return the filtered spot coordinates
    return filter_coord_spots 


def spot_trace_extractor(file_path_D, file_path_A, coord_spots, matrix_align, 
                         drift_correct = 0, 
                         method_align = 'Optical Flow',
                         sigma = 3, 
                         DA_is = 'odd', 
                         back_corr = 'median_TV_min',
                         TV_lam = 5):
    """
    Function used to extract and concatenate the subtraces from DD, DA and AA channels from all submovies

    Parameters
    ----------
    coord_spots : Numpy array
        List of coordinates of the distance-filtered detected spots.
    img_stack_D : Numpy array
        First donor submovie as Numpy array.
    file_path_D : List of str
        List of paths of all donor submovies
    img_stack_A : Numpy array
        First acceptor submovie as Numpy array.
    file_path_A : List of str
        List of paths of all acceptor submovies
    matrix_align : list of arrays
        List of the optimized transformation matrix to use for the chromatic aberrations correction.
    sigma : int, optional
        radius in pixels of the area around a spot where to look for the intensity of the spot (so any slight drift or innacurate detection of the spot is handled).\
            Needs to be less than the minimal interspots distance. Default is 2.
    DA_is : str, optional
        'odd' if DA frames are the odd ones, anything else if there are the even ones

    Returns
    -------
    traces_data : dict
        Dictionnary containing all the traces for individual filtered spots, including DD, DA and AA traces; the keys of this dictionnary correspond to the individual traces ID\
        Each entry also contains the xy coordinates of the related spot. It doesn't include the frames for each timepoints, as we consider there are no skipped frames in the movies

    """
    
    print('Start of the traces extraction and concatenation step...')
    
    
    Open_dataset = Dataset(
    title="My FRET Experiment",
    traces=[],
    description="FRET data of protein folding",
    experiment_type="2-Color FRET",
    authors=["John Doe", "Jane Smith"],
    institution="University X",
    date=date(2026, 6, 22),
    metadata=Metadata({"experiment_id": "20240101_JD_JS_1", "movie_file": "20240101_CoolExperiment.TIF", "ALEX": 'yes'}),
    sample_details={"buffer_conditions": "Phosphate buffer", "other_details": Metadata({"ph": 7.4})}, #Example of nested metadata
    instrument_details={"microscope": "Olympus IX83", "other_details": Metadata({"objective": "60x oil 1.5 NA"})}, #Example of nested metadata
)
    
    frame_sum = 0
    
    coord_spots_corr_DD = np.zeros(coord_spots.shape)
    
    # define if DA channel is odd or even frames, default is odd, anything else is considered even
    if DA_is == 'odd':
        DA_start = 1
        AA_start = 0
    else:
        DA_start = 0
        AA_start = 1
    
    # after extraction of the first submovie, loop over the next submovies until all have been processed, concatenate the subtraces to their corresponding traces ID inside the dict\
        # really similar to the previous section, with addition of a chromatic correction of the donor submovie
    for j in range(len(file_path_D)):
        print('Extracting traces from movie ' + str(j+1) + '...')
        
        # correct chromatic aberrations donor submovie
        
        if method_align == 'Optical Flow':
            img_stack_D = Warp_OpticalFlow(load_submovies(file_path_D[j]), matrix_align[0], matrix_align[1])
        else:
            img_stack_D = generate_chrom_ab_corr_movie(load_submovies(file_path_D[j]), matrix_align)
            
        img_stack_A = load_submovies(file_path_A[j])
        
        img_stack_D_proj = np.sum(img_stack_D, axis = 2)
        
        if j == 0:
            drift_j = [0, 0]
        else:
            drift_j = [drift_correct['y'].values[j-1], drift_correct['x'].values[j-1]]
        frame_sum = frame_sum + img_stack_D.shape[2]
        
        # loop of detected spots
        for k in range(len(coord_spots)):
            #search area corners management
            coord_spot_A = np.round(coord_spots[k,:] + drift_j).astype(int)
            xl = int(np.max([coord_spot_A[0]-sigma, 0]))
            xr = int(np.min([coord_spot_A[0]+ sigma, img_stack_A.shape[0]]))
            yu = int(np.max([coord_spot_A[1]-sigma, 0]))
            yd = int(np.min([coord_spot_A[1]+ sigma, img_stack_A.shape[1]]))
            
            xl_back = int(np.max([coord_spot_A[0]-2*sigma, 0]))
            xr_back = int(np.min([coord_spot_A[0]+ 2*sigma, img_stack_D.shape[0]]))
            yu_back = int(np.max([coord_spot_A[1]-2*sigma, 0]))
            yd_back = int(np.min([coord_spot_A[1]+ 2*sigma, img_stack_D.shape[1]]))
            
            xl_back_D_det = int(np.max([coord_spot_A[0]-20, 0]))
            xr_back_D_det = int(np.min([coord_spot_A[0]+ 20, img_stack_D.shape[0]]))
            yu_back_D_det = int(np.max([coord_spot_A[1]-20, 0]))
            yd_back_D_det = int(np.min([coord_spot_A[1]+ 20, img_stack_D.shape[1]]))
            
            if j == 0:
                img_spot = img_stack_D_proj[xl_back_D_det:xr_back_D_det,yu_back_D_det:yd_back_D_det]
                try:
                    coord_detec_D = detection.detect_spots(img_spot, log_kernel_size=2, minimum_distance=2)
                    coord_spots_corr_DD[k,:] = coord_detec_D[np.argmin(np.sum((coord_spot_A-np.array([xl_back_D_det, yu_back_D_det])-coord_detec_D)**2, axis = 1))] + [xl_back_D_det, yu_back_D_det]
                except ValueError:
                    coord_spots_corr_DD[k,:] = coord_spot_A
            
            coord_spot_D = np.round(coord_spots_corr_DD[k,:] + drift_j).astype(int)
            
            xl_D = int(np.max([coord_spot_D[0]-sigma, 0]))
            xr_D = int(np.min([coord_spot_D[0]+ sigma, img_stack_D.shape[0]]))
            yu_D = int(np.max([coord_spot_D[1]-sigma, 0]))
            yd_D = int(np.min([coord_spot_D[1]+ sigma, img_stack_D.shape[1]]))
            
            xl_back_D = int(np.max([coord_spot_D[0]-2*sigma, 0]))
            xr_back_D = int(np.min([coord_spot_D[0]+ 2*sigma, img_stack_D.shape[0]]))
            yu_back_D = int(np.max([coord_spot_D[1]-2*sigma, 0]))
            yd_back_D = int(np.min([coord_spot_D[1]+ 2*sigma, img_stack_D.shape[1]]))
            
            # extraction of subtraces
            mask_spot_signal = generate_mask_background(coord_spot_A, xl, xr, yu, yd, img_stack_D.shape[2], mask_radius = sigma)
            area_mask =  np.sum(np.sum(mask_spot_signal, axis = 0), axis = 0)
            
            mask_spot_signal_D = generate_mask_background(coord_spot_D, xl_D, xr_D, yu_D, yd_D, img_stack_D.shape[2], mask_radius = sigma)
            area_mask_D =  np.sum(np.sum(mask_spot_signal_D, axis = 0), axis = 0)
            trace_DD = np.sum(np.sum(img_stack_D[xl_D:xr_D,yu_D:yd_D,:] * mask_spot_signal_D, axis = 0),axis=0)  # example
            trace_DA = np.sum(np.sum(img_stack_A[xl:xr,yu:yd,[i for i in range(DA_start,img_stack_A.shape[2],2)]] * mask_spot_signal, axis = 0),axis=0)
            trace_AA = np.sum(np.sum(img_stack_A[xl:xr,yu:yd,[i for i in range(AA_start,img_stack_A.shape[2],2)]] * mask_spot_signal, axis = 0),axis=0)
            
            # crop the search area around the spot if it's at the corner of the image; xl, xr, yl and yr define the search aera corners coordinates
            
            # measure the background noise level around the spot using the mean inside the search area with masked spot
            mask_spot_back = generate_mask_background(coord_spot_A, xl_back, xr_back, yu_back, yd_back, img_stack_D.shape[2], mask_radius = 2*sigma)
            mask_spot_back_center_spot = generate_mask_background(coord_spot_A, xl_back, xr_back, yu_back, yd_back, img_stack_D.shape[2], mask_radius = sigma)
            mask_spot_back = mask_spot_back ^ mask_spot_back_center_spot # remove the spot signal part from the background mask, for more accurate background estimation
            area_mask_back =  np.sum(np.sum(mask_spot_back, axis = 0), axis = 0)
            
            
            mask_spot_back_D = generate_mask_background(coord_spot_D, xl_back_D, xr_back_D, yu_back_D, yd_back_D, img_stack_D.shape[2], mask_radius = 2*sigma)
            mask_spot_back_center_spot_D = generate_mask_background(coord_spot_D, xl_back_D, xr_back_D, yu_back_D, yd_back_D, img_stack_D.shape[2], mask_radius = sigma)
            mask_spot_back_D = mask_spot_back_D ^ mask_spot_back_center_spot_D # remove the spot signal part from the background mask, for more accurate background estimation
            area_mask_back_D =  np.sum(np.sum(mask_spot_back_D, axis = 0), axis = 0)
            dim2_back_mask = img_stack_D[xl_back:xr_back,yu_back:yd_back,:].shape[2]
            
            img_stack_D_spot = img_stack_D[xl_back_D:xr_back_D,yu_back_D:yd_back_D,:]
            img_stack_DA_spot = img_stack_A[xl_back:xr_back,yu_back:yd_back,[i for i in range(DA_start,img_stack_A.shape[2],2)]]
            img_stack_AA_spot = img_stack_A[xl_back:xr_back,yu_back:yd_back,[i for i in range(AA_start,img_stack_A.shape[2],2)]]
            
            back_trace_DD = np.median(img_stack_D_spot[mask_spot_back_D].reshape((area_mask_back_D[0], dim2_back_mask)), axis = 0)
            back_trace_DA = np.median(img_stack_DA_spot[mask_spot_back].reshape((area_mask_back[0], dim2_back_mask)), axis = 0)
            back_trace_AA = np.median(img_stack_AA_spot[mask_spot_back].reshape((area_mask_back[0], dim2_back_mask)), axis = 0)
            
            
            if j == 0:
                channel_DD = Channel("DD", list(np.round(trace_DD,0).astype(np.int64)))
                channel_DA = Channel("DA", list(np.round(trace_DA,0).astype(np.int64)))
                channel_AA = Channel("AA", list(np.round(trace_AA,0).astype(np.int64)))
                
                channel_back_DD = Channel("back_DD", list(np.round(back_trace_DD * area_mask_D,0).astype(np.int64)))
                channel_back_DA = Channel("back_DA", list(np.round(back_trace_DA * area_mask,0).astype(np.int64)))
                channel_back_AA = Channel("back_AA", list(np.round(back_trace_AA * area_mask,0).astype(np.int64)))
                
                trace1_i = Trace([channel_DD, channel_DA, channel_AA, channel_back_DD, channel_back_DA, channel_back_AA],
                                 metadata=Metadata({"molecule_id": str(k), 'x_coord': coord_spots[k,0], 'y_coord': coord_spots[k,1], "UUID_v7": str(uuid7())}))
                
                Open_dataset.traces.append(trace1_i)
                
            else:
            # concatenation of subtraces to their cooresponding trace ID
                
                Open_dataset.traces[k].channels[0].data = Open_dataset.traces[k].channels[0].data + list(np.round(trace_DD,0).astype(np.int64))
                Open_dataset.traces[k].channels[1].data = Open_dataset.traces[k].channels[1].data + list(np.round(trace_DA,0).astype(np.int64))
                Open_dataset.traces[k].channels[2].data = Open_dataset.traces[k].channels[2].data + list(np.round(trace_AA,0).astype(np.int64))
                Open_dataset.traces[k].channels[3].data = Open_dataset.traces[k].channels[3].data + list(np.round(back_trace_DD * area_mask_D,0).astype(np.int64))
                Open_dataset.traces[k].channels[4].data = Open_dataset.traces[k].channels[4].data + list(np.round(back_trace_DA * area_mask,0).astype(np.int64))
                Open_dataset.traces[k].channels[5].data = Open_dataset.traces[k].channels[5].data + list(np.round(back_trace_AA * area_mask,0).astype(np.int64))
                
            
    print('Traces extraction and concatenation completed!')
    return Open_dataset

def spot_trace_extractor_no_ALEX(file_path_D, file_path_A, coord_spots, matrix_align, 
                         drift_correct = 0, 
                         method_align = 'Optical Flow',
                         sigma = 3,
                         back_corr = 'median_TV_min',
                         TV_lam = 5):
    """
    Function used to extract and concatenate the subtraces from DD, DA and AA channels from all submovies

    Parameters
    ----------
    coord_spots : Numpy array
        List of coordinates of the distance-filtered detected spots.
    img_stack_D : Numpy array
        First donor submovie as Numpy array.
    file_path_D : List of str
        List of paths of all donor submovies
    img_stack_A : Numpy array
        First acceptor submovie as Numpy array.
    file_path_A : List of str
        List of paths of all acceptor submovies
    matrix_align : list of arrays
        List of the optimized transformation matrix to use for the chromatic aberrations correction.
    sigma : int, optional
        radius in pixels of the area around a spot where to look for the intensity of the spot (so any slight drift or innacurate detection of the spot is handled).\
            Needs to be less than the minimal interspots distance. Default is 2.
    DA_is : str, optional
        'odd' if DA frames are the odd ones, anything else if there are the even ones

    Returns
    -------
    traces_data : dict
        Dictionnary containing all the traces for individual filtered spots, including DD, DA and AA traces; the keys of this dictionnary correspond to the individual traces ID\
        Each entry also contains the xy coordinates of the related spot. It doesn't include the frames for each timepoints, as we consider there are no skipped frames in the movies

    """
    
    print('Start of the traces extraction and concatenation step...')
    
    
    Open_dataset = Dataset(
    title="My FRET Experiment",
    traces=[],
    description="FRET data of protein folding",
    experiment_type="2-Color FRET",
    authors=["John Doe", "Jane Smith"],
    institution="University X",
    date=date(2026, 6, 22),
    metadata=Metadata({"experiment_id": "20240101_JD_JS_1", "movie_file": "20240101_CoolExperiment.TIF", "ALEX": 'no'}),
    sample_details={"buffer_conditions": "Phosphate buffer", "other_details": Metadata({"ph": 7.4})}, #Example of nested metadata
    instrument_details={"microscope": "Olympus IX83", "other_details": Metadata({"objective": "60x oil 1.5 NA"})}, #Example of nested metadata
)
    
    frame_sum = 0
    
    coord_spots_corr_DD = np.zeros(coord_spots.shape)
    
    # define if DA channel is odd or even frames, default is odd, anything else is considered even

    
    # after extraction of the first submovie, loop over the next submovies until all have been processed, concatenate the subtraces to their corresponding traces ID inside the dict\
        # really similar to the previous section, with addition of a chromatic correction of the donor submovie
    for j in range(len(file_path_D)):
        print('Extracting traces from movie ' + str(j+1) + '...')
        
        # correct chromatic aberrations donor submovie
        
        if method_align == 'Optical Flow':
            img_stack_D = Warp_OpticalFlow(load_submovies(file_path_D[j]), matrix_align[0], matrix_align[1])
        else:
            img_stack_D = generate_chrom_ab_corr_movie(load_submovies(file_path_D[j]), matrix_align)
            
        img_stack_A = load_submovies(file_path_A[j])
        
        img_stack_D_proj = np.sum(img_stack_D, axis = 2)
        
        if j == 0:
            drift_j = [0, 0]
        else:
            drift_j = [drift_correct['y'].values[j-1], drift_correct['x'].values[j-1]]
        frame_sum = frame_sum + img_stack_D.shape[2]
        
        # loop of detected spots
        for k in range(len(coord_spots)):
            #search area corners management
            coord_spot_A = np.round(coord_spots[k,:] + drift_j).astype(int)
            xl = int(np.max([coord_spot_A[0]-sigma, 0]))
            xr = int(np.min([coord_spot_A[0]+ sigma, img_stack_A.shape[0]]))
            yu = int(np.max([coord_spot_A[1]-sigma, 0]))
            yd = int(np.min([coord_spot_A[1]+ sigma, img_stack_A.shape[1]]))
            
            xl_back = int(np.max([coord_spot_A[0]-2*sigma, 0]))
            xr_back = int(np.min([coord_spot_A[0]+ 2*sigma, img_stack_D.shape[0]]))
            yu_back = int(np.max([coord_spot_A[1]-2*sigma, 0]))
            yd_back = int(np.min([coord_spot_A[1]+ 2*sigma, img_stack_D.shape[1]]))
            
            xl_back_D_det = int(np.max([coord_spot_A[0]-20, 0]))
            xr_back_D_det = int(np.min([coord_spot_A[0]+ 20, img_stack_D.shape[0]]))
            yu_back_D_det = int(np.max([coord_spot_A[1]-20, 0]))
            yd_back_D_det = int(np.min([coord_spot_A[1]+ 20, img_stack_D.shape[1]]))
            
            if j == 0:
                img_spot = img_stack_D_proj[xl_back_D_det:xr_back_D_det,yu_back_D_det:yd_back_D_det]
                try:
                    coord_detec_D = detection.detect_spots(img_spot, log_kernel_size=2, minimum_distance=2)
                    coord_spots_corr_DD[k,:] = coord_detec_D[np.argmin(np.sum((coord_spot_A-np.array([xl_back_D_det, yu_back_D_det])-coord_detec_D)**2, axis = 1))] + [xl_back_D_det, yu_back_D_det]
                except ValueError:
                    coord_spots_corr_DD[k,:] = coord_spot_A
            
            coord_spot_D = np.round(coord_spots_corr_DD[k,:] + drift_j).astype(int)
            
            xl_D = int(np.max([coord_spot_D[0]-sigma, 0]))
            xr_D = int(np.min([coord_spot_D[0]+ sigma, img_stack_D.shape[0]]))
            yu_D = int(np.max([coord_spot_D[1]-sigma, 0]))
            yd_D = int(np.min([coord_spot_D[1]+ sigma, img_stack_D.shape[1]]))
            
            xl_back_D = int(np.max([coord_spot_D[0]-2*sigma, 0]))
            xr_back_D = int(np.min([coord_spot_D[0]+ 2*sigma, img_stack_D.shape[0]]))
            yu_back_D = int(np.max([coord_spot_D[1]-2*sigma, 0]))
            yd_back_D = int(np.min([coord_spot_D[1]+ 2*sigma, img_stack_D.shape[1]]))
            
            # extraction of subtraces
            mask_spot_signal = generate_mask_background(coord_spot_A, xl, xr, yu, yd, img_stack_D.shape[2], mask_radius = sigma)
            area_mask =  np.sum(np.sum(mask_spot_signal, axis = 0), axis = 0)
            
            mask_spot_signal_D = generate_mask_background(coord_spot_D, xl_D, xr_D, yu_D, yd_D, img_stack_D.shape[2], mask_radius = sigma)
            area_mask_D =  np.sum(np.sum(mask_spot_signal_D, axis = 0), axis = 0)
            trace_DD = np.sum(np.sum(img_stack_D[xl_D:xr_D,yu_D:yd_D,:] * mask_spot_signal_D, axis = 0),axis=0)  # example
            trace_DA = np.sum(np.sum(img_stack_A[xl:xr,yu:yd,:] * mask_spot_signal, axis = 0),axis=0)
            trace_AA = np.zeros(len(trace_DA))
            
            # crop the search area around the spot if it's at the corner of the image; xl, xr, yl and yr define the search aera corners coordinates
            
            # measure the background noise level around the spot using the mean inside the search area with masked spot
            mask_spot_back = generate_mask_background(coord_spot_A, xl_back, xr_back, yu_back, yd_back, img_stack_D.shape[2], mask_radius = 2*sigma)
            mask_spot_back_center_spot = generate_mask_background(coord_spot_A, xl_back, xr_back, yu_back, yd_back, img_stack_D.shape[2], mask_radius = sigma)
            mask_spot_back = mask_spot_back ^ mask_spot_back_center_spot # remove the spot signal part from the background mask, for more accurate background estimation
            area_mask_back =  np.sum(np.sum(mask_spot_back, axis = 0), axis = 0)
            
            
            mask_spot_back_D = generate_mask_background(coord_spot_D, xl_back_D, xr_back_D, yu_back_D, yd_back_D, img_stack_D.shape[2], mask_radius = 2*sigma)
            mask_spot_back_center_spot_D = generate_mask_background(coord_spot_D, xl_back_D, xr_back_D, yu_back_D, yd_back_D, img_stack_D.shape[2], mask_radius = sigma)
            mask_spot_back_D = mask_spot_back_D ^ mask_spot_back_center_spot_D # remove the spot signal part from the background mask, for more accurate background estimation
            area_mask_back_D =  np.sum(np.sum(mask_spot_back_D, axis = 0), axis = 0)
            dim2_back_mask = img_stack_D[xl_back:xr_back,yu_back:yd_back,:].shape[2]
            
            img_stack_D_spot = img_stack_D[xl_back_D:xr_back_D,yu_back_D:yd_back_D,:]
            img_stack_DA_spot = img_stack_A[xl_back:xr_back,yu_back:yd_back,:]
            
            back_trace_DD = np.median(img_stack_D_spot[mask_spot_back_D].reshape((area_mask_back_D[0], dim2_back_mask)), axis = 0)
            back_trace_DA = np.median(img_stack_DA_spot[mask_spot_back].reshape((area_mask_back[0], dim2_back_mask)), axis = 0)
            
            
            if j == 0:
                
                channel_DD = Channel("DD", list(np.round(trace_DD,0).astype(np.int64)))
                channel_DA = Channel("DA", list(np.round(trace_DA,0).astype(np.int64)))
                
                channel_back_DD = Channel("back_DD", list(np.round(back_trace_DD * area_mask_D,0).astype(np.int64)))
                channel_back_DD = Channel("back_DA", list(np.round(back_trace_DA * area_mask,0).astype(np.int64)))
                
                trace1_i = Trace([channel_DD, channel_DA, channel_back_DD, channel_back_DA],
                                 metadata=Metadata({"molecule_id": str(k), 'x_coord': coord_spots[k,0], 'y_coord': coord_spots[k,1], "UUID_v7": str(uuid7())}))
                
                Open_dataset.traces.append(trace1_i)
                
                
            else:
            # concatenation of subtraces to their cooresponding trace ID
                
                
                Open_dataset.traces[k].channels[0].data = Open_dataset.traces[k].channels[0].data + list(np.round(trace_DD,0).astype(np.int64))
                Open_dataset.traces[k].channels[1].data = Open_dataset.traces[k].channels[1].data + list(np.round(trace_DA,0).astype(np.int64))
                Open_dataset.traces[k].channels[2].data = Open_dataset.traces[k].channels[2].data + list(np.round(back_trace_DD * area_mask_D,0).astype(np.int64))
                Open_dataset.traces[k].channels[3].data = Open_dataset.traces[k].channels[3].data + list(np.round(back_trace_DA * area_mask,0).astype(np.int64))
                
            
    print('Traces extraction and concatenation completed!')
    return Open_dataset

def generate_mask_background(coord_spot, xl, xr, yu, yd, lengh_trace, mask_radius = 2):
    a, b = coord_spot[0] - xl, coord_spot[1] - yu
    nx = xr - xl
    ny = yd - yu

    y,x = np.ogrid[-a:nx-a, -b:ny-b]
    mask = x*x + y*y <= mask_radius*mask_radius
    
    mask_3D = np.repeat(mask[:, :, np.newaxis], lengh_trace, axis=2)
    
    return mask_3D

def detect_bleaching_event(traces_data, percentage_value = 0.05, KCPD_model = 'linear'):
    """
    Function used to estimate the bleach point of each traces for both the donor and acceptor

    Parameters
    ----------
    traces_data : dict
        Dictionnary containing all the traces results.
    percentage_value : float, optional
        Percentage of the trace duration taken to estimate the median beginning / end signal levels. The default is 0.05.
    KCPD_model : str, optional
        model used by the Kernel change point detection from Ruptures library: "linear", "rbf or "cosine". Default is "rbf" (Gqussian model)

    Returns
    -------
    traces_data : Dict
        Updated dictionnary results with the donor and acceptor bleach point.

    """
    print('Start of the bleach point detection step...')
    
    # loop over each individual traces
    for i in range(len(traces_data)):
        print(str(i + 1) + ' / ' + str(len(traces_data)) + ' processed')

        # calculate the number of frames corresponding to X % of the total trace duration, that will be used to estimate the median level at beginning and end of the trace
        one_perc = int(len(traces_data[str(i)]['Intensity_DD'])*percentage_value)
        traceDD_start_value = np.median(traces_data[str(i)]['Intensity_DD'][0:one_perc])
        traceDD_final_value = np.median(traces_data[str(i)]['Intensity_DD'][-one_perc:])
        
        # calculate the std at the end of the trace
        stdDD_trace = np.std(traces_data[str(i)]['Intensity_DD'][-one_perc:])
        
        # if the beginning median level is not above the std of the ending section, then we consider that the trace doesn't show significant DD signal and\
        # the bleaching event is set to 0
        # Otherwise, the bleaching event is predicted by a Kernel change point detection method with a Gaussian model (from ruptures library)
        if (traceDD_start_value - traceDD_final_value) > stdDD_trace:
            algo_c = rpt.KernelCPD(kernel=KCPD_model, min_size=2).fit(traces_data[str(i)]['Intensity_DD'])
            
            # we want only one break point, as we only expect one bleaching event
            result = algo_c.predict(n_bkps=1)
            
            # the bleach event is added to the dict of traces as a new sub-entry
            traces_data[str(i)]['bleaching_event_DD'] = result[0]
        else:
            # if no significant DD signal, bleach event is set to 0
            traces_data[str(i)]['bleaching_event_DD'] = 0
        
        # Similar process for the AA channel, so we can distinguish the donor and acceptor bleaching events
        traceAA_start_value = np.median(traces_data[str(i)]['Intensity_AA'][0:one_perc])
        traceAA_final_value = np.median(traces_data[str(i)]['Intensity_AA'][-one_perc:])
        stdAA_trace = np.std(traces_data[str(i)]['Intensity_AA'][-one_perc:])
        if (traceAA_start_value - traceAA_final_value) > stdAA_trace:
            algo_c = rpt.KernelCPD(kernel=KCPD_model, min_size=2).fit(traces_data[str(i)]['Intensity_AA'])
            result = algo_c.predict(n_bkps=1)
            traces_data[str(i)]['bleaching_event_AA'] = result[0]
        else:
            traces_data[str(i)]['bleaching_event_AA'] = 0
            
    print('Bleach point detection step completed!')
    
    # return the updated dictionnary with all traces results
    return traces_data
        

def get_SNR(traces_data):
    """
    Function used to calculate the SNR of each traces for all channels (DD, DA and AA)

    Parameters
    ----------
    traces_data : dict
        Dictionnary containing all the traces results.

    Returns
    -------
    traces_data : dict
        Dictionnary containing all the traces results, including the SNR for each traces (DD, DA, AA)

    """
    print('Start of the SNR calculation step...')
    
    # loop over all traces
    for i in range(len(traces_data)):
        
        # If the DD channel doesn't show a bleach event, then we set its SNR to 0; we assume it's just noise, but this could also be a donor with long lifetime\
        # Then this part will need further refinement to take into account these long lifetime donor (same remark for AA channel)
        if traces_data[str(i)]['bleaching_event_DD'] == 0:
            traces_data[str(i)]['SNR_DD'] = 0
        else:
            # We calculate the SNR as the mean signal before the bleaching event minus the mean signal after bleaching (noise) and divided by the std after the bleach event
            # Same calculation for DA and AA channels
            bleach_event = traces_data[str(i)]['bleaching_event_DD']
            trace_i = traces_data[str(i)]['Intensity_DD']
            traces_data[str(i)]['SNR_DD'] = (np.mean(trace_i[0:bleach_event]) - np.mean(trace_i[bleach_event:])) / np.std(trace_i[bleach_event:])
            
        # if no bleaching from the acceptor, then SNR = 0 for both DA and AA channels
        if traces_data[str(i)]['bleaching_event_AA'] == 0:
            traces_data[str(i)]['SNR_AA'] = 0
            traces_data[str(i)]['SNR_DA'] = 0
        else:
            # if bleaching of the acceptor, first calculate SNR for AA
            bleach_event = traces_data[str(i)]['bleaching_event_AA']
            trace_i = traces_data[str(i)]['Intensity_AA']
            traces_data[str(i)]['SNR_AA'] = (np.mean(trace_i[0:bleach_event]) - np.mean(trace_i[bleach_event:])) / np.std(trace_i[bleach_event:])
            
            # Then check if the donor is also bleached, if not DA SNR = 0
            if traces_data[str(i)]['bleaching_event_DD'] == 0:
                traces_data[str(i)]['SNR_DA'] = 0
            else:
                bleach_event = int(np.min((traces_data[str(i)]['bleaching_event_DD'], traces_data[str(i)]['bleaching_event_AA'])))
                trace_i = traces_data[str(i)]['Intensity_DA']
                traces_data[str(i)]['SNR_DA'] = (np.mean(trace_i[0:bleach_event]) - np.mean(trace_i[bleach_event:])) / np.std(trace_i[bleach_event:])
   
    print('SNR calculation step completed!')
    return traces_data

def get_high_signal(traces_data):
    for i in range(len(traces_data)):
        max_trace_i = np.max(traces_data[str(i)]['Intensity_DA'])
        traces_data[str(i)]['SNR_DA'] = max_trace_i
        traces_data[str(i)]['SNR_DD'] = max_trace_i
        traces_data[str(i)]['SNR_AA'] = max_trace_i
    return traces_data

def filter_traces_on_SNR(traces_data, SNR_thresh = 8):
    """
    Function used to filter the traces with high SNRs for both DD and AA channels

    Parameters
    ----------
    traces_data : dict
        Dictionnary containing all the traces results.
    SNR_thresh : float or int, optional
        Threshold level for high SNR. The default is 8.

    Returns
    -------
    traces_high_SNR : dict
        Dictionnary containing all the high SNR traces
    list_index_high_SNR : list of str
        List of traces index (str) with high SNRs for both DD and AA channels
    
    """
    print('Start of the high SNR traces filter step...')
    
    # new dictionnary that will contain all results from traces with high SNR for both DD and AA channels
    traces_high_SNR = {}
    
    traces_data_keys = traces_data.keys()
    
    # loop over all traces
    for i in traces_data_keys:
        # select the traces with high SNR in both DD and AA channels
        if i != 'chromatic_aberration_corr_matrix':
            if (traces_data[i]['SNR_AA'] >= SNR_thresh) and (traces_data[i]['SNR_DD'] >= SNR_thresh):
                traces_high_SNR[i] = traces_data[i]
    # get the list of traces index with high SNRs
    list_index_high_SNR = list(traces_high_SNR.keys())
    
    print('High SNR traces filter step completed!')
    
    return traces_high_SNR, list_index_high_SNR



def old_spot_trace_extractor_No_ALEX(file_path_D, file_path_A, coord_spots, matrix_align, list_frames_drift_AA, drift_correct = 0, sigma = 3):
    """
    Function used to extract and concatenate the subtraces from DD, DA and AA channels from all submovies

    Parameters
    ----------
    coord_spots : Numpy array
        List of coordinates of the distance-filtered detected spots.
    img_stack_D : Numpy array
        First donor submovie as Numpy array.
    file_path_D : List of str
        List of paths of all donor submovies
    img_stack_A : Numpy array
        First acceptor submovie as Numpy array.
    file_path_A : List of str
        List of paths of all acceptor submovies
    matrix_align : list of arrays
        List of the optimized transformation matrix to use for the chromatic aberrations correction.
    sigma : int, optional
        radius in pixels of the area around a spot where to look for the intensity of the spot (so any slight drift or innacurate detection of the spot is handled).\
            Needs to be less than the minimal interspots distance. Default is 2.
    DA_is : str, optional
        'odd' if DA frames are the odd ones, anything else if there are the even ones

    Returns
    -------
    traces_data : dict
        Dictionnary containing all the traces for individual filtered spots, including DD, DA and AA traces; the keys of this dictionnary correspond to the individual traces ID\
        Each entry also contains the xy coordinates of the related spot. It doesn't include the frames for each timepoints, as we consider there are no skipped frames in the movies

    """
    
    print('Start of the traces extraction and concatenation step...')
    
    # create a dict that will gather all traces and related analysis
    traces_data = {}
    
    frame_sum = 0
    
    coord_spots_corr_DD = np.zeros(coord_spots.shape)
    # after extraction of the first submovie, loop over the next submovies until all have been processed, concatenate the subtraces to their corresponding traces ID inside the dict\
        # really similar to the previous section, with addition of a chromatic correction of the donor submovie
    for j in range(len(file_path_D)):
        print('Extracting traces from movie ' + str(j+1) + '...')
        
        # correct chromatic aberrations donor submovie
        img_stack_D = generate_chrom_ab_corr_movie(load_submovies(file_path_D[j]), matrix_align)
        img_stack_A = load_submovies(file_path_A[j])
        
        img_stack_D_proj = np.sum(img_stack_D, axis = 2)
        
        if j == 0:
            drift_j = [0, 0]
        else:
            drift_j = [drift_correct['y'].values[j-1], drift_correct['x'].values[j-1]]
        frame_sum = frame_sum + img_stack_D.shape[2]
        
        # loop of detected spots
        for k in range(len(coord_spots)):
            #search area corners management
            coord_spot_A = np.round(coord_spots[k,:] + drift_j).astype(int)
            xl = int(np.max([coord_spot_A[0]-sigma, 0]))
            xr = int(np.min([coord_spot_A[0]+ sigma, img_stack_A.shape[0]]))
            yu = int(np.max([coord_spot_A[1]-sigma, 0]))
            yd = int(np.min([coord_spot_A[1]+ sigma, img_stack_A.shape[1]]))
            
            xl_back = int(np.max([coord_spot_A[0]-2*sigma, 0]))
            xr_back = int(np.min([coord_spot_A[0]+ 2*sigma, img_stack_D.shape[0]]))
            yu_back = int(np.max([coord_spot_A[1]-2*sigma, 0]))
            yd_back = int(np.min([coord_spot_A[1]+ 2*sigma, img_stack_D.shape[1]]))
            
            xl_back_D_det = int(np.max([coord_spot_A[0]-20, 0]))
            xr_back_D_det = int(np.min([coord_spot_A[0]+ 20, img_stack_D.shape[0]]))
            yu_back_D_det = int(np.max([coord_spot_A[1]-20, 0]))
            yd_back_D_det = int(np.min([coord_spot_A[1]+ 20, img_stack_D.shape[1]]))
            
            if j == 0:
                img_spot = img_stack_D_proj[xl_back_D_det:xr_back_D_det,yu_back_D_det:yd_back_D_det]
                try:
                    coord_detec_D = detection.detect_spots(img_spot, log_kernel_size=2, minimum_distance=2)
                    coord_spots_corr_DD[k,:] = coord_detec_D[np.argmin(np.sum((coord_spot_A-np.array([xl_back_D_det, yu_back_D_det])-coord_detec_D)**2, axis = 1))] + [xl_back_D_det, yu_back_D_det]
                except ValueError:
                    coord_spots_corr_DD[k,:] = coord_spot_A
            
            coord_spot_D =coord_spots_corr_DD[k,:] + drift_j
            
            xl_D = int(np.max([coord_spot_D[0]-sigma, 0]))
            xr_D = int(np.min([coord_spot_D[0]+ sigma, img_stack_D.shape[0]]))
            yu_D = int(np.max([coord_spot_D[1]-sigma, 0]))
            yd_D = int(np.min([coord_spot_D[1]+ sigma, img_stack_D.shape[1]]))
            
            xl_back_D = int(np.max([coord_spot_D[0]-2*sigma, 0]))
            xr_back_D = int(np.min([coord_spot_D[0]+ 2*sigma, img_stack_D.shape[0]]))
            yu_back_D = int(np.max([coord_spot_D[1]-2*sigma, 0]))
            yd_back_D = int(np.min([coord_spot_D[1]+ 2*sigma, img_stack_D.shape[1]]))
            
            # extraction of subtraces
            mask_spot_signal = generate_mask_background(coord_spot_A, xl, xr, yu, yd, img_stack_D.shape[2], mask_radius = sigma)
            area_mask =  np.sum(np.sum(mask_spot_signal, axis = 0), axis = 0)
            
            mask_spot_signal_D = generate_mask_background(coord_spot_D, xl_D, xr_D, yu_D, yd_D, img_stack_D.shape[2], mask_radius = sigma)
            area_mask_D =  np.sum(np.sum(mask_spot_signal_D, axis = 0), axis = 0)
            trace_DD = np.sum(np.sum(img_stack_D[xl_D:xr_D,yu_D:yd_D,:] * mask_spot_signal_D, axis = 0),axis=0)  # example
            trace_DA = np.sum(np.sum(img_stack_A[xl:xr,yu:yd,:] * mask_spot_signal, axis = 0),axis=0)
            
            # crop the search area around the spot if it's at the corner of the image; xl, xr, yl and yr define the search aera corners coordinates
            
            # measure the background noise level around the spot using the mean inside the search area with masked spot
            mask_spot_back = generate_mask_background(coord_spot_A, xl_back, xr_back, yu_back, yd_back, img_stack_D.shape[2], mask_radius = 2*sigma)
            mask_spot_back_center_spot = generate_mask_background(coord_spot_A, xl_back, xr_back, yu_back, yd_back, img_stack_D.shape[2], mask_radius = sigma)
            mask_spot_back = mask_spot_back ^ mask_spot_back_center_spot # remove the spot signal part from the background mask, for more accurate background estimation
            area_mask_back =  np.sum(np.sum(mask_spot_back, axis = 0), axis = 0)
            
            
            mask_spot_back_D = generate_mask_background(coord_spot_D, xl_back_D, xr_back_D, yu_back_D, yd_back_D, img_stack_D.shape[2], mask_radius = 2*sigma)
            mask_spot_back_center_spot_D = generate_mask_background(coord_spot_D, xl_back_D, xr_back_D, yu_back_D, yd_back_D, img_stack_D.shape[2], mask_radius = sigma)
            mask_spot_back_D = mask_spot_back_D ^ mask_spot_back_center_spot_D # remove the spot signal part from the background mask, for more accurate background estimation
            area_mask_back_D =  np.sum(np.sum(mask_spot_back_D, axis = 0), axis = 0)

            dim2_back_mask = img_stack_D[xl_back:xr_back,yu_back:yd_back,:].shape[2]
            
            img_stack_D_spot = img_stack_D[xl_back_D:xr_back_D,yu_back_D:yd_back_D,:]
            img_stack_DA_spot = img_stack_A[xl_back:xr_back,yu_back:yd_back,:]
            
            back_trace_DD = np.median(img_stack_D_spot[mask_spot_back_D].reshape((area_mask_back_D[0], dim2_back_mask)), axis = 0)
            back_trace_DA = np.median(img_stack_DA_spot[mask_spot_back].reshape((area_mask_back[0], dim2_back_mask)), axis = 0)
            
            # remove the background noise level from the signals
            corr_trace_DD = trace_DD - back_trace_DD * area_mask_D #/ area_mask_back
            corr_trace_DA = trace_DA - back_trace_DA * area_mask #/ area_mask_back
            
            if j == 0:
                traces_data[str(k)] = {'x_coord': coord_spots[k,0], 'y_coord': coord_spots[k,1]}
                traces_data[str(k)]['Intensity_DD'] = np.round(corr_trace_DD,0).astype(np.int64)
                traces_data[str(k)]['Intensity_DA'] = np.round(corr_trace_DA,0).astype(np.int64)
                traces_data[str(k)]['Intensity_AA'] = np.zeros(corr_trace_DA.shape)
            else:
            # concatenation of subtraces to their cooresponding trace ID
                traces_data[str(k)]['Intensity_DD'] = np.concat((traces_data[str(k)]['Intensity_DD'], np.round(corr_trace_DD,0).astype(np.int64)))
                traces_data[str(k)]['Intensity_DA'] = np.concat((traces_data[str(k)]['Intensity_DA'], np.round(corr_trace_DA,0).astype(np.int64)))
                traces_data[str(k)]['Intensity_AA'] = np.concat((traces_data[str(k)]['Intensity_AA'], np.zeros(corr_trace_DA.shape)))
            
    print('Traces extraction and concatenation completed!')
    
    return traces_data

def shift_peak(img, shift_YX):
    M = np.float32([
        [1.0, 0.0, -shift_YX[1]],  # x' = x + 10.5
        [0.0, 1.0, -shift_YX[0]]   # y' = y + 20.3
        ])
    H, W = img.shape[0:2]
    
    img_i = np.zeros(img.shape)
    
    length_img = img.shape[2]
    
    for m in range(1 + length_img // 512):
        img_i[:,:,m*512:np.min((512*(m+1), length_img))] = cv2.warpAffine(img[:,:,m*512:np.min((512*(m+1), length_img))], M, (W, H)) #, flags=cv2.INTER_LINEAR)
    
    return img_i

def spot_trace_extractor_syn_movie(movie_syn, coord_spots, drift_correct, sigma = 3):
    """
    Function used to extract and concatenate the subtraces from DD, DA and AA channels from all submovies

    Parameters
    ----------
    coord_spots : Numpy array
        List of coordinates of the distance-filtered detected spots.
    img_stack_D : Numpy array
        First donor submovie as Numpy array.
    file_path_D : List of str
        List of paths of all donor submovies
    img_stack_A : Numpy array
        First acceptor submovie as Numpy array.
    file_path_A : List of str
        List of paths of all acceptor submovies
    matrix_align : list of arrays
        List of the optimized transformation matrix to use for the chromatic aberrations correction.
    sigma : int, optional
        radius in pixels of the area around a spot where to look for the intensity of the spot (so any slight drift or innacurate detection of the spot is handled).\
            Needs to be less than the minimal interspots distance. Default is 2.
    DA_is : str, optional
        'odd' if DA frames are the odd ones, anything else if there are the even ones

    Returns
    -------
    traces_data : dict
        Dictionnary containing all the traces for individual filtered spots, including DD, DA and AA traces; the keys of this dictionnary correspond to the individual traces ID\
        Each entry also contains the xy coordinates of the related spot. It doesn't include the frames for each timepoints, as we consider there are no skipped frames in the movies

    """
    
    print('Start of the traces extraction and concatenation step...')
    
    # create a dict that will gather all traces and related analysis
    traces_data = {}
    
    # define if DA channel is odd or even frames, default is odd, anything else is considered even
    
    # after extraction of the first submovie, loop over the next submovies until all have been processed, concatenate the subtraces to their corresponding traces ID inside the dict\
        # really similar to the previous section, with addition of a chromatic correction of the donor submovie
    for j in range(movie_syn.shape[2]):
        print('Extracting traces from movie ' + str(j+1) + '...')
        if j == 0:
            drift_j = [0, 0]
        else:
            drift_j = [drift_correct['y'].values[j-1], drift_correct['x'].values[j-1]]
        
        # loop of detected spots
        for k in range(len(coord_spots)):
            #search area corners management
            coord_spot_A = np.round(coord_spots[k,:] + drift_j).astype(int)  # should be +
            if coord_spot_A[0] < 0:
                coord_spot_A[0] = 0
            if coord_spot_A[1] < 0:
                coord_spot_A[1] = 0    
            xl = int(np.max([coord_spot_A[0]-sigma, 0]))
            xr = int(np.min([coord_spot_A[0]+ sigma, movie_syn.shape[0]]))
            yu = int(np.max([coord_spot_A[1]-sigma, 0]))
            yd = int(np.min([coord_spot_A[1]+ sigma, movie_syn.shape[1]]))
            
            xl_back = int(np.max([coord_spot_A[0]-2*sigma, 0]))
            xr_back = int(np.min([coord_spot_A[0]+ 2*sigma, movie_syn.shape[0]]))
            yu_back = int(np.max([coord_spot_A[1]-2*sigma, 0]))
            yd_back = int(np.min([coord_spot_A[1]+ 2*sigma, movie_syn.shape[1]]))
            
            
            # extraction of subtraces
            mask_spot_signal = generate_mask_background(coord_spot_A, xl, xr, yu, yd, 1, mask_radius = sigma)
            mask_spot_signal = mask_spot_signal[:,:,0]
            area_mask =  np.sum(np.sum(mask_spot_signal, axis = 0), axis = 0)
            trace_syn = np.sum(np.sum(movie_syn[xl:xr,yu:yd,j] * mask_spot_signal, axis = 0),axis=0)  # example
            
            # crop the search area around the spot if it's at the corner of the image; xl, xr, yl and yr define the search aera corners coordinates
            
            # measure the background noise level around the spot using the mean inside the search area with masked spot
            mask_spot_back = generate_mask_background(coord_spot_A, xl_back, xr_back, yu_back, yd_back, 1, mask_radius = 2*sigma)
            mask_spot_back = mask_spot_back[:,:,0]
            mask_spot_back_center_spot = generate_mask_background(coord_spot_A, xl_back, xr_back, yu_back, yd_back, 1, mask_radius = sigma)
            mask_spot_back_center_spot = mask_spot_back_center_spot[:,:,0]
            mask_spot_back = mask_spot_back ^ mask_spot_back_center_spot # remove the spot signal part from the background mask, for more accurate background estimation
            area_mask_back =  np.sum(np.sum(mask_spot_back, axis = 0), axis = 0)
            
            
            dim2_back_mask = 1
            
            img_stack_syn_spot = movie_syn[xl_back:xr_back,yu_back:yd_back,j]
            
            back_trace_syn = np.median(img_stack_syn_spot[mask_spot_back])
            
            # remove the background noise level from the signals
            corr_trace_syn = [trace_syn - back_trace_syn * area_mask] #/ area_mask_back
            
            if j == 0:
                traces_data[str(k)] = {'x_coord': coord_spots[k,0], 'y_coord': coord_spots[k,1]}
                traces_data[str(k)]['Intensity_DD'] = np.round(corr_trace_syn,0).astype(np.int64)
                traces_data[str(k)]['Intensity_DA'] = np.round(corr_trace_syn,0).astype(np.int64)
                traces_data[str(k)]['Intensity_AA'] = np.round(corr_trace_syn,0).astype(np.int64)
            else:
            # concatenation of subtraces to their cooresponding trace ID
                traces_data[str(k)]['Intensity_DD'] = np.concat((traces_data[str(k)]['Intensity_DD'], np.round(corr_trace_syn,0).astype(np.int64)))
                traces_data[str(k)]['Intensity_DA'] = np.concat((traces_data[str(k)]['Intensity_DA'], np.round(corr_trace_syn,0).astype(np.int64)))
                traces_data[str(k)]['Intensity_AA'] = np.concat((traces_data[str(k)]['Intensity_AA'], np.round(corr_trace_syn,0).astype(np.int64)))
            
    print('Traces extraction and concatenation completed!')
    
    return traces_data

def extract_traces_from_DD_spots(file_path_D, file_path_A, coord_spots, matrix_align,
                                 drift_correct = 0,
                                 sigma = 3,
                                 DA_is = 'odd',
                                 method_align = 'Optical Flow'):
    
    
    Open_dataset = Dataset(
    title="My FRET Experiment",
    traces=[],
    description="FRET data of protein folding",
    experiment_type="2-Color FRET",
    authors=["John Doe", "Jane Smith"],
    institution="University X",
    date=date(2026, 6, 22),
    metadata=Metadata({"experiment_id": "20240101_JD_JS_1", "movie_file": "20240101_CoolExperiment.TIF", "ALEX": 'yes'}),
    sample_details={"buffer_conditions": "Phosphate buffer", "other_details": Metadata({"ph": 7.4})}, #Example of nested metadata
    instrument_details={"microscope": "Olympus IX83", "other_details": Metadata({"objective": "60x oil 1.5 NA"})}, #Example of nested metadata
)
    
    frame_sum = 0
    
    if method_align == 'Optical Flow':
        shapeY_img, shapeX_img = matrix_align[0].shape
        for i in range(len(coord_spots)):
            coord_spots[i,0] = coord_spots[i,0] + matrix_align[0][np.min((int(coord_spots[i,0]),shapeY_img-1)),np.min((int(coord_spots[i,1]),shapeX_img-1))]
            coord_spots[i,1] = coord_spots[i,1] + matrix_align[1][np.min((int(coord_spots[i,0]),shapeY_img-1)),np.min((int(coord_spots[i,1]),shapeX_img-1))]
        
    else:
        coord_3D = np.stack((coord_spots[:,1], coord_spots[:,0], np.ones(len(coord_spots))))
        
        img_temp = Image.open(file_path_D[0])
        
        H, W = img_temp.height, img_temp.width
        
        omega = generate_trans_matrix(matrix_align, H, W)
        
        coord_3D_corr = omega @ coord_3D
        
        coord_spots = np.array([coord_3D_corr[1,:], coord_3D_corr[0,:]]).T
        
        mask = ~np.any(coord_spots <= 0, axis=1)
        coord_spots = coord_spots[mask]
    
    coord_spots_corr_A = np.zeros(coord_spots.shape)
    
    # define if DA channel is odd or even frames, default is odd, anything else is considered even
    if DA_is == 'odd':
        DA_start = 1
        AA_start = 0
    else:
        DA_start = 0
        AA_start = 1
    
    # after extraction of the first submovie, loop over the next submovies until all have been processed, concatenate the subtraces to their corresponding traces ID inside the dict\
        # really similar to the previous section, with addition of a chromatic correction of the donor submovie
    for j in range(len(file_path_D)):
        print('Extracting traces from movie ' + str(j+1) + '...')
        
        
        if method_align == 'Optical Flow':
            img_stack_D = Warp_OpticalFlow(load_submovies(file_path_D[j]), matrix_align[0], matrix_align[1])
        else:
            img_stack_D = generate_chrom_ab_corr_movie(load_submovies(file_path_D[j]), matrix_align)
        
        # correct chromatic aberrations donor submovie
        img_stack_A = load_submovies(file_path_A[j])
        
        img_stack_A_proj = np.sum(img_stack_A[:,:,[i for i in range(DA_start,img_stack_A.shape[2],2)]], axis = 2)
                               
        if j == 0:
            drift_j = [0, 0]
        else:
            drift_j = [drift_correct['y'].values[j-1], drift_correct['x'].values[j-1]]
        frame_sum = frame_sum + img_stack_D.shape[2]
        
        # loop of detected spots
        for k in range(len(coord_spots)):
            #search area corners management
            coord_spot_D = np.round(coord_spots[k,:] + drift_j).astype(int)
            xl = int(np.max([coord_spot_D[0]-sigma, 0]))
            xr = int(np.min([coord_spot_D[0]+ sigma, img_stack_D.shape[0]]))
            yu = int(np.max([coord_spot_D[1]-sigma, 0]))
            yd = int(np.min([coord_spot_D[1]+ sigma, img_stack_D.shape[1]]))
            
            xl_back = int(np.max([coord_spot_D[0]-2*sigma, 0]))
            xr_back = int(np.min([coord_spot_D[0]+ 2*sigma, img_stack_D.shape[0]]))
            yu_back = int(np.max([coord_spot_D[1]-2*sigma, 0]))
            yd_back = int(np.min([coord_spot_D[1]+ 2*sigma, img_stack_D.shape[1]]))
            
            xl_back_A_det = int(np.max([coord_spot_D[0]-20, 0]))
            xr_back_A_det = int(np.min([coord_spot_D[0]+ 20, img_stack_D.shape[0]]))
            yu_back_A_det = int(np.max([coord_spot_D[1]-20, 0]))
            yd_back_A_det = int(np.min([coord_spot_D[1]+ 20, img_stack_D.shape[1]]))
            
            if j == 0:
                img_spot = img_stack_A_proj[xl_back_A_det:xr_back_A_det,yu_back_A_det:yd_back_A_det]
                try:
                    coord_detec_A = detection.detect_spots(img_spot, log_kernel_size=2, minimum_distance=2)
                    coord_spots_corr_A[k,:] = coord_detec_A[np.argmin(np.sum((coord_spot_D-np.array([xl_back_A_det, yu_back_A_det])-coord_detec_A)**2, axis = 1))] + [xl_back_A_det, yu_back_A_det]
                except ValueError:
                    coord_spots_corr_A[k,:] = coord_spot_D
            
            coord_spot_A =coord_spots_corr_A[k,:] + drift_j
            
            xl_A = int(np.max([coord_spot_A[0]-sigma, 0]))
            xr_A = int(np.min([coord_spot_A[0]+ sigma, img_stack_A.shape[0]]))
            yu_A = int(np.max([coord_spot_A[1]-sigma, 0]))
            yd_A = int(np.min([coord_spot_A[1]+ sigma, img_stack_A.shape[1]]))
            
            xl_back_A = int(np.max([coord_spot_A[0]-2*sigma, 0]))
            xr_back_A = int(np.min([coord_spot_A[0]+ 2*sigma, img_stack_A.shape[0]]))
            yu_back_A = int(np.max([coord_spot_A[1]-2*sigma, 0]))
            yd_back_A = int(np.min([coord_spot_A[1]+ 2*sigma, img_stack_A.shape[1]]))
            
            # extraction of subtraces
            mask_spot_signal = generate_mask_background(coord_spot_D, xl, xr, yu, yd, img_stack_D.shape[2], mask_radius = sigma)
            area_mask =  np.sum(np.sum(mask_spot_signal, axis = 0), axis = 0)
            
            mask_spot_signal_A = generate_mask_background(coord_spot_A, xl_A, xr_A, yu_A, yd_A, img_stack_D.shape[2], mask_radius = sigma)
            area_mask_A =  np.sum(np.sum(mask_spot_signal_A, axis = 0), axis = 0)
            trace_DD = np.sum(np.sum(img_stack_D[xl:xr,yu:yd,:] * mask_spot_signal, axis = 0),axis=0)  # example
            trace_DA = np.sum(np.sum(img_stack_A[xl_A:xr_A,yu_A:yd_A,[i for i in range(DA_start,img_stack_A.shape[2],2)]] * mask_spot_signal_A, axis = 0),axis=0)
            trace_AA = np.sum(np.sum(img_stack_A[xl_A:xr_A,yu_A:yd_A,[i for i in range(AA_start,img_stack_A.shape[2],2)]] * mask_spot_signal_A, axis = 0),axis=0)
            
            # crop the search area around the spot if it's at the corner of the image; xl, xr, yl and yr define the search aera corners coordinates
            
            # measure the background noise level around the spot using the mean inside the search area with masked spot
            mask_spot_back = generate_mask_background(coord_spot_D, xl_back, xr_back, yu_back, yd_back, img_stack_D.shape[2], mask_radius = 2*sigma)
            mask_spot_back_center_spot = generate_mask_background(coord_spot_D, xl_back, xr_back, yu_back, yd_back, img_stack_D.shape[2], mask_radius = sigma)
            mask_spot_back = mask_spot_back ^ mask_spot_back_center_spot # remove the spot signal part from the background mask, for more accurate background estimation
            area_mask_back =  np.sum(np.sum(mask_spot_back, axis = 0), axis = 0)
            
            
            mask_spot_back_A = generate_mask_background(coord_spot_A, xl_back_A, xr_back_A, yu_back_A, yd_back_A, img_stack_D.shape[2], mask_radius = 2*sigma)
            mask_spot_back_center_spot_A = generate_mask_background(coord_spot_A, xl_back_A, xr_back_A, yu_back_A, yd_back_A, img_stack_D.shape[2], mask_radius = sigma)
            mask_spot_back_A = mask_spot_back_A ^ mask_spot_back_center_spot_A # remove the spot signal part from the background mask, for more accurate background estimation
            area_mask_back_A =  np.sum(np.sum(mask_spot_back_A, axis = 0), axis = 0)

            dim2_back_mask = img_stack_D[xl_back:xr_back,yu_back:yd_back,:].shape[2]
            
            img_stack_D_spot = img_stack_D[xl_back:xr_back,yu_back:yd_back,:]
            img_stack_DA_spot = img_stack_A[xl_back_A:xr_back_A,yu_back_A:yd_back_A,[i for i in range(DA_start,img_stack_A.shape[2],2)]]
            img_stack_AA_spot = img_stack_A[xl_back_A:xr_back_A,yu_back_A:yd_back_A,[i for i in range(AA_start,img_stack_A.shape[2],2)]]
            
            back_trace_DD = np.median(img_stack_D_spot[mask_spot_back].reshape((area_mask_back[0], dim2_back_mask)), axis = 0)
            back_trace_DA = np.median(img_stack_DA_spot[mask_spot_back_A].reshape((area_mask_back_A[0], dim2_back_mask)), axis = 0)
            back_trace_AA = np.median(img_stack_AA_spot[mask_spot_back_A].reshape((area_mask_back_A[0], dim2_back_mask)), axis = 0)
            
            if j == 0:

                channel_DD = Channel("DD", list(np.round(trace_DD,0).astype(np.int64)))

                channel_DA = Channel("DA", list(np.round(trace_DA,0).astype(np.int64)))

                channel_AA = Channel("AA", list(np.round(trace_AA,0).astype(np.int64)))
                

                channel_back_DD = Channel("back_DD", list(np.round(back_trace_DD * area_mask,0).astype(np.int64)))

                channel_back_DA = Channel("back_DA", list(np.round(back_trace_DA * area_mask_A,0).astype(np.int64)))

                channel_back_AA = Channel("back_AA", list(np.round(back_trace_AA * area_mask_A,0).astype(np.int64)))
                
                trace1_i = Trace([channel_DD, channel_DA, channel_AA, channel_back_DD, channel_back_DA, channel_back_AA],
                                 metadata=Metadata({"molecule_id": str(k), 'x_coord': coord_spots[k,0], 'y_coord': coord_spots[k,1], "UUID_v7": str(uuid7())}))
                
                Open_dataset.traces.append(trace1_i)
                
            else:
            # concatenation of subtraces to their cooresponding trace ID

                
                Open_dataset.traces[k].channels[0].data = Open_dataset.traces[k].channels[0].data + list(np.round(trace_DD,0).astype(np.int64))
                Open_dataset.traces[k].channels[1].data = Open_dataset.traces[k].channels[1].data + list(np.round(trace_DA,0).astype(np.int64))
                Open_dataset.traces[k].channels[2].data = Open_dataset.traces[k].channels[2].data + list(np.round(trace_AA,0).astype(np.int64))
                Open_dataset.traces[k].channels[3].data = Open_dataset.traces[k].channels[3].data + list(np.round(back_trace_DD * area_mask,0).astype(np.int64))
                Open_dataset.traces[k].channels[4].data = Open_dataset.traces[k].channels[4].data + list(np.round(back_trace_DA * area_mask_A,0).astype(np.int64))
                Open_dataset.traces[k].channels[5].data = Open_dataset.traces[k].channels[5].data + list(np.round(back_trace_AA * area_mask_A,0).astype(np.int64))
            
    print('Traces extraction and concatenation completed!')
    
    return Open_dataset

def calculate_dist_peaks(img_stack_acceptor, img_stack_donor, V = 0, U = 0, matrix_align = 0, k = 0, 
                         kernel_size = 2, min_distance = 1, method = 'Optical FLow'):
    
    peaks_A = detection.detect_spots(img_stack_acceptor[:,:,k], log_kernel_size=kernel_size, minimum_distance=min_distance)
    
    if method == 'Optical FLow':
        img_warp_donor = Warp_OpticalFlow(img_stack_donor, V, U)
    else:
        img_warp_donor = generate_chrom_ab_corr_movie(img_stack_donor, matrix_align)
    
    peaks_D = detection.detect_spots(img_warp_donor[:,:,k], log_kernel_size=kernel_size, minimum_distance=min_distance)
    
    test_V1, test_V2 = np.meshgrid(peaks_A[:,0], peaks_D[:,0])
    test_U1, test_U2 = np.meshgrid(peaks_A[:,1], peaks_D[:,1])
    
    dist_mat = ((test_V1 - test_V2)**2+(test_U1 - test_U2)**2)**0.5
    
    len_peaks_A = len(peaks_A)
    
    len_peaks_D = len(peaks_D)
    
    if len_peaks_A > len_peaks_D:
        min_dist_vec = np.min(dist_mat, axis = 0)
    else:
        min_dist_vec = np.min(dist_mat, axis = 1)
        
    data = np.sort(min_dist_vec)
    cdf = np.arange(1, len(min_dist_vec) + 1) / len(min_dist_vec)
    plt.figure()
    plt.plot(data, cdf, color='blue')
    plt.xlabel('inter-peaks distance (px)')
    plt.ylabel('CDF')
    plt.title('Cumulative distribution of calibration inter-peaks distance')
    plt.grid()
    plt.show()
