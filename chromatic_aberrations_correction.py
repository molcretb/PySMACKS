# -*- coding: utf-8 -*-
"""
Created on Thu Dec  4 08:45:46 2025

@author: Bastien Molcrette, Schmid lab, Chemistry department, University of Basel
"""
import numpy as np
from tkinter.filedialog import askopenfilenames
from tkinter import Tk
import matplotlib.pyplot as plt
from PIL import Image
from scipy.optimize import minimize
from cv2 import warpPerspective

def generate_trans_matrix(x, H, W):
    """
    Function that return the transformation matrix for chromatic aberrations, based on see https://doi.org/10.1038/s41467-021-26466-7

    Parameters
    ----------
    x : list of float
        [sx, sy, dx, dy, cx, cy, theta] parameters for the chromatic correction matrix
    H : int
        Height of the image to be corrected.
    W : int
        Width of the image to be corrected.

    Returns
    -------
    omega : Numpy matrix
        Transformation matrix to be used for chromatic aberrations correction.

    """
    
    # matA: scaling factors matrix; matB: translation matrix; matC: rotation matrix, with center of rotation given by matD \
        # see https://doi.org/10.1038/s41467-021-26466-7 for more details
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
    # we combine all transformation matrix into a single one
    omega = np.matmul(matA, np.matmul(matB, np.matmul(matC, matD)))
    
    return omega

def add_trans_matrix_to_traces_dict(traces_dict, matrix_align):
    """
    Function that add the chromatic aberrations correction matrix to the traces dictionary

    Parameters
    ----------
    traces_dict : dict
        Dictionnary containing all the traces results.
    matrix_align : list of float
        [sx, sy, dx, dy, cx, cy, theta] parameters for the optimized chromatic correction matrix

    Returns
    -------
    traces_dict : dict
        Dictionnary containing all the traces results with updated optimized chromatic correction matrix.

    """
    
    traces_dict['chromatic_aberration_corr_matrix'] = matrix_align
    
    return traces_dict
    

def logP_chrom_corr(x, img1, img2):
    """
    Function that return the negative log of the sum of square difference between the shifted image compared to the reference one

    Parameters
    ----------
    x : list of float
        [sx, sy, dx, dy, cx, cy, theta] parameters for the chromatic correction matrix
    img1 : Numpy array
        image to be corrected.
    img2 : Numpy array
        reference image.

    Returns
    -------
    neglogP : float
        negative log of the sum of square difference between the shifted image compared to the reference one, to be minimized

    """
    # Height and width of the image
    H, W = img1.shape[0:2]
    # generate the transformation matrix
    omega = generate_trans_matrix(x, H, W)
    
    # transform image1 (the one to be corrected) given the transformation matrix omega, using CV2
    trans_img1 = warpPerspective(img1,omega, (W, H))
    
    # calculate the negative log of the squared differences between the corrected image1 and the reference image2
    neglogP = 0.5 * W * H * np.log(np.sum((trans_img1-img2)**2))
    
    return neglogP

def minimize_logP(img1, img2, method_min='Nelder-Mead', maxfev_value = 5000):
    """
    Function that minimize the negative log squared differences between the corrected image1 and reference image 2

    Parameters
    ----------
    img1 : Numpy array
        Image to be corrected.
    img2 : Numpy array
        Reference image.
    method_min : Str, optional
        minimization method to be employed, default is 'Nelder-Mead', see Scipy doc for more details
    maxfev_value : int, optional
        Maximum allowed number of function evaluations, default is set at 5000, see Scipy doc for more details
        

    Returns
    -------
    res : OptimizeResult
        Results from Scipy minimize function, the transformation matrix coefficients are stored in res.x

    """
    # Height and width of the image
    H, W = img1.shape[0:2]
    # we want to minimize logP_chrom_corr with regard to the transformation matrix coefficients, with initial guess, using the minimze function from Scipy
    res = minimize(logP_chrom_corr, [1, 1, 1, 1, int(W/2), int(H/2), 0.01], method=method_min, args=(img1, img2), options={'maxfev': maxfev_value})
    
    return res

def generate_chrom_ab_corr_movie(img1, x):
    """
    Function that generate a chromatic aberrations corrected image

    Parameters
    ----------
    img1 : Numpy array
        Image to be corrected.
    x : list of float
        [sx, sy, dx, dy, cx, cy, theta] parameters for the chromatic correction matrix

    Returns
    -------
    trans_img1 : Numpy array
        Chromatic aberrations corrected image.

    """
    # Height and width of the image
    H, W = img1.shape[0:2]
    
    # generate the transformation matrix
    omega = generate_trans_matrix(x, H, W)
    
    # transform image1 (the one to be corrected) given the optimized transformation matrix previously calculated, using CV2
    trans_img1 = warpPerspective(img1,omega, (W, H))
   
    return trans_img1

def show_align_res(matrix_align, img_donor, img_acceptor, k = 0):
    """
    Function that generates a chromatic aberrations corrected image in RGB with red channel as the corrected donor image and green channel as the donor image (blue channel is 0)

    Parameters
    ----------
    matrix_align : list of float
        [sx, sy, dx, dy, cx, cy, theta] parameters for the optimized chromatic correction matrix
    img_donor : Numpy array
        Donor image stack.
    img_acceptor : Numpy array
        Acceptor image stack.
    k : int, optional
        the frame from the donor/acceptor stacks we want to display. The default is 0.

    Returns
    -------
    None.

    """
    # generate the corrected donor image (frame k for timelapse movie)
    trans_img_donor = generate_chrom_ab_corr_movie(img_donor[:,:,k], matrix_align)
    
    # Create the red and green channels with corrected donor and reference acceptor images, and normalized them with their respective maximum values
    R_channel = trans_img_donor/np.max(trans_img_donor)
    G_channel = img_acceptor[:,:,k]/np.max(img_acceptor[:,:,k])
    B_channel = np.zeros(trans_img_donor.shape)
    
    # Show the merge channels
    plt.figure('Chromatic abberation correction, frame = ' + str(k))
    plt.imshow(np.stack((R_channel, G_channel, B_channel), axis=-1))
    return

def pipeline_chrom_ab_correction(nb_ref_frames = 5, method_min_opt='Nelder-Mead', maxfev_value_opt = 5000):
    """
    Main function to run the chromatic aberrations correction

    Parameters
    ----------
    nb_ref_frames : int, optional
        DESCRIPTION. Number of frames to be used for the optimization process; accuracy of the optimization should increase with this number, nut also takes more time. The default is 5.
    method_min_opt : Str, optional
        minimization method to be employed, default is 'Nelder-Mead', see Scipy doc for more details
    maxfev_value_opt : int, optional
        Maximum allowed number of function evaluations, default is set at 5000, see Scipy doc for more details

    Returns
    -------
    matrix_align : list of arrays
        List of the optimized transformation matrix to use for the chromatic aberrations correction.

    """
    
    print('Start of the chromatic aberrations correction step')
    # Select first the donor calibration movie and second the acceptor calibration movie \
        # (we used the donor channel as reference for the chromatic aberrations correction)
    root = Tk(className='Open calibration channels', )
    file_path_donor = askopenfilenames(title="Select the Donor calibration channel movie")
    file_path_acceptor = askopenfilenames(title="Select the Acceptor calibration channel movie")
    root.destroy()
    # try to open the movies, end the process if the files don't match the requirements as file type (likely TIFF multi-planes files)
    try:
        img_donor = Image.open(file_path_donor[0])
        print(f"File '{file_path_donor}' loaded.")
        img_acceptor = Image.open(file_path_acceptor[0])
        print(f"File '{file_path_acceptor}' loaded.")
    except FileNotFoundError:
        print(f"File '{file_path_donor}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
    
    # get the length of the movie
    nb_frame = img_donor.n_frames
    
    # pick a given number of frames to be used as references for the correction optimization, constant interval between each picked frame\
        # so it spans the entire length of the movie
    select_ref_frame_list = np.linspace(0,nb_frame-1,nb_ref_frames).astype(int)
    
    # Numpy arrays that are used to stored the selected frames
    img_stack_donor = np.zeros((img_donor.height, img_donor.width, nb_ref_frames))
    img_stack_acceptor = np.zeros((img_acceptor.height, img_acceptor.width, nb_ref_frames))
    count_i = 0
    
    # we use the seek method from PIL to directly pick the selected frame without needing to load the other ones
    for i in select_ref_frame_list:
        img_donor.seek(i)
        img_stack_donor[:,:,count_i] = np.array(img_donor)
        img_acceptor.seek(i)
        img_stack_acceptor[:,:,count_i] = np.array(img_acceptor)
        count_i = count_i + 1
    
    # we run the optimization algorithm with the selected stacks of frames from the donor and acceptor movies\
        # the optional parameters method_min and maxfev_value can be provided to finetune the optimization process, if needed
    print('Optimization of the transformation matrix ongoing...')
    res_align = minimize_logP(img_stack_donor, img_stack_acceptor, method_min=method_min_opt, maxfev_value = maxfev_value_opt)
    print('Optimization of the transformation matrix done!')
    
    # we extract the optimized transformation matrix coefficients
    matrix_align = res_align.x
    
    # Show the result on the first frame
    show_align_res(matrix_align, img_stack_donor, img_stack_acceptor, k = 0)
    
    print('Chromatic aberrations correction step completed!')
    
    return matrix_align