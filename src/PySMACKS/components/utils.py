# -*- coding: utf-8 -*-
"""
Created on Thu Dec  4 10:31:28 2025

@author: Bastien Molcrette, Schmid lab, Chemistry department, University of Basel
"""
import numpy as np
from tkinter.filedialog import askopenfilenames
from tkinter import Tk
from PIL import Image
from tkinter.filedialog import asksaveasfile
from json import dump, JSONEncoder, load
import matplotlib.pyplot as plt
from pandas import DataFrame as DF
from seaborn import jointplot as jointplot
from scipy.optimize import curve_fit
import scipy.optimize as opt
from scipy.optimize import least_squares
from TVDCondat2013 import tvd_2013
import shutil
import os
import copy
from openfret import read_data, write_data

def load_data():
    """
    Function used to select the set of TIFF files for a given experiment, and return the first submovie as Numpy array and list of paths
    
    Parameters
    ----------
    None

    Returns
    -------
    img_stack : Numpy array
        First submovie from the list of selected TIFF files.
    file_path : list of str
        List of the absolute paths of all selected TIFF files.

    """
    # Select the TIFF movies you want to load
    root = Tk(className='Open TIFF movie', )
    file_path = askopenfilenames(title="Select the TIFF movie")
    root.destroy()
    
    # try to open the first movie of the list of movies, end the process if the file doesn't match the requirements as file type (likely TIFF multi-planes files)
    try:
        img = Image.open(file_path[0])
        print(f"File '{file_path}' loaded.")
    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

    # get the length of the movie
    nb_frame = img.n_frames
    
    # Numpy array that is used to stored the movie
    img_stack = np.zeros((img.height, img.width, nb_frame))
    
    # we use the seek method from PIL to pick the frames one by one
    for i in range(nb_frame):
        img.seek(i)
        img_stack[:,:,i] = np.array(img)
        
    # return the 1st movie as Numpy array and the list of path of selected movies
    return img_stack, file_path

def load_submovies(filename):
    """
    Function used to load a submovie from its path

    Parameters
    ----------
    filename : str
        Absolute path of the TIFF movie to be loaded.

    Returns
    -------
    img_stack : Numpy array
        Submovie as Numpy array.

    """
    
    # Open the TIFF file with PIL
    img = Image.open(filename)
    
    # get the length of the movie
    nb_frame = img.n_frames
    
    # Numpy array that is used to stored the movie
    img_stack = np.zeros((img.height, img.width, nb_frame))
    
    # we use the seek method from PIL to pick the frames one by one
    for i in range(nb_frame):
        img.seek(i)
        img_stack[:,:,i] = np.array(img)
        
    # return the movie as Numpy array
    return img_stack

def select_donor_and_acceptor_movies():
    
    root = Tk(className='Open TIFF movies', )
    file_path_A = askopenfilenames(title="Select the acceptor TIFF movies")
    file_path_D = askopenfilenames(title="Select the donor TIFF movies")
    root.destroy()
    
    return file_path_D, file_path_A

def deinterleave_acceptor_channel(img_A, DA_is = 'odd'):
    """
    Function used to deinterleave the acceptor channel when doing ALEX FRET, even frames are supposed to be the AA, odd ones are DA

    Parameters
    ----------
    img_A : Numpy array
        Acceptor movie as Numpy array stack with interleaved DA and AA channels.
    DA_is : str, optional
        'odd' if DA frames are the odd ones, anything else if there are the even ones
        

    Returns
    -------
    img_DA : Numpy array
        DA channel as Numpy array.
    img_AA : Numpy array
        AA channel as Numpy array.

    """
    
    # if the DA frames are the odd ones
    if DA_is == 'odd':
        img_AA = img_A[:,:,[i for i in range(0,img_A.shape[2],2)]]
        img_DA = img_A[:,:,[i for i in range(1,img_A.shape[2],2)]]
    else:
        img_AA = img_A[:,:,[i for i in range(1,img_A.shape[2],2)]]
        img_DA = img_A[:,:,[i for i in range(0,img_A.shape[2],2)]]
    
    #return the deinterleaved DA and AA channels as separate Numpy arrays
    return img_DA, img_AA

class NpEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)


def save_traces_dict(traces_dict):
    """
    This function is used to save the traces results dictionary as a JSON file.

    Parameters
    ----------
    traces_dict : dict
        Dictionnary containing all the traces results.

    Returns
    -------
    None.

    """
    
    # A Tkinter window first appears and asks the user to input the name and path of the JSON file to save the traces dataset
    root = Tk(className='Save traces datasets', )
    file_save = asksaveasfile(title="Save traces dataset", defaultextension=".json", initialfile="NAME_of_your_DATASET.json")
    file_save.close()
    root.destroy()
    
    # The traces dictionary is then saved at the path location provided by the user
    with open(file_save.name, "w") as fp:
        dump(traces_dict, fp, cls=NpEncoder, separators=(',', ':')) 
    print('Traces dataset saved as JSON file to path: ' + file_save.name)
    return

def select_autosave_file():
    root = Tk(className='Create an autosave file...', )
    file_save = asksaveasfile(title="Create an autosave file", defaultextension=".json", initialfile="NAME_of_your_temporary_file.json")
    file_save.close()
    root.destroy()
    return file_save.name

def autosave_intermed_res(autosave_file, dict_auto_save):
    
    with open(autosave_file, "w") as fp:
        dump(dict_auto_save, fp, cls=NpEncoder, separators=(',', ':')) 
    print('Autosave done!' + autosave_file)

def convert_calib_dict2Var(calib_dict):
    align_method = calib_dict['method']
    
    if align_method == 'Optical Flow':
        V = np.array(calib_dict['matrix_align'][0])
        U = np.array(calib_dict['matrix_align'][1])
        matrix_align = (V, U)
    else:
        matrix_align = np.array(calib_dict['matrix_align'])
    
    return matrix_align, align_method

def load_JSON_traces_data():
    """
    Function used to load a previously saved traces dataset JSON file
    
    Parameters
    ----------
    None.

    Returns
    -------
    traces_dict : dict
        Dictionnary containing all the traces results.

    """
    
    # Select the JSON dataset file you want to load
    root = Tk(className='Open JSON dataset', )
    file_path = askopenfilenames(title="Open the JSON dataset")
    root.destroy()
    
    # Read the JSON file and convert it to a dictionary
    with open(file_path[0]) as json_file:
        traces_dict = load(json_file)
        
    print('Traces dataset loaded!')
    
    # Return the traces dataset as dictionary
    return traces_dict

def plot_multiple_traces(traces_data, traces_ID = 0, DD_on = 1, DA_on = 1, AA_on = 1, DD_DA_on = 1, nb_plot = 5):
    """
    Function used to plot multiple traces on subplots, can include or exclude each channels

    Parameters
    ----------
    traces_data : dict
        Dictionnary containing all the traces results.
    traces_ID : int, optional
        First trace ID we want to plot, the other nb_plot - 1 traces are simply the next ones in the traces IDs order.
    nb_plot : int, optional
        Number of different traces we want to subplot. The default is 5.
    DD_on : bool, optional
        1 if DD trace should be included, 0 to exclude
    DA_on : bool, optional
        1 if DA trace should be included, 0 to exclude
    AA_on : bool, optional
        1 if AA trace should be included, 0 to exclude
    DD_DA_on : bool, optional
        1 if DD + DA trace should be included, 0 to exclude

    Returns
    -------
    None.

    """
    list_traces_IDs = list(traces_data.keys())
    
    fig, axs = plt.subplots(nb_plot)
    fig.suptitle('Traces ' + list_traces_IDs[traces_ID] + ' to ' + list_traces_IDs[traces_ID + nb_plot - 1])
    for k in range(nb_plot):
        if DD_on == 1:
            axs[k].plot(traces_data[list_traces_IDs[traces_ID + k]]['Intensity_DD'],color='orange');
            axs[k].axvline(x = traces_data[list_traces_IDs[traces_ID + k]]['bleaching_event_DD'], color = 'orange', label = 'predicted breakpoint');
        if DA_on == 1:
            axs[k].plot(traces_data[list_traces_IDs[traces_ID + k]]['Intensity_DA'],color='red');
        if AA_on == 1:
            axs[k].plot(traces_data[list_traces_IDs[traces_ID + k]]['Intensity_AA'],color='gray');
            axs[k].axvline(x = traces_data[list_traces_IDs[traces_ID + k]]['bleaching_event_AA'], color = 'gray', label = 'predicted breakpoint');
        if DD_DA_on == 1:
            axs[k].plot([i+j for i,j in zip(traces_data[list_traces_IDs[traces_ID + k]]['Intensity_DD'],traces_data[list_traces_IDs[traces_ID + k]]['Intensity_DA'])],color='blue');
    fig.supxlabel('Absolute raw intensity')
    fig.supylabel('Frames')
    plt.show()
    return


def plot_Stoichio_FRETeff_2D_hist(traces_data, nb_bin = 100):
    list_Eff_Stoichio = []
    
    for k in list(traces_data.keys()):
        
        bleach_event_DD = traces_data[k]['bleaching_event_DD']
        bleach_event_AA = traces_data[k]['bleaching_event_AA']
        end_trace = np.min((bleach_event_DD, bleach_event_AA))
        
        if end_trace > 0:
        
            DD_trace = np.array(traces_data[k]['Intensity_DD'][:end_trace])
            DA_trace = np.array(traces_data[k]['Intensity_DA'][:end_trace])
            AA_trace = np.array(traces_data[k]['Intensity_AA'][:end_trace])
        
        
        
            FRET_eff = DA_trace / (DD_trace + DA_trace)
            FRET_stoichio = (DD_trace + DA_trace) / (DD_trace + DA_trace + AA_trace)
            if len(list_Eff_Stoichio) == 0:
                list_Eff_Stoichio = np.stack((FRET_eff, FRET_stoichio), axis = 1)
            else:
                list_Eff_Stoichio = np.concat((list_Eff_Stoichio, np.stack((FRET_eff, FRET_stoichio), axis = 1)))
    
    df_data = DF(data={'Fret_eff': list_Eff_Stoichio[:,0], 'stoichio': list_Eff_Stoichio[:,1]})
    
    jointplot(data=df_data, x="Fret_eff", y="stoichio",kind="hist", 
              joint_kws={'cmap':'viridis'}, xlim=[0, 1], ylim=[0, 1],
              marginal_kws=dict(bins=nb_bin))
    
    return


def concat_submovies(file_path):
    nb_submovies = len(file_path)
    
    for i in range(nb_submovies):
        img_stack_i = load_submovies(file_path[i])
        
        if i == 0:
            img_stack_tot = img_stack_i
        else:
            img_stack_tot = np.concatenate((img_stack_tot, img_stack_i), axis=2)
    return img_stack_tot

def generate_syn_movie(sizeX, sizeY, movie_len, nb_spot, sigma_x = 1, sigma_y = 1):
    
    x_vec = np.linspace(0, sizeX - 1, sizeX)
    y_vec = np.linspace(0, sizeY - 1, sizeY)
    x, y = np.meshgrid(y_vec, x_vec)
    
    x = np.rollaxis(np.stack([x]*nb_spot),0,3)
    y = np.rollaxis(np.stack([y]*nb_spot),0,3)
    
    movie_syn = np.random.poisson(size = (sizeX, sizeY, movie_len))
    
    # Generate list of spot coordinates
    coord_spots = np.random.rand(nb_spot, 2)
    coord_spots[:,0] = coord_spots[:,0] * sizeY
    coord_spots[:,1] = coord_spots[:,1] * sizeX
    
    # Generate a drift
    
    drift = np.zeros((movie_len, 2))
    
    for i in range(movie_len):
        print('Processing frame ' + str(i))
        if i > 0:
            drift[i,:] = drift[i-1,:] + np.random.normal(0, 0.7, 2)
        movie_syn[:,:,i] = movie_syn[:,:,i] + 7*np.max(np.exp( - ((x-coord_spots[:,1]-drift[i,1])**2/(2*sigma_x**2) + (y-coord_spots[:,0]-drift[i,0])**2/(2*sigma_y**2))), axis = 2)
        
    return movie_syn, coord_spots, drift

def merge_trace_datasets():
    root = Tk(className='Open TIFF movie', )
    file_path = askopenfilenames(title="Select the TIFF movie")
    root.destroy()
    
    traces_dict = {}
    
    for i in file_path:
        with open(i) as json_file:
            traces_dict_i = load(json_file)
        list_keys_i = traces_dict_i.keys()
        list_keys = traces_dict.keys()
        for key in list_keys_i:
            if key in list_keys:
                compt = 0
                while key + '_' + str(compt) in list_keys:
                    compt = compt + 1
                traces_dict[key + '_' + str(compt)] = traces_dict_i[key]
            else:
                traces_dict[key] = traces_dict_i[key]
    
    save_traces_dict(traces_dict)
                
    return


def FRET_corr_pipeline(nb_gauss = 2, nb_bins = 100):
    
    root = Tk(className='Open JSON Donor-Acceptor dataset', )
    file_path_DA = askopenfilenames(title="Select JSON Donor-Acceptor dataset")
    root.destroy()
    root = Tk(className='Open JSON Donor-Only dataset', )
    file_path_DO = askopenfilenames(title="Select JSON Donor-Only dataset")
    root.destroy()
    root = Tk(className='Open JSON Acceptor-Only dataset', )
    file_path_AO = askopenfilenames(title="Select JSON Acceptor-Only dataset")
    root.destroy()
    
    # get alpha correction factor from DD data
    
    print('Calculation of alpha correction factor...')
    with open(file_path_DO[0]) as json_file:
        traces_dict_DO = load(json_file)
    traces_dict_DO = calculate_FRET_Eff(traces_dict_DO)
    alpha = calculate_alpha_corr_factor(traces_dict_DO)
    plot_corr_histo_FRET(traces_dict_DO, alpha, 0, 1, 1)
    
    # get delta correction factor from AA data
    
    print('Calculation of delta correction factor...')
    with open(file_path_AO[0]) as json_file:
        traces_dict_AO = load(json_file)
    traces_dict_AO = calculate_FRET_Eff(traces_dict_AO)
    delta = calculate_delta_corr_factor(traces_dict_AO)
    plot_corr_histo_FRET(traces_dict_AO, 0, delta, 1, 1)
    
    # get beta and gamma correction factors from DA data
    
    print('Calculation of beta and gamma correction factors...')
    with open(file_path_DA[0]) as json_file:
        traces_dict_DA = load(json_file)
    traces_dict_DA = calculate_FRET_Eff(traces_dict_DA)
    traces_dict_DA = plot_corr_histo_FRET(traces_dict_DA, alpha, delta, 1, 1, plot_hist = 0, return_dict = 1, bin_size = 0.01, ax_lim_eps = 0.2)
    beta, gamma = calculate_beta_gamma_corr_factor(traces_dict_DA, nb_gauss = nb_gauss, bin_size = 0.01, ax_lim_eps = 0.2)
    
    # get FRET corrected 2D histogram
    
    corr_histo = plot_corr_histo_FRET(traces_dict_DA, alpha, delta, beta, gamma, bin_size = 0.01, ax_lim_eps = 0.2)
    
    corr_FRET_eff_hist = np.sum(corr_histo, axis = 0)[1:-1]
    
    x_fret_eff = np.linspace(1/nb_bins, (nb_bins - 2)/nb_bins, nb_bins-2)
    
    hist_1D = np.stack((x_fret_eff, corr_FRET_eff_hist), axis = 0)
    
    df_FRET_peaks = multi_fit_Gauss_spot(hist_1D, nb_gauss = nb_gauss)
    
    return corr_histo, hist_1D, [alpha, delta, beta, gamma], df_FRET_peaks

def calculate_FRET_Eff(traces_dict):
    
    ID_traces = list(traces_dict.keys())
    
    FRET_conc = np.array([])
    
    FRET_S_conc = np.array([])
    
    for i in ID_traces:
        trace_DD = np.array(traces_dict[i]['Intensity_DD'])
        trace_DA = np.array(traces_dict[i]['Intensity_DA'])
        trace_AA = np.array(traces_dict[i]['Intensity_AA'])
        
        FRET_eff_i = trace_DA / (trace_DD + trace_DA)
        
        FRET_eff_i_list = np.ma.masked_invalid(FRET_eff_i).mask
        
        FRET_stoi_i = (trace_DA + trace_DD) / (trace_DD + trace_DA + trace_AA)
        
        FRET_stoi_i_list = np.ma.masked_invalid(FRET_stoi_i).mask
        
        list_FRET_OK = ~np.logical_or(FRET_eff_i_list, FRET_stoi_i_list)
        
        traces_dict[i]['FRET_eff'] = FRET_eff_i[list_FRET_OK]
        
        traces_dict[i]['FRET_stoi'] = FRET_stoi_i[list_FRET_OK]
        
        FRET_conc = np.concat((FRET_conc, traces_dict[i]['FRET_eff']))
        
        FRET_S_conc = np.concat((FRET_S_conc, traces_dict[i]['FRET_stoi']))
    
    return traces_dict

def calculate_alpha_corr_factor(traces_dict):
    
    ID_traces = list(traces_dict.keys())
    
    FRET_conc = np.array([])
    
    for i in ID_traces:
        FRET_conc = np.concat((FRET_conc, traces_dict[i]['FRET_eff']))
    fig, ax = plt.subplots()
    counts, bins, bars = ax.hist(FRET_conc, 100, range=(-1,1))
    plt.close(fig)
    bins_center = 0.5*(bins[1:]+bins[0:-1])
    popt, pcov = curve_fit(Gauss_data, bins_center[1:-1], counts[1:-1], p0 = [np.max(counts[1:-1]),np.argmax(counts[1:-1]),np.argmax(counts[1:-1])])
    alpha = popt[2]/(1-popt[2])
    
    x = np.linspace(-1, 1, 1000)
    
    y_gauss = Gauss_data(x, popt[0], popt[1], popt[2])
    
    fig, ax = plt.subplots()
    ax.plot(bins_center[1:-1], counts[1:-1],'kx', label='data')
    label_fit = r'$\alpha$ = '+str(np.round(alpha, 4))
    ax.plot(x, y_gauss, label=label_fit)
    plt.title(r'$\alpha$ correction from Donor-Only traces')
    plt.xlabel('FRET efficiency')
    plt.ylabel('Counts')
    ax.legend()
    
    return alpha

def calculate_delta_corr_factor(traces_dict):
    
    ID_traces = list(traces_dict.keys())
    
    FRET_S_conc = np.array([])
    
    for i in ID_traces:
        FRET_S_conc = np.concat((FRET_S_conc, traces_dict[i]['FRET_stoi']))
    fig, ax = plt.subplots()
    counts, bins, bars = ax.hist(FRET_S_conc, 100, range=(-1,1))
    plt.close(fig)
    bins_center = 0.5*(bins[1:]+bins[0:-1])
    popt, pcov = curve_fit(Gauss_data, bins_center[1:-1], counts[1:-1], p0 = [np.max(counts[1:-1]),np.argmax(counts[1:-1]),np.argmax(counts[1:-1])])
    delta = popt[2]/(1-popt[2])
    
    x = np.linspace(-1, 1, 1000)
    
    y_gauss = Gauss_data(x, popt[0], popt[1], popt[2])
    
    fig, ax = plt.subplots()
    ax.plot(bins_center[1:-1], counts[1:-1],'kx', label='data')
    label_fit = r'$\delta$ = '+str(np.round(delta, 4))
    ax.plot(x, y_gauss, label=label_fit)
    plt.title(r'$\delta$ correction from Acceptor-Only traces')
    plt.xlabel('Stoichiometry')
    plt.ylabel('Counts')
    ax.legend()
    
    return delta

def plot_corr_histo_FRET(traces_dict, alpha, delta, beta, gamma, plot_hist = 1, return_dict = 0, return_fret_conc = 0, bin_size = 0.01, ax_lim_eps = 0.2):
    
    ID_traces = list(traces_dict.keys())
    
    FRET_conc = np.array([])
    
    FRET_S_conc = np.array([])
    
    for i in ID_traces:
        trace_DD = np.array(traces_dict[i]['Intensity_DD'])
        trace_DA = np.array(traces_dict[i]['Intensity_DA'])
        trace_AA = np.array(traces_dict[i]['Intensity_AA'])
        
        FRET_eff_i = (trace_DA - alpha * trace_DD - delta * trace_AA) / (gamma * trace_DD + trace_DA - alpha * trace_DD - delta * trace_AA)
        
        FRET_stoi_i = (trace_DA - alpha * trace_DD + gamma * trace_DD - delta * trace_AA) / (gamma * trace_DD + trace_DA - alpha * trace_DD + trace_AA / beta - delta * trace_AA)
        
        FRET_eff_i_list = np.ma.masked_invalid(FRET_eff_i).mask
        
        FRET_stoi_i_list = np.ma.masked_invalid(FRET_stoi_i).mask
        
        list_FRET_OK = ~np.logical_or(FRET_eff_i_list, FRET_stoi_i_list)
        
        traces_dict[i]['FRET_eff_corr'] = FRET_eff_i[list_FRET_OK]
        
        traces_dict[i]['FRET_stoi_corr'] = FRET_stoi_i[list_FRET_OK]
        
        FRET_conc = np.concat((FRET_conc, traces_dict[i]['FRET_eff_corr']))
        
        FRET_S_conc = np.concat((FRET_S_conc, traces_dict[i]['FRET_stoi_corr']))
        
    hist2D_map = generate_2D_histogram(FRET_conc, FRET_S_conc, bin_size = bin_size, ax_lim_eps = ax_lim_eps)
    
    if plot_hist == 1:
        plt.figure()
        plt.imshow(hist2D_map, interpolation='none', origin='lower', extent=[-ax_lim_eps,1+ax_lim_eps,-ax_lim_eps,1+ax_lim_eps])
        plt.xlabel('FRET efficiency')
        plt.ylabel('Stoichiometry')
        plt.title(r'$\alpha$-$\delta$-$\beta$-$\gamma$ corrected FRET histogram')
    if return_dict == 1:
        if return_fret_conc == 1:
            return traces_dict, FRET_conc, FRET_S_conc
        else:
            return traces_dict
    else:
        if return_fret_conc == 1:
            return hist2D_map, FRET_conc, FRET_S_conc
        else:
            return hist2D_map

def Gauss_data(x, A, B, C):
    return A * np.exp(-(x-C)**2/B**2)

def generate_2D_histogram(FRET_conc, FRET_S_conc, bin_size = 0.01, ax_lim_eps = 0.2):
    
    nb_bins = int((1+2*ax_lim_eps)/bin_size)
    
    hist2D_map = np.zeros((nb_bins, nb_bins))
    
    compt_nan = 0
    
    
    extra_bins = int(ax_lim_eps // bin_size)
    
    FRET_hist_bin = [FRET_S_conc // bin_size + extra_bins, FRET_conc // bin_size + extra_bins]
    
    for i in range(len(FRET_conc)):
        try:
            x_i = FRET_hist_bin[0][i].astype(int)
            y_i = FRET_hist_bin[1][i].astype(int)
            if (x_i >= 0) and (y_i >=0):
                hist2D_map[x_i, y_i] += 1
        except IndexError:
            compt_nan += 1
    
    return hist2D_map
    
def calculate_beta_gamma_corr_factor(traces_dict, nb_gauss = 2, bin_size = 0.01, ax_lim_eps = 0.2):
    
    ID_traces = list(traces_dict.keys())
    
    FRET_conc = np.array([])
    FRET_S_conc = np.array([])
    
    for i in ID_traces:
        FRET_conc = np.concat((FRET_conc, traces_dict[i]['FRET_eff_corr']))
        FRET_S_conc = np.concat((FRET_S_conc, traces_dict[i]['FRET_stoi_corr']))
        
    hist2D_map = generate_2D_histogram(FRET_conc, FRET_S_conc, bin_size = bin_size, ax_lim_eps = ax_lim_eps)
    
    peak_coord = fit_Gauss_spot(hist2D_map, nb_gauss = nb_gauss, bin_size = bin_size, ax_lim_eps = ax_lim_eps)
    
    nb_bins = int((1+2*ax_lim_eps)/bin_size)
    
    x_ax = np.linspace(-ax_lim_eps,1+ax_lim_eps, nb_bins)
    
    b, a = np.polyfit(x_ax[peak_coord[:,0].astype(int)], 1/x_ax[peak_coord[:,1].astype(int)], 1)
    
    beta = a + b -1
    gamma = (a - 1) / (a + b - 1)
    
    return beta, gamma

def twoD_Gaussian(xy, gauss_param, nb_gauss):
    y, x = xy
    g = gauss_param[0] + np.zeros(x.shape)
    for i in range(nb_gauss):
        g = g + gauss_param[5*i+1]*np.exp( - ((x-float(gauss_param[5*i+2]))**2/(2*gauss_param[5*i+4]**2) + (y-float(gauss_param[5*i+3]))**2/(2*gauss_param[5*i+5]**2)))
    return g.ravel()

def error_func(gauss_param, nb_gauss, x, y):
    return twoD_Gaussian(x, gauss_param, nb_gauss) - y

def fit_Gauss_spot(img_spot, nb_gauss = 2, bin_size = 0.01, ax_lim_eps = 0.2):
    x_len = img_spot.shape[0]
    y_len = img_spot.shape[1]
    x_vec = np.linspace(0, x_len-1, x_len)
    y_vec = np.linspace(0, y_len-1, y_len)
    x, y = np.meshgrid(y_vec, x_vec)
    max_img = np.max(img_spot)
    noise_img = np.median(img_spot)
    gauss_param_estim = [noise_img]
    Fret_eff_estim = np.linspace(0.3, 0.7, nb_gauss)
    for i in range(nb_gauss):
        gauss_param_estim.append(max_img-noise_img)
        gauss_param_estim.append(x_len*0.5)
        gauss_param_estim.append(x_len*Fret_eff_estim[i])
        gauss_param_estim.append(10)
        gauss_param_estim.append(10)
    initial_guess = tuple(gauss_param_estim)
    try:
        XX = (x, y)
        output = least_squares(error_func, x0=initial_guess, jac='2-point', bounds=(0, np.inf), method='trf', ftol=1e-08,
                       xtol=1e-08, gtol=1e-08, x_scale=1.0, loss='linear', f_scale=1.0, diff_step=None,
                       tr_solver=None,
                       tr_options={}, jac_sparsity=None, max_nfev=None, verbose=0, args=(nb_gauss, XX, img_spot.ravel()))
        
    except RuntimeError:
        print('No peak found')
        return
    
    peak_coord = np.zeros((nb_gauss, 2))
    
    for j in range(nb_gauss):
        peak_coord[j,0] = output.x[5*j+3]  # FRET eff
        peak_coord[j,1] = output.x[5*j+2]  # FRET stoich
    
    
    nb_bins = int((1+2*ax_lim_eps)/bin_size)
    
    x_ax = np.linspace(-ax_lim_eps,1+ax_lim_eps, nb_bins)
        
    plt.figure()
    plt.imshow(img_spot,origin='lower', interpolation='none', extent=[-ax_lim_eps,1+ax_lim_eps,-ax_lim_eps,1+ax_lim_eps])
    plt.plot(x_ax[peak_coord[:,0].astype(int)], x_ax[peak_coord[:,1].astype(int)], 'ro')
    plt.xlabel('FRET efficiency')
    plt.ylabel('FRET stoichiometry')
    plt.title(r'$\alpha$-$\delta$ corrected FRET histogram')
    
    return peak_coord


def multi_Gaussian(x, gauss_param, nb_gauss):
    g = np.zeros(len(x))
    for i in range(nb_gauss):
        g = g + gauss_param[3*i]*np.exp( - ((x-float(gauss_param[3*i+1]))**2/(2*gauss_param[3*i+2]**2)))
    return g

def multi_error_func(gauss_param, nb_gauss, x, y):
    return multi_Gaussian(x, gauss_param, nb_gauss) - y

def multi_fit_Gauss_spot(hist_1D, nb_gauss = 2, ax_lim_eps = 0.2):
    x_vec = hist_1D[0,:]
    y_vec = hist_1D[1,:]
    gauss_param_estim = []
    Fret_eff_estim = np.linspace(0.3, 0.7, nb_gauss)
    for i in range(nb_gauss):
        gauss_param_estim.append(np.max(y_vec))
        gauss_param_estim.append(Fret_eff_estim[i])
        gauss_param_estim.append(0.2)
    initial_guess = tuple(gauss_param_estim)
    try:
        output = least_squares(multi_error_func, x0=initial_guess, jac='2-point', bounds=(0, np.inf), method='trf', ftol=1e-08,
                       xtol=1e-08, gtol=1e-08, x_scale=1.0, loss='linear', f_scale=1.0, diff_step=None,
                       tr_solver=None,
                       tr_options={}, jac_sparsity=None, max_nfev=None, verbose=0, args=(nb_gauss, x_vec, y_vec))
        
    except RuntimeError:
        print('No peak found')
        return
    
    x_fit = np.linspace(-ax_lim_eps, 1+ax_lim_eps, 1000)

    multi_gauss_fit = multi_Gaussian(x_fit, output.x, nb_gauss)
    
    fig, ax = plt.subplots()
    ax.plot(x_vec, y_vec, 'kx', label='data')
    ax.plot(x_fit, multi_gauss_fit, 'b-', label='Mixed Gaussian model')
    
    individual_fit = []
    for k in range(nb_gauss):
        individual_fit.append(multi_Gaussian(x_fit, output.x[3*k:3*k+3], 1))
        label_k = r'$\mu$ = ' + str(np.round(output.x[3*k+1],3)) + r', $\sigma$ = '+ str(np.round(output.x[3*k+2], 3))
        ax.plot(x_fit, individual_fit[k], '--', label=label_k)
    plt.title(str(nb_gauss)+' Gaussian(s) model')
    plt.xlabel('FRET corrected Efficiency')
    plt.ylabel('Counts')
    ax.legend()
    
    FRET_peak_amp = []
    FRET_peak_pos = []
    FRET_peak_std = []
    
    for j in range(nb_gauss):
        FRET_peak_amp.append(output.x[3*j])
        FRET_peak_pos.append(output.x[3*j+1])
        FRET_peak_std.append(output.x[3*j+2])
    
    df_FRET_peaks = DF(data={'FRET_peak_amplitude': FRET_peak_amp, 'FRET_peak_position': FRET_peak_pos, 'FRET_peak_STD': FRET_peak_std})
    
    return df_FRET_peaks


def plot_detection_modif_BM(
        image,
        spots,
        shape="circle",
        factor_con_max = 0.5,
        factor_con_min = 1,
        radius=3,
        color="red",
        linewidth=1,
        fill=False,
        rescale=False,
        contrast=False,
        title=None,
        framesize=(15, 10),
        remove_frame=True,
        path_output=None,
        ext="png",
        show=True):
    """Plot detected spots and foci on a 2-d image.

    Parameters
    ----------
    image : np.ndarray
        A 2-d image with shape (y, x).
    spots : list or np.ndarray
        Array with coordinates and shape (nb_spots, 3) or (nb_spots, 2). To
        plot different kind of detected spots with different symbols, use a
        list of arrays.
    shape : list or str, default='circle'
        List of symbols used to localized the detected spots in the image,
        among `circle`, `square` or `polygon`. One symbol per array in `spots`.
        If `shape` is a string, the same symbol is used for every elements of
        'spots'.
    radius : list or int or float, default=3
        List of yx radii of the detected spots, in pixel. One radius per array
        in `spots`. If `radius` is a scalar, the same value is applied for
        every elements of `spots`.
    color : list or str, default='red'
        List of colors of the detected spots. One color per array in `spots`.
        If `color` is a string, the same color is applied for every elements
        of `spots`.
    linewidth : list or int, default=1
        List of widths or width of the border symbol. One integer per array
        in `spots`. If `linewidth` is an integer, the same width is applied
        for every elements of `spots`.
    fill : list or bool, default=False
        List of boolean to fill the symbol of the detected spots. If `fill` is
        a boolean, it is applied for every symbols.
    rescale : bool, default=False
        Rescale pixel values of the image (made by default in matplotlib).
    contrast : bool, default=False
        Contrast image.
    title : str, optional
        Title of the image.
    framesize : tuple, default=(15, 10)
        Size of the frame used to plot with ``plt.figure(figsize=framesize)``.
    remove_frame : bool, default=True
        Remove axes and frame.
    path_output : str, optional
        Path to save the image (without extension).
    ext : str or list, default='png'
        Extension used to save the plot. If it is a list of strings, the plot
        will be saved several times.
    show : bool, default=True
        Show the figure or not.

    """
    # enlist and format parameters
    if not isinstance(spots, list):
        spots = [spots]

    # plot
    fig, ax = plt.subplots(1, 1, figsize=framesize)

    ax.imshow(image, vmin=np.min(image)*factor_con_min, vmax=np.max(image)*factor_con_max)

    for i, coordinates in enumerate(spots):

        # get 2-d coordinates
        if coordinates.shape[1] == 3:
            coordinates_2d = coordinates[:, 1:]
        else:
            coordinates_2d = coordinates

        # plot symbols
        for y, x in coordinates_2d:
            x = _define_patch(
                x, y, shape, radius, color, linewidth, fill)
            ax.add_patch(x)

    # titles and frames
    if title is not None:
        ax.set_title("Detection results", fontweight="bold", fontsize=10)
    if remove_frame:
        ax.axis("off")
    plt.tight_layout()

    plt.show()

def _define_patch(x, y, shape, radius, color, linewidth, fill):
    """Define a matplotlib.patches to plot.

    Parameters
    ----------
    x : int or float
        Coordinate x for the patch center.
    y : int or float
        Coordinate y for the patch center.
    shape : str
        Shape of the patch to define (among `circle`, `square` or `polygon`)
    radius : int or float
        Radius of the patch.
    color : str
        Color of the patch.
    linewidth : int
        Width of the patch border.
    fill : bool
        Make the patch shape empty or not.

    Returns
    -------
    x : matplotlib.patches object
        Geometric form to add to a plot.

    """
    # circle
    x = plt.Circle(
        (x, y),
        radius,
        color=color,
        linewidth=linewidth,
        fill=fill)
    return x

def plot_results_cluster_filter(image, coord_spot_before, coord_spot_after,
                                factor_con_max = 0.5,
                                factor_con_min = 1,
                                shape="circle",
                                radius=3,
                                color="red",
                                clusterDisp=1,
                                colorCluster="black",
                                linewidth=1,
                                fill=False,
                                rescale=False,
                                contrast=False,
                                title=1,
                                framesize=(15, 10),
                                remove_frame=True,
                                path_output=None,
                                ext="png",
                                show=True):
    
    fig, ax = plt.subplots(1, 2, figsize=framesize)
    ax[0].imshow(image, vmin=np.min(image)*factor_con_min, vmax=np.max(image)*factor_con_max)
    ax[1].imshow(image, vmin=np.min(image)*factor_con_min, vmax=np.max(image)*factor_con_max)
    
    
    for coord_i in coord_spot_before:
        if coord_i in coord_spot_after:
            y, x = coord_i
            x1 = _define_patch(
                x, y, shape, radius, color, linewidth, fill)
            ax[0].add_patch(x1)
            x2 = _define_patch(
                x, y, shape, radius, color, linewidth, fill)
            ax[1].add_patch(x2)
        else:
            y, x = coord_i
            x = _define_patch(
                x, y, shape, radius, colorCluster, linewidth, fill)
            ax[0].add_patch(x)


    # titles and frames
    if title is not None:
        ax[0].set_title("Before clusters filtering", fontweight="bold", fontsize=10)
        ax[1].set_title("After clusters filtering", fontweight="bold", fontsize=10)
    if remove_frame:
        ax[0].axis("off")
        ax[1].axis("off")
    plt.tight_layout()

    plt.show()
    
def export_OpenFRET_to_CSV_traces(traces_dict, select_key, path):
    
    nb_traces = len(traces_dict.traces)
    
    #list_keys = []
    
    if not os.path.exists(path+'\\Exported_traces_CSV\\' + select_key):
        os.makedirs(path+'\\Exported_traces_CSV\\' + select_key)
    
        
    for i in range(nb_traces):
        if (select_key in list(traces_dict.traces[i].metadata.keys())) or (select_key == 'Full'):
            
            if select_key == 'Full':
                indmin = 0
                indmax = len(traces_dict.traces[i].channels[0].data) - 1
            else:
                indmin = traces_dict.traces[i].metadata[select_key]['indmin']
                indmax = traces_dict.traces[i].metadata[select_key]['indmax']
            
            method_back_corr = traces_dict.traces[i].metadata['background_correction']
            
            DD_trace = np.array(traces_dict.traces[i].channels[0].data[indmin:indmax])
            DA_trace = np.array(traces_dict.traces[i].channels[1].data[indmin:indmax])
            if traces_dict.metadata['ALEX'] == 'yes':
                AA_trace = np.array(traces_dict.traces[i].channels[2].data[indmin:indmax])
                back_DD_trace = np.array(traces_dict.traces[i].channels[3].data[indmin:indmax])
                back_DA_trace = np.array(traces_dict.traces[i].channels[4].data[indmin:indmax])
                back_AA_trace = np.array(traces_dict.traces[i].channels[5].data[indmin:indmax])
            else:
                back_DD_trace = np.array(traces_dict.traces[i].channels[2].data[indmin:indmax])
                back_DA_trace = np.array(traces_dict.traces[i].channels[3].data[indmin:indmax])
            
            match method_back_corr:
                case 'None':
                    back_DD_trace = np.zeros(len(DD_trace))
                    back_DA_trace = np.zeros(len(DA_trace))
                    if traces_dict.metadata['ALEX'] == 'yes':
                        back_AA_trace = np.zeros(len(AA_trace))
                    
                case 'Total variation':
                    max_back_DD = np.max(back_DD_trace)
                    norm_back_trace_DD = np.array(back_DD_trace) / max_back_DD
                    TV_back_DD = tvd_2013(norm_back_trace_DD.astype('float'), traces_dict.traces[i].metadata['TV_lambda'])
                    back_DD_trace = TV_back_DD * max_back_DD
                    
                    max_back_DA = np.max(back_DA_trace)
                    norm_back_trace_DA = np.array(back_DA_trace) / max_back_DA
                    TV_back_DA = tvd_2013(norm_back_trace_DA.astype('float'), traces_dict.traces[i].metadata['TV_lambda'])
                    back_DA_trace = TV_back_DA * max_back_DA
                    
                    if traces_dict.metadata['ALEX'] == 'yes':
                        max_back_AA = np.max(back_AA_trace)
                        norm_back_trace_AA = np.array(back_AA_trace) / max_back_AA
                        TV_back_AA = tvd_2013(norm_back_trace_AA.astype('float'), traces_dict.traces[i].metadata['TV_lambda'])
                        back_AA_trace = TV_back_AA * max_back_AA
                    
                case 'Min. of TV':
                    max_back_DD = np.max(back_DD_trace)
                    norm_back_trace_DD = np.array(back_DD_trace) / max_back_DD
                    TV_back_DD = tvd_2013(norm_back_trace_DD.astype('float'), traces_dict.traces[i].metadata['TV_lambda'])
                    back_DD_trace = np.min(TV_back_DD) * max_back_DD
                    
                    max_back_DA = np.max(back_DA_trace)
                    norm_back_trace_DA = np.array(back_DA_trace) / max_back_DA
                    TV_back_DA = tvd_2013(norm_back_trace_DA.astype('float'), traces_dict.traces[i].metadata['TV_lambda'])
                    back_DA_trace = np.min(TV_back_DA) * max_back_DA
                    
                    if traces_dict.metadata['ALEX'] == 'yes':
                        max_back_AA = np.max(back_AA_trace)
                        norm_back_trace_AA = np.array(back_AA_trace) / max_back_AA
                        TV_back_AA = tvd_2013(norm_back_trace_AA.astype('float'), traces_dict.traces[i].metadata['TV_lambda'])
                        back_AA_trace = np.min(TV_back_AA) * max_back_AA
            corr_DD = DD_trace - back_DD_trace
            corr_DA = DA_trace - back_DA_trace
            if traces_dict.metadata['ALEX'] == 'yes':
                corr_AA = AA_trace - back_AA_trace
            
            if traces_dict.metadata['ALEX'] == 'yes':
                df_i = np.round(DF(data = {'DD': corr_DD, 'DA': corr_DA, 'AA':corr_AA}), 3)
            else:
                df_i = np.round(DF(data = {'DD': corr_DD, 'DA': corr_DA}), 3)
            
            uuid_trace_i = traces_dict.traces[i].metadata['UUID_v7']
            
            print(uuid_trace_i)
            
            df_i.to_csv(path+'\\Exported_traces_CSV\\' + select_key + '\\' + uuid_trace_i + '.csv', index = False)
    
    shutil.make_archive(path+'\\Exported_traces_CSV\\' + select_key, 'zip', path+'\\Exported_traces_CSV\\' + select_key)
    if os.path.exists(path+'\\Exported_traces_CSV\\' + select_key):
        shutil.rmtree(path+'\\Exported_traces_CSV\\' + select_key)
    
    print('Traces successfully exported to: ' + path+'\\Exported_traces_CSV\\' + select_key + '.zip')
    
def merge_OpenFRET_datasets(files_path):
    
    nb_datasets = len(files_path)
    
    if nb_datasets == 1:
        print('Only one dataset selected')
        return
    else:
        merge_datasets =  copy.deepcopy(read_data(files_path[0]))
        for i in range(1,nb_datasets):
            dataset_i = read_data(files_path[i])
            merge_datasets.traces = merge_datasets.traces + dataset_i.traces
        
        merge_folder_path = os.path.dirname(files_path[0]) + '\\Merge_dataset'
        if not os.path.exists(merge_folder_path):
            os.makedirs(merge_folder_path)
            
        write_data(merge_datasets, merge_folder_path + '\\merged_datasets.json', compress=True)
        
        print('Datasets successfully merged and saved to: ' + merge_folder_path + '\\merged_datasets.json.zip')
        return
    
    
    
    
    
    
    