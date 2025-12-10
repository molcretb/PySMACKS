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
    file_save = asksaveasfile(title="Save traces dataset", defaultextension=".json")
    file_save.close()
    root.destroy()
    
    # The traces dictionary is then saved at the path location provided by the user
    with open(file_save.name, "w") as fp:
        dump(traces_dict, fp, cls=NpEncoder, separators=(',', ':')) 
    print('Traces dataset saved as JSON file to path: ' + file_save.name)
    return

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

