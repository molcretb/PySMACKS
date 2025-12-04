# -*- coding: utf-8 -*-
"""
Created on Thu Dec  4 15:18:43 2025

@author: Bastien Molcrette, Schmid lab, Chemistry department, University of Basel
"""
from chromatic_aberrations_correction import *
from utils import *
from traces_extractor import *


def main():
    # get the transformation matrix for the chromatic aberrations correction
    matrix_align = pipeline_chrom_ab_correction(nb_ref_frames = 5)
    # load donor stack
    img_stack_D, file_path_D = load_data()
    # generate the calibrated donor movie
    img_DD = generate_chrom_ab_corr_movie(img_stack_D, matrix_align)
    # load acceptor stack
    img_stack_A, file_path_A = load_data()
    # desinterleave the acceptor raw channel
    img_DA, img_AA = deinterleave_acceptor_channel(img_stack_A)
    # detect the spots
    coord_spots = detect_spot(img_DD, img_DA)
    #filter the clustered spots
    coord_spots_filter = filter_close_prox_spots(coord_spots, min_dist = 7)
    # extract DD and DA traces
    traces_data = spot_trace_extractor(coord_spots_filter, img_DD, file_path_D, img_stack_A, file_path_A, matrix_align)
    # detect bleaching event
    traces_data = detect_bleaching_event(traces_data)
    # calculate SNR for each traces
    traces_data = get_SNR(traces_data)
    # filter only high SNR traces
    traces_high_SNR, list_highSNR_IDs = filter_traces_on_SNR(traces_data, SNR_thresh = 8)
    return traces_high_SNR, list_highSNR_IDs, traces_data, matrix_align
    

if __name__ == "__main__":
    traces_high_SNR, list_highSNR_IDs, traces_data, matrix_align = main()