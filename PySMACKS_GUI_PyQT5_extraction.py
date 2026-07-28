# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 11:04:58 2026

@author: molcre0000
"""

import sys
from PyQt5.QtWidgets import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PyQt5.QtGui import QFont
from PyQt5.QtCore import Qt, pyqtSignal
from chromatic_aberrations_correction import *
from utils import *
from traces_extractor import *
from drift_correction import *

class QCDriftWindow(QMainWindow):
    submitClicked = pyqtSignal(object)
    def __init__(self, drift_dict):
        super().__init__()
        self.setWindowTitle("QC drift")
        self.resize(1200, 900)
        
        drift_channels = list(drift_dict.keys())
        
        central_widget = QWidget()
        layout1 = QVBoxLayout(central_widget)
        layout_radius = QHBoxLayout()
        label = QLabel("Channel")
        
        channel_drift = QComboBox()
        #self.channel_drift.addItems(drift_channels)
        if len(drift_channels) == 2:
            drift_channels.append('Average')
        drift_channels.append('None')
        channel_drift.addItems(drift_channels)
        layout1.addWidget(channel_drift)
        
        
        
        spin_box_drift_channel = QSpinBox()
        spin_box_drift_channel.setMinimum(1)  # Set minimum value
        spin_box_drift_channel.setMaximum(len(drift_dict[drift_channels[0]])-1)   # Set maximum value
        spin_box_drift_channel.setValue(np.min([10,len(drift_dict[drift_channels[0]])-1]))
        button_plot = QPushButton("Plot")
        button_plot.clicked.connect(lambda: self.drift_plot(drift_dict, channel_drift.currentText(), spin_box_drift_channel.value()))
        layout_radius.addWidget(label)
        layout_radius.addWidget(spin_box_drift_channel)
        layout_radius.addWidget(button_plot)
        layout1.addLayout(layout_radius)
        self.figure_drift_channel = plt.figure()
        self.canvas_drift_channel = FigureCanvas(self.figure_drift_channel)
        toolbar = NavigationToolbar(self.canvas_drift_channel, self)
        layout1.addWidget(toolbar)
        layout1.addWidget(self.canvas_drift_channel)
        #self.central_widget.addLayout(layout1)
        self.setCentralWidget(central_widget)
    
        
        
        # layout1 = QVBoxLayout()
        # self.figure = plt.figure()
        # self.canvas = FigureCanvas(self.figure)
        # self.toolbar = NavigationToolbar(self.canvas, self)
        # layout1.addWidget(self.canvas)
        
        #self.QC_regis_plot(img_stack_acceptor, img_stack_donor, matrix_align, method)
        # label = QLabel("This is the second window!", self)
        # layout = QVBoxLayout()
        # layout.addWidget(label)
        # self.setLayout(layout)
        
        
    def drift_plot(self, df_drift, channel2plot, savgol_value):
        
        if channel2plot == 'Average':
            drift_channels = list(df_drift.keys())
            mean_drift_X = np.mean([df_drift[drift_channels[0]]['x'], df_drift[drift_channels[1]]['x']], axis = 0)
            mean_drift_Y = np.mean([df_drift[drift_channels[0]]['y'], df_drift[drift_channels[1]]['y']], axis = 0)
            self.smooth_drift = pd.DataFrame(data={'y': savgol_filter(mean_drift_Y, savgol_value, 3), 'x': savgol_filter(mean_drift_X, savgol_value, 3)})
        
        elif channel2plot == 'None':
            drift_channels = list(df_drift.keys())
            self.smooth_drift = pd.DataFrame(data={'y': list(np.zeros(len(df_drift[drift_channels[0]]['y']))), 'x': list(np.zeros(len(df_drift[drift_channels[0]]['x'])))})
        
        else:
            self.smooth_drift = pd.DataFrame(data={'y': df_drift[channel2plot]['y'].tolist(), 'x': df_drift[channel2plot]['x'].tolist()})
        
        
        self.figure_drift_channel.clear()
        ax = self.figure_drift_channel.add_subplot(111)
        
        # Show the merge channels
        #plt.figure('Chromatic abberation correction, frame = ' + str(k))
        #plt.figure()
        ax.plot(self.smooth_drift['y'], color='blue')
        ax.plot(self.smooth_drift['x'], color='orange')
        
        ax.set_xlabel('Submovies #')
        ax.set_ylabel('Drift (px)')
        ax.set_title('Stage drift')
        ax.grid()
        #plt.show()
        
        self.canvas_drift_channel.draw()
        
    def emit_drift(self):

        return self.smooth_drift

def create_extraction_tab(self):
    # Second tab
    # tab2 = QWidget()
    # layout2 = QVBoxLayout()
    # label2 = QLabel("This is Tab 2")
    # button2 = QPushButton("Button in Tab 2")
    # layout2.addWidget(label2)
    # layout2.addWidget(button2)
    # tab2.setLayout(layout2)
    
    # self.tabs.addTab(tab2, "Extraction")
    
    self.isSpotQC = 0
    self.isDriftDD = 0
    self.isDriftAA = 0
    self.loadDonor = 0
    self.loadAccept = 0
    
    self.drift_window = None
    self.traces_data = None
    
    
    tab1 = QWidget()
    layout1 = QVBoxLayout()
    
    self.ALEX_checkbox = QCheckBox("ALEX")
    self.ALEX_checkbox.setToolTip('Tick if your data are ALEX-smFRET')
    self.ALEX_checkbox.setChecked(True)
    self.ALEX_checkbox.stateChanged.connect(lambda: toggle_ALEX_buttons(self))
    layout1.addWidget(self.ALEX_checkbox)
    
    layout_acceptor_calib = QHBoxLayout()
    label_acceptor_calib = QLabel("Acceptor raw movies")
    self.button_acceptor_raw_data = QPushButton("1. Browse file")
    self.button_acceptor_raw_data.setToolTip('Select all the acceptor TIFF submovies')
    self.button_acceptor_raw_data.clicked.connect(lambda: open_file_dialog_acceptor_data(self))
    self.button_acceptor_raw_data.setStyleSheet("background-color : yellow")
    layout_acceptor_calib.addWidget(label_acceptor_calib)
    layout_acceptor_calib.addWidget(self.button_acceptor_raw_data)
    layout_donor_calib = QHBoxLayout()
    label_donor_calib = QLabel("Donor raw movies")
    self.button_donor_raw_data = QPushButton("2. Browse file")
    self.button_donor_raw_data.setToolTip('Select all the donor TIFF submovies')
    self.button_donor_raw_data.clicked.connect(lambda: open_file_dialog_donor_data(self))
    self.button_donor_raw_data.setStyleSheet("background-color : yellow")
    layout_donor_calib.addWidget(label_donor_calib)
    layout_donor_calib.addWidget(self.button_donor_raw_data)
    layout1.addLayout(layout_acceptor_calib)
    layout1.addLayout(layout_donor_calib)
    
    # widget to select spot radius
    
    layout_frames_donor_rad = QHBoxLayout()
    label = QLabel("Donor-Donor spot radius")
    layout_frames_donor_rad.addWidget(label)
    spin_box_donor_rad = QSpinBox()
    spin_box_donor_rad.setToolTip('Select the approximate DD spot radius (in pixels)')
    spin_box_donor_rad.setMinimum(1)  # Set minimum value
    spin_box_donor_rad.setMaximum(1000)   # Set maximum value
    spin_box_donor_rad.setValue(2)
    layout_frames_donor_rad.addWidget(spin_box_donor_rad)
    self.button_runQC_donor = QPushButton("Run QC DD")
    self.button_runQC_donor.setToolTip('Run a quick spot detection to visualize the detected DD spots and confirm the spot radius')
    self.button_runQC_donor.clicked.connect(lambda: Run_QC_spot(self, self.filename_donor_data, spin_box_donor_rad.value(), isALEX_DA = 0))
    self.button_runQC_donor.setStyleSheet("background-color : yellow")
    layout_frames_donor_rad.addWidget(self.button_runQC_donor)
    
    self.button_run_drift_DD = QPushButton("3. Run DD Drift")
    self.button_run_drift_DD.setToolTip('Detect all the spots in the DD channel overtime and calculate the drift of the stage based on DD spots displacements')
    self.button_run_drift_DD.clicked.connect(lambda: Run_DD_drift_calc(self, self.filename_donor_data, kernel_size=spin_box_donor_rad.value()))
    self.button_run_drift_DD.setStyleSheet("background-color : yellow")
    layout_frames_donor_rad.addWidget(self.button_run_drift_DD)
    
    layout1.addLayout(layout_frames_donor_rad)
    
    self.frame_AA = QFrame()
    layout_frames_accept_rad = QHBoxLayout()
    self.frame_AA.setLayout(layout_frames_accept_rad)
    label = QLabel("Acceptor-Acceptor spot radius")
    layout_frames_accept_rad.addWidget(label)
    spin_box_accept_rad = QSpinBox()
    spin_box_accept_rad.setToolTip('Select the approximate AA spot radius (in pixels)')
    spin_box_accept_rad.setMinimum(1)  # Set minimum value
    spin_box_accept_rad.setMaximum(1000)   # Set maximum value
    spin_box_accept_rad.setValue(2)
    layout_frames_accept_rad.addWidget(spin_box_accept_rad)
    self.button_runQC_accept = QPushButton("Run QC AA")
    self.button_runQC_accept.setToolTip('Run a quick spot detection to visualize the detected AA spots and confirm the spot radius')
    self.button_runQC_accept.clicked.connect(lambda: Run_QC_spot(self, self.filename_acceptor_data, spin_box_accept_rad.value(), isALEX_DA = 1))
    self.button_runQC_accept.setStyleSheet("background-color : yellow")
    layout_frames_accept_rad.addWidget(self.button_runQC_accept)
    
    self.button_run_drift_AA = QPushButton("4. Run AA Drift")
    self.button_run_drift_AA.setToolTip('Detect all the spots in the AA channel overtime and calculate the drift of the stage based on AA spots displacements')
    self.button_run_drift_AA.clicked.connect(lambda: Run_AA_drift_calc(self, self.filename_acceptor_data, kernel_size=spin_box_accept_rad.value()))
    self.button_run_drift_AA.setStyleSheet("background-color : yellow")
    layout_frames_accept_rad.addWidget(self.button_run_drift_AA)
    
    layout1.addWidget(self.frame_AA)
    #layout1.addLayout(layout_frames_accept_rad)
    
    button_plot_graph_drift = QPushButton("5. Plot drift")
    button_plot_graph_drift.setToolTip('Open a new window to visualize and select the stage drift input data')
    button_plot_graph_drift.clicked.connect(lambda: plot_stage_drift(self))
    layout1.addWidget(button_plot_graph_drift)
    
    layout_frames_DA_rad = QHBoxLayout()
    label = QLabel("Donor-Acceptor spot radius")
    layout_frames_DA_rad.addWidget(label)
    spin_box_DA_rad = QSpinBox()
    spin_box_DA_rad.setToolTip('Select the approximate AA spot radius (in pixels)')
    spin_box_DA_rad.setMinimum(1)  # Set minimum value
    spin_box_DA_rad.setMaximum(1000)   # Set maximum value
    spin_box_DA_rad.setValue(2)
    layout_frames_DA_rad.addWidget(spin_box_DA_rad)
    self.button_runQC_DA = QPushButton("Run QC DA")
    self.button_runQC_DA.setToolTip('Run a quick spot detection to visualize the detected DA spots and confirm the spot radius')
    self.button_runQC_DA.clicked.connect(lambda: Run_QC_spot(self, self.filename_acceptor_data, spin_box_DA_rad.value(), isALEX_DA = 2))
    self.button_runQC_DA.setStyleSheet("background-color : yellow")
    layout_frames_DA_rad.addWidget(self.button_runQC_DA)
    
    self.button_run_drift_DA = QPushButton("6. Run DA spot coordinates extraction")
    self.button_run_drift_DA.setToolTip('Detect all the spots in the DA channel overtime with stage drift correction')
    self.button_run_drift_DA.clicked.connect(lambda: Run_DA_drift_calc(self, self.filename_acceptor_data, spin_box_DA_rad.value()))
    self.button_run_drift_DA.setStyleSheet("background-color : yellow")
    layout_frames_DA_rad.addWidget(self.button_run_drift_DA)
    
    layout1.addLayout(layout_frames_DA_rad)
    
    layout_remove_clusters = QHBoxLayout()
    label = QLabel("Inter-spots min. distance (px)")
    layout_remove_clusters.addWidget(label)
    spin_box_remove_spots = QSpinBox()
    spin_box_remove_spots.setToolTip('Minimal inter-spots distance (in pixels) for cluster filtering')
    spin_box_remove_spots.setMinimum(1)  # Set minimum value
    spin_box_remove_spots.setMaximum(1000)   # Set maximum value
    spin_box_remove_spots.setValue(9)
    layout_remove_clusters.addWidget(spin_box_remove_spots)
    self.button_remove_clusters = QPushButton("7. Remove clusters")
    self.button_remove_clusters.setToolTip('Filter the spots by removing the clusters based on a minimal inter-spots distance (in pixels)')
    self.button_remove_clusters.clicked.connect(lambda: filter_clusters(self, spin_box_remove_spots.value()))
    self.button_remove_clusters.setStyleSheet("background-color : yellow")
    layout_remove_clusters.addWidget(self.button_remove_clusters)
    layout1.addLayout(layout_remove_clusters)
    
    
    layout_traces_extract = QHBoxLayout()
    label = QLabel("sigma (px)")
    layout_traces_extract.addWidget(label)
    spin_box_traces_extract = QSpinBox()
    spin_box_traces_extract.setToolTip('Radius (in pixels) for extraction of signal (1 sigma) and background (2 sigma)')
    spin_box_traces_extract.setMinimum(1)  # Set minimum value
    spin_box_traces_extract.setMaximum(1000)   # Set maximum value
    spin_box_traces_extract.setValue(3)
    layout_traces_extract.addWidget(spin_box_traces_extract)
    self.button_run_traces_extract = QPushButton("8. Run DA traces extraction")
    self.button_run_traces_extract.setToolTip('Run the donor-acceptor traces extraction using the filter spots and the sigma value')
    self.button_run_traces_extract.clicked.connect(lambda: run_traces_extraction(self, spin_box_traces_extract.value()))
    self.button_run_traces_extract.setStyleSheet("background-color : yellow")
    layout_traces_extract.addWidget(self.button_run_traces_extract)
    layout1.addLayout(layout_traces_extract)
    
    layout_extract_DO = QHBoxLayout()
    button_cluster_DO = QPushButton("Remove DO clusters")
    button_cluster_DO.setToolTip('Filter the clusters on donor movies using the fitler spots value')
    button_cluster_DO.clicked.connect(lambda: run_DO_cluster(self, spin_box_remove_spots.value()))
    layout_extract_DO.addWidget(button_cluster_DO)
    self.button_run_DO = QPushButton("Run DO extraction")
    self.button_run_DO.setToolTip('Run traces extraction on spots detected on the donor movie (for Donor-Only traces and alpha FRET correction)')
    self.button_run_DO.clicked.connect(lambda: run_DO_extraction(self, spin_box_traces_extract.value()))
    layout_extract_DO.addWidget(self.button_run_DO)
    layout1.addLayout(layout_extract_DO)
    
    layout_extract_AO = QHBoxLayout()
    button_cluster_AO = QPushButton("Remove AO clusters")
    button_cluster_AO.setToolTip('Filter the clusters on AA movies using the fitler spots value')
    button_cluster_AO.clicked.connect(lambda: run_AO_cluster(self, spin_box_remove_spots.value()))
    layout_extract_AO.addWidget(button_cluster_AO)
    self.button_run_AO = QPushButton("Run AO extraction")
    self.button_run_AO.setToolTip('Run traces extraction on spots detected on the AA movie (for Acceptor-Only traces and delta FRET correction)')
    self.button_run_AO.clicked.connect(lambda: run_AO_extraction(self, spin_box_traces_extract.value()))
    layout_extract_AO.addWidget(self.button_run_AO)
    layout1.addLayout(layout_extract_AO)
    
    self.figure_extract = plt.figure()
    self.canvas_extract = FigureCanvas(self.figure_extract)
    toolbar = NavigationToolbar(self.canvas_extract, self)
    layout1.addWidget(toolbar)
    
    layout_calib_slider1 = QHBoxLayout()
    layout_calib_slider1.addWidget(self.canvas_extract)
    self.sliders = {}
    vbox = QVBoxLayout()
    label = QLabel("Contrast")
    slider = QSlider(Qt.Vertical)
    slider.setMinimum(0)   # Corresponds to 0.5
    slider.setMaximum(150)  # Corresponds to 1.5
    slider.setValue(100)    # Default 1.0
    slider.setTickPosition(QSlider.TicksBelow)
    slider.setTickInterval(10)
    #slider.valueChanged.connect(self.update_contrast)
    slider.valueChanged.connect(lambda: slider_QC_contrast(self))
    self.sliders['0'] = slider

        #vbox.addWidget(label)
    vbox.addWidget(label)
    vbox.addWidget(slider)
    layout_calib_slider1.addLayout(vbox)
    layout1.addLayout(layout_calib_slider1)
    
    tab1.setLayout(layout1)
    self.tabs.addTab(tab1, "Extraction")
    
def toggle_ALEX_buttons(self):
    if self.ALEX_checkbox.isChecked():  # Checked
        self.frame_AA.show()
    else:  # Unchecked
        self.frame_AA.hide()
    
def filter_clusters(self, radius):
    
    self.coord_spots_DA_filter = filter_close_prox_spots(self.coord_spots_track_DA, min_dist = radius)
    
    print(len(self.coord_spots_DA_filter))
    
    self.button_remove_clusters.setStyleSheet("background-color : green")
    
def run_DO_cluster(self, radius):
    self.coord_spots_DO_filter = filter_close_prox_spots(self.coord_spots_track_DD, min_dist = radius)
    
    print(len(self.coord_spots_DO_filter))
    
def run_AO_cluster(self, radius):
    self.coord_spots_AO_filter = filter_close_prox_spots(self.coord_spots_track_AA, min_dist = radius)
    
    print(len(self.coord_spots_AO_filter))
    
    
def run_DO_extraction(self, radius):
    self.traces_DO_data = extract_traces_from_DD_spots(self.filename_donor_data, 
                                        self.filename_acceptor_data, 
                                        self.coord_spots_DO_filter, 
                                        self.matrix_align,
                                        drift_correct = self.drift,
                                        method_align = self.method_calib.currentText(), 
                                        sigma = radius)
    
    print(len(self.traces_DO_data.traces))
    
    self.button_run_DO.setStyleSheet("background-color : green")
    
def run_AO_extraction(self, radius):
    self.traces_AO_data = spot_trace_extractor(self.filename_donor_data, 
                                        self.filename_acceptor_data, 
                                        self.coord_spots_AO_filter, 
                                        self.matrix_align,
                                        drift_correct = self.drift,
                                        method_align = self.method_calib.currentText(), 
                                        sigma = radius)
    
    print(len(self.traces_AO_data.traces))
    
    self.button_run_AO.setStyleSheet("background-color : green")
    
def run_traces_extraction(self, radius):
    
    if self.ALEX_checkbox.isChecked():
        self.traces_data = spot_trace_extractor(self.filename_donor_data, 
                                            self.filename_acceptor_data, 
                                            self.coord_spots_DA_filter, 
                                            self.matrix_align,
                                            drift_correct = self.drift,
                                            method_align = self.method_calib.currentText(), 
                                            sigma = radius)
    else:
        self.traces_data = spot_trace_extractor_no_ALEX(self.filename_donor_data, 
                                            self.filename_acceptor_data, 
                                            self.coord_spots_DA_filter, 
                                            self.matrix_align,
                                            drift_correct = self.drift,
                                            method_align = self.method_calib.currentText(), 
                                            sigma = radius)
    
    #self.traces_data = get_high_signal(self.traces_data)
    
    print(len(self.traces_data.traces))
    
    self.trace_IDs = [str(i) for i in range(len(self.traces_data.traces))] #list(self.traces_data.keys())
    self.choose_trace_ID.addItems(self.trace_IDs)
    
    self.button_run_traces_extract.setStyleSheet("background-color : green")
    
def plot_stage_drift(self):
    
    drift_dict = {}
    
    if self.isDriftDD == 1:
        drift_dict['DD'] = self.drift_DD
    if self.isDriftAA == 1:
        drift_dict['AA'] = self.drift_AA
    
    self.drift_window = QCDriftWindow(drift_dict)
    self.drift_window.show()

    
def Run_DD_drift_calc(self, file_path, kernel_size=2, min_distance=2, ref_AA = 0, perc_frames = 0.05, get_coord = 1, Z_project = 1):
    
    # root = Tk(className='Open TIFF movies', )
    # file_path = askopenfilenames(title="Select the donor TIFF movies")
    # root.destroy()
    
    time_frames = 0
    
    comp_frame_drift = 0
    
    self.button_run_drift_DD.setStyleSheet("background-color : yellow")
    
    for k in range(len(file_path)):
        print('Processing movie ' + str(k))
        img = Image.open(file_path[k])
        nb_frame = img.n_frames
        if k < len(file_path) - 1:
            nb_frame_old = img.n_frames
        if Z_project == 0:
            list_frames = [int(i) for i in np.floor(np.linspace(0,nb_frame-1,int(nb_frame*perc_frames)))]
            for j in list_frames:
                #print('Processing frame ' + str(j))
                img.seek(j)
                img_raw = np.array(img)
            #img_back = median_filter(img_raw, size = size_median_filt)
            #img_med_filter = img_raw - img_back
                coord_spots = detection.detect_spots(img_raw, log_kernel_size=kernel_size, minimum_distance=min_distance)
                if (k+j) == 0:
                    df_coords = pd.DataFrame(data={'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots)), 'frame':np.zeros(len(coord_spots))})
                    time_frames = nb_frame_old*k + j
                    comp_frame_drift = comp_frame_drift + 1
                else:
                    #df_coords = df_coords._append({'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame':np.zeros(len(coord_spots))+time_frames}, ignore_index=True)
                    
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
        #time_frames = time_frames + nb_frame
        
    df_linked = tp.link(df_coords, 2, memory=time_frames)  # 5
    
    
    list_frames_drift = df_linked['frame_real'].unique()
    self.drift_DD = tp.compute_drift(df_linked)
    
    spot_IDs = df_linked['particle'].unique()
    
    self.coord_spots_track_DD = np.zeros((len(spot_IDs), 2))
    
    if get_coord == 1:
    
        for j in spot_IDs:
            df_i = df_linked[df_linked['particle']==j]
            if len(df_i) >  2:
                min_frame_i = np.min(df_i['frame'])
                if min_frame_i == 0:
                    drift_corr_i = [0, 0]
                else:
                    #drift_corr_i = np.round([drift['y'][min_frame_i], drift['x'][min_frame_i]])
                    drift_corr_i = [self.drift_DD['y'][min_frame_i], self.drift_DD['x'][min_frame_i]]
                    
                y_spot = int(df_i[df_i['frame']==min_frame_i]['y'].tolist()[0] - drift_corr_i[0])
                x_spot = int(df_i[df_i['frame']==min_frame_i]['x'].tolist()[0] - drift_corr_i[1])
                
                self.coord_spots_track_DD[j, 0] = y_spot
                self.coord_spots_track_DD[j, 1] = x_spot
            
        mask = ~np.any(self.coord_spots_track_DD <= 0, axis=1)
        self.coord_spots_track_DD = self.coord_spots_track_DD[mask]
        
        self.button_run_drift_DD.setStyleSheet("background-color : green")
        
        self.isDriftDD = 1
    

    
def Run_AA_drift_calc(self, file_path, kernel_size=2, min_distance=2, ref_AA = 0, perc_frames = 0.05, get_coord = 1, Z_project = 1):
    
    time_frames = 0
    
    comp_frame_drift = 0
    
    self.button_run_drift_AA.setStyleSheet("background-color : yellow")
    
    for k in range(len(file_path)):
        print('Processing movie ' + str(k))
        img = Image.open(file_path[k])
        nb_frame = img.n_frames
        if k < len(file_path) - 1:
            nb_frame_old = img.n_frames
        if Z_project == 0:
            list_frames = [int(i)*2+ref_AA for i in np.floor(np.linspace(0,int(nb_frame/2)-1,int(nb_frame/2*perc_frames)))]
            for j in list_frames:
                #print('Processing frame ' + str(j))
                img.seek(j)
                img_raw = np.array(img)
            #img_back = median_filter(img_raw, size = size_median_filt)
            #img_med_filter = img_raw - img_back
                coord_spots = detection.detect_spots(img_raw, log_kernel_size=kernel_size, minimum_distance=min_distance)
                if (k+j) == 0:
                    df_coords = pd.DataFrame(data={'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots)), 'frame':np.zeros(len(coord_spots))})
                    time_frames = nb_frame_old*k + j
                    comp_frame_drift = comp_frame_drift + 1
                else:
                    #df_coords = df_coords._append({'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame':np.zeros(len(coord_spots))+time_frames}, ignore_index=True)
                    
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
                
        #time_frames = time_frames + nb_frame
        
    df_linked = tp.link(df_coords, 2, memory=time_frames)  # 5
    
    
    list_frames_drift = df_linked['frame_real'].unique()
    self.drift_AA = tp.compute_drift(df_linked)
    
    spot_IDs = df_linked['particle'].unique()
    
    self.coord_spots_track_AA = np.zeros((len(spot_IDs), 2))
    
    
    if get_coord == 1:
        for j in spot_IDs:
            df_i = df_linked[df_linked['particle']==j]
            if len(df_i) > 2:
                min_frame_i = np.min(df_i['frame'])
                if min_frame_i == 0:
                    drift_corr_i = [0, 0]
                else:
                    drift_corr_i = [self.drift_AA['y'][min_frame_i], self.drift_AA['x'][min_frame_i]]
                    
                y_spot = int(df_i[df_i['frame']==min_frame_i]['y'].tolist()[0] - drift_corr_i[0])
                x_spot = int(df_i[df_i['frame']==min_frame_i]['x'].tolist()[0] - drift_corr_i[1])
                
                self.coord_spots_track_AA[j, 0] = y_spot
                self.coord_spots_track_AA[j, 1] = x_spot
            
        mask = ~np.any(self.coord_spots_track_AA <= 0, axis=1)
        self.coord_spots_track_AA = self.coord_spots_track_AA[mask]
        
    self.button_run_drift_AA.setStyleSheet("background-color : green")
    
    self.isDriftAA = 1
        
def Run_DA_drift_calc(self, file_path, kernel_size, min_distance=2, ref_AA = 1, perc_frames = 0.05, Z_project = 1): # drift
    
    self.drift = self.drift_window.emit_drift()
    
    
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
                #print('Processing frame ' + str(j))
                img.seek(j)
                img_raw = np.array(img)
            #img_back = median_filter(img_raw, size = size_median_filt)
            #img_med_filter = img_raw - img_back
                coord_spots = detection.detect_spots(img_raw, log_kernel_size=kernel_size, minimum_distance=min_distance)
                # coord_spots = filter_close_prox_spots(coord_spots, min_dist = 7)
                #print(len(coord_spots))
                if (k+j-ref_AA) == 0:
                    df_coords = pd.DataFrame(data={'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots)), 'frame':np.zeros(len(coord_spots))})
                    time_frames = nb_frame_old*k + j
                    comp_frame_drift = comp_frame_drift + 1
                else:
                    #df_coords = df_coords._append({'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame':np.zeros(len(coord_spots))+time_frames}, ignore_index=True)
                    
                    time_frames = nb_frame_old*k + j
                    new_coords = {'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots))+time_frames, 'frame':np.zeros(len(coord_spots)) + comp_frame_drift}
                    df_coords = pd.concat([df_coords, pd.DataFrame(data=new_coords)]).reset_index(drop=True)
                    comp_frame_drift = comp_frame_drift + 1
        else:
            if self.ALEX_checkbox.isChecked():
                img.seek(1)
                img0 = np.array(img)
                dim_img = img0.shape
                img_raw = np.zeros(((dim_img[0], dim_img[1], int(nb_frame/2))))
                img_raw[:,:,0] = img0
                for j in range(1,int(nb_frame/2)):
                    img.seek(2*j+1)
                    img_raw[:,:,j] = np.array(img)
            else:
                img.seek(0)
                img0 = np.array(img)
                dim_img = img0.shape
                img_raw = np.zeros((dim_img[0], dim_img[1], int(nb_frame)))
                img_raw[:,:,0] = img0
                for j in range(1,int(nb_frame)):
                    img.seek(j)
                    img_raw[:,:,j] = np.array(img)
                
            img_raw_projZ = np.sum(img_raw, axis = 2)
            coord_spots = detection.detect_spots(img_raw_projZ, log_kernel_size=kernel_size, minimum_distance=min_distance)
            
            if k ==0:
                df_coords = pd.DataFrame(data={'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots)), 'frame':np.zeros(len(coord_spots))})
                if self.ALEX_checkbox.isChecked():
                    time_frames = time_frames + int(nb_frame/2)
                else:
                    time_frames = time_frames + int(nb_frame)
            else:
                new_coords = {'y': coord_spots[:,0].tolist(), 'x': coord_spots[:,1].tolist(), 'frame_real':np.zeros(len(coord_spots))+time_frames, 'frame':np.zeros(len(coord_spots)) + k}
                df_coords = pd.concat([df_coords, pd.DataFrame(data=new_coords)]).reset_index(drop=True)
                if self.ALEX_checkbox.isChecked():
                    time_frames = time_frames + int(nb_frame/2)
                else:
                    time_frames = time_frames + int(nb_frame)
                
        #time_frames = time_frames + nb_frame
        
    df_linked = tp.link(df_coords, 3, memory=time_frames, adaptive_stop = 0.5, adaptive_step = 0.95)  # 5
    
    
    list_frames_drift = df_linked['frame_real'].unique()
    
    spot_IDs = df_linked['particle'].unique()
    
    self.coord_spots_track_DA = np.zeros((len(spot_IDs), 2))
    
    for j in spot_IDs:
        df_i = df_linked[df_linked['particle']==j]
        if len(df_i) > 2:
            min_frame_i = np.min(df_i['frame'])
            if min_frame_i == 0:
                drift_corr_i = [0, 0]
            else:
                #drift_corr_i = np.round([drift['y'][min_frame_i], drift['x'][min_frame_i]])
                drift_corr_i = [self.drift['y'][min_frame_i-1], self.drift['x'][min_frame_i-1]]
                
            y_spot = int(df_i[df_i['frame']==min_frame_i]['y'].tolist()[0] - drift_corr_i[0])
            x_spot = int(df_i[df_i['frame']==min_frame_i]['x'].tolist()[0] - drift_corr_i[1])
            
            self.coord_spots_track_DA[j, 0] = y_spot
            self.coord_spots_track_DA[j, 1] = x_spot
        
    mask = ~np.any(self.coord_spots_track_DA <= 0, axis=1)
    self.coord_spots_track_DA = self.coord_spots_track_DA[mask]
    
    self.button_run_drift_DA.setStyleSheet("background-color : green")
    
def slider_QC_contrast(self):
    
    # img_copy = np.array(list(self.image))
    # max_val = np.max(img_copy)
    # img_copy = img_copy / (self.sliders['0'].value()/100)
    # img_copy[img_copy > max_val] = max_val
  
    # self.img_spot_QC.set_data(img_copy)
    
    # self.canvas_extract.draw()
    
    self.img_spot_QC.set_clim(vmin=0, vmax=self.image.max()*self.sliders['0'].value()/100)
    
    self.canvas_extract.draw_idle()
    

def open_file_dialog_acceptor_data(self):
    # Open file dialog and get selected file path
    options = QFileDialog.Options()
    # You can set options like QFileDialog.DontUseNativeDialog if needed
    file_name, _ = QFileDialog.getOpenFileNames(self, "Select the Acceptor raw movies", "", "All Files (*);;Text Files (*.txt)", options=options)
    if file_name:
        self.filename_acceptor_data = file_name
        self.loadAccept = 1
    if (self.loadDonor == 1):
        if len(self.filename_acceptor_data) - len(self.filename_donor_data) != 0:
            self.button_donor_raw_data.setStyleSheet("background-color : red")
            self.button_acceptor_raw_data.setStyleSheet("background-color : red")
        else:
            self.button_acceptor_raw_data.setStyleSheet("background-color : green")
            self.button_donor_raw_data.setStyleSheet("background-color : green")
    else:
        self.button_acceptor_raw_data.setStyleSheet("background-color : green")
    
        
def open_file_dialog_donor_data(self):
    # Open file dialog and get selected file path
    options = QFileDialog.Options()
    # You can set options like QFileDialog.DontUseNativeDialog if needed
    file_name, _ = QFileDialog.getOpenFileNames(self, "Select the Donor raw movies", "", "All Files (*);;Text Files (*.txt)", options=options)
    if file_name:
        self.filename_donor_data = file_name
        self.loadDonor = 1
    if (self.loadAccept == 1):
        if len(self.filename_acceptor_data) - len(self.filename_donor_data) != 0:
            self.button_donor_raw_data.setStyleSheet("background-color : red")
            self.button_acceptor_raw_data.setStyleSheet("background-color : red")
        else:
            self.button_donor_raw_data.setStyleSheet("background-color : green")
            self.button_acceptor_raw_data.setStyleSheet("background-color : green")
    else:
        self.button_donor_raw_data.setStyleSheet("background-color : green")
        
def Run_QC_spot(self, filenames, radius, isALEX_DA = 0):
    img_stack = load_submovies(filenames[0])
    if isALEX_DA == 1: # AA channel
        img_Z_proj = np.sum(img_stack[:,:,0::2], axis = 2)
    elif isALEX_DA == 2:
        img_Z_proj = np.sum(img_stack[:,:,1::2], axis = 2)
    else:
        img_Z_proj = np.sum(img_stack, axis = 2)
    coord_spots_QC = detection.detect_spots(img_Z_proj, log_kernel_size=radius, minimum_distance=2)
    plot_QC_detection(self, img_Z_proj, coord_spots_QC)
    
    if isALEX_DA == 1: # AA channel
        self.button_runQC_accept.setStyleSheet("background-color : green")
    elif isALEX_DA == 2:
        self.button_runQC_DA.setStyleSheet("background-color : green")
    else:
        self.button_runQC_donor.setStyleSheet("background-color : green")
    
    

def plot_QC_detection(
        self,
        image,
        spots,
        shape="circle",
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
    self.image = image
    
    img_copy = np.array(list(image))
    max_val = np.max(img_copy)
    img_copy = img_copy / (self.sliders['0'].value()/100)
    img_copy[img_copy > max_val] = max_val
    
    if self.isSpotQC == 0:
        self.figure_extract.clear()
        self.ax_patch_spots = self.figure_extract.add_subplot(111)
        self.img_spot_QC = self.ax_patch_spots.imshow(img_copy, vmin=np.min(img_copy), vmax=np.max(img_copy))
        self.isSpotQC = 1
    else:
        [p.remove() for p in reversed(self.ax_patch_spots.patches)]
        self.img_spot_QC.set_data(img_copy)
    #self.img_spot_QC.set_data(image)
    #self.canvas.draw()
    
    #fig, ax = plt.subplots(1, 1, figsize=framesize)

    #ax.imshow(image, vmin=np.min(image)*factor_con_min, vmax=np.max(image)*factor_con_max)

    for i, coordinates in enumerate(spots):

        # get 2-d coordinates
        if coordinates.shape[1] == 3:
            coordinates_2d = coordinates[:, 1:]
        else:
            coordinates_2d = coordinates

        # plot symbols
        for y, x in coordinates_2d:
            x = _define_patch_QC_GUI(
                x, y, shape, radius, color, linewidth, fill)
            #ax.add_patch(x)
            self.ax_patch_spots.add_patch(x)
    self.canvas_extract.draw()

    # titles and frames
    # if title is not None:
    #     ax.set_title("Detection results", fontweight="bold", fontsize=10)
    # if remove_frame:
    #     ax.axis("off")
    #plt.tight_layout()

    #plt.show()

def _define_patch_QC_GUI(x, y, shape, radius, color, linewidth, fill):
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