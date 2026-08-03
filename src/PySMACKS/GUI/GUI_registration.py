# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 10:23:32 2026

@author: molcre0000
"""
import sys
from PyQt5.QtWidgets import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PyQt5.QtGui import QFont
from PyQt5.QtCore import Qt
from PySMACKS.components.chromatic_aberrations_correction import *
from PySMACKS.components.utils import *
from PySMACKS.components.traces_extractor import *
from PySMACKS.components.drift_correction import *

class QCRegistrationWindow(QMainWindow):
    def __init__(self, img_stack_acceptor, img_stack_donor, matrix_align, method, count_windows):
        super().__init__()
        self.setWindowTitle("QC registration #" + str(count_windows))
        self.resize(700, 700)
        
        
        central_widget = QWidget()
        layout1 = QVBoxLayout(central_widget)
        layout_radius = QHBoxLayout()
        label = QLabel("Spot radius")
        spin_box_spot_rad = QSpinBox()
        spin_box_spot_rad.setMinimum(1)  # Set minimum value
        spin_box_spot_rad.setMaximum(20)   # Set maximum value
        spin_box_spot_rad.setValue(3)
        button_plot = QPushButton("Plot")
        button_plot.clicked.connect(lambda: self.QC_regis_plot(img_stack_acceptor, img_stack_donor, matrix_align, method, spin_box_spot_rad.value()))
        layout_radius.addWidget(label)
        layout_radius.addWidget(spin_box_spot_rad)
        layout_radius.addWidget(button_plot)
        layout1.addLayout(layout_radius)
        self.figure_QC_calib = plt.figure()
        self.canvas_QC_calib = FigureCanvas(self.figure_QC_calib)
        toolbar_QC_calib = NavigationToolbar(self.canvas_QC_calib, self)
        layout1.addWidget(toolbar_QC_calib)
        layout1.addWidget(self.canvas_QC_calib)
        self.setCentralWidget(central_widget)
        
        self.QC_regis_plot(img_stack_acceptor, img_stack_donor, matrix_align, method, spin_box_spot_rad.value())

        
    def QC_regis_plot(self, img_stack_acceptor, img_stack_donor, matrix_align, method, radius):

        peaks_A = detection.detect_spots(img_stack_acceptor[:,:,0], log_kernel_size=radius, minimum_distance=1)
        
        if method == 'Optical Flow':
            img_warp_donor = Warp_OpticalFlow(img_stack_donor, matrix_align[0], matrix_align[1])
        else:
            img_warp_donor = generate_chrom_ab_corr_movie(img_stack_donor, matrix_align)
        
        peaks_D = detection.detect_spots(img_warp_donor[:,:,0], log_kernel_size=radius, minimum_distance=1)
        
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
        
        
        self.figure_QC_calib.clear()
        ax = self.figure_QC_calib.add_subplot(111)
        
        ax.plot(data, cdf, color='blue')
        
        ax.set_xlabel('inter-peaks distance (px)')
        ax.set_ylabel('CDF')
        ax.set_title('Cumulative distribution of calibration inter-peaks distance')
        ax.grid()
        
        self.canvas_QC_calib.draw_idle()




def run_registration(self, nb_frame_OF):
    nb_ref_frames = nb_frame_OF
    method = self.method_calib.currentText()
    radius_OF = 10
    method_min_opt='Nelder-Mead'
    maxfev_value_opt = 5000
    path_acceptor_calib = self.widget_path_acceptor_calib.text()
    path_donor_calib = self.widget_path_donor_calib.text()
    print(path_acceptor_calib)
    
    self.button_run_calib.setStyleSheet("background-color : yellow")
    self.pbar.setValue(0)
    try:
        img_acceptor = Image.open(path_acceptor_calib)
        print(f"File '{path_acceptor_calib}' loaded.")
        img_donor = Image.open(path_donor_calib)
        print(f"File '{path_donor_calib}' loaded.")
    except FileNotFoundError:
        print(f"File '{path_donor_calib}' not found.")
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
    
    self.img_acceptor_0 = np.array(list(img_stack_acceptor))
    self.img_donor_0 = np.array(list(img_stack_donor))
    
    if method == 'Optical Flow':
        
        print('Optical Flow registration started...')
        V_mean, U_mean = GUI_OpticalFlow_registration(self, img_stack_acceptor, img_stack_donor, radius_OF)
        print('Optical Flow registration finished!')
        print('Warping of donor movie...')
        img_warp_donor = Warp_OpticalFlow(img_stack_donor, V_mean, U_mean)
        
        self.pbar.setValue(100)
        
        print('Registration done!')
        self.button_run_calib.setStyleSheet("background-color : green")
        
        
        GUI_plot_OF_align(self, img_stack_acceptor, img_warp_donor, k = 0)
        
        self.matrix_align = (V_mean, U_mean)
        
    elif method == 'Affine':
    
        # we run the optimization algorithm with the selected stacks of frames from the donor and acceptor movies\
            # the optional parameters method_min and maxfev_value can be provided to finetune the optimization process, if needed
        print('Optimization of the transformation matrix ongoing...')
        res_align = minimize_logP(img_stack_donor, img_stack_acceptor, method_min=method_min_opt, maxfev_value = maxfev_value_opt)
        print('Optimization of the transformation matrix done!')
        
        # we extract the optimized transformation matrix coefficients
        self.matrix_align = res_align.x
        
        self.pbar.setValue(100)
        self.button_run_calib.setStyleSheet("background-color : green")
        
        print('Chromatic aberrations correction step completed!')
    
    self.setWindowModified(True)
    
def create_registration_tab(self):
    # First tab
    self.matrix_align = None
    
    tab1 = QWidget()
    layout1 = QVBoxLayout()
    
    layout_load_calib = QHBoxLayout()
    self.button_load_regist_matrix = QPushButton("load registration matrix")
    self.button_load_regist_matrix.clicked.connect(lambda: load_registration_matrix(self))
    layout_load_calib.addWidget(self.button_load_regist_matrix)
    button_show_regist_matrix = QPushButton("Show registration results")
    button_show_regist_matrix.clicked.connect(lambda: show_loaded_matrix(self))
    layout_load_calib.addWidget(button_show_regist_matrix)
    layout1.addLayout(layout_load_calib)
    
    layout_acceptor_calib = QHBoxLayout()
    label_acceptor_calib = QLabel("Acceptor calibration movie")
    button_acceptor_calib = QPushButton("Browse file")
    button_acceptor_calib.clicked.connect(lambda: open_file_dialog_acceptor_regis(self))
    self.widget_path_acceptor_calib = QLineEdit()
    layout_acceptor_calib.addWidget(label_acceptor_calib)
    layout_acceptor_calib.addWidget(button_acceptor_calib)
    layout_acceptor_calib.addWidget(self.widget_path_acceptor_calib)
    layout_donor_calib = QHBoxLayout()
    label_donor_calib = QLabel("Donor calibration movie")
    button_donor_calib = QPushButton("Browse file")
    button_donor_calib.clicked.connect(lambda: open_file_dialog_donor_regis(self))
    self.widget_path_donor_calib = QLineEdit()
    layout_donor_calib.addWidget(label_donor_calib)
    layout_donor_calib.addWidget(button_donor_calib)
    layout_donor_calib.addWidget(self.widget_path_donor_calib)
    layout1.addLayout(layout_acceptor_calib)
    layout1.addLayout(layout_donor_calib)
    
    layout_frames_calib = QHBoxLayout()
    label = QLabel("# of frames")
    layout_frames_calib.addWidget(label)

    # Create a QSpinBox widget
    spin_box_nb_frame = QSpinBox()
    spin_box_nb_frame.setMinimum(1)  # Set minimum value
    spin_box_nb_frame.setMaximum(1000000)   # Set maximum value
    spin_box_nb_frame.setValue(10)
    layout_frames_calib.addWidget(spin_box_nb_frame)
    layout1.addLayout(layout_frames_calib)
    
    self.method_calib = QComboBox()
    self.method_calib.addItems(["Optical Flow", "Affine"])
    layout1.addWidget(self.method_calib)
    
    layout_progress = QHBoxLayout()
    self.button_run_calib = QPushButton("Run registration")
    self.button_run_calib.clicked.connect(lambda: run_registration(self, spin_box_nb_frame.value()))
    self.button_run_calib.setStyleSheet("background-color : yellow")
    layout_progress.addWidget(self.button_run_calib)
    self.pbar = QProgressBar(self)
    layout_progress.addWidget(self.pbar)
    layout1.addLayout(layout_progress)
    
    layout_calib_vizu_slider = QVBoxLayout()
    self.figure_calib_res = plt.figure()
    self.canvas_calib_res = FigureCanvas(self.figure_calib_res)
    toolbar_calib_res = NavigationToolbar(self.canvas_calib_res, self)
    layout_calib_vizu_slider.addWidget(toolbar_calib_res)
    layout1.addLayout(layout_calib_vizu_slider)
    layout_calib_slider1 = QHBoxLayout()
    layout_calib_slider1.addWidget(self.canvas_calib_res)
    sliders_layout = QHBoxLayout()

    self.sliders_calib = {}
    self.checkboxs_calib = {}
    channels = ['Acceptor', 'Donor']
    for idx, channel in enumerate(channels):
        vbox = QVBoxLayout()
        checkbox = QCheckBox(f"{channel}", self)
        checkbox.setChecked(True)
        checkbox.stateChanged.connect(lambda: checkbox_state_changed(self))
        slider = QSlider(Qt.Vertical)
        slider.setMinimum(0)   # Corresponds to 0.5
        slider.setMaximum(150)  # Corresponds to 1.5
        slider.setValue(100)    # Default 1.0
        slider.setTickPosition(QSlider.TicksBelow)
        slider.setTickInterval(10)
        slider.valueChanged.connect(lambda: checkbox_state_changed(self))
        self.sliders_calib[channel] = slider
        self.checkboxs_calib[channel] = checkbox

        vbox.addWidget(checkbox)
        vbox.addWidget(slider)
        sliders_layout.addLayout(vbox)
    
    
    
    
    
    layout_calib_slider1.addLayout(sliders_layout)
    layout1.addLayout(layout_calib_slider1)
    
    
    tab1.setLayout(layout1)
    
    self.tabs.addTab(tab1, "Registration")
    
def load_registration_matrix(self):
    options = QFileDialog.Options()
    # You can set options like QFileDialog.DontUseNativeDialog if needed
    file_name, _ = QFileDialog.getOpenFileName(self, "Select the registration matrix file", "", "All Files (*);;Numpy Files (*.npy)", options=options)
    if file_name:
        self.matrix_align = np.load(file_name)
        print('Registration matrix loaded')
        self.button_load_regist_matrix.setStyleSheet("background-color : green")
        
def show_loaded_matrix(self):
    path_acceptor_calib = self.widget_path_acceptor_calib.text()
    path_donor_calib = self.widget_path_donor_calib.text()
    try:
        img_acceptor = Image.open(path_acceptor_calib)
        print(f"File '{path_acceptor_calib}' loaded.")
        img_donor = Image.open(path_donor_calib)
        print(f"File '{path_donor_calib}' loaded.")
    except FileNotFoundError:
        print(f"File '{path_donor_calib}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
        
    self.img_donor_0 = np.zeros((img_donor.height, img_donor.width, 10))
    self.img_acceptor_0 = np.zeros((img_acceptor.height, img_acceptor.width, 10))
        
    img_donor.seek(0)
    img_acceptor.seek(0)
    self.img_donor_0[:,:,0] = np.array(img_donor)
    self.img_acceptor_0[:,:,0] = np.array(img_acceptor)
    
    method = self.method_calib.currentText()

    if method == 'Optical Flow':
        
        img_warp_donor = Warp_OpticalFlow(self.img_donor_0, self.matrix_align[0], self.matrix_align[1])
        
        GUI_plot_OF_align(self, self.img_acceptor_0, img_warp_donor, k = 0)
        
    elif method == 'Affine':
    
        print('Development ongoing...')
        
def checkbox_state_changed(self):
    
    if self.checkboxs_calib['Acceptor'].isChecked():
        R_corr = np.array(list(self.R_channel))
        R_corr = R_corr / (self.sliders_calib['Acceptor'].value()/100)
        R_corr[R_corr > 1] = 1
    else:
        R_corr = np.array(list(self.B_channel))  # copy the zero intensity Blue channel (independent copy)
        
    
    if self.checkboxs_calib['Donor'].isChecked():
        G_corr = np.array(list(self.G_channel))
        G_corr = G_corr / (self.sliders_calib['Donor'].value()/100)
        G_corr[G_corr > 1] = 1
    else:
        G_corr = np.array(list(self.B_channel)) # copy the zero intensity Blue channel (independent copy)
        
    self.img_calib.set_data(np.stack((R_corr, G_corr, self.B_channel), axis=-1))
    
    self.canvas_calib_res.draw_idle()
    
def update_contrast(self):
    
    R_corr = np.array(list(self.R_channel))
    R_corr = R_corr / (self.sliders_calib['Acceptor'].value()/100)
    R_corr[R_corr > 1] = 1
    
    G_corr = np.array(list(self.G_channel))
    G_corr = G_corr / (self.sliders_calib['Donor'].value()/100)
    G_corr[G_corr > 1] = 1
    
    self.img_calib.set_data(np.stack((R_corr, G_corr, self.B_channel), axis=-1))
    
    self.canvas_calib_res.draw_idle()
    
def GUI_plot_OF_align(self, img_stack_acceptor, img_stack_donor, k = 0):
    
    self.R_channel = (img_stack_acceptor[:,:,k]-np.min(img_stack_acceptor[:,:,k]))/(np.max(img_stack_acceptor[:,:,k])-np.min(img_stack_acceptor[:,:,k]))
    self.G_channel = (img_stack_donor[:,:,k]-np.min(img_stack_donor[:,:,k]))/(np.max(img_stack_donor[:,:,k])-np.min(img_stack_donor[:,:,k]))
    self.B_channel = np.zeros(img_stack_donor[:,:,0].shape)
    
    self.figure_calib_res.clear()
    ax = self.figure_calib_res.add_subplot(111)
    
    # Show the merge channels
    self.img_calib = ax.imshow(np.stack((self.R_channel, self.G_channel, self.B_channel), axis=-1))
    
    self.canvas_calib_res.draw_idle()
    
def regist_QC_plot(self):
    self.count_windows = self.count_windows + 1
    second_window = QCRegistrationWindow(self.img_acceptor_0, self.img_donor_0, self.matrix_align, self.method_calib.currentText(), self.count_windows)
    second_window.show()
    self.plot_windows.append(second_window)
        
def GUI_OpticalFlow_registration(self, img_stack_acceptor, img_stack_donor, radius_OF):
    
    nb_frame = img_stack_acceptor.shape[2]
    
    v_stack = np.zeros(img_stack_acceptor.shape)
    u_stack = np.zeros(img_stack_acceptor.shape)
    
    for k in range(nb_frame):
        self.pbar.setValue(int((k/nb_frame*100)))
        v_stack[:,:,k], u_stack[:,:,k] = optical_flow_ilk(img_stack_acceptor[:,:,k], img_stack_donor[:,:,k], 
                                                          radius=radius_OF, 
                                                          prefilter = False,
                                                          gaussian = False)

    
    V_mean = np.mean(v_stack, axis = 2)
    U_mean = np.mean(u_stack, axis = 2)
    
    return V_mean, U_mean


def open_file_dialog_acceptor_regis(self):
    # Open file dialog and get selected file path
    options = QFileDialog.Options()
    # You can set options like QFileDialog.DontUseNativeDialog if needed
    file_name, _ = QFileDialog.getOpenFileName(self, "Select the Acceptor calibration movie", "", "All Files (*);;Text Files (*.txt)", options=options)
    if file_name:
        self.widget_path_acceptor_calib.setText(file_name)
        
def open_file_dialog_donor_regis(self):
    # Open file dialog and get selected file path
    options = QFileDialog.Options()
    # You can set options like QFileDialog.DontUseNativeDialog if needed
    file_name, _ = QFileDialog.getOpenFileName(self, "Select the Donor calibration movie", "", "All Files (*);;Text Files (*.txt)", options=options)
    if file_name:
        self.widget_path_donor_calib.setText(file_name)