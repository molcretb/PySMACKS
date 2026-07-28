# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 09:26:31 2026

@author: molcre0000
"""

import sys
from PyQt5.QtWidgets import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.widgets import SpanSelector
from PyQt5.QtGui import QFont
from PyQt5.QtCore import Qt, pyqtSignal, QSignalBlocker
from openfret import read_data
from chromatic_aberrations_correction import *
from utils import *
from traces_extractor import *
from drift_correction import *
import copy


class PlotFRETWindow(QMainWindow):
    def __init__(self, SE_matrix, bin_size, ax_lim_eps, nb_gauss):
        super().__init__()
        self.setWindowTitle("Corrected FRET histogram")
        self.resize(1400, 700)
        
        self.x = np.linspace(-ax_lim_eps, 1+ax_lim_eps, len(SE_matrix)) # -1
        
        self.lines_fit = [[] for i in range(nb_gauss)]
        
        self.y = np.sum(SE_matrix, 0)
        self.nb_gauss = nb_gauss
        self.bin_size = bin_size
        self.ax_lim_eps = ax_lim_eps
        
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout1 = QHBoxLayout(central_widget)
        
        
        left = QVBoxLayout()
        self.figure_FRET_hist = plt.figure()
        self.canvas_FRET_hist = FigureCanvas(self.figure_FRET_hist)
        toolbar = NavigationToolbar(self.canvas_FRET_hist, self)
        left.addWidget(toolbar)
        left.addWidget(self.canvas_FRET_hist)
        
        self.ax_FRET_hist = self.figure_FRET_hist.add_subplot(111)
        self.ax_FRET_hist.plot(self.x, self.y, 'kx')
        plt.title('FRET corrected histogram')
        plt.xlabel('FRET efficiency')
        plt.ylabel('Counts')
        
        for i in range(nb_gauss):
            self.lines_fit[i], = self.ax_FRET_hist.plot([],[])
        
        self.canvas_FRET_hist.draw_idle()
        
        
        right = QVBoxLayout()

        right.addWidget(QLabel("Fit parameters"))

        self.table = QTableWidget(3, nb_gauss + 1)
        
        list_fit_param_names = ["Parameter"]
        
        for j in range(nb_gauss):
            list_fit_param_names.append("Peak " + str(j))

        self.table.setHorizontalHeaderLabels(
            list_fit_param_names
        )

        #names = ["Amplitude", "Mean", "Sigma"]
        
        self.table.setItem(0, 0, QTableWidgetItem("Amplitude"))
        self.table.setItem(1, 0, QTableWidgetItem("Mean"))
        self.table.setItem(2, 0, QTableWidgetItem("Sigma"))

        for i in range(nb_gauss):
            self.table.setItem(0, i+1, QTableWidgetItem(str(np.max(self.y))))
            self.table.setItem(1, i+1, QTableWidgetItem(str(0.5)))
            self.table.setItem(2, i+1, QTableWidgetItem(str(0.5)))

        right.addWidget(self.table)

        self.fit_button = QPushButton("Fit")

        self.fit_button.clicked.connect(self.fit_curve)

        right.addWidget(self.fit_button)
        
        
        self.table_area = QTableWidget(nb_gauss + 1, 2)
        self.table_area.setHorizontalHeaderLabels(['Peak', 'Fraction (%)'])
        
        for j in range(nb_gauss):
            self.table_area.setItem(j, 0, QTableWidgetItem(str(j)))
            self.table_area.setItem(j, 1, QTableWidgetItem("NA"))
        
        
        right.addWidget(self.table_area)
        
        
        
        right.addStretch()

        layout1.addLayout(left, 3)
        layout1.addLayout(right, 1)
        

        #self.plot()
        
        
        
        
        # ax_FRET_hist = figure_FRET_hist.add_subplot(111)
        # ax_FRET_hist.plot(x, y, 'kx')
        # plt.title('FRET corrected histogram')
        # plt.xlabel('FRET efficiency')
        # plt.ylabel('Counts')
        # canvas_FRET_hist.draw_idle()
        
        # fit_res = self.multi_fit_Gauss_spot(x, y, nb_gauss = nb_gauss, bin_size = bin_size, ax_lim_eps = ax_lim_eps)
        
        # y_sum = np.zeros(len(y))
        
        # for j in range(nb_gauss):
        #     y_j = fit_res[0] + self.generate_Gauss_data(x, fit_res[3*j+1], fit_res[3*j+2], fit_res[3*j+3])
        #     ax_FRET_hist.plot(x, y_j, '-')
        #     y_sum = y_sum + y_j
        # ax_FRET_hist.plot(x, y_sum - (nb_gauss - 1) * fit_res[0], '--')
        
    def plot(self):

        #self.figure.clear()
        
        sum_area = np.zeros(self.nb_gauss)
        
        for i in range(self.nb_gauss):
            self.lines_fit[i].set_data(self.x, self.generate_Gauss_data(self.x, self.fit_res[3*i], self.fit_res[3*i+1], self.fit_res[3*i+2]))
            sum_area[i] = self.fit_res[3*i]*self.fit_res[3*i+2]
        
        for i in range(self.nb_gauss):
            self.table_area.item(i, 1).setText(str(np.round(sum_area[i]/np.sum(sum_area), 3)))

        self.canvas_FRET_hist.draw_idle()
        
        
    def get_parameters(self):

        pars = []

        for j in range(1,self.nb_gauss + 1):
            pars.append(float(self.table.item(0, j).text()))
            pars.append(float(self.table.item(1, j).text()))
            pars.append(float(self.table.item(2, j).text()))

        return pars

    # -------------------------------------------------------------
    def fit_curve(self):

        p0 = self.get_parameters()

        self.fit_res = self.multi_fit_Gauss_spot(self.x, self.y, p0, nb_gauss = self.nb_gauss, bin_size = self.bin_size, ax_lim_eps = self.ax_lim_eps)

        count = 0
        for j in range(1,self.nb_gauss + 1):
            self.table.item(0, j).setText(str(np.round(self.fit_res[count], 5)))
            self.table.item(1, j).setText(str(np.round(self.fit_res[count+1], 5)))
            self.table.item(2, j).setText(str(np.round(self.fit_res[count+2], 5)))
            count += 3

        self.plot()
        
    def FRET_hist_multi_Gauss(self, x, gauss_param, nb_gauss):
        g = np.zeros(len(x))
        for i in range(nb_gauss):
            g = g + gauss_param[3*i]*np.exp( - (x-float(gauss_param[3*i+1]))**2/(2*gauss_param[3*i+2]**2))
        return g.ravel()
    
    def multi_Gauss_error_func(self, gauss_param, nb_gauss, x, y):
        return self.FRET_hist_multi_Gauss(x, gauss_param, nb_gauss) - y
    
    def multi_fit_Gauss_spot(self, x, y, gauss_param, nb_gauss = 2, bin_size = 0.01, ax_lim_eps = 0.2):
        initial_guess = tuple(gauss_param)
        # custom_gaussian = lambda x, offset, mu: twoD_Gaussian(x, offset, mu, nb_gauss)
        try:
            # popt, pcov = opt.curve_fit(custom_gaussian, (x, y), img_spot.ravel(), p0=initial_guess)
            output = least_squares(self.multi_Gauss_error_func, x0=initial_guess, jac='2-point', bounds=(0, np.inf), method='trf', ftol=1e-08,
                           xtol=1e-08, gtol=1e-08, x_scale=1.0, loss='linear', f_scale=1.0, diff_step=None,
                           tr_solver=None,
                           tr_options={}, jac_sparsity=None, max_nfev=None, verbose=0, args=(nb_gauss, x, y.ravel()))
            
        except RuntimeError:
            print('No peak found')
            return
        
        
        for j in range(nb_gauss):
            print('Peak ' + str(j) + ': ' + str(output.x[3*j+1]))
        
        return output.x
        
    def generate_Gauss_data(self, x, A, B, C):
        return A * np.exp(-(x-B)**2/C**2)
        

class AlphaCorrWindow(QMainWindow):
    def __init__(self, popt, bins_center, counts, alpha, FRET_conc, FRET_S_conc, bin_size, ax_lim_eps):
        super().__init__()
        self.setWindowTitle("Alpha factor correction")
        self.resize(1400, 700)
        
        
        toolbar_container = QWidget()
        toolbar_stack = QStackedLayout(toolbar_container)
        
        central_widget = QWidget()
        layout1 = QVBoxLayout(central_widget)
        figure_QC_alpha = plt.figure()
        self.canvas_QC_alpha = FigureCanvas(figure_QC_alpha)
        self.toolbar_QC = NavigationToolbar(self.canvas_QC_alpha, self)
        toolbar_stack.addWidget(self.toolbar_QC)
        
        layout1.addWidget(toolbar_container)
        layout1.addWidget(self.canvas_QC_alpha, 2)
        #self.central_widget.addLayout(layout1)
        self.setCentralWidget(central_widget)
        
        x = np.linspace(0, 1, 1000) # -1
        
        y_gauss = GUI_Gauss_data(x, popt[0], popt[1], popt[2])
        
        
        ax_QC_alpha = figure_QC_alpha.add_subplot(111)
        
        ax_QC_alpha.plot(bins_center[1:-1], counts[1:-1],'kx', label='data')
        label_fit = r'$\alpha$ = '+str(np.round(alpha, 4))
        ax_QC_alpha.plot(x, y_gauss, label=label_fit)
        plt.title(r'$\alpha$ correction from Donor-Only traces')
        plt.xlabel('FRET efficiency')
        plt.ylabel('Counts')
        ax_QC_alpha.legend()
        self.canvas_QC_alpha.draw_idle()
        
        
        self.figure_alpha_map = plt.figure()
        self.canvas_alpha_map = FigureCanvas(self.figure_alpha_map)
        self.toolbar_alpha_map = NavigationToolbar(self.canvas_alpha_map, self)
        toolbar_stack.addWidget(self.toolbar_alpha_map)
        
        self.toolbar_stack = toolbar_stack
        
        self.ax_alpha_map = self.figure_alpha_map.add_subplot(111)
        
        self.img_alpha_map = self.ax_alpha_map.imshow(np.zeros((1, 1)), origin='lower')
        plt.xlabel('FRET efficiency')
        plt.ylabel('Stoichiometry')
        self.canvas_alpha_map.draw_idle()
        
        layout1.addWidget(self.canvas_alpha_map, 4)
        
        self.GUI_generate_2D_histogram(FRET_conc, FRET_S_conc, bin_size, ax_lim_eps)
        
        self.canvas_QC_alpha.mpl_connect(
            "button_press_event",
            lambda event: self.activate_canvas(self.canvas_QC_alpha))
            
        self.canvas_alpha_map.mpl_connect(
            "button_press_event",
            lambda event: self.activate_canvas(self.canvas_alpha_map)
            )
        
    def activate_canvas(self, canvas):
        if canvas is self.canvas_QC_alpha:
            self.toolbar_stack.setCurrentWidget(self.toolbar_QC)
        else:
            self.toolbar_stack.setCurrentWidget(self.toolbar_alpha_map)
        
    def GUI_generate_2D_histogram(self, FRET_conc, FRET_S_conc, bin_size, ax_lim_eps):
        
        nb_bins = int((1+2*ax_lim_eps)/bin_size)
        
        self.hist2D_map = np.zeros((nb_bins, nb_bins))
        
        compt_nan = 0
        
        extra_bins = int(ax_lim_eps // bin_size)
        
        FRET_hist_bin = [FRET_S_conc // bin_size + extra_bins, FRET_conc // bin_size + extra_bins]
        
        for i in range(len(FRET_conc)):
            try:
                x_i = FRET_hist_bin[0][i].astype(int)
                y_i = FRET_hist_bin[1][i].astype(int)
                if (x_i >= 0) and (y_i >=0):
                    self.hist2D_map[x_i, y_i] += 1
            except IndexError:
                compt_nan += 1
        
        self.img_alpha_map.set_data(self.hist2D_map)
        self.img_alpha_map.set_clim(vmin=0, vmax=self.hist2D_map.max())
        self.img_alpha_map.set_extent([-ax_lim_eps,1+ax_lim_eps,-ax_lim_eps,1+ax_lim_eps])
         #, interpolation='none', origin='lower', extent=[-ax_lim_eps,1+ax_lim_eps,-ax_lim_eps,1+ax_lim_eps])
        self.ax_alpha_map.set_ylim([-ax_lim_eps,1+ax_lim_eps])
        self.ax_alpha_map.set_xlim([-ax_lim_eps,1+ax_lim_eps])
        
        self.ax_alpha_map.autoscale_view()
        self.canvas_alpha_map.draw_idle()
        
        
class DeltaCorrWindow(QMainWindow):
    def __init__(self, popt, bins_center, counts, delta, FRET_conc, FRET_S_conc, bin_size, ax_lim_eps):
        super().__init__()
        self.setWindowTitle("Delta factor correction")
        self.resize(1400, 700)
        
        
        toolbar_container = QWidget()
        toolbar_stack = QStackedLayout(toolbar_container)
        
        central_widget = QWidget()
        layout1 = QVBoxLayout(central_widget)
        figure_QC_delta = plt.figure()
        self.canvas_QC_delta = FigureCanvas(figure_QC_delta)
        self.toolbar_QC = NavigationToolbar(self.canvas_QC_delta, self)
        toolbar_stack.addWidget(self.toolbar_QC)
        
        layout1.addWidget(toolbar_container)
        layout1.addWidget(self.canvas_QC_delta, 2)
        #self.central_widget.addLayout(layout1)
        self.setCentralWidget(central_widget)
        
        x = np.linspace(0, 1, 1000) # -1
        
        y_gauss = GUI_Gauss_data(x, popt[0], popt[1], popt[2])
        
        
        ax_QC_delta = figure_QC_delta.add_subplot(111)
        
        ax_QC_delta.plot(bins_center[1:-1], counts[1:-1],'kx', label='data')
        label_fit = r'$\delta$ = '+str(np.round(delta, 4))
        ax_QC_delta.plot(x, y_gauss, label=label_fit)
        plt.title(r'$\delta$ correction from Acceptor-Only traces')
        plt.xlabel('FRET efficiency')
        plt.ylabel('Counts')
        ax_QC_delta.legend()
        self.canvas_QC_delta.draw_idle()
        
        
        self.figure_delta_map = plt.figure()
        self.canvas_delta_map = FigureCanvas(self.figure_delta_map)
        self.toolbar_delta_map = NavigationToolbar(self.canvas_delta_map, self)
        toolbar_stack.addWidget(self.toolbar_delta_map)
        
        self.toolbar_stack = toolbar_stack
        
        self.ax_delta_map = self.figure_delta_map.add_subplot(111)
        
        self.img_delta_map = self.ax_delta_map.imshow(np.zeros((1, 1)), origin='lower')
        plt.xlabel('FRET efficiency')
        plt.ylabel('Stoichiometry')
        self.canvas_delta_map.draw_idle()
        
        layout1.addWidget(self.canvas_delta_map, 4)
        
        self.GUI_generate_2D_histogram(FRET_conc, FRET_S_conc, bin_size, ax_lim_eps)
        
        self.canvas_QC_delta.mpl_connect(
            "button_press_event",
            lambda event: self.activate_canvas(self.canvas_QC_delta))
            
        self.canvas_delta_map.mpl_connect(
            "button_press_event",
            lambda event: self.activate_canvas(self.canvas_delta_map)
            )
        
    def activate_canvas(self, canvas):
        if canvas is self.canvas_QC_delta:
            self.toolbar_stack.setCurrentWidget(self.toolbar_QC)
        else:
            self.toolbar_stack.setCurrentWidget(self.toolbar_delta_map)
        
    def GUI_generate_2D_histogram(self, FRET_conc, FRET_S_conc, bin_size, ax_lim_eps):
        
        nb_bins = int((1+2*ax_lim_eps)/bin_size)
        
        self.hist2D_map = np.zeros((nb_bins, nb_bins))
        
        compt_nan = 0
        
        extra_bins = int(ax_lim_eps // bin_size)
        
        FRET_hist_bin = [FRET_S_conc // bin_size + extra_bins, FRET_conc // bin_size + extra_bins]
        
        for i in range(len(FRET_conc)):
            try:
                x_i = FRET_hist_bin[0][i].astype(int)
                y_i = FRET_hist_bin[1][i].astype(int)
                if (x_i >= 0) and (y_i >=0):
                    self.hist2D_map[x_i, y_i] += 1
            except IndexError:
                compt_nan += 1
        
        self.img_delta_map.set_data(self.hist2D_map)
        self.img_delta_map.set_clim(vmin=0, vmax=self.hist2D_map.max())
        self.img_delta_map.set_extent([-ax_lim_eps,1+ax_lim_eps,-ax_lim_eps,1+ax_lim_eps])
         #, interpolation='none', origin='lower', extent=[-ax_lim_eps,1+ax_lim_eps,-ax_lim_eps,1+ax_lim_eps])
        self.ax_delta_map.set_ylim([-ax_lim_eps,1+ax_lim_eps])
        self.ax_delta_map.set_xlim([-ax_lim_eps,1+ax_lim_eps])
        
        self.ax_delta_map.autoscale_view()
        self.canvas_delta_map.draw_idle()
        
# class DeltaCorrWindow(QMainWindow):
#     def __init__(self, popt, bins_center, counts, delta, FRET_conc, FRET_S_conc, bin_size, ax_lim_eps):
#         super().__init__()
#         self.setWindowTitle("Delta factor correction")
#         self.resize(1400, 700)
        
        
#         central_widget = QWidget()
#         layout1 = QVBoxLayout(central_widget)
#         figure_QC_delta = plt.figure()
#         canvas_QC_delta = FigureCanvas(figure_QC_delta)
#         toolbar_QC_delta = NavigationToolbar(canvas_QC_delta, self)
#         layout1.addWidget(toolbar_QC_delta)
#         layout1.addWidget(canvas_QC_delta)
#         #self.central_widget.addLayout(layout1)
#         self.setCentralWidget(central_widget)
        
#         x = np.linspace(0, 1, 1000)
        
#         y_gauss = GUI_Gauss_data(x, popt[0], popt[1], popt[2])
        
        
#         ax_QC_delta = figure_QC_delta.add_subplot(111)
        
#         ax_QC_delta.plot(bins_center[1:-1], counts[1:-1],'kx', label='data')
#         label_fit = r'$\delta$ = '+str(np.round(delta, 4))
#         ax_QC_delta.plot(x, y_gauss, label=label_fit)
#         plt.title(r'$\delta$ correction from Acceptor-Only traces')
#         plt.xlabel('FRET stoichiometry')
#         plt.ylabel('Counts')
#         ax_QC_delta.legend()
#         canvas_QC_delta.draw_idle()
        
#         self.figure_delta_map = plt.figure()
#         self.canvas_delta_map = FigureCanvas(self.figure_delta_map)
#         toolbar = NavigationToolbar(self.canvas_delta_map, self)
#         layout1.addWidget(toolbar)
        
#         self.ax_delta_map = self.figure_delta_map.add_subplot(111)
        
#         self.img_delta_map = self.ax_delta_map.imshow(np.zeros((1, 1)), origin='lower')
#         plt.xlabel('FRET efficiency')
#         plt.ylabel('Stoichiometry')
#         self.canvas_delta_map.draw_idle()
        
#         layout1.addWidget(self.canvas_delta_map)
        
#         self.GUI_generate_2D_histogram(FRET_conc, FRET_S_conc, bin_size, ax_lim_eps)
        
#     def GUI_generate_2D_histogram(self, FRET_conc, FRET_S_conc, bin_size, ax_lim_eps):
        
#         nb_bins = int((1+2*ax_lim_eps)/bin_size)
        
#         self.hist2D_map = np.zeros((nb_bins, nb_bins))
        
#         compt_nan = 0
        
#         # ID_traces = list(traces_dict.keys())
        
#         # FRET_conc = np.array([])
        
#         # FRET_S_conc = np.array([])
        
#         # for i in ID_traces:
#         #     FRET_conc = np.concat((FRET_conc, traces_dict[i]['FRET_eff']))
#         #     FRET_S_conc = np.concat((FRET_S_conc, traces_dict[i]['FRET_stoi']))
#         # x_ax = np.linspace(-ax_lim_eps,1+ax_lim_eps, nb_bins)
        
#         extra_bins = int(ax_lim_eps // bin_size)
        
#         #FRET_conc = FRET_conc * (FRET_conc >= 0)
        
#         #FRET_S_conc = FRET_S_conc * (FRET_S_conc >= 0)
        
#         FRET_hist_bin = [FRET_S_conc // bin_size + extra_bins, FRET_conc // bin_size + extra_bins]
        
#         for i in range(len(FRET_conc)):
#             try:
#                 x_i = FRET_hist_bin[0][i].astype(int)
#                 y_i = FRET_hist_bin[1][i].astype(int)
#                 if (x_i >= 0) and (y_i >=0):
#                     self.hist2D_map[x_i, y_i] += 1
#             except IndexError:
#                 compt_nan += 1
        
#         # hist2D_map[0,:] = 0
#         # hist2D_map[-1,:] = 0
#         # hist2D_map[:,0] = 0
#         # hist2D_map[:,-1] = 0
        
#         self.img_delta_map.set_data(self.hist2D_map)
#         self.img_delta_map.set_clim(vmin=0, vmax=self.hist2D_map.max())
#         self.img_delta_map.set_extent([-ax_lim_eps,1+ax_lim_eps,-ax_lim_eps,1+ax_lim_eps])
#          #, interpolation='none', origin='lower', extent=[-ax_lim_eps,1+ax_lim_eps,-ax_lim_eps,1+ax_lim_eps])
#         self.ax_delta_map.set_ylim([-ax_lim_eps,1+ax_lim_eps])
#         self.ax_delta_map.set_xlim([-ax_lim_eps,1+ax_lim_eps])
        
#         self.ax_delta_map.autoscale_view()
#         self.canvas_delta_map.draw_idle()
        
        

def create_SEhist_tab(self):
    
    self.filt_DA_data = None
    self.filt_DO_data = None
    self.filt_AO_data = None
    self.FRET_corr_factors = None
    
    tab1 = QWidget()
    layout1 = QVBoxLayout()
    
    self.multiclass_corr_checkbox = QCheckBox("Multiclass labels")
    self.multiclass_corr_checkbox.setToolTip('Tick if your dataset contains both Donor-Acceptor, Donor-Only and Acceptor-Only labels; otherwise you need to load the datasets separately')
    self.multiclass_corr_checkbox.setChecked(True)
    self.multiclass_corr_checkbox.stateChanged.connect(lambda: toggle_multiclass_corr_buttons(self))
    layout1.addWidget(self.multiclass_corr_checkbox)
    
    self.frame_multi_datasets = QFrame()
    
    layout_loading_filter_multiclassdata = QHBoxLayout()
    self.frame_multi_datasets.setLayout(layout_loading_filter_multiclassdata)
    self.button_load_multi_data = QPushButton("Load filtered dataset")
    self.button_load_multi_data.setToolTip('Select the filtered multilabelled traces dataset')
    self.button_load_multi_data.clicked.connect(lambda: load_filt_multi_data(self))
    layout_loading_filter_multiclassdata.addWidget(self.button_load_multi_data)
    
    layout1.addWidget(self.frame_multi_datasets)
    
    ##################### HLayout for displaying individual channels
    
    self.frame_DA_DO_AO_datasets = QFrame()
    
    layout_loading_filter_data = QHBoxLayout()
    
    self.frame_DA_DO_AO_datasets.setLayout(layout_loading_filter_data)
    
    self.button_load_DA_data = QPushButton("Load DA dataset")
    self.button_load_DA_data.setToolTip('Select the filtered DA traces dataset')
    self.button_load_DA_data.clicked.connect(lambda: load_filt_DA_data(self))
    layout_loading_filter_data.addWidget(self.button_load_DA_data)
    
    self.button_load_DO_data = QPushButton("Load DO dataset")
    self.button_load_DO_data.setToolTip('Select the filtered Donor-Only traces dataset')
    self.button_load_DO_data.clicked.connect(lambda: load_filt_DO_data(self))
    layout_loading_filter_data.addWidget(self.button_load_DO_data)
    
    self.button_load_AO_data = QPushButton("Load AO dataset")
    self.button_load_AO_data.setToolTip('Select the filtered Acceptor-Only traces dataset')
    self.button_load_AO_data.clicked.connect(lambda: load_filt_AO_data(self))
    layout_loading_filter_data.addWidget(self.button_load_AO_data)
    
    layout1.addWidget(self.frame_DA_DO_AO_datasets)
    
    self.frame_DA_DO_AO_datasets.hide()
    
    #layout1.addLayout(layout_loading_filter_data)
    
    layout_map_param = QHBoxLayout()
    label = QLabel("Bin size")
    layout_map_param.addWidget(label)
    self.spin_box_bin_map = QDoubleSpinBox()
    self.spin_box_bin_map.setDecimals(3)
    self.spin_box_bin_map.setToolTip('Bin size of the 2D SE histogram')
    self.spin_box_bin_map.setRange(0.001, 0.1)  # Set minimum value
    self.spin_box_bin_map.setSingleStep(0.001)
    self.spin_box_bin_map.setValue(0.01)
    layout_map_param.addWidget(self.spin_box_bin_map)
    
    label = QLabel("Axis extra range")
    layout_map_param.addWidget(label)
    self.spin_box_extra_range = QDoubleSpinBox()
    self.spin_box_extra_range.setDecimals(2)
    self.spin_box_extra_range.setToolTip('Extra range for axis')
    self.spin_box_extra_range.setRange(0, 1)  # Set minimum value
    self.spin_box_extra_range.setSingleStep(0.05)
    self.spin_box_extra_range.setValue(0.2)
    layout_map_param.addWidget(self.spin_box_extra_range)
    
    label = QLabel("Nb. states")
    layout_map_param.addWidget(label)
    self.spin_box_nb_states = QSpinBox()
    self.spin_box_nb_states.setToolTip('Select the expected number of states (for beta-gamma correction)')
    self.spin_box_nb_states.setMinimum(1)  # Set minimum value
    self.spin_box_nb_states.setMaximum(100)   # Set maximum value
    self.spin_box_nb_states.setValue(2)
    layout_map_param.addWidget(self.spin_box_nb_states)
    
    layout1.addLayout(layout_map_param)
    
    
    layout_corr_factor = QHBoxLayout()
    label = QLabel("Correction factors:")
    layout_corr_factor.addWidget(label)
    
    self.alpha_corr_checkbox = QCheckBox('\u03B1')
    self.alpha_corr_checkbox.setToolTip('\u03B1 correction factor')
    self.alpha_corr_checkbox.setChecked(False)
    self.alpha_corr_checkbox.stateChanged.connect(lambda: GUI_calculate_alpha(self, self.filt_DO_data))
    layout_corr_factor.addWidget(self.alpha_corr_checkbox)
    self.widget_alpha_value = QLineEdit()
    self.widget_alpha_value.setToolTip('\u03B1 value')
    self.widget_alpha_value.setText('0')
    layout_corr_factor.addWidget(self.widget_alpha_value)
    
    self.delta_corr_checkbox = QCheckBox('\u03B4')
    self.delta_corr_checkbox.setToolTip('\u03B4 correction factor')
    self.delta_corr_checkbox.setChecked(False)
    self.delta_corr_checkbox.stateChanged.connect(lambda: GUI_calculate_delta(self, self.filt_AO_data))
    layout_corr_factor.addWidget(self.delta_corr_checkbox)
    self.widget_delta_value = QLineEdit()
    self.widget_delta_value.setToolTip('\u03B4 value')
    self.widget_delta_value.setText('0')
    layout_corr_factor.addWidget(self.widget_delta_value)
    
    self.beta_corr_checkbox = QCheckBox('\u03B2')
    self.beta_corr_checkbox.setToolTip('\u03B2 correction factor')
    self.beta_corr_checkbox.setChecked(False)
    self.beta_corr_checkbox.stateChanged.connect(lambda: GUI_beta_corr(self))
    layout_corr_factor.addWidget(self.beta_corr_checkbox)
    self.widget_beta_value = QLineEdit()
    self.widget_beta_value.setToolTip('\u03B2 value')
    self.widget_beta_value.setText('1')
    layout_corr_factor.addWidget(self.widget_beta_value)
    
    self.gamma_corr_checkbox = QCheckBox('\u03B3')
    self.gamma_corr_checkbox.setToolTip('\u03B3 correction factor')
    self.gamma_corr_checkbox.setChecked(False)
    self.gamma_corr_checkbox.stateChanged.connect(lambda: GUI_gamma_corr(self))
    layout_corr_factor.addWidget(self.gamma_corr_checkbox)
    self.widget_gamma_value = QLineEdit()
    self.widget_gamma_value.setToolTip('\u03B3 value')
    self.widget_gamma_value.setText('1')
    layout_corr_factor.addWidget(self.widget_gamma_value)
    
    layout1.addLayout(layout_corr_factor)
    
    
    layout_run_corr = QHBoxLayout()
    
    button_calc_FRET_eff = QPushButton("Calculate FRET eff.")
    button_calc_FRET_eff.setToolTip('Calculate the FRET efficiency')
    button_calc_FRET_eff.clicked.connect(lambda: GUI_calculate_FRET_eff(self, self.filt_DA_data))
    layout_run_corr.addWidget(button_calc_FRET_eff)
    
    layout1.addLayout(layout_run_corr)
    
    
    layout_project_FRET_eff = QHBoxLayout()
    button_project_FRET_eff = QPushButton("Plot FRET histogram")
    button_project_FRET_eff.setToolTip('Plot the corrected FRET Efficiency histogram')
    button_project_FRET_eff.clicked.connect(lambda: GUI_plot_FRET_eff(self))
    layout_project_FRET_eff.addWidget(button_project_FRET_eff)
    layout1.addLayout(layout_project_FRET_eff)
    
    self.figure_SEhist = plt.figure()
    self.canvas_SEhist = FigureCanvas(self.figure_SEhist)
    toolbar = NavigationToolbar(self.canvas_SEhist, self)
    layout1.addWidget(toolbar)
    
    self.ax_SEhist = self.figure_SEhist.add_subplot(111)
    
    self.img_SEhist = self.ax_SEhist.imshow(np.zeros((1, 1)), origin='lower')
    plt.xlabel('FRET efficiency')
    plt.ylabel('Stoichiometry')
    self.canvas_SEhist.draw_idle()
    
    #layout1.addWidget(self.canvas_SEhist)
    
    
    layout_calib_slider1 = QHBoxLayout()
    layout_calib_slider1.addWidget(self.canvas_SEhist)
    self.sliders_SEhist = {}
    vbox = QVBoxLayout()
    label = QLabel("Contrast")
    slider = QSlider(Qt.Vertical)
    slider.setMinimum(0)   # Corresponds to 0.5
    slider.setMaximum(100)  # Corresponds to 1.5
    slider.setValue(100)    # Default 1.0
    slider.setTickPosition(QSlider.TicksBelow)
    slider.setTickInterval(10)
    #slider.valueChanged.connect(self.update_contrast)
    slider.valueChanged.connect(lambda: slider_SEhist_contrast(self))
    self.sliders_SEhist['0'] = slider

        #vbox.addWidget(label)
    vbox.addWidget(label)
    vbox.addWidget(slider)
    layout_calib_slider1.addLayout(vbox)
    layout1.addLayout(layout_calib_slider1, 4)
    
    tab1.setLayout(layout1)
    self.tabs.addTab(tab1, "FRET hist.")
    

def GUI_beta_corr(self):
    print('Beta running')
    if self.beta_corr_checkbox.isChecked():
        with QSignalBlocker(self.gamma_corr_checkbox):
            self.gamma_corr_checkbox.setChecked(True)
        GUI_calculate_beta_gamma(self, self.filt_DA_data)
    else:
        with QSignalBlocker(self.gamma_corr_checkbox):
            self.gamma_corr_checkbox.setChecked(False)
        self.widget_beta_value.setText('1')
        self.widget_gamma_value.setText('1')
        
def toggle_multiclass_corr_buttons(self):
    if self.multiclass_corr_checkbox.isChecked():  # Checked
        self.frame_multi_datasets.show()
        self.frame_DA_DO_AO_datasets.hide()
    else:  # Unchecked
        self.frame_multi_datasets.hide()
        self.frame_DA_DO_AO_datasets.show()
        

def GUI_gamma_corr(self):
    print('Gamma running')
    if self.gamma_corr_checkbox.isChecked():
        with QSignalBlocker(self.beta_corr_checkbox):
            self.beta_corr_checkbox.setChecked(True)
        GUI_calculate_beta_gamma(self, self.filt_DA_data)
    else:
        with QSignalBlocker(self.beta_corr_checkbox):
            self.beta_corr_checkbox.setChecked(False)
        self.widget_beta_value.setText('1')
        self.widget_gamma_value.setText('1')

def slider_SEhist_contrast(self):
    
    self.img_SEhist.set_clim(vmin=0, vmax=self.hist2D_map.max()*self.sliders_SEhist['0'].value()/100)
    
    self.canvas_SEhist.draw_idle()
    
def load_filt_multi_data(self):
    options = QFileDialog.Options()
    file_name, _ = QFileDialog.getOpenFileName(self, "Select the multilabelled filtered dataset", "", "All Files (*);;Text Files (*.txt)", options=options)
    
    if file_name:
        # with open(file_name) as json_file:
        #     self.filt_DA_data = load(json_file)
        self.filt_multi_data = read_data(file_name)
        
        self.filt_DA_data = extract_labelled_traces(self.filt_multi_data, 'FRET')
        self.filt_DO_data = extract_labelled_traces(self.filt_multi_data, 'DO')
        self.filt_AO_data = extract_labelled_traces(self.filt_multi_data, 'AO')
        
        self.trace_IDs = []
        for j in range(len(self.filt_DA_data.traces)):
            self.trace_IDs = self.trace_IDs + [j]
        
        self.trace_IDs_DO = []
        for j in range(len(self.filt_DO_data.traces)):
            self.trace_IDs_DO = self.trace_IDs_DO + [j]
                
        self.trace_IDs_AO = []
        for j in range(len(self.filt_AO_data.traces)):
            self.trace_IDs_AO = self.trace_IDs_AO + [j]
        
        
        print('Multi labelled traces dataset loaded!')
        self.button_load_multi_data.setStyleSheet("background-color : green")
        
        
def extract_labelled_traces(traces_dataset, select_key):
    
    traces_dict =  copy.deepcopy(traces_dataset)
    
    nb_traces = len(traces_dict.traces)
    
    toBeDeleted = []
        
    for i in range(nb_traces):
        if select_key in list(traces_dict.traces[i].metadata.keys()):
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
                    
                # case 'Median':  # case 'median' --> keep background traces as it is
                    # if traces_dict.metadata['ALEX'] == 'yes':
                    #     back_DD_trace = np.array(traces_dict.traces[i].channels[3].data[indmin:indmax])
                    # else:
                    #     back_DD_trace = np.array(traces_dict.traces[i].channels[2].data[indmin:indmax])
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
            #self.ax_trace_plot.plot(np.array(self.traces_data[trace_ID].get('Intensity_DD')), 'orange')
            #self.line_DD.set_data(time_data, np.array(self.traces_data[trace_ID].get('Intensity_DD')))
            traces_dict.traces[i].channels[0].data = DD_trace - back_DD_trace
            traces_dict.traces[i].channels[1].data = DA_trace - back_DA_trace
            if traces_dict.metadata['ALEX'] == 'yes':
                traces_dict.traces[i].channels[2].data = AA_trace - back_AA_trace
            
        else:
            toBeDeleted.append(i)
    count = 0        
    for j in toBeDeleted:
        del traces_dict.traces[j-count]
        count += 1
    return traces_dict
            

def load_filt_DA_data(self):
    options = QFileDialog.Options()
    file_name, _ = QFileDialog.getOpenFileName(self, "Select the donor-acceptor filtered dataset", "", "All Files (*);;Text Files (*.txt)", options=options)
    
    if file_name:
        # with open(file_name) as json_file:
        #     self.filt_DA_data = load(json_file)
        self.filt_DA_data = read_data(file_name)
        
        self.trace_IDs = []
        for j in range(len(self.filt_DA_data.traces)):
            if 'indmin' in list(self.filt_DA_data.traces[j].metadata.keys()):
                self.trace_IDs = self.trace_IDs + [j]
        
        
        print('DA traces dataset loaded!')
        self.button_load_DA_data.setStyleSheet("background-color : green")

def load_filt_DO_data(self):
    options = QFileDialog.Options()
    file_name, _ = QFileDialog.getOpenFileName(self, "Select the donor-only filtered dataset", "", "All Files (*);;Text Files (*.txt)", options=options)
    
    if file_name:
        # with open(file_name) as json_file:
        #     self.filt_DO_data = load(json_file)
        self.filt_DO_data = read_data(file_name)
        
        self.trace_IDs_DO = []
        for j in range(len(self.filt_DO_data.traces)):
            if 'indmin' in list(self.filt_DO_data.traces[j].metadata.keys()):
                self.trace_IDs_DO = self.trace_IDs_DO + [j]
        
        
        print('DO traces dataset loaded!')
        self.button_load_DO_data.setStyleSheet("background-color : green")

def load_filt_AO_data(self):
    options = QFileDialog.Options()
    file_name, _ = QFileDialog.getOpenFileName(self, "Select the acceptor-only filtered dataset", "", "All Files (*);;Text Files (*.txt)", options=options)
    
    if file_name:
        # with open(file_name) as json_file:
        #     self.filt_AO_data = load(json_file)
        self.filt_AO_data = read_data(file_name)
        
        self.trace_IDs_AO = []
        for j in range(len(self.filt_AO_data.traces)):
            if 'indmin' in list(self.filt_AO_data.traces[j].metadata.keys()):
                self.trace_IDs_AO = self.trace_IDs_AO + [j]
        
        print('AO traces dataset loaded!')
        self.button_load_AO_data.setStyleSheet("background-color : green")
        
def GUI_plot_FRET_eff(self):
    
    plot_FRET_window = PlotFRETWindow(self.hist2D_map, self.spin_box_bin_map.value(), self.spin_box_extra_range.value(), self.spin_box_nb_states.value())
    plot_FRET_window.show()
    self.plot_windows.append(plot_FRET_window)
    

def GUI_calculate_FRET_eff(self, traces_dict):
    # self.trace_IDs = []
    # for j in range(len(traces_data.traces)):
    #     if 'indmin' in list(self.traces_data.traces[j].metadata.keys()):
    #         self.trace_IDs = self.trace_IDs + j
    #ID_traces = list(traces_dict.keys())
    
    #self.trace_IDs = [str(i) for i in range(len(self.traces_data.traces))]
    
    FRET_conc = np.array([])
    
    FRET_S_conc = np.array([])
    
    alpha = float(self.widget_alpha_value.text())
    delta = float(self.widget_delta_value.text())
    beta = float(self.widget_beta_value.text())
    gamma = float(self.widget_gamma_value.text())
    
    self.FRET_corr_factors = [alpha, delta, beta, gamma]
    
    
    for i in self.trace_IDs:
        # trace_DD = np.array(traces_dict[i]['Intensity_DD'])
        # trace_DA = np.array(traces_dict[i]['Intensity_DA'])
        # trace_AA = np.array(traces_dict[i]['Intensity_AA'])
        
        # trace_DD = np.array(traces_dict.traces[i].channels[0].data)
        # trace_DA = np.array(traces_dict.traces[i].channels[1].data)
        # trace_AA = np.array(traces_dict.traces[i].channels[2].data)
        
        corr_DD_trace = np.array(traces_dict.traces[i].channels[0].data)
        corr_DA_trace = np.array(traces_dict.traces[i].channels[1].data)
        corr_AA_trace = np.array(traces_dict.traces[i].channels[2].data)
        
        # back_method_i = traces_dict.traces[i].metadata['background_correction']
        
        # match back_method_i:
            
        #     case 'None':
        #         back_DD_trace = np.zeros(len(trace_DD))
        #         back_DA_trace = np.zeros(len(trace_DA))
        #         back_AA_trace = np.zeros(len(trace_AA))
                
        #     case 'Median':
        #         back_DD_trace = np.array(traces_dict.traces[i].channels[3].data)
        #         back_DA_trace = np.array(traces_dict.traces[i].channels[4].data)
        #         back_AA_trace = np.array(traces_dict.traces[i].channels[5].data)
                
        #     case 'Total variation':
        #         raw_back_DD = np.array(traces_dict.traces[i].channels[3].data)
        #         max_back_DD = np.max(raw_back_DD)
        #         norm_back_trace_DD = np.array(raw_back_DD) / max_back_DD
        #         TV_back_DD = tvd_2013(norm_back_trace_DD.astype('float'), traces_dict.traces[i].metadata['TV_lambda'])
        #         back_DD_trace = TV_back_DD * max_back_DD
                
        #         raw_back_DA = np.array(traces_dict.traces[i].channels[4].data)
        #         max_back_DA = np.max(raw_back_DA)
        #         norm_back_trace_DA = np.array(raw_back_DA) / max_back_DA
        #         TV_back_DA = tvd_2013(norm_back_trace_DA.astype('float'), traces_dict.traces[i].metadata['TV_lambda'])
        #         back_DA_trace = TV_back_DA * max_back_DA
                
        #         raw_back_AA = np.array(traces_dict.traces[i].channels[5].data)
        #         max_back_AA = np.max(raw_back_AA)
        #         norm_back_trace_AA = np.array(raw_back_AA) / max_back_AA
        #         TV_back_AA = tvd_2013(norm_back_trace_AA.astype('float'), traces_dict.traces[i].metadata['TV_lambda'])
        #         back_AA_trace = TV_back_AA * max_back_AA
                
        #     case 'Min. of TV':
        #         raw_back_DD = np.array(traces_dict.traces[i].channels[3].data)
        #         max_back_DD = np.max(raw_back_DD)
        #         norm_back_trace_DD = np.array(raw_back_DD) / max_back_DD
        #         TV_back_DD = tvd_2013(norm_back_trace_DD.astype('float'), traces_dict.traces[i].metadata['TV_lambda'])
        #         back_DD_trace = np.min(TV_back_DD) * max_back_DD
                
        #         raw_back_DA = np.array(traces_dict.traces[i].channels[4].data)
        #         max_back_DA = np.max(raw_back_DA)
        #         norm_back_trace_DA = np.array(raw_back_DA) / max_back_DA
        #         TV_back_DA = tvd_2013(norm_back_trace_DA.astype('float'), traces_dict.traces[i].metadata['TV_lambda'])
        #         back_DA_trace = np.min(TV_back_DA) * max_back_DA
                
        #         raw_back_AA = np.array(traces_dict.traces[i].channels[5].data)
        #         max_back_AA = np.max(raw_back_AA)
        #         norm_back_trace_AA = np.array(raw_back_AA) / max_back_AA
        #         TV_back_AA = tvd_2013(norm_back_trace_AA.astype('float'), traces_dict.traces[i].metadata['TV_lambda'])
        #         back_AA_trace = np.min(TV_back_AA) * max_back_AA
                
        # corr_DD_trace = trace_DD - back_DD_trace
        # corr_DA_trace = trace_DA - back_DA_trace
        # corr_AA_trace = trace_AA - back_AA_trace
        
        # list_all_positive = ~np.logical_or((trace_DD < 0), (trace_DA < 0), (trace_AA < 0))
        
        # trace_DD = trace_DD [list_all_positive]
        # trace_DA = trace_DA[list_all_positive]
        # trace_AA = trace_AA [list_all_positive]
        

        FRET_eff_i = (corr_DA_trace - alpha * corr_DD_trace - delta * corr_AA_trace) / (gamma * corr_DD_trace + corr_DA_trace - alpha * corr_DD_trace - delta * corr_AA_trace)
        
        FRET_stoi_i = (corr_DA_trace - alpha * corr_DD_trace + gamma * corr_DD_trace - delta * corr_AA_trace) / (gamma * corr_DD_trace + corr_DA_trace - alpha * corr_DD_trace + corr_AA_trace / beta - delta * corr_AA_trace)
        
        FRET_eff_i_list = np.ma.masked_invalid(FRET_eff_i).mask
        
        FRET_stoi_i_list = np.ma.masked_invalid(FRET_stoi_i).mask
        
        list_FRET_OK = ~np.logical_or(FRET_eff_i_list, FRET_stoi_i_list)
        
        #traces_dict[i]['FRET_eff'] = FRET_eff_i[list_FRET_OK]
        
        #traces_dict[i]['FRET_stoi'] = FRET_stoi_i[list_FRET_OK]
        
        # FRET_conc = np.concat((FRET_conc, traces_dict[i]['FRET_eff']))
        
        # FRET_S_conc = np.concat((FRET_S_conc, traces_dict[i]['FRET_stoi']))
        
        FRET_conc = np.concat((FRET_conc, FRET_eff_i[list_FRET_OK]))
        
        FRET_S_conc = np.concat((FRET_S_conc, FRET_stoi_i[list_FRET_OK]))
        
        
    GUI_generate_2D_histogram(self, FRET_conc, FRET_S_conc)
        
        
        # FRET_conc = np.concat((FRET_conc, traces_dict[i]['FRET_eff']))
        
        # FRET_S_conc = np.concat((FRET_S_conc, traces_dict[i]['FRET_stoi']))


def GUI_generate_2D_histogram(self, FRET_conc, FRET_S_conc):
    
    bin_size = self.spin_box_bin_map.value()
    ax_lim_eps = self.spin_box_extra_range.value()
    
    nb_bins = int((1+2*ax_lim_eps)/bin_size)
    
    self.hist2D_map = np.zeros((nb_bins, nb_bins))
    
    compt_nan = 0
    
    # ID_traces = list(traces_dict.keys())
    
    # FRET_conc = np.array([])
    
    # FRET_S_conc = np.array([])
    
    # for i in ID_traces:
    #     FRET_conc = np.concat((FRET_conc, traces_dict[i]['FRET_eff']))
    #     FRET_S_conc = np.concat((FRET_S_conc, traces_dict[i]['FRET_stoi']))
    # x_ax = np.linspace(-ax_lim_eps,1+ax_lim_eps, nb_bins)
    
    extra_bins = int(ax_lim_eps // bin_size)
    
    #FRET_conc = FRET_conc * (FRET_conc >= 0)
    
    #FRET_S_conc = FRET_S_conc * (FRET_S_conc >= 0)
    
    FRET_hist_bin = [FRET_S_conc // bin_size + extra_bins, FRET_conc // bin_size + extra_bins]
    
    for i in range(len(FRET_conc)):
        try:
            x_i = FRET_hist_bin[0][i].astype(int)
            y_i = FRET_hist_bin[1][i].astype(int)
            if (x_i >= 0) and (y_i >=0):
                self.hist2D_map[x_i, y_i] += 1
        except IndexError:
            compt_nan += 1
    
    # hist2D_map[0,:] = 0
    # hist2D_map[-1,:] = 0
    # hist2D_map[:,0] = 0
    # hist2D_map[:,-1] = 0
    
    self.img_SEhist.set_data(self.hist2D_map)
    self.img_SEhist.set_clim(vmin=0, vmax=self.hist2D_map.max())
    self.img_SEhist.set_extent([-ax_lim_eps,1+ax_lim_eps,-ax_lim_eps,1+ax_lim_eps])
     #, interpolation='none', origin='lower', extent=[-ax_lim_eps,1+ax_lim_eps,-ax_lim_eps,1+ax_lim_eps])
    self.ax_SEhist.set_ylim([-ax_lim_eps,1+ax_lim_eps])
    self.ax_SEhist.set_xlim([-ax_lim_eps,1+ax_lim_eps])
    
    self.ax_SEhist.autoscale_view()
    self.canvas_SEhist.draw_idle()

def GUI_calculate_alpha(self, traces_dict):
    
    if self.alpha_corr_checkbox.isChecked():
        
        #ID_traces = list(traces_dict.keys())
        
        FRET_conc = np.array([])
        FRET_S_conc = np.array([])
        
        for i in self.trace_IDs_DO:
            # trace_DD = np.array(traces_dict[i]['Intensity_DD'])
            # trace_DA = np.array(traces_dict[i]['Intensity_DA'])
            # trace_AA = np.array(traces_dict[i]['Intensity_AA'])
            
            # list_all_positive = ~np.logical_or((trace_DD < 0), (trace_DA < 0), (trace_AA < 0))
            
            # trace_DD = trace_DD [list_all_positive]
            # trace_DA = trace_DA[list_all_positive]
            # trace_AA = trace_AA [list_all_positive]
            
            # trace_DD = np.array(traces_dict.traces[i].channels[0].data)
            # trace_DA = np.array(traces_dict.traces[i].channels[1].data)
            # trace_AA = np.array(traces_dict.traces[i].channels[2].data)
            
            corr_DD_trace = np.array(traces_dict.traces[i].channels[0].data)
            corr_DA_trace = np.array(traces_dict.traces[i].channels[1].data)
            corr_AA_trace = np.array(traces_dict.traces[i].channels[2].data)
            
            # back_method_i = traces_dict.traces[i].metadata['background_correction']
            
            # match back_method_i:
                
            #     case 'None':
            #         back_DD_trace = np.zeros(len(trace_DD))
            #         back_DA_trace = np.zeros(len(trace_DA))
            #         back_AA_trace = np.zeros(len(trace_AA))
                    
            #     case 'Median':
            #         back_DD_trace = np.array(traces_dict.traces[i].channels[3].data)
            #         back_DA_trace = np.array(traces_dict.traces[i].channels[4].data)
            #         back_AA_trace = np.array(traces_dict.traces[i].channels[5].data)
                    
            #     case 'Total variation':
            #         raw_back_DD = np.array(traces_dict.traces[i].channels[3].data)
            #         max_back_DD = np.max(raw_back_DD)
            #         norm_back_trace_DD = np.array(raw_back_DD) / max_back_DD
            #         TV_back_DD = tvd_2013(norm_back_trace_DD.astype('float'), traces_dict.traces[i].metadata['TV_lambda'])
            #         back_DD_trace = TV_back_DD * max_back_DD
                    
            #         raw_back_DA = np.array(traces_dict.traces[i].channels[4].data)
            #         max_back_DA = np.max(raw_back_DA)
            #         norm_back_trace_DA = np.array(raw_back_DA) / max_back_DA
            #         TV_back_DA = tvd_2013(norm_back_trace_DA.astype('float'), traces_dict.traces[i].metadata['TV_lambda'])
            #         back_DA_trace = TV_back_DA * max_back_DA
                    
            #         raw_back_AA = np.array(traces_dict.traces[i].channels[5].data)
            #         max_back_AA = np.max(raw_back_AA)
            #         norm_back_trace_AA = np.array(raw_back_AA) / max_back_AA
            #         TV_back_AA = tvd_2013(norm_back_trace_AA.astype('float'), traces_dict.traces[i].metadata['TV_lambda'])
            #         back_AA_trace = TV_back_AA * max_back_AA
                    
            #     case 'Min. of TV':
            #         raw_back_DD = np.array(traces_dict.traces[i].channels[3].data)
            #         max_back_DD = np.max(raw_back_DD)
            #         norm_back_trace_DD = np.array(raw_back_DD) / max_back_DD
            #         TV_back_DD = tvd_2013(norm_back_trace_DD.astype('float'), traces_dict.traces[i].metadata['TV_lambda'])
            #         back_DD_trace = np.min(TV_back_DD) * max_back_DD
                    
            #         raw_back_DA = np.array(traces_dict.traces[i].channels[4].data)
            #         max_back_DA = np.max(raw_back_DA)
            #         norm_back_trace_DA = np.array(raw_back_DA) / max_back_DA
            #         TV_back_DA = tvd_2013(norm_back_trace_DA.astype('float'), traces_dict.traces[i].metadata['TV_lambda'])
            #         back_DA_trace = np.min(TV_back_DA) * max_back_DA
                    
            #         raw_back_AA = np.array(traces_dict.traces[i].channels[5].data)
            #         max_back_AA = np.max(raw_back_AA)
            #         norm_back_trace_AA = np.array(raw_back_AA) / max_back_AA
            #         TV_back_AA = tvd_2013(norm_back_trace_AA.astype('float'), traces_dict.traces[i].metadata['TV_lambda'])
            #         back_AA_trace = np.min(TV_back_AA) * max_back_AA
                    
            # corr_DD_trace = trace_DD - back_DD_trace
            # corr_DA_trace = trace_DA - back_DA_trace
            # corr_AA_trace = trace_AA - back_AA_trace
            
            FRET_eff_i = corr_DA_trace / (corr_DD_trace + corr_DA_trace)
            
            FRET_eff_i_list = np.ma.masked_invalid(FRET_eff_i).mask
            
            FRET_stoi_i = (corr_DA_trace + corr_DD_trace) / (corr_DD_trace + corr_DA_trace + corr_AA_trace)
            
            FRET_stoi_i_list = np.ma.masked_invalid(FRET_stoi_i).mask
            
            list_FRET_OK = ~np.logical_or(FRET_eff_i_list, FRET_stoi_i_list)
            
            #traces_dict[i]['FRET_eff'] = FRET_eff_i[list_FRET_OK]
        
            FRET_conc = np.concat((FRET_conc, FRET_eff_i[list_FRET_OK]))
            
            FRET_S_conc = np.concat((FRET_S_conc, FRET_stoi_i[list_FRET_OK]))
            
        fig, ax = plt.subplots()
        counts, bins, bars = ax.hist(FRET_conc, 100, range=(0,1))  # range=(-1,1))
        plt.close(fig)
        bins_center = 0.5*(bins[1:]+bins[0:-1])
        popt, pcov = curve_fit(GUI_Gauss_data, bins_center[1:-1], counts[1:-1], p0 = [np.max(counts[1:-1]),np.argmax(counts[1:-1]),np.argmax(counts[1:-1])])
        alpha = popt[2]/(1-popt[2])
        
        self.widget_alpha_value.setText(str(alpha))
        
        alpha_QC_window = AlphaCorrWindow(popt, bins_center, counts, alpha, FRET_conc, FRET_S_conc, self.spin_box_bin_map.value(), self.spin_box_extra_range.value())
        alpha_QC_window.show()
        self.plot_windows.append(alpha_QC_window)
    else:
        self.widget_alpha_value.setText('0')
    
    # x = np.linspace(-1, 1, 1000)
    
    # y_gauss = Gauss_data(x, popt[0], popt[1], popt[2])
    
    # fig, ax = plt.subplots()
    # ax.plot(bins_center[1:-1], counts[1:-1],'kx', label='data')
    # label_fit = r'$\alpha$ = '+str(np.round(alpha, 4))
    # ax.plot(x, y_gauss, label=label_fit)
    # plt.title(r'$\alpha$ correction from Donor-Only traces')
    # plt.xlabel('FRET efficiency')
    # plt.ylabel('Counts')
    # ax.legend()

def GUI_calculate_delta(self, traces_dict):
    
    if self.delta_corr_checkbox.isChecked():
        
        #ID_traces = list(traces_dict.keys())
        
        FRET_conc = np.array([])
        FRET_S_conc = np.array([])
        
        for i in self.trace_IDs_AO:
            
            # trace_DD = np.array(traces_dict.traces[i].channels[0].data)
            # trace_DA = np.array(traces_dict.traces[i].channels[1].data)
            # trace_AA = np.array(traces_dict.traces[i].channels[2].data)
            
            corr_DD_trace = np.array(traces_dict.traces[i].channels[0].data)
            corr_DA_trace = np.array(traces_dict.traces[i].channels[1].data)
            corr_AA_trace = np.array(traces_dict.traces[i].channels[2].data)
            
            # back_method_i = traces_dict.traces[i].metadata['background_correction']
            
            # match back_method_i:
                
            #     case 'None':
            #         back_DD_trace = np.zeros(len(trace_DD))
            #         back_DA_trace = np.zeros(len(trace_DA))
            #         back_AA_trace = np.zeros(len(trace_AA))
                    
            #     case 'Median':
            #         back_DD_trace = np.array(traces_dict.traces[i].channels[3].data)
            #         back_DA_trace = np.array(traces_dict.traces[i].channels[4].data)
            #         back_AA_trace = np.array(traces_dict.traces[i].channels[5].data)
                    
            #     case 'Total variation':
            #         raw_back_DD = np.array(traces_dict.traces[i].channels[3].data)
            #         max_back_DD = np.max(raw_back_DD)
            #         norm_back_trace_DD = np.array(raw_back_DD) / max_back_DD
            #         TV_back_DD = tvd_2013(norm_back_trace_DD.astype('float'), traces_dict.traces[i].metadata['TV_lambda'])
            #         back_DD_trace = TV_back_DD * max_back_DD
                    
            #         raw_back_DA = np.array(traces_dict.traces[i].channels[4].data)
            #         max_back_DA = np.max(raw_back_DA)
            #         norm_back_trace_DA = np.array(raw_back_DA) / max_back_DA
            #         TV_back_DA = tvd_2013(norm_back_trace_DA.astype('float'), traces_dict.traces[i].metadata['TV_lambda'])
            #         back_DA_trace = TV_back_DA * max_back_DA
                    
            #         raw_back_AA = np.array(traces_dict.traces[i].channels[5].data)
            #         max_back_AA = np.max(raw_back_AA)
            #         norm_back_trace_AA = np.array(raw_back_AA) / max_back_AA
            #         TV_back_AA = tvd_2013(norm_back_trace_AA.astype('float'), traces_dict.traces[i].metadata['TV_lambda'])
            #         back_AA_trace = TV_back_AA * max_back_AA
                    
            #     case 'Min. of TV':
            #         raw_back_DD = np.array(traces_dict.traces[i].channels[3].data)
            #         max_back_DD = np.max(raw_back_DD)
            #         norm_back_trace_DD = np.array(raw_back_DD) / max_back_DD
            #         TV_back_DD = tvd_2013(norm_back_trace_DD.astype('float'), traces_dict.traces[i].metadata['TV_lambda'])
            #         back_DD_trace = np.min(TV_back_DD) * max_back_DD
                    
            #         raw_back_DA = np.array(traces_dict.traces[i].channels[4].data)
            #         max_back_DA = np.max(raw_back_DA)
            #         norm_back_trace_DA = np.array(raw_back_DA) / max_back_DA
            #         TV_back_DA = tvd_2013(norm_back_trace_DA.astype('float'), traces_dict.traces[i].metadata['TV_lambda'])
            #         back_DA_trace = np.min(TV_back_DA) * max_back_DA
                    
            #         raw_back_AA = np.array(traces_dict.traces[i].channels[5].data)
            #         max_back_AA = np.max(raw_back_AA)
            #         norm_back_trace_AA = np.array(raw_back_AA) / max_back_AA
            #         TV_back_AA = tvd_2013(norm_back_trace_AA.astype('float'), traces_dict.traces[i].metadata['TV_lambda'])
            #         back_AA_trace = np.min(TV_back_AA) * max_back_AA
                    
            # corr_DD_trace = trace_DD - back_DD_trace
            # corr_DA_trace = trace_DA - back_DA_trace
            # corr_AA_trace = trace_AA - back_AA_trace
            # trace_DD = np.array(traces_dict[i]['Intensity_DD'])
            # trace_DA = np.array(traces_dict[i]['Intensity_DA'])
            # trace_AA = np.array(traces_dict[i]['Intensity_AA'])
            
            # list_all_positive = ~np.logical_or((trace_DD < 0), (trace_DA < 0), (trace_AA < 0))
            
            # trace_DD = trace_DD [list_all_positive]
            # trace_DA = trace_DA[list_all_positive]
            # trace_AA = trace_AA [list_all_positive]
            
            FRET_eff_i = corr_DA_trace / (corr_DD_trace + corr_DA_trace)
            
            FRET_eff_i_list = np.ma.masked_invalid(FRET_eff_i).mask
            
            FRET_stoi_i = (corr_DA_trace + corr_DD_trace) / (corr_DD_trace + corr_DA_trace + corr_AA_trace)
            
            FRET_stoi_i_list = np.ma.masked_invalid(FRET_stoi_i).mask
            
            list_FRET_OK = ~np.logical_or(FRET_eff_i_list, FRET_stoi_i_list)
            
            #traces_dict[i]['FRET_stoi'] = FRET_stoi_i[list_FRET_OK]
            
            FRET_conc = np.concat((FRET_conc, FRET_eff_i[list_FRET_OK]))
        
            FRET_S_conc = np.concat((FRET_S_conc, FRET_stoi_i[list_FRET_OK]))
            
        fig, ax = plt.subplots()
        counts, bins, bars = ax.hist(FRET_S_conc, 100, range=(0,1))
        plt.close(fig)
        bins_center = 0.5*(bins[1:]+bins[0:-1])
        popt, pcov = curve_fit(GUI_Gauss_data, bins_center[1:-1], counts[1:-1], p0 = [np.max(counts[1:-1]),np.argmax(counts[1:-1]),np.argmax(counts[1:-1])])
        delta = popt[2]/(1-popt[2])
        
        self.widget_delta_value.setText(str(delta))
        
        delta_QC_window = DeltaCorrWindow(popt, bins_center, counts, delta, FRET_conc, FRET_S_conc, self.spin_box_bin_map.value(), self.spin_box_extra_range.value())
        delta_QC_window.show()
        self.plot_windows.append(delta_QC_window)
    else:
        self.widget_delta_value.setText('0')
        
def GUI_calculate_beta_gamma(self, traces_dict):
    peak_coord = GUI_fit_Gauss_spot(self.hist2D_map, nb_gauss = self.spin_box_nb_states.value(), bin_size = self.spin_box_bin_map.value(), ax_lim_eps = self.spin_box_extra_range.value())
    
    nb_bins = int((1+2*self.spin_box_extra_range.value())/self.spin_box_bin_map.value())
    
    x_ax = np.linspace(-self.spin_box_extra_range.value(),1+self.spin_box_extra_range.value(), nb_bins)
    
    b, a = np.polyfit(x_ax[peak_coord[:,0].astype(int)], 1/x_ax[peak_coord[:,1].astype(int)], 1)
    
    beta = a + b -1
    gamma = (a - 1) / (a + b - 1)
    
    self.widget_beta_value.setText(str(beta))
    self.widget_gamma_value.setText(str(gamma))

def GUI_twoD_Gaussian(xy, gauss_param, nb_gauss):
    y, x = xy
    g = gauss_param[0] + np.zeros(x.shape)
    for i in range(nb_gauss):
        g = g + gauss_param[5*i+1]*np.exp( - ((x-float(gauss_param[5*i+2]))**2/(2*gauss_param[5*i+4]**2) + (y-float(gauss_param[5*i+3]))**2/(2*gauss_param[5*i+5]**2)))
    return g.ravel()

def GUI_error_func(gauss_param, nb_gauss, x, y):
    return GUI_twoD_Gaussian(x, gauss_param, nb_gauss) - y

def GUI_fit_Gauss_spot(img_spot, nb_gauss = 2, bin_size = 0.01, ax_lim_eps = 0.2):
    x_len = img_spot.shape[0]
    y_len = img_spot.shape[1]
    x_vec = np.linspace(0, x_len-1, x_len)
    y_vec = np.linspace(0, y_len-1, y_len)
    x, y = np.meshgrid(y_vec, x_vec)
    max_img = np.max(img_spot)
    noise_img = np.median(img_spot)
    # initial_guess = (max_img-noise_img, x_len*0.5, x_len*0.6, 10, 10, max_img-noise_img, x_len*0.5, x_len*0.2, 10, 10, noise_img)
    gauss_param_estim = [noise_img]
    Fret_eff_estim = np.linspace(0.3, 0.7, nb_gauss)
    for i in range(nb_gauss):
        gauss_param_estim.append(max_img-noise_img)
        gauss_param_estim.append(x_len*0.5)
        gauss_param_estim.append(x_len*Fret_eff_estim[i])
        gauss_param_estim.append(10)
        gauss_param_estim.append(10)
    initial_guess = tuple(gauss_param_estim)
    # custom_gaussian = lambda x, offset, mu: twoD_Gaussian(x, offset, mu, nb_gauss)
    try:
        # popt, pcov = opt.curve_fit(custom_gaussian, (x, y), img_spot.ravel(), p0=initial_guess)
        XX = (x, y)
        output = least_squares(GUI_error_func, x0=initial_guess, jac='2-point', bounds=(0, np.inf), method='trf', ftol=1e-08,
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
        
    
    
    # nb_bins = int((1+2*ax_lim_eps)/bin_size)
    
    # x_ax = np.linspace(-ax_lim_eps,1+ax_lim_eps, nb_bins)
        
    # plt.figure()
    # plt.imshow(img_spot,origin='lower', interpolation='none', extent=[-ax_lim_eps,1+ax_lim_eps,-ax_lim_eps,1+ax_lim_eps])
    # #plt.plot(peak_coord[:,0]/len(img_spot), peak_coord[:,1]/len(img_spot), 'ro')
    # plt.plot(x_ax[peak_coord[:,0].astype(int)], x_ax[peak_coord[:,1].astype(int)], 'ro')
    # plt.xlabel('FRET efficiency')
    # plt.ylabel('FRET stoichiometry')
    # plt.title(r'$\alpha$-$\delta$ corrected FRET histogram')
    
    return peak_coord

def GUI_Gauss_data(x, A, B, C):
    return A * np.exp(-(x-C)**2/B**2)
