# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 12:01:07 2026

@author: molcre0000
"""

import sys
from PyQt5.QtWidgets import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.widgets import SpanSelector
from PyQt5.QtGui import QFont, QPainter
from PyQt5.QtCore import Qt, pyqtSignal, QSignalBlocker, QAbstractTableModel, QRectF
from openfret import read_data
from chromatic_aberrations_correction import *
from utils import *
from traces_extractor import *
from drift_correction import *
from GUI_HMM_utils import *
import copy
from PySMACKS_GUI_PyQT5_SEhist import extract_labelled_traces


class TbT_ratesWindow(QMainWindow):
    def __init__(self, stack_rate_mat, framerate):
        super().__init__()
        self.setWindowTitle("TbT kinetic rates")
        self.resize(700, 700)
        
        
        central_widget = QWidget()
        layout1 = QVBoxLayout(central_widget)
        figure_TbT_rates = plt.figure()
        canvas_TbT_rates = FigureCanvas(figure_TbT_rates)
        toolbar_TbT_rates = NavigationToolbar(canvas_TbT_rates, self)
        layout1.addWidget(toolbar_TbT_rates)
        layout1.addWidget(canvas_TbT_rates)
        #self.central_widget.addLayout(layout1)
        self.setCentralWidget(central_widget)
        
        nb_hidden_states = stack_rate_mat.shape[1]
        
        # low FRET is state 0, high FRET is state 1
        
        compt = 0
        rate_dict = {}
        
        ax_TbT_rates = figure_TbT_rates.add_subplot(111)
        
        for i in range(nb_hidden_states):
            for j in range(nb_hidden_states):
                if i != j:
                    label_i = r'$k_{' + str(i) + str(j) + r'}$'
                    rate_dict[label_i] = stack_rate_mat[:,i,j]/framerate
                    compt += 1
        ax_TbT_rates.boxplot(rate_dict.values(), labels=list(rate_dict.keys()))
        plt.ylabel(r'Kinetic rate ($s^{-1}$)')
        #ax_TbT_rates.legend()
        canvas_TbT_rates.draw_idle()
        
class EnsHMM_ratesWindow(QMainWindow):
    def __init__(self, stack_rate_mat, Ens_trans_mat, framerate):
        super().__init__()
        self.setWindowTitle("EnsHMM kinetic rates")
        self.resize(700, 700)
        
        
        central_widget = QWidget()
        layout1 = QVBoxLayout(central_widget)
        figure_EnsHMM_rates = plt.figure()
        canvas_EnsHMM_rates = FigureCanvas(figure_EnsHMM_rates)
        toolbar_EnsHMM_rates = NavigationToolbar(canvas_EnsHMM_rates, self)
        layout1.addWidget(toolbar_EnsHMM_rates)
        layout1.addWidget(canvas_EnsHMM_rates)
        #self.central_widget.addLayout(layout1)
        self.setCentralWidget(central_widget)
        
        nb_hidden_states = stack_rate_mat.shape[1]
        
        # low FRET is state 0, high FRET is state 1
        
        compt = 0
        rate_dict = {}
        
        ax_EnsHMM_rates = figure_EnsHMM_rates.add_subplot(111)
        
        for i in range(nb_hidden_states):
            for j in range(nb_hidden_states):
                if i != j:
                    label_i = r'$k_{' + str(i) + str(j) + r'}$'
                    rate_dict[label_i] = stack_rate_mat[:,i,j]/framerate
                    compt += 1
        ax_EnsHMM_rates.boxplot(rate_dict.values(), labels=list(rate_dict.keys()))
        
        compt_EnsHMM = 1
        for i in range(nb_hidden_states):
            for j in range(nb_hidden_states):
                if i != j:
                        if compt_EnsHMM == 1:
                            ax_EnsHMM_rates.scatter(compt_EnsHMM, Ens_trans_mat[i,j]/framerate, color='red', marker='x', zorder=3, label='Ens. HMM')
                            compt_EnsHMM += 1
                        else:
                            ax_EnsHMM_rates.scatter(compt_EnsHMM, Ens_trans_mat[i,j]/framerate, color='red', marker='x', zorder=3)
                            compt_EnsHMM += 1
        plt.ylabel(r'Kinetic rate ($s^{-1}$)')
        ax_EnsHMM_rates.legend()
        canvas_EnsHMM_rates.draw_idle()

def create_kinetic_tab(self):
    
    self.isHMMDictLoaded = 0
    #self.fill_green = None
    #self.fill_red = None
    #self.fill_redFRETEff = None
    #self.fill_greenFRETEff = None
    
    self.fill_color = None
    self.fill_FRETEff = None
    
    
    self.stack_rate_mat = None
    self.Ens_trans_mat = None
    self.HMM_colors = ['lightgreen', 'lightcoral', 'lightblue', 'wheat', 'lavender', 'lightpink', 'lightyellow', 'mistyrose', 'powderblue', 'paleturquoise']
    
    
    tab1 = QWidget()
    layout1 = QVBoxLayout()
    
    layout_load_HMM_data = QHBoxLayout()
    
    self.button_load_HMM_data = QPushButton("Load DA dataset")
    self.button_load_HMM_data.setToolTip('Select the filtered DA traces dataset')
    self.button_load_HMM_data.clicked.connect(lambda: load_HMM_data(self))
    layout_load_HMM_data.addWidget(self.button_load_HMM_data)
    
    label = QLabel("Nb. states")
    layout_load_HMM_data.addWidget(label)
    self.spin_box_HMM_nb_states = QSpinBox()
    self.spin_box_HMM_nb_states.setToolTip('Select the number of states')
    self.spin_box_HMM_nb_states.setMinimum(2)  # Set minimum value
    self.spin_box_HMM_nb_states.setMaximum(100)   # Set maximum value
    self.spin_box_HMM_nb_states.setValue(2)
    layout_load_HMM_data.addWidget(self.spin_box_HMM_nb_states)
    
    label = QLabel("Delay (s)")
    layout_load_HMM_data.addWidget(label)
    self.widget_frame_delay = QLineEdit()
    self.widget_frame_delay.setToolTip('Delay between successive frames (in seconds)')
    self.widget_frame_delay.setText('0.1')
    layout_load_HMM_data.addWidget(self.widget_frame_delay)
    
    self.HMM_display_checkbox = QCheckBox("Pretty HMM")
    self.HMM_display_checkbox.setToolTip('Switch HMM display')
    self.HMM_display_checkbox.setChecked(True)
    self.HMM_display_checkbox.stateChanged.connect(lambda: display_HMM_plot(self))
    layout_load_HMM_data.addWidget(self.HMM_display_checkbox)
    
    layout1.addLayout(layout_load_HMM_data)
    
    layout_run_TbT_HMM = QHBoxLayout()
    self.button_run_TbT_HMM = QPushButton("Run TbT HMM")
    self.button_run_TbT_HMM.setToolTip('Run the Trace-by-Trace HMM analysis')
    self.button_run_TbT_HMM.clicked.connect(lambda: GUI_run_TbT_HMM(self))
    layout_run_TbT_HMM.addWidget(self.button_run_TbT_HMM)
    
    button_plot_TbT_rates = QPushButton("Plot TbT rates")
    button_plot_TbT_rates.setToolTip('Plot Trace-by-Trace rates as boxplots')
    button_plot_TbT_rates.clicked.connect(lambda: plot_TbT_rates_boxplot(self))
    layout_run_TbT_HMM.addWidget(button_plot_TbT_rates)
    
    layout1.addLayout(layout_run_TbT_HMM)
    
    layout_run_EnsHMM = QHBoxLayout()
    self.button_run_EnsHMM = QPushButton("Run Ensemble HMM")
    self.button_run_EnsHMM.setToolTip('Run the Ensemble HMM analysis')
    self.button_run_EnsHMM.clicked.connect(lambda: GUI_run_EnsHMM_HMM(self))
    layout_run_EnsHMM.addWidget(self.button_run_EnsHMM)
    
    button_plot_EnsHMM_rates = QPushButton("Plot EnsHMM rates")
    button_plot_EnsHMM_rates.setToolTip('Plot Ensemble HMM rates as boxplots')
    button_plot_EnsHMM_rates.clicked.connect(lambda: plot_EnsHMM_rates_boxplot(self))
    layout_run_EnsHMM.addWidget(button_plot_EnsHMM_rates)
    
    layout1.addLayout(layout_run_EnsHMM)
    
    ##################### HLayout for displaying individual channels
    
    layout_channel_checkbox = QHBoxLayout()
    
    
    # Checkbox for DD channel
    self.DD_HMM_checkbox = QCheckBox("DD")
    self.DD_HMM_checkbox.setToolTip('Display Donor-Donor channel')
    self.DD_HMM_checkbox.setChecked(True)
    self.DD_HMM_checkbox.stateChanged.connect(lambda: display_HMM_plot(self))
    layout_channel_checkbox.addWidget(self.DD_HMM_checkbox)
    
    # Checkbox for DA channel
    self.DA_HMM_checkbox = QCheckBox("DA")
    self.DA_HMM_checkbox.setToolTip('Display Donor-Acceptor channel')
    self.DA_HMM_checkbox.setChecked(True)
    self.DA_HMM_checkbox.stateChanged.connect(lambda: display_HMM_plot(self))
    layout_channel_checkbox.addWidget(self.DA_HMM_checkbox)
    
    # Checkbox for AA channel
    self.AA_HMM_checkbox = QCheckBox("AA")
    self.AA_HMM_checkbox.setToolTip('Display Acceptor-Acceptor channel')
    self.AA_HMM_checkbox.setChecked(False)
    self.AA_HMM_checkbox.stateChanged.connect(lambda: display_HMM_plot(self))
    layout_channel_checkbox.addWidget(self.AA_HMM_checkbox)
    
    label = QLabel("Trace ID")
    layout_channel_checkbox.addWidget(label)
    self.choose_trace_ID_HMM = QComboBox()
    layout_channel_checkbox.addWidget(self.choose_trace_ID_HMM)
    button_plot_TbT_HMM = QPushButton("Plot")
    button_plot_TbT_HMM.setToolTip('Plot the TbT HMM result for the selected trace ID')
    #button_plot_TbT_HMM.clicked.connect(lambda: plot_HMM_button_function(self))
    button_plot_TbT_HMM.clicked.connect(lambda: display_HMM_plot(self))
    layout_channel_checkbox.addWidget(button_plot_TbT_HMM)
    
    button_left_arrow = QPushButton()
    button_left_arrow.setToolTip('Plot previous trace')
    button_left_arrow.clicked.connect(lambda: move_left_HMM(self))
    arrow_left_icon = QApplication.style().standardIcon(QStyle.SP_ArrowLeft)
    button_left_arrow.setIcon(arrow_left_icon)
    layout_channel_checkbox.addWidget(button_left_arrow)
    button_right_arrow = QPushButton()
    button_right_arrow.setToolTip('Plot next trace')
    button_right_arrow.clicked.connect(lambda: move_right_HMM(self))
    arrow_right_icon = QApplication.style().standardIcon(QStyle.SP_ArrowRight)
    button_right_arrow.setIcon(arrow_right_icon)
    layout_channel_checkbox.addWidget(button_right_arrow)
    
    layout1.addLayout(layout_channel_checkbox)
    
    self.figure_TbT_HMM = plt.figure()
    self.canvas_TbT_HMM = FigureCanvas(self.figure_TbT_HMM)
    self.toolbar_HMM = NavigationToolbar(self.canvas_TbT_HMM, self)
    layout1.addWidget(self.toolbar_HMM)
    self.ax_TbT_HMM = self.figure_TbT_HMM.add_subplot(211)
    self.line_TbT_HMM_DD, = self.ax_TbT_HMM.plot([],[], 'orange', label='DD')
    self.line_TbT_HMM_DA, = self.ax_TbT_HMM.plot([],[], 'red', label='DA')
    self.line_TbT_HMM_AA, = self.ax_TbT_HMM.plot([],[], 'gray', label='AA')
    self.line_TbT_HMM_TbT, = self.ax_TbT_HMM.plot([],[], 'black', label='HMM')
    plt.xlabel('Time (s)')
    plt.ylabel('Intensity (A.U.)')
    plt.legend(loc='upper right')
    
    self.ax_TbT_HMM_FRETEff = self.figure_TbT_HMM.add_subplot(212, sharex=self.ax_TbT_HMM)
    self.line_TbT_HMM_FRETEff, = self.ax_TbT_HMM_FRETEff.plot([],[], 'blue')
    self.line_TbT_HMM_FRETEff_predict, = self.ax_TbT_HMM_FRETEff.plot([],[], 'black')
    plt.xlabel('Time (s)')
    plt.ylabel('FRET eff.')
    
    self.canvas_TbT_HMM.draw_idle()
    layout1.addWidget(self.canvas_TbT_HMM, 5)
    
    
    tab1.setLayout(layout1)
    self.tabs.addTab(tab1, "Kinetic")
    

def load_HMM_data(self):
    self.trace_IDs_HMM = []
    options = QFileDialog.Options()
    file_name, _ = QFileDialog.getOpenFileName(self, "Select the donor-acceptor filtered dataset", "", "All Files (*);;Text Files (*.json.zip)", options=options)
    
    if file_name:
        # with open(file_name) as json_file:
        #     self.filt_HMM_data = load(json_file)
        self.traces_data_raw =  read_data(file_name)
        print('DA traces dataset loaded!')
        
        # for i in range(len(self.traces_data_raw.traces)):
        #     if 'FRET' in list(self.traces_data_raw.traces[i].metadata.keys()):
        #         self.trace_IDs_HMM = self.trace_IDs_HMM + [str(i)]
        
        #self.trace_IDs_HMM = list(self.filt_HMM_data.keys())
        # self.choose_trace_ID_HMM.addItems(self.trace_IDs_HMM)
        # self.button_load_HMM_data.setStyleSheet("background-color : green")
        
        self.traces_data = extract_labelled_traces(self.traces_data_raw, 'FRET')
        
        self.choose_trace_ID_HMM.addItems([str(i) for i in range(len(self.traces_data.traces))])
        self.button_load_HMM_data.setStyleSheet("background-color : green")
        

def GUI_run_TbT_HMM(self):
    #self.TbT_dataset, self.stack_rate_mat, self.stack_prob_init = TbT_HMM_pipeline(self.filt_HMM_data, self.trace_IDs_HMM, nb_hidden_states = self.spin_box_HMM_nb_states.value())
    self.TbT_dataset, self.stack_rate_mat, self.stack_prob_init = TbT_HMM_pipeline(self.traces_data, 
                                                                                   nb_hidden_states = self.spin_box_HMM_nb_states.value())
    
    self.button_run_TbT_HMM.setStyleSheet("background-color : green")
    self.isHMMDictLoaded = 1
    
    display_HMM_plot(self)
    
def display_HMM_plot(self):
    
    trace_ID = self.choose_trace_ID_HMM.currentText()
    
    # time_data = np.linspace(0, len(np.array(self.TbT_dataset[trace_ID].get('Intensity_DD')))-1, 
    #                         len(np.array(self.TbT_dataset[trace_ID].get('Intensity_DD'))))*float(self.widget_frame_delay.text())
    
    time_data = np.linspace(0,
                            (len(np.array(self.TbT_dataset.traces[int(trace_ID)].channels[0].data))-1)*float(self.widget_frame_delay.text()),
                            len(np.array(self.TbT_dataset.traces[int(trace_ID)].channels[0].data)))
    
    # scale_factor_HMM = np.percentile(np.concat((np.array(self.TbT_dataset[trace_ID].get('Intensity_DD')), 
    #                                             np.array(self.TbT_dataset[trace_ID].get('Intensity_DA')))), 95)
    
    scale_factor_HMM = np.percentile(np.concat((np.array(self.TbT_dataset.traces[int(trace_ID)].channels[0].data), 
                                                np.array(self.TbT_dataset.traces[int(trace_ID)].channels[1].data))), 95)
    
    if self.DD_HMM_checkbox.isChecked():
        #self.line_TbT_HMM_DD.set_data(time_data, np.array(self.TbT_dataset[trace_ID].get('Intensity_DD')))
        self.line_TbT_HMM_DD.set_data(time_data, np.array(self.TbT_dataset.traces[int(trace_ID)].channels[0].data))
        self.line_TbT_HMM_DD.set_label('DD')
    else:
        self.line_TbT_HMM_DD.set_data([],[])
        self.line_TbT_HMM_DD.set_label('_nolegend_')
    
    if self.DA_HMM_checkbox.isChecked():
        self.line_TbT_HMM_DA.set_data(time_data, np.array(self.TbT_dataset.traces[int(trace_ID)].channels[1].data))
        self.line_TbT_HMM_DA.set_label('DA')
    else:
        self.line_TbT_HMM_DA.set_data([],[])
        self.line_TbT_HMM_DA.set_label('_nolegend_')
    
    if self.AA_HMM_checkbox.isChecked():
        self.line_TbT_HMM_AA.set_data(time_data, np.array(self.TbT_dataset.traces[int(trace_ID)].channels[2].data))
        self.line_TbT_HMM_AA.set_label('AA')
    else:
        self.line_TbT_HMM_AA.set_data([],[])
        self.line_TbT_HMM_AA.set_label('_nolegend_')
        
    # self.line_TbT_HMM_FRETEff.set_data(time_data, 
    #                                    np.array(self.TbT_dataset[trace_ID].get('Intensity_DA')) / (np.array(self.TbT_dataset[trace_ID].get('Intensity_DA')) + np.array(self.TbT_dataset[trace_ID].get('Intensity_DD'))))
    
    self.line_TbT_HMM_FRETEff.set_data(time_data, 
                                       np.array(self.TbT_dataset.traces[int(trace_ID)].channels[1].data) / (np.array(self.TbT_dataset.traces[int(trace_ID)].channels[1].data) + np.array(self.TbT_dataset.traces[int(trace_ID)].channels[0].data)))
    
    #self.ax_TbT_HMM_FRETEff.set_xlim(int(-len(np.array(self.TbT_dataset[trace_ID].get('Intensity_DD')))*0.05*float(self.widget_frame_delay.text())), int(len(np.array(self.TbT_dataset[trace_ID].get('Intensity_DD')))*1.05*float(self.widget_frame_delay.text())))
    
    # if self.HMM_display_checkbox.isChecked():
    
    
    #     if self.fill_green is not None:
    #         self.fill_green.remove()
    #         self.fill_red.remove()
    #         self.fill_greenFRETEff.remove()
    #         self.fill_redFRETEff.remove()
    #     self.line_TbT_HMM_TbT.set_data([],[])
    #     self.line_TbT_HMM_FRETEff_predict.set_data([],[])
    #     self.fill_green = self.ax_TbT_HMM.fill_between(
    #         time_data, 0, 1,            #self.ax_TbT_HMM.get_ylim()[0], self.ax_TbT_HMM.get_ylim()[1]
    #         where=np.array(self.TbT_dataset[trace_ID]['HMM_param'].get('predict')) == 0,
    #         color='lightgreen',
    #         alpha=0.3,
    #         interpolate=False,
    #         transform=self.ax_TbT_HMM.get_xaxis_transform()
    #     )
    
    #     self.fill_red = self.ax_TbT_HMM.fill_between(
    #         time_data, 0, 1,
    #         where=np.array(self.TbT_dataset[trace_ID]['HMM_param'].get('predict')) == 1,
    #         color='lightcoral',
    #         alpha=0.3,
    #         interpolate=False,
    #         transform=self.ax_TbT_HMM.get_xaxis_transform()
    #     )
        
    #     self.fill_greenFRETEff = self.ax_TbT_HMM_FRETEff.fill_between(
    #         time_data, 0, 1,            #self.ax_TbT_HMM.get_ylim()[0], self.ax_TbT_HMM.get_ylim()[1]
    #         where=np.array(self.TbT_dataset[trace_ID]['HMM_param'].get('predict')) == 0,
    #         color='lightgreen',
    #         alpha=0.3,
    #         interpolate=False,
    #         transform=self.ax_TbT_HMM_FRETEff.get_xaxis_transform()
    #     )
    
    #     self.fill_redFRETEff = self.ax_TbT_HMM_FRETEff.fill_between(
    #         time_data, 0, 1,
    #         where=np.array(self.TbT_dataset[trace_ID]['HMM_param'].get('predict')) == 1,
    #         color='lightcoral',
    #         alpha=0.3,
    #         interpolate=False,
    #         transform=self.ax_TbT_HMM_FRETEff.get_xaxis_transform()
    #     )
    
    # else:
    #     if self.fill_green is not None:
    #         self.fill_green.remove()
    #         self.fill_red.remove()
    #         self.fill_greenFRETEff.remove()
    #         self.fill_redFRETEff.remove()
    #         self.fill_green = None
    #         self.fill_red = None
    #         self.fill_greenFRETEff = None
    #         self.fill_redFRETEff = None
    #     self.line_TbT_HMM_TbT.set_data(time_data, scale_factor_HMM * np.array(self.TbT_dataset[trace_ID]['HMM_param'].get('predict')))
    #     self.line_TbT_HMM_FRETEff_predict.set_data(time_data, np.array(self.TbT_dataset[trace_ID]['HMM_param'].get('predict')) / (self.spin_box_HMM_nb_states.value() - 1))
    
    if self.HMM_display_checkbox.isChecked():


        if self.fill_color is not None:
            for j in range(len(self.fill_color)):
                self.fill_color[j].remove()
                self.fill_FRETEff[j].remove()
        self.fill_color = [None for i in range(self.stack_rate_mat.shape[1])]
        self.fill_FRETEff = [None for i in range(self.stack_rate_mat.shape[1])]
        for j in range(len(self.fill_color)):
            self.fill_color[j] = self.ax_TbT_HMM.fill_between(
                time_data, 0, 1,            #self.ax_TbT_HMM.get_ylim()[0], self.ax_TbT_HMM.get_ylim()[1]
                #where=np.array(self.TbT_dataset[trace_ID]['HMM_param'].get('predict')) == j,
                where=np.array(self.TbT_dataset.traces[int(trace_ID)].metadata['HMM_param']['predict']) == j,
                color=self.HMM_colors[j],
                alpha=0.5,
                interpolate=False,
                transform=self.ax_TbT_HMM.get_xaxis_transform()
            )
            
            self.fill_FRETEff[j] = self.ax_TbT_HMM_FRETEff.fill_between(
                time_data, 0, 1,            #self.ax_TbT_HMM.get_ylim()[0], self.ax_TbT_HMM.get_ylim()[1]
                #where=np.array(self.TbT_dataset[trace_ID]['HMM_param'].get('predict')) == j,
                where=np.array(self.TbT_dataset.traces[int(trace_ID)].metadata['HMM_param']['predict']) == j,
                color=self.HMM_colors[j],
                alpha=0.5,
                interpolate=False,
                transform=self.ax_TbT_HMM_FRETEff.get_xaxis_transform()
            )
        
        self.line_TbT_HMM_TbT.set_data([],[])
        self.line_TbT_HMM_FRETEff_predict.set_data([],[])

    else:
        if self.fill_color is not None:
            for j in range(len(self.fill_color)):
                self.fill_color[j].remove()
                self.fill_FRETEff[j].remove()
        self.fill_color = None
        self.fill_FRETEff = None
        #self.line_TbT_HMM_TbT.set_data(time_data, scale_factor_HMM * np.array(self.TbT_dataset[trace_ID]['HMM_param'].get('predict')))
        self.line_TbT_HMM_TbT.set_data(time_data, scale_factor_HMM * np.array(self.TbT_dataset.traces[int(trace_ID)].metadata['HMM_param']['predict']))
        #self.line_TbT_HMM_FRETEff_predict.set_data(time_data, np.array(self.TbT_dataset[trace_ID]['HMM_param'].get('predict')) / (self.spin_box_HMM_nb_states.value() - 1))
        self.line_TbT_HMM_FRETEff_predict.set_data(time_data, np.array(self.TbT_dataset.traces[int(trace_ID)].metadata['HMM_param']['predict']) / (self.spin_box_HMM_nb_states.value() - 1))
    
    
    self.ax_TbT_HMM.legend(loc='upper right')
    self.ax_TbT_HMM.relim()
    self.ax_TbT_HMM.autoscale_view(scalex=False, scaley=True)
    self.ax_TbT_HMM.set_xlim(np.min(time_data), np.max(time_data))
    
    #self.ax_TbT_HMM_FRETEff.autoscale_view(scalex=True, scaley=False)
    self.ax_TbT_HMM_FRETEff.set_ylim(-0.1, 1.1)
    self.ax_TbT_HMM_FRETEff.set_xlim(np.min(time_data), np.max(time_data))
    
    self.canvas_TbT_HMM.draw_idle()
    
    self.toolbar_HMM.update()


def plot_HMM_button_function(self):
    trace_ID = self.choose_trace_ID_HMM.currentText()
    
    time_data = np.linspace(0, len(np.array(self.TbT_dataset[trace_ID].get('Intensity_DD')))-1, len(np.array(self.TbT_dataset[trace_ID].get('Intensity_DD'))))*float(self.widget_frame_delay.text())
    
    scale_factor_HMM = np.percentile(np.concat((np.array(self.TbT_dataset[trace_ID].get('Intensity_DD')), np.array(self.TbT_dataset[trace_ID].get('Intensity_DA')))), 95)
    
    if self.DD_HMM_checkbox.isChecked():
        self.line_TbT_HMM_DD.set_data(time_data, np.array(self.TbT_dataset[trace_ID].get('Intensity_DD')))
        self.line_TbT_HMM_DD.set_label('DD')
    else:
        self.line_TbT_HMM_DD.set_data([],[])
        self.line_TbT_HMM_DD.set_label('_nolegend_')
    
    if self.DA_HMM_checkbox.isChecked():
        self.line_TbT_HMM_DA.set_data(time_data, np.array(self.TbT_dataset[trace_ID].get('Intensity_DA')))
        self.line_TbT_HMM_DA.set_label('DA')
    else:
        self.line_TbT_HMM_DA.set_data([],[])
        self.line_TbT_HMM_DA.set_label('_nolegend_')
    
    if self.AA_HMM_checkbox.isChecked():
        self.line_TbT_HMM_AA.set_data(time_data, np.array(self.TbT_dataset[trace_ID].get('Intensity_AA')))
        self.line_TbT_HMM_AA.set_label('AA')
    else:
        self.line_TbT_HMM_AA.set_data([],[])
        self.line_TbT_HMM_AA.set_label('_nolegend_')
        
    self.line_TbT_HMM_FRETEff.set_data(time_data, np.array(self.TbT_dataset[trace_ID].get('Intensity_DA')) / (np.array(self.TbT_dataset[trace_ID].get('Intensity_DA')) + np.array(self.TbT_dataset[trace_ID].get('Intensity_DD'))))
    
    
    self.ax_TbT_HMM.legend(loc='upper right')
    
    self.ax_TbT_HMM.set_xlim(int(-len(np.array(self.TbT_dataset[trace_ID].get('Intensity_DD')))*0.05*float(self.widget_frame_delay.text())), int(len(np.array(self.TbT_dataset[trace_ID].get('Intensity_DD')))*1.05*float(self.widget_frame_delay.text())))
    self.ax_TbT_HMM.set_ylim(np.min(np.concat((np.array(self.TbT_dataset[trace_ID].get('Intensity_DD')), np.array(self.TbT_dataset[trace_ID].get('Intensity_DA'))))), np.max(np.concat((np.array(self.TbT_dataset[trace_ID].get('Intensity_DD')), np.array(self.TbT_dataset[trace_ID].get('Intensity_DA'))))))
    
    
    self.ax_TbT_HMM_FRETEff.autoscale_view(scalex=True, scaley=False)
    self.ax_TbT_HMM_FRETEff.set_ylim(-1.5, 1.5)
    self.ax_TbT_HMM_FRETEff.set_xlim(int(-len(np.array(self.TbT_dataset[trace_ID].get('Intensity_DD')))*0.05*float(self.widget_frame_delay.text())), int(len(np.array(self.TbT_dataset[trace_ID].get('Intensity_DD')))*1.05*float(self.widget_frame_delay.text())))
    
    # if self.HMM_display_checkbox.isChecked():
    #     if self.fill_green is not None:
    #         self.fill_green.remove()
    #         self.fill_red.remove()
    #         self.fill_greenFRETEff.remove()
    #         self.fill_redFRETEff.remove()
    #     self.line_TbT_HMM_TbT.set_data([],[])
    #     self.line_TbT_HMM_FRETEff_predict.set_data([],[])
    #     self.fill_green = self.ax_TbT_HMM.fill_between(
    #         time_data, 0, 1,            #self.ax_TbT_HMM.get_ylim()[0], self.ax_TbT_HMM.get_ylim()[1]
    #         where=np.array(self.TbT_dataset[trace_ID]['HMM_param'].get('predict')) == 0,
    #         color='lightgreen',
    #         alpha=0.3,
    #         interpolate=False,
    #         transform=self.ax_TbT_HMM.get_xaxis_transform()
    #     )
    
    #     self.fill_red = self.ax_TbT_HMM.fill_between(
    #         time_data, 0, 1,
    #         where=np.array(self.TbT_dataset[trace_ID]['HMM_param'].get('predict')) == 1,
    #         color='lightcoral',
    #         alpha=0.3,
    #         interpolate=False,
    #         transform=self.ax_TbT_HMM.get_xaxis_transform()
    #     )
        
    #     self.fill_greenFRETEff = self.ax_TbT_HMM_FRETEff.fill_between(
    #         time_data, 0, 1,            #self.ax_TbT_HMM.get_ylim()[0], self.ax_TbT_HMM.get_ylim()[1]
    #         where=np.array(self.TbT_dataset[trace_ID]['HMM_param'].get('predict')) == 0,
    #         color='lightgreen',
    #         alpha=0.3,
    #         interpolate=False,
    #         transform=self.ax_TbT_HMM_FRETEff.get_xaxis_transform()
    #     )
    
    #     self.fill_redFRETEff = self.ax_TbT_HMM_FRETEff.fill_between(
    #         time_data, 0, 1,
    #         where=np.array(self.TbT_dataset[trace_ID]['HMM_param'].get('predict')) == 1,
    #         color='lightcoral',
    #         alpha=0.3,
    #         interpolate=False,
    #         transform=self.ax_TbT_HMM_FRETEff.get_xaxis_transform()
    #     )
    
    # else:
    #     if self.fill_green is not None:
    #         self.fill_green.remove()
    #         self.fill_red.remove()
    #         self.fill_greenFRETEff.remove()
    #         self.fill_redFRETEff.remove()
    #         self.fill_green = None
    #         self.fill_red = None
    #         self.fill_greenFRETEff = None
    #         self.fill_redFRETEff = None
    #     self.line_TbT_HMM_TbT.set_data(time_data, scale_factor_HMM * np.array(self.TbT_dataset[trace_ID]['HMM_param'].get('predict')))
    #     self.line_TbT_HMM_FRETEff_predict.set_data(time_data, np.array(self.TbT_dataset[trace_ID]['HMM_param'].get('predict')) / (self.spin_box_HMM_nb_states.value() - 1))
    
    if self.HMM_display_checkbox.isChecked():


        if self.fill_color is not None:
            for j in range(len(self.fill_color)):
                self.fill_color[j].remove()
                self.fill_FRETEff[j].remove()
        self.fill_color = [None for i in range(self.stack_rate_mat.shape[1])]
        self.fill_FRETEff = [None for i in range(self.stack_rate_mat.shape[1])]
        for j in range(len(self.fill_color)):
            self.fill_color[j] = self.ax_TbT_HMM.fill_between(
                time_data, 0, 1,            #self.ax_TbT_HMM.get_ylim()[0], self.ax_TbT_HMM.get_ylim()[1]
                where=np.array(self.TbT_dataset[trace_ID]['HMM_param'].get('predict')) == j,
                color=self.HMM_colors[j],
                alpha=0.5,
                interpolate=False,
                transform=self.ax_TbT_HMM.get_xaxis_transform()
            )
            
            self.fill_FRETEff[j] = self.ax_TbT_HMM_FRETEff.fill_between(
                time_data, 0, 1,            #self.ax_TbT_HMM.get_ylim()[0], self.ax_TbT_HMM.get_ylim()[1]
                where=np.array(self.TbT_dataset[trace_ID]['HMM_param'].get('predict')) == j,
                color=self.HMM_colors[j],
                alpha=0.5,
                interpolate=False,
                transform=self.ax_TbT_HMM_FRETEff.get_xaxis_transform()
            )
        
        self.line_TbT_HMM_TbT.set_data([],[])
        self.line_TbT_HMM_FRETEff_predict.set_data([],[])

    else:
        if self.fill_color is not None:
            for j in range(len(self.fill_color)):
                self.fill_color[j].remove()
                self.fill_FRETEff[j].remove()
        self.fill_color = None
        self.fill_FRETEff = None
        self.line_TbT_HMM_TbT.set_data(time_data, scale_factor_HMM * np.array(self.TbT_dataset[trace_ID]['HMM_param'].get('predict')))
        self.line_TbT_HMM_FRETEff_predict.set_data(time_data, np.array(self.TbT_dataset[trace_ID]['HMM_param'].get('predict')) / (self.spin_box_HMM_nb_states.value() - 1))
    
    
    
    self.canvas_TbT_HMM.draw_idle()
    
    self.toolbar_HMM.update()
    
def move_left_HMM(self):
    if self.isHMMDictLoaded == 0:
        print('Missing HMM analysis')
        return
    current_trace_ID = self.choose_trace_ID_HMM.currentIndex()
    nb_traces = len(self.traces_data.traces)
    if current_trace_ID - 1 >= 0:
        self.choose_trace_ID_HMM.setCurrentIndex(current_trace_ID - 1)
        #plot_HMM_button_function(self)
    display_HMM_plot(self)
        
def move_right_HMM(self):
    if self.isHMMDictLoaded == 0:
        print('Missing HMM analysis')
        return
    current_trace_ID = self.choose_trace_ID_HMM.currentIndex()
    nb_traces = len(self.traces_data.traces)
    if current_trace_ID + 1 < nb_traces:
        self.choose_trace_ID_HMM.setCurrentIndex(current_trace_ID + 1)
    display_HMM_plot(self)
        
def plot_TbT_rates_boxplot(self):
    if self.stack_rate_mat is None:
        print('Need to run TbT analysis first')
        return
    TbT_window = TbT_ratesWindow(self.stack_rate_mat, float(self.widget_frame_delay.text()))
    TbT_window.show()
    self.plot_windows.append(TbT_window)
    
    TbT_HMM_table = Table_rate_Window(np.mean(self.stack_rate_mat, axis = 0)/float(self.widget_frame_delay.text()), 'TbT HMM mean rates matrix')
    TbT_HMM_table.show()
    self.plot_windows.append(TbT_HMM_table)
    
def GUI_run_EnsHMM_HMM(self):
    if self.stack_rate_mat is None:
        print('Need to run TbT analysis first')
        return
    trans_mat_0 = np.mean(self.stack_rate_mat, axis = 0)
    prob_init_0 = np.mean(self.stack_prob_init, axis = 0)
    
    print(trans_mat_0)
    print(prob_init_0)
    
    self.Ens_trans_mat, self.Ens_start_prob = Ens_HMM(self.TbT_dataset, trans_mat_0, prob_init_0)
    
    nb_hidden_states = self.stack_rate_mat.shape[1]
    
    print('Rates estimated from the Ensemble HMM:')
    for i in range(nb_hidden_states):
            for j in range(nb_hidden_states):
                if i != j:
                    label_i = r'$k_{' + str(i) + str(j) + '} = ' + str(np.round(self.Ens_trans_mat[i,j]/float(self.widget_frame_delay.text()), 4)) + ' s^{-1}$'
                    print(label_i)

def plot_EnsHMM_rates_boxplot(self):
    if self.Ens_trans_mat is None:
        print('Need to run EnsHMM analysis first')
        return
    EnsHMM_window = EnsHMM_ratesWindow(self.stack_rate_mat, self.Ens_trans_mat, float(self.widget_frame_delay.text()))
    EnsHMM_window.show()
    self.plot_windows.append(EnsHMM_window)
    
    EnsHMM_table = Table_rate_Window(self.Ens_trans_mat/float(self.widget_frame_delay.text()), 'EnsHMM rates matrix')
    EnsHMM_table.show()
    self.plot_windows.append(EnsHMM_table)


class MatrixModel(QAbstractTableModel):
    def __init__(self, matrix):
        super().__init__()
        self.matrix = matrix

    def rowCount(self, parent=None):
        return self.matrix.shape[0]

    def columnCount(self, parent=None):
        return self.matrix.shape[1]

    def data(self, index, role):
        if role == Qt.DisplayRole:
            return f"{self.matrix[index.row(), index.column()]:.3f}"
        return None

    def headerData(self, section, orientation, role):
        if role != Qt.DisplayRole:
            return None

        if orientation == Qt.Horizontal:
            return f"State {section}"

        if orientation == Qt.Vertical:
            return f"State {section}"

        return None


class VerticalLabel(QWidget):
    def __init__(self, text="", parent=None):
        super().__init__(parent)
        self._text = text

        self.setMinimumWidth(40)
        self.setSizePolicy(
            QSizePolicy.Fixed,
            QSizePolicy.Expanding
        )

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.TextAntialiasing)

        # Move origin to bottom-left and rotate
        painter.translate(0, self.height())
        painter.rotate(-90)

        # After rotation:
        # width  <-> height
        rect = QRectF(
            0,
            0,
            self.height(),
            self.width()
        )

        painter.drawText(
            rect,
            Qt.AlignCenter,
            self._text
        )


class Table_rate_Window(QWidget):
    def __init__(self, matrix, name):
        super().__init__()

        self.setWindowTitle(name)

        #matrix = np.array([[1, 2, 3],[4, 5, 6],[7, 8, 9]]) # np.random.rand(20, 6)

        table = QTableView()
        table.setModel(MatrixModel(matrix))

        table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeToContents
        )

        title = QLabel("Final")
        title.setAlignment(Qt.AlignCenter)

        y_label = VerticalLabel("Initial")

        center_layout = QHBoxLayout()
        center_layout.addWidget(y_label)
        center_layout.addWidget(table)

        layout = QVBoxLayout(self)
        layout.addWidget(title)
        layout.addLayout(center_layout)

        self.resize(800, 500)
