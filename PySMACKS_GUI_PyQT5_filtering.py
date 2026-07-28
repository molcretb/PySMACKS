# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 08:54:48 2026

@author: molcre0000
"""

import sys
from PyQt5.QtWidgets import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.widgets import SpanSelector
from PyQt5.QtGui import QFont
from PyQt5.QtCore import Qt, pyqtSignal
from chromatic_aberrations_correction import *
from utils import *
from traces_extractor import *
from drift_correction import *



def create_filtering_tab(self):
    
    self.isTraceDictLoaded = 0
    
    self.axvspan = {}
    
    #self.axvspan_DO = None
    
    #self.axvspan_AO = None
    
    self.span = None
    
    self.selected_traces = None
    
    tab1 = QWidget()
    layout1 = QVBoxLayout()
    
    ##################### HLayout for displaying individual channels
    
    layout_channel_checkbox = QHBoxLayout()
    
    
    # Checkbox for DD channel
    self.DD_trace_checkbox = QCheckBox("DD")
    self.DD_trace_checkbox.setToolTip('Display Donor-Donor channel')
    self.DD_trace_checkbox.setChecked(True)
    self.DD_trace_checkbox.stateChanged.connect(lambda: display_traces_plot(self))
    layout_channel_checkbox.addWidget(self.DD_trace_checkbox)
    
    # Checkbox for DA channel
    self.DA_trace_checkbox = QCheckBox("DA")
    self.DA_trace_checkbox.setToolTip('Display Donor-Acceptor channel')
    self.DA_trace_checkbox.setChecked(True)
    self.DA_trace_checkbox.stateChanged.connect(lambda: display_traces_plot(self))
    layout_channel_checkbox.addWidget(self.DA_trace_checkbox)
    
    # Checkbox for AA channel
    self.AA_trace_checkbox = QCheckBox("AA")
    self.AA_trace_checkbox.setToolTip('Display Acceptor-Acceptor channel')
    self.AA_trace_checkbox.setChecked(False)
    self.AA_trace_checkbox.stateChanged.connect(lambda: display_traces_plot(self))
    layout_channel_checkbox.addWidget(self.AA_trace_checkbox)
    
    # Checkbox for DD+DA channel
    self.DD_DA_trace_checkbox = QCheckBox("DD+DA")
    self.DD_DA_trace_checkbox.setToolTip('Display DD + DA sum')
    self.DD_DA_trace_checkbox.setChecked(False)
    self.DD_DA_trace_checkbox.stateChanged.connect(lambda: display_traces_plot(self))
    layout_channel_checkbox.addWidget(self.DD_DA_trace_checkbox)
    
    layout1.addLayout(layout_channel_checkbox)
    
    ####################### HLayout for trace ID selection and plot
    
    layout_choose_ID = QHBoxLayout()
    label = QLabel("Trace ID")
    layout_choose_ID.addWidget(label)
    self.choose_trace_ID = QComboBox()
    layout_choose_ID.addWidget(self.choose_trace_ID)
    self.plot_trace = QPushButton("Plot trace")
    self.plot_trace.clicked.connect(lambda: display_traces_plot(self))
    layout_choose_ID.addWidget(self.plot_trace)
    button_left_arrow = QPushButton()
    button_left_arrow.setToolTip('Plot previous trace')
    button_left_arrow.clicked.connect(lambda: move_left(self))
    arrow_left_icon = QApplication.style().standardIcon(QStyle.SP_ArrowLeft)
    button_left_arrow.setIcon(arrow_left_icon)
    layout_choose_ID.addWidget(button_left_arrow)
    button_right_arrow = QPushButton()
    button_right_arrow.setToolTip('Plot next trace')
    button_right_arrow.clicked.connect(lambda: move_right(self))
    arrow_right_icon = QApplication.style().standardIcon(QStyle.SP_ArrowRight)
    button_right_arrow.setIcon(arrow_right_icon)
    layout_choose_ID.addWidget(button_right_arrow)
    
    remove_trace = QPushButton("Clear")
    remove_trace.setToolTip('Clear this trace from the filtered dataset')
    remove_trace.clicked.connect(lambda: clear_saved_traces(self))
    layout_choose_ID.addWidget(remove_trace)
    
    layout1.addLayout(layout_choose_ID)
    
    
    ##################### HLayout for background correction method
    
    layout_choose_back_method = QHBoxLayout()
    label = QLabel("Background correction method")
    layout_choose_back_method.addWidget(label)
    self.back_method = QComboBox()
    self.back_method.addItems(['None', 'Median', 'Total variation', 'Min. of TV'])
    self.back_method.setCurrentText("Min. of TV")
    self.back_method.currentTextChanged.connect(lambda: display_TV_layout(self))
    layout_choose_back_method.addWidget(self.back_method)
    layout1.addLayout(layout_choose_back_method)
    self.TV_param_frame = QFrame()
    layout_select_TV_param = QHBoxLayout()
    self.TV_param_frame.setLayout(layout_select_TV_param)
    label = QLabel("Total variation smoothing")
    layout_select_TV_param.addWidget(label)
    
    self.spin_box_TV_param = QSpinBox()
    self.spin_box_TV_param.setToolTip('Select the smoothing value for the Total Variation denoising of background')
    self.spin_box_TV_param.setMinimum(1)  # Set minimum value
    self.spin_box_TV_param.setMaximum(1000)   # Set maximum value
    self.spin_box_TV_param.setValue(5)
    layout_select_TV_param.addWidget(self.spin_box_TV_param)
    
    layout1.addWidget(self.TV_param_frame)
    
    
    layout_choose_trace_label = QHBoxLayout()
    label = QLabel("Label class")
    layout_choose_trace_label.addWidget(label)
    self.label_choose = QComboBox()
    self.label_choose.setToolTip('Select the class for the manual labelling: FRET, Donor-Only, Acceptor-Only')
    #self.label_choose.addItems(['FRET', 'DO', 'AO'])
    self.label_choose.addItems(self.list_label_classes)
    self.label_choose.setCurrentText(self.list_label_classes[0])
    layout_choose_trace_label.addWidget(self.label_choose)
    layout1.addLayout(layout_choose_trace_label)
    
    
    ##################### Canva for the trace plot
    
    self.figure_filtering = plt.figure()
    self.canvas_filtering = FigureCanvas(self.figure_filtering)
    self.canvas_filtering.setFocus()
    toolbar = NavigationToolbar(self.canvas_filtering, self)
    layout1.addWidget(toolbar)
    
    layout1.addWidget(self.canvas_filtering)
    
    tab1.setLayout(layout1)
    self.tabs.addTab(tab1, "Filtering")


def clear_saved_traces(self):
    trace_ID = self.choose_trace_ID.currentText()
    change = 0
    #del self.selected_traces[trace_ID]
    # if 'indmin' in list(self.traces_data.traces[int(trace_ID)].metadata.keys()):
    #     del self.traces_data.traces[int(trace_ID)].metadata['indmin']
    #     del self.traces_data.traces[int(trace_ID)].metadata['indmax']
    #     change = 1
    # if 'indmin_DO' in list(self.traces_data.traces[int(trace_ID)].metadata.keys()):
    #     del self.traces_data.traces[int(trace_ID)].metadata['indmin_DO']
    #     del self.traces_data.traces[int(trace_ID)].metadata['indmax_DO']
    #     change = 1
    # if 'indmin_AO' in list(self.traces_data.traces[int(trace_ID)].metadata.keys()):
    #     del self.traces_data.traces[int(trace_ID)].metadata['indmin_AO']
    #     del self.traces_data.traces[int(trace_ID)].metadata['indmax_AO']
    #     change = 1
    # if change == 1:
    #     display_traces_plot(self)
    # else:
    #     return
    
    for label_i in self.list_label_classes:
        if label_i in list(self.traces_data.traces[int(trace_ID)].metadata.keys()):
            del self.traces_data.traces[int(trace_ID)].metadata[label_i]
            #del self.traces_data.traces[int(trace_ID)].metadata['FRET']['indmax']
            change = 1
    # if 'DO' in list(self.traces_data.traces[int(trace_ID)].metadata.keys()):
    #     del self.traces_data.traces[int(trace_ID)].metadata['DO']
    #     #del self.traces_data.traces[int(trace_ID)].metadata['DO']['indmax']
    #     change = 1
    # if 'AO' in list(self.traces_data.traces[int(trace_ID)].metadata.keys()):
    #     del self.traces_data.traces[int(trace_ID)].metadata['AO']
    #     #del self.traces_data.traces[int(trace_ID)].metadata['AO']['indmax']
    #     change = 1
    if change == 1:
        display_traces_plot(self)
    else:
        return

def move_left(self):
    if self.isTraceDictLoaded == 0:
        print('No loaded dataset')
        return
    current_trace_ID = self.choose_trace_ID.currentIndex()
    nb_traces = len(self.trace_IDs)
    if current_trace_ID - 1 >= 0:
        self.choose_trace_ID.setCurrentIndex(current_trace_ID - 1)
        display_traces_plot(self)
    else:
        display_traces_plot(self)
        
def move_right(self):
    if self.isTraceDictLoaded == 0:
        print('No loaded dataset')
        return
    current_trace_ID = self.choose_trace_ID.currentIndex()
    nb_traces = len(self.trace_IDs)
    if current_trace_ID + 1 < nb_traces:
        self.choose_trace_ID.setCurrentIndex(current_trace_ID + 1)
        display_traces_plot(self)
    else:
        display_traces_plot(self)
        
def keyPressEvent(self, event):
    if event.key() == Qt.Key_Right:
        current_trace_ID = self.choose_trace_ID.currentIndex()
        nb_traces = len(self.trace_IDs)
        if current_trace_ID + 1 < nb_traces:
            self.choose_trace_ID.setCurrentIndex(current_trace_ID + 1)
            display_traces_plot(self)
        else:
            display_traces_plot(self)
    elif event.key() == Qt.Key_Left:
        current_trace_ID = self.choose_trace_ID.currentIndex()
        nb_traces = len(self.trace_IDs)
        if current_trace_ID - 1 > 0:
            self.choose_trace_ID.setCurrentIndex(current_trace_ID - 1)
            display_traces_plot(self)
        else:
            display_traces_plot(self)
    
def display_traces_plot(self):
    if self.isTraceDictLoaded == 0:
        print('No loaded dataset')
        return
    trace_ID = self.choose_trace_ID.currentText()
    
    method_back_corr = self.back_method.currentText()
    
    #time_data = np.linspace(0, len(np.array(self.traces_data[trace_ID].get('Intensity_DD')))-1, len(np.array(self.traces_data[trace_ID].get('Intensity_DD'))))
    
    time_data = np.linspace(0, len(np.array(self.traces_data.traces[int(trace_ID)].channels[0].data))-1, len(np.array(self.traces_data.traces[int(trace_ID)].channels[0].data)))
    
    #self.figure_filtering.clear()
    
    #self.ax_trace_plot = self.figure_filtering.add_subplot(111)
    
    if self.DD_trace_checkbox.isChecked():
        
        DD_trace = np.array(self.traces_data.traces[int(trace_ID)].channels[0].data)
        
        match method_back_corr:
            case 'None':
                back_DD_trace = np.zeros(len(DD_trace))
            case 'Median':
                if self.traces_data.metadata['ALEX'] == 'yes':
                    back_DD_trace = np.array(self.traces_data.traces[int(trace_ID)].channels[3].data)
                else:
                    back_DD_trace = np.array(self.traces_data.traces[int(trace_ID)].channels[2].data)
            case 'Total variation':
                if self.traces_data.metadata['ALEX'] == 'yes':
                    raw_back_DD = np.array(self.traces_data.traces[int(trace_ID)].channels[3].data)
                else:
                    raw_back_DD = np.array(self.traces_data.traces[int(trace_ID)].channels[2].data)
                max_back_DD = np.max(raw_back_DD)
                norm_back_trace_DD = np.array(raw_back_DD) / max_back_DD
                TV_back_DD = tvd_2013(norm_back_trace_DD.astype('float'), self.spin_box_TV_param.value())
                back_DD_trace = TV_back_DD * max_back_DD
                
            case 'Min. of TV':
                if self.traces_data.metadata['ALEX'] == 'yes':
                    raw_back_DD = np.array(self.traces_data.traces[int(trace_ID)].channels[3].data)
                else:
                    raw_back_DD = np.array(self.traces_data.traces[int(trace_ID)].channels[2].data)
                max_back_DD = np.max(raw_back_DD)
                norm_back_trace_DD = np.array(raw_back_DD) / max_back_DD
                #print('TV lambda:' + str(self.spin_box_TV_param.value()))
                TV_back_DD = tvd_2013(norm_back_trace_DD.astype('float'), self.spin_box_TV_param.value())
                back_DD_trace = np.min(TV_back_DD) * max_back_DD
        #self.ax_trace_plot.plot(np.array(self.traces_data[trace_ID].get('Intensity_DD')), 'orange')
        #self.line_DD.set_data(time_data, np.array(self.traces_data[trace_ID].get('Intensity_DD')))
        self.line_DD.set_data(time_data, DD_trace - back_DD_trace)
    else:
        self.line_DD.set_data([],[])
    if self.DA_trace_checkbox.isChecked():
        #self.ax_trace_plot.plot(np.array(self.traces_data[trace_ID].get('Intensity_DA')),'red')
        #self.line_DA.set_data(time_data, np.array(self.traces_data[trace_ID].get('Intensity_DA')))
        #self.line_DA.set_data(time_data, np.array(self.traces_data.traces[int(trace_ID)].channels[1].data))
        DA_trace = np.array(self.traces_data.traces[int(trace_ID)].channels[1].data)
        
        match method_back_corr:
            case 'None':
                back_DA_trace = np.zeros(len(DA_trace))
            case 'Median':
                if self.traces_data.metadata['ALEX'] == 'yes':
                    back_DA_trace = np.array(self.traces_data.traces[int(trace_ID)].channels[4].data)
                else:
                    back_DA_trace = np.array(self.traces_data.traces[int(trace_ID)].channels[3].data)
            case 'Total variation':
                if self.traces_data.metadata['ALEX'] == 'yes':
                    raw_back_DA = np.array(self.traces_data.traces[int(trace_ID)].channels[4].data)
                else:
                    raw_back_DA = np.array(self.traces_data.traces[int(trace_ID)].channels[3].data)
                max_back_DA = np.max(raw_back_DA)
                norm_back_trace_DA = np.array(raw_back_DA) / max_back_DA
                TV_back_DA = tvd_2013(norm_back_trace_DA.astype('float'), self.spin_box_TV_param.value())
                back_DA_trace = TV_back_DA * max_back_DA
                
            case 'Min. of TV':
                if self.traces_data.metadata['ALEX'] == 'yes':
                    raw_back_DA = np.array(self.traces_data.traces[int(trace_ID)].channels[4].data)
                else:
                    raw_back_DA = np.array(self.traces_data.traces[int(trace_ID)].channels[3].data)
                max_back_DA = np.max(raw_back_DA)
                norm_back_trace_DA = np.array(raw_back_DA) / max_back_DA
                TV_back_DA = tvd_2013(norm_back_trace_DA.astype('float'), self.spin_box_TV_param.value())
                back_DA_trace = np.min(TV_back_DA) * max_back_DA
        self.line_DA.set_data(time_data, DA_trace - back_DA_trace)
    else:
        self.line_DA.set_data([],[])
    if self.AA_trace_checkbox.isChecked():
        #self.ax_trace_plot.plot(np.array(self.traces_data[trace_ID].get('Intensity_AA')),'gray')
        #self.line_AA.set_data(time_data, np.array(self.traces_data[trace_ID].get('Intensity_AA')))
        #self.line_AA.set_data(time_data, np.array(self.traces_data.traces[int(trace_ID)].channels[2].data))
        
        if self.traces_data.metadata['ALEX'] == 'yes':
            AA_trace = np.array(self.traces_data.traces[int(trace_ID)].channels[2].data)
            
            match method_back_corr:
                case 'None':
                    back_AA_trace = np.zeros(len(AA_trace))
                case 'Median':
                    back_AA_trace = np.array(self.traces_data.traces[int(trace_ID)].channels[5].data)
                case 'Total variation':
                    raw_back_AA = np.array(self.traces_data.traces[int(trace_ID)].channels[5].data)
                    max_back_AA = np.max(raw_back_AA)
                    norm_back_trace_AA = np.array(raw_back_AA) / max_back_AA
                    TV_back_AA = tvd_2013(norm_back_trace_AA.astype('float'), self.spin_box_TV_param.value())
                    back_AA_trace = TV_back_AA * max_back_AA
                    
                case 'Min. of TV':
                    raw_back_AA = np.array(self.traces_data.traces[int(trace_ID)].channels[5].data)
                    max_back_AA = np.max(raw_back_AA)
                    norm_back_trace_AA = np.array(raw_back_AA) / max_back_AA
                    TV_back_AA = tvd_2013(norm_back_trace_AA.astype('float'), self.spin_box_TV_param.value())
                    back_AA_trace = np.min(TV_back_AA) * max_back_AA
            self.line_AA.set_data(time_data, AA_trace - back_AA_trace)
        else:
            self.line_AA.set_data([],[])
            print('No ALEX data')
    else:
        self.line_AA.set_data([],[])
    if self.DD_DA_trace_checkbox.isChecked():
        #self.ax_trace_plot.plot(np.array(self.traces_data[trace_ID].get('Intensity_DD')) + np.array(self.traces_data[trace_ID].get('Intensity_DA')),'blue')
        #self.line_DD_DA.set_data(time_data, np.array(self.traces_data[trace_ID].get('Intensity_DD')) + np.array(self.traces_data[trace_ID].get('Intensity_DA')))
        #self.line_DD_DA.set_data(time_data, np.array(self.traces_data.traces[int(trace_ID)].channels[0].data) + np.array(self.traces_data.traces[int(trace_ID)].channels[1].data))
        
        DD_trace = np.array(self.traces_data.traces[int(trace_ID)].channels[0].data)
        DA_trace = np.array(self.traces_data.traces[int(trace_ID)].channels[1].data)
        
        match method_back_corr:
            case 'None':
                back_DD_trace = np.zeros(len(DD_trace))
                back_DA_trace = np.zeros(len(DA_trace))
            case 'Median':
                if self.traces_data.metadata['ALEX'] == 'yes':
                    back_DD_trace = np.array(self.traces_data.traces[int(trace_ID)].channels[3].data)
                    back_DA_trace = np.array(self.traces_data.traces[int(trace_ID)].channels[4].data)
                else:
                    back_DD_trace = np.array(self.traces_data.traces[int(trace_ID)].channels[2].data)
                    back_DA_trace = np.array(self.traces_data.traces[int(trace_ID)].channels[3].data)
            case 'Total variation':
                if self.traces_data.metadata['ALEX'] == 'yes':
                    raw_back_DD = np.array(self.traces_data.traces[int(trace_ID)].channels[3].data)
                    raw_back_DA = np.array(self.traces_data.traces[int(trace_ID)].channels[4].data)
                else:
                    raw_back_DD = np.array(self.traces_data.traces[int(trace_ID)].channels[2].data)
                    raw_back_DA = np.array(self.traces_data.traces[int(trace_ID)].channels[3].data)
                max_back_DD = np.max(raw_back_DD)
                max_back_DA = np.max(raw_back_DA)
                norm_back_trace_DD = np.array(raw_back_DD) / max_back_DD
                norm_back_trace_DA = np.array(raw_back_DA) / max_back_DA
                TV_back_DD = tvd_2013(norm_back_trace_DD.astype('float'), self.spin_box_TV_param.value())
                TV_back_DA = tvd_2013(norm_back_trace_DA.astype('float'), self.spin_box_TV_param.value())
                back_DD_trace = TV_back_DD * max_back_DD
                back_DA_trace = TV_back_DA * max_back_DA
                
            case 'Min. of TV':
                if self.traces_data.metadata['ALEX'] == 'yes':
                    raw_back_DD = np.array(self.traces_data.traces[int(trace_ID)].channels[3].data)
                    raw_back_DA = np.array(self.traces_data.traces[int(trace_ID)].channels[4].data)
                else:
                    raw_back_DD = np.array(self.traces_data.traces[int(trace_ID)].channels[2].data)
                    raw_back_DA = np.array(self.traces_data.traces[int(trace_ID)].channels[3].data)
                max_back_DD = np.max(raw_back_DD)
                max_back_DA = np.max(raw_back_DA)
                norm_back_trace_DD = np.array(raw_back_DD) / max_back_DD
                norm_back_trace_DA = np.array(raw_back_DA) / max_back_DA
                TV_back_DD = tvd_2013(norm_back_trace_DD.astype('float'), self.spin_box_TV_param.value())
                TV_back_DA = tvd_2013(norm_back_trace_DA.astype('float'), self.spin_box_TV_param.value())
                back_DD_trace = np.min(TV_back_DD) * max_back_DD
                back_DA_trace = np.min(TV_back_DA) * max_back_DA
        self.line_DD_DA.set_data(time_data, DD_trace + DA_trace - back_DD_trace - back_DA_trace)
        
        
    else:
        self.line_DD_DA.set_data([],[])
        
    self.ax_trace_plot.relim()
    self.ax_trace_plot.autoscale_view(scalex=False, scaley=True)
    self.ax_trace_plot.set_xlim(int(-len(time_data)*0.05), int(len(time_data)*1.05))
    
    for i in list(self.axvspan.keys()):
        if self.axvspan[i] is not None:
            self.axvspan[i].remove()
            self.axvspan[i] = None
    
    # if self.axvspan is not None:
    #     self.axvspan.remove()
    #     self.axvspan = None
        
    # if self.axvspan_DO is not None:
    #     self.axvspan_DO.remove()
    #     self.axvspan_DO = None
        
    # if self.axvspan_AO is not None:
    #     self.axvspan_AO.remove()
    #     self.axvspan_AO = None
        
    
    if self.span is not None:
        self.span.set_visible(False)
        #self.span.disconnect_events()
        #self.span = None
    
    # if trace_ID in list(self.selected_traces.keys()):
    #     self.axvspan = self.ax_trace_plot.axvspan(self.selected_traces[trace_ID]['indmin'], self.selected_traces[trace_ID]['indmax'], alpha=0.5, facecolor="tab:green")
    count_color = 0
    for label_i in self.list_label_classes:
    
        if label_i in list(self.traces_data.traces[int(trace_ID)].metadata.keys()):
            self.axvspan[label_i] = self.ax_trace_plot.axvspan(self.traces_data.traces[int(trace_ID)].metadata[label_i]['indmin'], self.traces_data.traces[int(trace_ID)].metadata[label_i]['indmax'], alpha=0.5, facecolor=self.HMM_colors[count_color])
        count_color += 1
    # if 'DO' in list(self.traces_data.traces[int(trace_ID)].metadata.keys()):
    #     self.axvspan_DO = self.ax_trace_plot.axvspan(self.traces_data.traces[int(trace_ID)].metadata['DO']['indmin'], self.traces_data.traces[int(trace_ID)].metadata['DO']['indmax'], alpha=0.5, facecolor="lightcoral")
    
    # if 'AO' in list(self.traces_data.traces[int(trace_ID)].metadata.keys()):
    #     self.axvspan_AO = self.ax_trace_plot.axvspan(self.traces_data.traces[int(trace_ID)].metadata['AO']['indmin'], self.traces_data.traces[int(trace_ID)].metadata['AO']['indmax'], alpha=0.5, facecolor="lightblue")
    
    
    
    self.canvas_filtering.draw_idle()
    
    # self.span = SpanSelector(
    #     self.ax_trace_plot,
    #     self.drag_select,
    #     "horizontal",
    #     useblit=True,
    #     props=dict(alpha=0.5, facecolor="tab:blue"),
    #     interactive=True,
    #     drag_from_anywhere=True
    # )
    
    
def display_TV_layout(self):
    back_method = self.back_method.currentText()
    if (back_method == 'Total variation') or (back_method == 'Min. of TV'):  # Checked
        self.TV_param_frame.show()
    else:  # Unchecked
        self.TV_param_frame.hide()
    
    
    
    
    