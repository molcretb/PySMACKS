# -*- coding: utf-8 -*-
"""
Created on Wed May 27 15:24:47 2026

@author: molcre0000
"""

import sys
import os
from PyQt5.QtWidgets import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PyQt5.QtGui import QFont, QKeySequence
from PyQt5.QtCore import Qt, pyqtSignal
from matplotlib.widgets import SpanSelector
from openfret import read_data, write_data
from PySMACKS.components.chromatic_aberrations_correction import *
from PySMACKS.components.utils import *
from PySMACKS.components.traces_extractor import *
from PySMACKS.components.drift_correction import *
from PySMACKS.GUI.GUI_registration import *
from PySMACKS.GUI.GUI_extraction import *
from PySMACKS.GUI.GUI_filtering import *
from PySMACKS.GUI.GUI_SEhist import *
from PySMACKS.GUI.GUI_kinetic import *


class Save_Window(QMainWindow):
    def __init__(self, registration_matrix, traces_data, traces_data_DO, traces_data_AO):
        super().__init__()
        self.setWindowTitle("Save results")
        self.resize(700, 700)
        
        
        central_widget = QWidget()
        layout1 = QVBoxLayout(central_widget)
        
        label = QLabel("Select the data you want to save (as separate files):")
        layout1.addWidget(label)
        self.save_regist_matrix = QCheckBox("Registration")
        self.save_regist_matrix.setToolTip('Save registration matrix')
        self.save_regist_matrix.setChecked(False)
        layout1.addWidget(self.save_regist_matrix)
        self.save_extract_traces = QCheckBox("Extraction")
        self.save_extract_traces.setToolTip('Save extracted traces')
        self.save_extract_traces.setChecked(False)
        layout1.addWidget(self.save_extract_traces)
        self.save_filter_traces = QCheckBox("Donor-Only")
        self.save_filter_traces.setToolTip('Save DO traces')
        self.save_filter_traces.setChecked(False)
        layout1.addWidget(self.save_filter_traces)
        self.save_FRET_corr_factors = QCheckBox("Acceptor-Only")
        self.save_FRET_corr_factors.setToolTip('Save AO traces')
        self.save_FRET_corr_factors.setChecked(False)
        layout1.addWidget(self.save_FRET_corr_factors)
        
        button_save = QPushButton("Save data")
        button_save.setToolTip('Save results from the analysis pipeline')
        button_save.clicked.connect(lambda: self.save_function(registration_matrix, traces_data, traces_data_DO, traces_data_AO))
        layout1.addWidget(button_save)

        self.setCentralWidget(central_widget)
        
        
        
    def save_function(self, registration_matrix, traces_data, traces_data_DO, traces_data_AO):

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self, 
            "Save File", os.getcwd() + "\\BASENAME_of_your_data", "All Files (*)", options = options)
        if fileName:
            if self.save_regist_matrix.isChecked():
                if registration_matrix is None:
                    print('No registration matrix')
                else:
                    self.save_registration(registration_matrix, fileName)
            if self.save_extract_traces.isChecked():
                if traces_data is None:
                    print('No traces dataset')
                else:
                    self.save_data_JSON(traces_data, fileName, 'traces_dataset')
            if self.save_filter_traces.isChecked():
                if traces_data_DO is None:
                    print('No DO traces dataset')
                else:
                    self.save_data_JSON(traces_data_DO, fileName, 'DO_dataset')
            if self.save_FRET_corr_factors.isChecked():
                if traces_data_AO is None:
                    print('No AO traces dataset')
                else:
                    self.save_data_JSON(traces_data_AO, fileName, 'AO_dataset')
        
        
    def save_data_JSON(self, traces_dict, file_save, type_data):
        
        
        name_complete = file_save + '_' + type_data + '.json'
        
        write_data(traces_dict, name_complete, compress=True)
        
        print(type_data + ' Traces dataset saved as compressed openFRET JSON file to path: ' + name_complete)
        
    def save_registration(self, matrix_align, file_save):
        name_complete = file_save + '_registration_matrix.npy'
        
        np.save(name_complete, matrix_align)
        
        print('Registration matrix saved as numpy file to path: ' + name_complete)
        
class Help_shortcuts(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Shortcuts")
        self.setGeometry(100, 100, 600, 600)

        # Create a table with 2 columns and 4 rows
        table = QTableWidget(4, 2)

        # Set column headers
        table.setHorizontalHeaderLabels(["Keys", "Functionality"])

        # Fill the table with text
        data = [
            ("1-9", "Select a class for trace labelling"),
            ("Left-Right arrows", "Select previous/next trace"),
            ("Return", "Plot the selected trace"),
            ("Backspace", "Clear the labelled selection"),
        ]

        for row, (name, occupation) in enumerate(data):
            table.setItem(row, 0, QTableWidgetItem(name))
            table.setItem(row, 1, QTableWidgetItem(occupation))
        table.verticalHeader().setVisible(False)
        table.resizeColumnsToContents()
        self.setCentralWidget(table)
        
class metadata_Viewer(QMainWindow):
    def __init__(self, dataset):
        super().__init__()

        self.setWindowTitle("Metadata Viewer")
        self.resize(500, 300)
        
        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)
        
        list_metadata_fields = [dataset.metadata, dataset.sample_details, dataset.instrument_details]
        
        list_field_name = ['dataset', 'sample', 'instrument']
        
        count_tab = 0
        for metadata_item in list_metadata_fields:
            tab = QWidget()
            layout = QVBoxLayout()

            table = QTableWidget()
            table.setColumnCount(2)
            table.setHorizontalHeaderLabels(["Key", "Value"])
            if list_field_name[count_tab] == 'dataset':
                table.setRowCount(len(metadata_item)+6)
            else:
                table.setRowCount(len(metadata_item))
            count_row = 0
            for row, (key, value) in enumerate(metadata_item.items()):
                table.setItem(row, 0, QTableWidgetItem(str(key)))
                table.setItem(row, 1, QTableWidgetItem(str(value)))
                count_row  += 1
            if list_field_name[count_tab] == 'dataset':
                metadata_item = dataset
                table.setItem(count_row, 0, QTableWidgetItem('title'))
                table.setItem(count_row, 1, QTableWidgetItem(str(metadata_item.title)))
                table.setItem(count_row+1, 0, QTableWidgetItem('description'))
                table.setItem(count_row+1, 1, QTableWidgetItem(str(metadata_item.description)))
                table.setItem(count_row+2, 0, QTableWidgetItem('experiment_type'))
                table.setItem(count_row+2, 1, QTableWidgetItem(str(metadata_item.experiment_type)))
                table.setItem(count_row+3, 0, QTableWidgetItem('authors'))
                table.setItem(count_row+3, 1, QTableWidgetItem(str(metadata_item.authors)))
                table.setItem(count_row+4, 0, QTableWidgetItem('institution'))
                table.setItem(count_row+4, 1, QTableWidgetItem(str(metadata_item.institution)))
                table.setItem(count_row+5, 0, QTableWidgetItem('date'))
                table.setItem(count_row+5, 1, QTableWidgetItem(str(metadata_item.date.isoformat())))
            table.resizeColumnsToContents()
            table.horizontalHeader().setStretchLastSection(True)
            table.verticalHeader().setVisible(False)

            layout.addWidget(table)

            
            tab.setLayout(layout)
            
            self.tabs.addTab(tab, list_field_name[count_tab])
            count_tab += 1
        
        tab_traces = QWidget()
        self.layout_traces = QVBoxLayout()
        
        layout_choose_ID = QHBoxLayout()
        label = QLabel("Trace ID")
        layout_choose_ID.addWidget(label)
        self.choose_trace_ID = QComboBox()
        trace_IDs = [str(i) for i in range(len(dataset.traces))]
        self.choose_trace_ID.addItems(trace_IDs)
        layout_choose_ID.addWidget(self.choose_trace_ID)
        self.plot_metadata = QPushButton("Plot trace")
        self.plot_metadata.clicked.connect(lambda: display_traces_metadata(self, dataset.traces[int(self.choose_trace_ID.currentText())]))
        layout_choose_ID.addWidget(self.plot_metadata)
        
        self.layout_traces.addLayout(layout_choose_ID)
        
        self.table = QTableWidget()
        self.table.setColumnCount(2)
        self.table.setHorizontalHeaderLabels(["Key", "Value"])
        self.table.setRowCount(1)
        self.layout_traces.addWidget(self.table)
        
        tab_traces.setLayout(self.layout_traces)
        
        self.tabs.addTab(tab_traces, 'Traces')
        
        def display_traces_metadata(self, dataset_trace_i):
            
            self.table.clearContents()
            self.table.setRowCount(4 + len(dataset_trace_i.metadata.items()))
            self.table.setItem(0, 0, QTableWidgetItem('exposure_time'))
            self.table.setItem(0, 1, QTableWidgetItem(str(dataset_trace_i.channels[0].exposure_time)))
            self.table.setItem(1, 0, QTableWidgetItem('excitation_wavelength'))
            self.table.setItem(1, 1, QTableWidgetItem(str(dataset_trace_i.channels[0].excitation_wavelength)))
            self.table.setItem(2, 0, QTableWidgetItem('emission_wavelength'))
            self.table.setItem(2, 1, QTableWidgetItem(str(dataset_trace_i.channels[0].emission_wavelength)))
            self.table.setItem(3, 0, QTableWidgetItem('channel_metadata'))
            self.table.setItem(3, 1, QTableWidgetItem(str(dataset_trace_i.channels[0].metadata)))

            for row, (key, value) in enumerate(dataset_trace_i.metadata.items()):
                self.table.setItem(row+4, 0, QTableWidgetItem(str(key)))
                self.table.setItem(row+4, 1, QTableWidgetItem(str(value)))
            
            self.table.resizeColumnsToContents()
            self.table.horizontalHeader().setStretchLastSection(True)
            self.table.verticalHeader().setVisible(False)
            
            self.layout_traces.addWidget(self.table)
        
        
class addLabelWindow(QWidget):
    addRequested = pyqtSignal(str)
    removeRequested = pyqtSignal(str)
    def __init__(self, list_labels):
        super().__init__()

        self.setWindowTitle("Labelling classes")

        # List widget
        self.list_widget = QListWidget()
        self.list_widget.addItems(list_labels)

        # Input field
        self.entry = QLineEdit()
        self.entry.setPlaceholderText("Enter a new item...")

        # Buttons
        self.add_button = QPushButton("Add")
        self.remove_button = QPushButton("Remove")

        # Top row: input + add
        input_layout = QHBoxLayout()
        input_layout.addWidget(self.entry)
        input_layout.addWidget(self.add_button)

        # Bottom row: remove button
        button_layout = QHBoxLayout()
        button_layout.addStretch()
        button_layout.addWidget(self.remove_button)
        
        # Bottom row: update button

        # Main layout
        layout = QVBoxLayout(self)
        layout.addWidget(self.list_widget)
        layout.addLayout(input_layout)
        layout.addLayout(button_layout)

        # Signals
        self.add_button.clicked.connect(self.add_item)
        self.entry.returnPressed.connect(self.add_item)
        self.remove_button.clicked.connect(self.remove_item)

    def add_item(self):
        text = self.entry.text().strip()
        if text:
            self.list_widget.addItem(text)
            self.addRequested.emit(text)
            self.entry.clear()
            self.entry.setFocus()

    def remove_item(self):
        item = self.list_widget.currentItem()
        if item is not None:
            row = self.list_widget.row(item)
            self.list_widget.takeItem(row)
            self.removeRequested.emit(item.text())
            
            
class export_traces(QMainWindow):
    def __init__(self, traces_data):
        super().__init__()
        self.setWindowTitle("Export results")
        self.resize(700, 700)
        
        
        ### section to collect all the classes contained in the dataset
        
        nb_traces = len(traces_data.traces)
        
        list_keys = []
        
        for i in range(nb_traces):
            list_keys = list_keys + list(traces_data.traces[i].metadata.keys())
        seen = set()
        list_keys[:] = [item for item in list_keys if item not in seen and not seen.add(item)]
        
        list_keys.append('Full')
        
        not_classes = ['molecule_id', 'x_coord', 'y_coord', 'UUID_v7', 'background_correction', 'TV_lambda']
        
        for item in not_classes:
            if item in list_keys:
                list_keys.remove(item)
        
        ###
        
        list_export_format = ['CSV']
        
        central_widget = QWidget()
        layout1 = QVBoxLayout(central_widget)
        
        layout_export_class = QHBoxLayout()
        
        label = QLabel("Select the class of data you want to save:")
        layout_export_class.addWidget(label)
        
        self.export_class = QComboBox()
        self.export_class.addItems(list_keys)
        layout_export_class.addWidget(self.export_class)
        
        layout1.addLayout(layout_export_class)
        
        layout_export_path = QHBoxLayout()
        
        button_export_path = QPushButton("Select export folder")
        button_export_path.clicked.connect(self.select_export_folder)
        layout_export_path.addWidget(button_export_path)
        self.path_export_folder = QLineEdit()
        layout_export_path.addWidget(self.path_export_folder)
        
        layout1.addLayout(layout_export_path)
        
        layout_export_format = QHBoxLayout()
        
        label = QLabel("Select the export data format:")
        layout_export_format.addWidget(label)
        
        self.export_format = QComboBox()
        self.export_format.addItems(list_export_format)
        layout_export_format.addWidget(self.export_format)
        
        layout1.addLayout(layout_export_format)
        
        button_export_data = QPushButton("Export data")
        button_export_data.clicked.connect(lambda: self.export_data(traces_data))
        layout1.addWidget(button_export_data)
        
        self.setCentralWidget(central_widget)
        
    def select_export_folder(self):
        # You can set options like QFileDialog.DontUseNativeDialog if needed
        folder_name = QFileDialog.getExistingDirectory(self, "Select the folder where to export the traces dataset")
        if folder_name:
            self.path_export_folder.setText(folder_name)
            
    def export_data(self, traces_data):
        if not self.path_export_folder.text().strip():
            print('No export folder selected, please select one.')
        else:
            export_format = self.export_format.currentText()
            export_class = self.export_class.currentText()
            export_path = self.path_export_folder.text()
            match export_format:
                case 'CSV':
                    export_OpenFRET_to_CSV_traces(traces_data, export_class, export_path)


class FileSelectorApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("PySMACKS GUI[*]")
        self.setGeometry(100, 100, 1600, 400)
        self.plot_windows = []
        
        self.auto_save_calib = {}
        self.fileName = None
        
        self.second_window = None
        self.count_windows = 0
        
        self.traces_data = None
        self.traces_DO_data = None
        self.traces_AO_data = None
        self.shortcuts = []
        
        self.list_label_classes = ["FRET", "DO", "AO"]
        

        # Initialize the central widget and layout
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)

        self.layout = QVBoxLayout()
        self.central_widget.setLayout(self.layout)

        # Label to display selected file path
        self.file_path_label = QLabel("No file selected")
        self.layout.addWidget(self.file_path_label)

        # Create menu bar and add 'File' menu
        menu_bar = self.menuBar()
        file_menu = menu_bar.addMenu("File")
        QC_menu = menu_bar.addMenu("QC")
        help_menu = menu_bar.addMenu("Help")
        label_menu = menu_bar.addMenu("Labels")
        
        # Create the tab widget
        self.tabs = QTabWidget()
        #self.tabs.setFocusPolicy(Qt.NoFocus)
        self.setCentralWidget(self.tabs)
        
        create_registration_tab(self)
        create_extraction_tab(self)
        create_filtering_tab(self)
        create_SEhist_tab(self)
        create_kinetic_tab(self)
        
        # shortcuts for selecting label classes
        for j in range(len(self.list_label_classes)):
            shortcut = QShortcut(QKeySequence(str(j+1)), self, activated=lambda j=j: self.label_choose.setCurrentIndex(j))
            self.shortcuts.append(shortcut)
        QShortcut(QKeySequence(Qt.Key_Left), self, activated=lambda: move_left(self))
        QShortcut(QKeySequence(Qt.Key_Right), self, activated=lambda: move_right(self))
        QShortcut(QKeySequence("Return"), self, activated=lambda: display_traces_plot(self))
        QShortcut(QKeySequence("Backspace"), self, activated=lambda: clear_saved_traces(self))
        
        font_tab = QFont("Consolas", 12, QFont.Bold)
        self.tabs.tabBar().setFont(font_tab)

        # Create 'Open' action
        open_action = QAction("Open", self)
        open_action.triggered.connect(lambda: open_file_dialog_acceptor_regis(self))
        file_menu.addAction(open_action)
        
        # Create 'Load traces dataset' action
        load_data_action = QAction("Load traces dataset", self)
        load_data_action.triggered.connect(self.load_traces_dataset)
        file_menu.addAction(load_data_action)
        
        # Create 'merge dataset' action
        merge_data_action = QAction("Merge traces dataset", self)
        merge_data_action.triggered.connect(self.merge_traces_dataset)
        file_menu.addAction(merge_data_action)
        
        # Create 'Save Registration' action
        saveAction = file_menu.addAction('Save')
        saveAction.triggered.connect(self.save)
        
        
        saveAsAction = file_menu.addAction('Save as...')
        saveAsAction.triggered.connect(self.saveAs)
        
        # Create 'export dataset' action
        export_data_action = QAction("Export traces dataset", self)
        export_data_action.triggered.connect(self.export_traces_dataset)
        file_menu.addAction(export_data_action)
        
        # Menu to access traces metadata from the viewer
        
        metadataView_Action = QC_menu.addAction('Metadata viewer')
        metadataView_Action.triggered.connect(self.open_metadata_viewer)
        
        # Menu button to generate the registration QC plot
        registQCAction = QC_menu.addAction('Registration QC')
        registQCAction.triggered.connect(lambda: regist_QC_plot(self))
        
        disp_help_window = help_menu.addAction('Shortcuts')
        disp_help_window.triggered.connect(self.display_help_window)
        
        edit_labels_window = label_menu.addAction('Edit')
        edit_labels_window.triggered.connect(self.display_label_window)
        
    def open_metadata_viewer(self):
        
        if self.isTraceDictLoaded == 0:
            print('No loaded dataset')
            return
        metadata_Window = metadata_Viewer(self.traces_data)
        metadata_Window.show()
        self.plot_windows.append(metadata_Window)
        
    def display_help_window(self):
        help_Window = Help_shortcuts()
        help_Window.show()
        self.plot_windows.append(help_Window)
        
    def export_traces_dataset(self):
        if self.isTraceDictLoaded == 0:
            print('No loaded dataset')
            return
        export_Window = export_traces(self.traces_data)
        export_Window.show()
        self.plot_windows.append(export_Window)
        
    def display_label_window(self):
        label_Window = addLabelWindow(self.list_label_classes)
        label_Window.addRequested.connect(self.add_label_item)
        label_Window.removeRequested.connect(self.remove_label_item)
        label_Window.show()
        self.plot_windows.append(label_Window)
        
    def add_label_item(self, text):
        
        self.list_label_classes = self.list_label_classes + [text]
        
        self.label_choose.clear()
        
        self.label_choose.addItems(self.list_label_classes)
        
        for sc in self.shortcuts:
            sc.deleteLater()
        self.shortcuts.clear()
        
        for j in range(len(self.list_label_classes)):
            shortcut = QShortcut(QKeySequence(str(j+1)), self, activated=lambda j=j: self.label_choose.setCurrentIndex(j))
            self.shortcuts.append(shortcut)
        
    def remove_label_item(self, text):
        
        self.list_label_classes.remove(text)
            
        self.label_choose.clear()
        
        self.label_choose.addItems(self.list_label_classes)
        
        for sc in self.shortcuts:
            sc.deleteLater()
        self.shortcuts.clear()
        
        for j in range(len(self.list_label_classes)):
            shortcut = QShortcut(QKeySequence(str(j+1)), self, activated=lambda j=j: self.label_choose.setCurrentIndex(j))
            self.shortcuts.append(shortcut)
        
    def closeEvent(self, event):
        reply = QMessageBox.question(
            self,
            "Exit Confirmation",
            "Are you sure you want to quit?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )

        if reply == QMessageBox.Yes:
            event.accept()  # Close the window
            QApplication.closeAllWindows()
        else:
            event.ignore()  # Keep the window open
            
    def save(self):
        if not self.isWindowModified():
            return
        if not self.fileName:
            self.saveAs()
        else:
            self.auto_save_calib['method'] = self.method_calib.currentText()
            self.auto_save_calib['matrix_align'] = self.matrix_align
            
            autosave_intermed_res(self.fileName, self.auto_save_calib)
        
        self.setWindowModified(False)
    
    def saveAs(self):

        saveWindow = Save_Window(self.matrix_align, self.traces_data, self.traces_DO_data, self.traces_AO_data)
        saveWindow.show()
        self.plot_windows.append(saveWindow)
        
    def merge_traces_dataset(self):
        
        files, _ = QFileDialog.getOpenFileNames(self, "Select the JSON traces datasets to be merged", "", "All Files (*);;Text Files (*.json)")
        
        if files:
            merge_OpenFRET_datasets(files)
        else:
            print('No files selected')
        
    
    def load_traces_dataset(self):
        
        self.selected_traces = {}
        
        options = QFileDialog.Options()
        # You can set options like QFileDialog.DontUseNativeDialog if needed
        file_name, _ = QFileDialog.getOpenFileName(self, "Select the JSON traces dataset", "", "All Files (*);;Text Files (*.json)", options=options)
        if file_name:
            self.traces_data =  read_data(file_name)#load(json_file)
            self.isTraceDictLoaded = 1
            self.trace_IDs = [str(i) for i in range(len(self.traces_data.traces))] #list(self.traces_data.keys())
            self.choose_trace_ID.addItems(self.trace_IDs)
            print('Traces dataset loaded!')
            
            self.figure_filtering.clear()
            self.ax_trace_plot = self.figure_filtering.add_subplot(111)
            self.ax_trace_plot.axhline(y=0, color ='k', linestyle = '--')

            if self.DD_trace_checkbox.isChecked():
                self.line_DD, = self.ax_trace_plot.plot(np.array(self.traces_data.traces[0].channels[0].data), 'orange')
            else:
                self.line_DD, = self.ax_trace_plot.plot([], [], 'orange')
            if self.DA_trace_checkbox.isChecked():
                self.line_DA, = self.ax_trace_plot.plot(np.array(self.traces_data.traces[0].channels[1].data), 'red')
            else:
                self.line_DA, = self.ax_trace_plot.plot([], [], 'red')
            if self.AA_trace_checkbox.isChecked():
                self.line_AA, = self.ax_trace_plot.plot(np.array(self.traces_data.traces[0].channels[2].data), 'gray')
            else:
                self.line_AA, = self.ax_trace_plot.plot([], [], 'gray')
            if self.DD_DA_trace_checkbox.isChecked():
                self.line_DD_DA, = self.ax_trace_plot.plot(np.array(self.traces_data.traces[0].channels[0].data) + np.array(self.traces_data.traces[0].channels[1].data), 'blue')
            else:
                self.line_DD_DA, = self.ax_trace_plot.plot([], [], 'blue')
                
            self.canvas_filtering.draw_idle()
        
        self.span = SpanSelector(
            self.ax_trace_plot,
            self.drag_select,
            "horizontal",
            useblit=True,
            props=dict(alpha=0.5, facecolor="tab:blue"),
            interactive=True,
            drag_from_anywhere=True
        )
    
    def drag_select(self, xmin, xmax):
        
        current_value = self.choose_trace_ID.currentText()
        
        len_x = len(self.traces_data.traces[int(current_value)].channels[0].data)
        
        x = np.linspace(0, len_x - 1, len_x)
        
        self.indmin, self.indmax = np.searchsorted(x, (xmin, xmax))
        self.indmax = min(len(x) - 1, self.indmax)
        print(self.indmin)
        print(self.indmax)
        
        
        class_label = self.label_choose.currentText()
        
        self.traces_data.traces[int(current_value)].metadata[class_label]= {'indmin': int(self.indmin), 'indmax': int(self.indmax)}
        
        
        self.traces_data.traces[int(current_value)].metadata['background_correction'] = self.back_method.currentText()
        if self.back_method.currentText() in ['Total variation', 'Min. of TV']:
            self.traces_data.traces[int(current_value)].metadata['TV_lambda'] = int(self.spin_box_TV_param.value())
        
def start_GUI():
    app = QApplication(sys.argv)
    window = FileSelectorApp()
    window.show()
    sys.exit(app.exec_())

# if __name__ == "__main__":
#     app = QApplication(sys.argv)
#     window = FileSelectorApp()
#     window.show()
#     sys.exit(app.exec_())