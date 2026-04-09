# -*- coding: utf-8 -*-
"""
Created on Fri Dec 12 12:20:59 2025

@author: molcre0000
"""

import tkinter as tk
from tkinter import ttk, messagebox
import matplotlib.pyplot as plt
from utils import *
from traces_extractor import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.widgets import RectangleSelector
from PIL import Image, ImageTk
from tifffile import TiffFile, imread
from GUI_HMM_utils import *

class TimeSeriesViewer(tk.Tk):
    def __init__(self):   # def __init__(self, data_dict):
        super().__init__()
        self.title("Traces viewer")
        self.state('zoomed')
        self.data_dict = {}
        self.HMM_predict_dict = {}
        self.current_slice = 0
        self.num_slices = 0
        
        self.HMM_predict = None
        
        self.data_dict_selected = {}
        
        
        self.movie_donor = []
        
        self.movie_DA = []
        
        self.movie_AA = []
        
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=0)
        self.columnconfigure(2, weight=0)
        self.columnconfigure(3, weight=0)
        self.rowconfigure(2, weight=1)
        # Create UI components
        self.create_menu()
        panel_1 = self.create_panel_1()
        panel_1.grid(column=0, row=0)#, sticky=tk.W)
        
        panel_2 = self.create_panel_2()
        panel_2.grid(column=1, row=0)#, sticky=tk.W)
        
        # panel_movie = self.create_panel_movie_viewer()
        # panel_movie.grid(column=2, row=0, sticky=tk.W)
        
        
        self.create_widgets()
        
        self.plot_frame = ttk.Frame(self)
        self.plot_frame.grid(row=2, column=0, sticky='NSWE', columnspan=5)

        # Placeholder for the matplotlib figure
        self.fig = None
        self.canvas = None
        
        self.protocol("WM_DELETE_WINDOW", self.on_closing)
        
        # rectangle selection
        # self.xlim_label = tk.Label(self, text="X Limits: ")
        # self.xlim_label.pack()
        
    def create_menu(self):
        self.menubar = tk.Menu(self)
        self.config(menu=self.menubar)
        
        # create the file_menu
        self.file_menu = tk.Menu(
            self.menubar,
            tearoff=0
        )
        
        # add menu items to the File menu
        self.file_menu.add_command(label='Load JSON traces dataset', command=self.load_JSON_file_in_viewer)
        self.file_menu.add_command(label='Save selected traces as JSON', command=self.save_selected_JSON)
        self.file_menu.add_command(label='Load TIFF movie', command=self.load_movie_in_viewer)
        self.file_menu.add_separator()

        # add Exit menu item
        self.file_menu.add_command(
            label='Exit',
            command=self.exit_viewer
        )

        # add the File menu to the menubar
        self.menubar.add_cascade(
            label="File",
            menu=self.file_menu
        )
    
    def create_panel_1(self):
        frame = ttk.Frame(self)

        # grid layout for the input frame
        frame.columnconfigure(0, weight=1)
        #frame.columnconfigure(0, weight=3)
        
        subtitle_panel1_label = ttk.Label(frame, text='Visualization', font=('Arial', 12, 'bold'))
        subtitle_panel1_label.grid(column=0, row=0, pady=3)
        
        # Find what
        ttk.Label(frame, text='Trace ID:').grid(column=0, row=1, sticky=tk.W, pady=3)
        self.series_var = tk.StringVar()
        self.series_combo = ttk.Combobox(frame, textvariable=self.series_var)
        self.series_combo['values'] = list(self.data_dict.keys())
        self.series_combo.state(['readonly'])
        self.series_combo.grid(column=1, row=1, sticky=tk.W)
        ttk.Button(frame, text='Plot', command=self.plot_series).grid(column=2, row=1)

        ttk.Label(frame, text='Show first %').grid(column=0, row=2, sticky=tk.W)
        self.firstPerc = tk.StringVar(value="10")
        self.entryfirstPerc = tk.Entry(frame, textvariable=self.firstPerc)
        self.entryfirstPerc.grid(row=2, column=1)
        firstPerc_button = tk.Button(frame, text="Show", command=self.plot_series)
        firstPerc_button.grid(row=2, column=2)
        
        # Select channel
        self.choiceDD = tk.IntVar()
        self.choiceDA = tk.IntVar()
        self.choiceAA = tk.IntVar()
        self.choiceDD_DA = tk.IntVar()
        
        chkbtn_DD = tk.Checkbutton(frame,text="DD", onvalue=1, command=self.plot_series, offvalue=0, variable=self.choiceDD)
        chkbtn_DD.grid(row=1, column=3)
        chkbtn_DA = tk.Checkbutton(frame,text="DA", onvalue=1, command=self.plot_series, offvalue=0, variable=self.choiceDA)
        chkbtn_DA.grid(row=1, column=4)
        chkbtn_AA = tk.Checkbutton(frame,text="AA", onvalue=1, command=self.plot_series, offvalue=0, variable=self.choiceAA)
        chkbtn_AA.grid(row=2, column=3)
        chkbtn_DD_DA = tk.Checkbutton(frame,text="DD+DA", onvalue=1, command=self.plot_series, offvalue=0, variable=self.choiceDD_DA)
        chkbtn_DD_DA.grid(row=2, column=4)
        
        self.choicebreak = tk.IntVar()
        chkbtn_break = tk.Checkbutton(frame, text="Show break events", onvalue=1, command=self.plot_series, offvalue=0, variable=self.choicebreak)
        chkbtn_break.grid(row=3, column=0)
        
        
        separator = ttk.Separator(frame, orient='vertical')
        separator.grid(row=0, column=5, sticky="ns", padx=10, rowspan=4)
        
        return frame
    
    def create_panel_2(self):
        frame = ttk.Frame(self)
        frame.columnconfigure(0, weight=1)
        subtitle_panel2_label = ttk.Label(frame, text='Filtering', font=('Arial', 12, 'bold'))
        subtitle_panel2_label.grid(column=0, row=0, pady=3)
        
        #ttk.Label(frame, text='Trace selected:').grid(column=0, row=1, sticky=tk.W)
        self.choice_trace = tk.IntVar()
        button_selected = tk.Checkbutton(frame, text='Trace selected', onvalue=1, command=self.select_trace_button, offvalue=0, variable=self.choice_trace)
        button_selected.grid(row=1, column=0)
        
        ttk.Label(frame, text='SNR threshold:').grid(column=1, row=1, sticky=tk.W)
        self.entry = tk.Entry(frame)
        self.entry.grid(row=1, column=2)
        filter_button = tk.Button(frame, text="Filter", command=self.filter_dictionary)
        filter_button.grid(row=1, column=3)
        HMM_button = tk.Button(frame, text="HMM", command=self.HMM_call)
        HMM_button.grid(row=2, column=3)
        
        separator = ttk.Separator(frame, orient='vertical')
        separator.grid(row=0, column=5, sticky="ns", padx=10, rowspan=4)
        
        
        return frame
    
    def create_panel_movie_viewer(self):
        frame = ttk.Frame(self)
        frame.columnconfigure(0, weight=1)
        subtitle_panel_movie_label = ttk.Label(frame, text='Movie spot viewer', font=('Arial', 12, 'bold'))
        subtitle_panel_movie_label.pack(fill=tk.BOTH, expand=True)
        
        self.movie_label = ttk.Label(frame)
        self.movie_label.pack(fill=tk.BOTH, expand=True)
        
        self.movie_slider = ttk.Scale(frame, from_=0, to=self.num_slices - 1, orient='horizontal', command=self.on_slider_move)
        self.movie_slider.pack(fill='x', padx=10, pady=10)

        # Label to show current slice index
        self.slice_label = ttk.Label(frame, text=f"Slice: {self.current_slice + 1}/{self.num_slices}")
        self.slice_label.pack()
        
        self.DD_movie = tk.IntVar()
        self.DD_movie.set(1)
        self.DA_movie = tk.IntVar()
        self.DA_movie.set(0)
        self.AA_movie = tk.IntVar()
        self.AA_movie.set(0)
        
        check_DD_movie = tk.Checkbutton(frame,text="DD", onvalue=1, command=self.reload_movie_chanel, offvalue=0, variable=self.DD_movie)
        check_DD_movie.pack()
        check_DA_movie = tk.Checkbutton(frame,text="DA", onvalue=1, command=self.reload_movie_chanel, offvalue=0, variable=self.DA_movie)
        check_DA_movie.pack()
        check_AA_movie = tk.Checkbutton(frame,text="AA", onvalue=1, command=self.reload_movie_chanel, offvalue=0, variable=self.AA_movie)
        check_AA_movie.pack()

        # Display the initial slice
        #self.display_slice(self.current_slice)
        
        return frame
    
    def display_slice(self, index):
        # Extract the slice
        self.XY_spot_coord()
        self.movie_stack_donor = self.movie_donor[:, self.xl:self.xr, self.yu:self.yd]
        self.movie_stack_DA= self.movie_DA[self.xl:self.xr, self.yu:self.yd, :]
        self.movie_stack_AA= self.movie_AA[self.xl:self.xr, self.yu:self.yd, :]
        
        if self.DD_movie.get() == 1:
            slice_img_DD = self.movie_stack_donor[index, :, :]
        else:
            slice_img_DD = np.zeros((self.xr - self.xl, self.yd - self.yu))
        
        if self.DA_movie.get() == 1:
            slice_img_DA = self.movie_stack_DA[:, :, index]
        else:
            slice_img_DA = np.zeros((self.xr - self.xl, self.yd - self.yu))
        
        if self.AA_movie.get() == 1:
            slice_img_AA = self.movie_stack_AA[:, :, index]
        else:
            slice_img_AA = np.zeros((self.xr - self.xl, self.yd - self.yu))
            
        RGB_img = np.stack(((slice_img_DD>255)*255+(slice_img_DD<=255)*slice_img_DD, (slice_img_DA>255)*255+(slice_img_DA<=255)*slice_img_DA, (slice_img_AA>255)*255+(slice_img_AA<=255)*slice_img_AA), axis = 2).astype(np.uint8)
        

        # Convert to PIL Image for Tkinter
        self.pil_img = Image.fromarray(RGB_img, mode='RGB')
        
        #self.img_copy = self.pil_img.copy()

        # Resize image if needed (optional)
        self.pil_img = self.pil_img.resize((int((self.winfo_height()+1)*0.25), int((self.winfo_height()+1)*0.25)), resample=0)

        # Convert to ImageTk
        self.tk_img = ImageTk.PhotoImage(self.pil_img)

        # Update label with image
        self.movie_label.config(image=self.tk_img)
        
        self.movie_label.pack(fill=tk.BOTH, expand=True)

        # Update slice label
        self.slice_label.config(text=f"Slice: {index + 1}/{self.num_slices}")
    
    def on_slider_move(self, val):
        index = int(float(val))
        self.current_slice = index
        self.display_slice(index)


    def create_widgets(self):
        
        self.bind('<Right>', self.select_next)
        self.bind('<Left>', self.select_previous)
        self.bind('<Return>', self.select_trace_enter)
        
    
    def load_JSON_file_in_viewer(self):
        traces_data = load_JSON_traces_data()
        self.data_dict_full = traces_data
        self.data_dict = traces_data
        
        self.series_combo['values'] = list(self.data_dict.keys())
        self.series_combo.set(self.series_combo['values'][0])
        self.series_combo.state(['readonly'])
        self.plot_series()
        
    def load_movie_in_viewer(self):
        root = Tk(className='Open TIFF movie', )
        file_path_donor = askopenfilenames(title="Select the donor first submovie")
        file_path_acceptor = askopenfilenames(title="Select the acceptor first submovie")
        root.destroy()
        with TiffFile(file_path_donor[0]) as tif:
            self.num_slices = len(tif.pages)
            self.X_len, self.Y_len = tif.pages[0].shape
            self.movie_donor = np.zeros((self.num_slices, self.X_len, self.Y_len))
            #self.movie_donor[k,:,:] = generate_chrom_ab_corr_movie(image, self.data_dict.get('chromatic_aberration_corr_matrix'))
            
            k = 0

            for page in tif.pages:
                image = page.asarray()
                self.movie_donor[k,:,:] = generate_chrom_ab_corr_movie(image, self.data_dict.get('chromatic_aberration_corr_matrix'))
                k = k + 1
        tif_acceptor = imread(file_path_acceptor[0])
        self.movie_DA, self.movie_AA = deinterleave_acceptor_channel(tif_acceptor.transpose(1, 2, 0), DA_is = 'odd')
        
        # with TiffFile(file_path_acceptor[0]) as tif:
        #     self.movie_DA = np.zeros((self.num_slices, self.X_len, self.Y_len))
        #     self.movie_AA = np.zeros((self.num_slices, self.X_len, self.Y_len))
        #     k = 0
        #     for page in tif.pages:
        #         image = page.asarray()
        #         self.movie_DA[k,:,:], self.movie_AA[k,:,:] = deinterleave_acceptor_channel(image, DA_is = 'odd')
        #         k = k + 1
        
        panel_movie = self.create_panel_movie_viewer()
        panel_movie.grid(column=4, row=0, sticky=tk.W)
        self.display_slice(0)
        
        
    def XY_spot_coord(self):
        x_coord = self.data_dict[self.series_var.get()].get('x_coord')
        y_coord = self.data_dict[self.series_var.get()].get('y_coord')
        self.xl = np.max((0, x_coord - 5))
        self.xr = np.min((self.X_len - 1, x_coord + 6))
        self.yu = np.max((0, y_coord - 5))
        self.yd = np.min((self.Y_len - 1, y_coord + 6))
    
    def save_selected_JSON(self):
        save_traces_dict(self.data_dict_selected)
    
    def exit_viewer(self):
        response=messagebox.askyesno('Exit','Are you sure you want to exit?')
        if response:
            self.destroy()

    def plot_series(self):
        series_name = self.series_var.get()
        if not series_name:
            messagebox.showwarning("Selection Error", "Please select a time series.")
            return
        data_DD = self.data_dict[series_name].get('Intensity_DD')
        data_DA = self.data_dict[series_name].get('Intensity_DA')
        data_AA = self.data_dict[series_name].get('Intensity_AA')
        if data_DD is None:
            messagebox.showerror("Data Error", f"No data found for {series_name}.")
            return

        # Clear previous plot if exists
        for widget in self.plot_frame.winfo_children():
            widget.destroy()

        # Create matplotlib figure
        factor = np.max((float(self.entryfirstPerc.get())/100, 0.01)) # we cannot display less than 1% of the traces
        
        self.fig, ax = plt.subplots(2)
        if self.choiceDD.get() == 1:
            ax[0].plot(data_DD[:int(len(data_DD)*factor)], 'orange', label='DD')
            ax[1].plot(data_DD, 'orange', label='DD')
        if self.choiceDA.get() == 1:
            ax[0].plot(data_DA[:int(len(data_DA)*factor)], 'red', label='DA')
            ax[1].plot(data_DA, 'red', label='DA')
        if self.choiceAA.get() == 1:
            ax[0].plot(data_AA[:int(len(data_AA)*factor)], 'gray', label='AA')
            ax[1].plot(data_AA, 'gray', label='AA')
        if self.choiceDD_DA.get() == 1:
            ax[0].plot((np.array(data_DD[:int(len(data_DD)*factor)])+np.array(data_DA[:int(len(data_DA)*factor)])), 'blue', label='DD + DA')
            ax[1].plot(np.array(data_DD)+np.array(data_DA), 'blue', label='DD + DA')
        if self.choicebreak.get() == 1:
            ax[1].axvline(x = self.data_dict[series_name].get('bleaching_event_DD'), color = 'orange', label = 'Donor bleaching', ls='--');
            ax[1].axvline(x = self.data_dict[series_name].get('bleaching_event_AA'), color = 'gray', label = 'Acceptor bleaching', ls='--');
        #if self.HMM_predict_dict[series_name] is not None:
        if series_name in list(self.HMM_predict_dict.keys()):
            HMM_seq = self.HMM_predict_dict[series_name]
            factor_HMM = np.max([data_DD, data_DA])
            ax[0].plot(factor_HMM*HMM_seq[:int(len(HMM_seq)*factor)], 'green', label='HMM predict')
            ax[1].plot(factor_HMM*HMM_seq, 'green', label='HMM predict')
            
        ax[0].set_title(f"Time Series: {series_name}")
        ax[1].set_xlabel("Time")
        ax[0].set_ylabel("Intensity")
        ax[1].set_ylabel("Intensity")
        ax[1].legend(loc="upper right")
        ax[0].grid(True)
        ax[1].grid(True)

        # Embed the plot in Tkinter
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2Tk(self.canvas,self.plot_frame)
        toolbar.update()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        plt.close(self.fig)
        self.is_trace_selected()
        
        if len(self.movie_donor) != 0 :
            self.display_slice(int(float(self.movie_slider.get())))
        
        # For the rectangle selection, not working right now
        
        # self.rs = RectangleSelector(ax, self.on_select,
        #                             useblit=True,
        #                             button=[1],  # Left mouse button
        #                             minspanx=5, minspany=5,
        #                             spancoords='pixels',
        #                             drag_from_anywhere = True,
        #                             interactive=True)
        
        # extent = self.rs.extents
        # print(extent[0])
        # self.xlim_label.config(text=f"X Limits: {extent[0]:.2f} to {extent[1]:.2f}")
        
        # self.bind("<Configure>", self.on_resize)
        # #self.update_xlim_label()
        
    def select_next(self, event):
        current_value = self.series_combo.get()
        values = self.series_combo['values']
        if not values:
            return
        try:
            current_index = values.index(current_value)
        except ValueError:
            # If current value not in list, select the first one
            self.series_combo.set(values[0])
            if len(self.movie_donor) != 0 :
                self.display_slice(int(float(self.movie_slider.get())))
            return
        
        next_index = (current_index + 1) % len(values)
        self.series_combo.set(values[next_index])
        self.plot_series()
        self.is_trace_selected()
        if len(self.movie_donor) != 0 :
            self.display_slice(int(float(self.movie_slider.get())))
            
        
    def select_previous(self, event):
        current_value = self.series_combo.get()
        values = self.series_combo['values']
        if not values:
            return
        try:
            current_index = values.index(current_value)
        except ValueError:
            # If current value not in list, select the last one
            self.series_combo.set(values[-1])
            if len(self.movie_donor) != 0 :
                self.display_slice(int(float(self.movie_slider.get())))
            return     
        previous_index = (current_index - 1) % len(values)
        self.series_combo.set(values[previous_index])
        self.plot_series()
        self.is_trace_selected()
        if len(self.movie_donor) != 0 :
            self.display_slice(int(float(self.movie_slider.get())))
    
    def reload_movie_chanel(self):
        self.display_slice(int(float(self.movie_slider.get())))    
    
    
    def select_trace_enter(self, event):
        current_value = self.series_var.get()
        self.data_dict_selected[current_value] = self.data_dict[current_value]
        self.choice_trace.set(1)
        
    def select_trace_button(self):
        current_value = self.series_var.get()
        self.data_dict_selected[current_value] = self.data_dict[current_value]
        self.choice_trace.set(1)
    
    def is_trace_selected(self):
        if self.series_var.get() in self.data_dict_selected.keys():
            self.choice_trace.set(1)
        else:
            self.choice_trace.set(0)
    def HMM_call(self):
        series_name = self.series_var.get()
        DD_channel = self.data_dict[series_name].get('Intensity_DD')
        DA_channel = self.data_dict[series_name].get('Intensity_DA')
        self.HMM_predict = TbT_init(DD_channel, DA_channel, nb_hidden_states = 3, dim = 2, HMM_iter = 10)
        self.HMM_predict_dict[series_name] = self.HMM_predict
        self.plot_series()
    
    def filter_dictionary(self):
        try:
            SNR_threshold = float(self.entry.get())
        except ValueError:
            messagebox.showerror("Invalid input", "Please enter a valid number.")
            return
        
        # Filter the dictionary based on the threshold
        self.data_dict, _ = filter_traces_on_SNR(self.data_dict_full, SNR_thresh = SNR_threshold)
        self.series_combo['values'] = list(self.data_dict.keys())
        self.series_combo.set(self.series_combo['values'][0])
        self.series_combo.state(['readonly'])
        self.plot_series()
        self.is_trace_selected()
        if len(self.movie_donor) != 0 :
            self.display_slice(int(float(self.movie_slider.get())))
   
    def on_closing(self):
        if messagebox.askokcancel('Exit','Are you sure you want to exit?'):
            self.destroy()
    
    # rectangle selection
    
    # def on_select(self, eclick, erelease):
    #     # eclick and erelease are matplotlib events
    #     x1, y1 = eclick.xdata, eclick.ydata
    #     x2, y2 = erelease.xdata, erelease.ydata

    #     # Create rectangle patch
    #     xmin, xmax = sorted([x1, x2])
    #     ymin, ymax = sorted([y1, y2])
    #     width = xmax - xmin
    #     height = ymax - ymin

    #     self.canvas.draw_idle()
    
    # def on_resize(self, event):
    #     # Resize the canvas when the window size changes
    #     self.canvas.draw_idle()
    
    # def update_xlim_label(self):
    #     try:
    #         extent = self.rs.extents
    #         print(extent[0])
    #         self.xlim_label.config(text=f"X Limits: {extent[0]:.2f} to {extent[1]:.2f}")
    #     except ValueError:
    #         print('No box selected')

if __name__ == "__main__":
    app = TimeSeriesViewer()
    app.mainloop()