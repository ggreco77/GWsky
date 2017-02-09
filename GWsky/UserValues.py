# -*- coding: utf-8 -*-

import sys
# checking the python version----------------------------------------#
#--------------------------------------------------------------------#
if sys.version[0] == "3":                                            #
    print('  ===================================================')   #
    print('  ***Please, run the script in  python 2.7.x***')         #
    print('  ===================================================')   #
    sys.exit()                                                       #
#--------------------------------------------------------------------#
    
import healpy as hp
import numpy as np

import json
import time
try:
   import cPickle as pickle
except:
   import pickle

from astropy.vo.samp import SAMPHubError # no work in py3
from astropy.time import Time
from astropy.io.votable import parse
from astropy.io.votable.tree import VOTableFile, Resource, Field, Param
from astropy.utils.data import download_file

from aladinSAMP import AladinViaSAMP, AladinScriptCommands 
samp = AladinViaSAMP()
aladin = AladinScriptCommands()

from Tkinter import *
import tkFileDialog
import tkFont
import tkMessageBox

class UserValuesGUI:
    """GUI with seven labelframe widgets:
        1. Enter skymap information
        2. MOC contour plot - altervative to json file
        3. Insert the Geodetic coordinates of the Observatory
        4. Insert the observation time (UTC)
        5. FoV shape
        6. Query
        7. Statistic window
        8. Open tiling Coverage."""
    
    
    def __init__(self, master):
        """Two dicts are initialized to create the config. files (GWsky_config and GWsky_coords)
          for working with the coverage module.
          Init. entry for the FoV shape: (1) box and (2) circle."""
        
        self.master = master

        self.GWsky_config = dict() # used to create the config. file GWsky_config
        self.GWsky_coords = dict() # used to create the config. file GWsky_coords
        self.shape = 0 # FoV shape: (1) box and (2) circle
        
        self.v = StringVar() # default statistic window
        self.v.set("A")  
        self.window_stat = self.v  # Active/Deactive statistic window
        self.GWsky_config.update({"GWsky_basic" : "A"})

        self.one_label_frame = LabelFrame(master,
                                          text=" 1. Enter skymap information: ")
        self.one_label_frame.grid(row=0, columnspan=7, sticky='W',
                 padx=15, pady=10, ipadx=5, ipady=5)

        self.two_label_frame = LabelFrame(master,
                                             text=" 2. MOC contour plot - altervative to json file: ")
        self.two_label_frame.grid(row=1, columnspan=7, sticky='W',
                 padx=15, pady=10, ipadx=5, ipady=5)

        self.three_label_frame = LabelFrame(master,
                                          text=" 3. Insert the Geodetic coordinates of the Observatory: ")
        self.three_label_frame.grid(row=2, columnspan=10, sticky='W',
                     padx=15, pady=10, ipadx=5, ipady=5)

        self.four_label_frame = LabelFrame(master,
                                            text=" 4. Insert the observation time (UTC): ")
        self.four_label_frame.grid(row=3, columnspan=7, sticky='W',
                       padx=15, pady=10, ipadx=5, ipady=5)

        self.five_label_frame = LabelFrame(master,
                                           text=" 5. FoV shape: ")
        self.five_label_frame.grid(row=4, columnspan=7, sticky='W',
                       padx=15, pady=10, ipadx=5, ipady=5)

        self.six_label_frame = LabelFrame(master,
                                           text=" 6. Query: ")
        self.six_label_frame.grid(row=5, columnspan=7, sticky='W',
                       padx=15, pady=10, ipadx=5, ipady=5)

        self.seven_label_frame = LabelFrame(master,
                                          text=" 7. Statistic window") 
        self.seven_label_frame.grid(row=6, columnspan=7, sticky='W',
                       padx=15, pady=10, ipadx=5, ipady=5)

        self.eight_label_frame = LabelFrame(master,
                                          text=" 8. Initialize tiling coverage") 
        self.eight_label_frame.grid(row=6, columnspan=7, sticky='E',
                       padx=15, pady=10, ipadx=5, ipady=5)


        # 1. Enter skymap information: healpix skymap
        self.entry_skymap_txt = Label(self.one_label_frame,
                                text="healpix skymap:")
        self.entry_skymap_txt.grid(row=0, column=0, padx=5, sticky="E", pady=2)
        
        self.entry_skymap = Entry(self.one_label_frame, width=40, justify=CENTER)
        self.entry_skymap.grid(row=0, column=1, columnspan=6, sticky='WE', padx=5, pady=2)

##        self.browse_1 = Button(self.one_label_frame,
##                               text="Browse ...", command=self.browsecsv)
##        self.browse_1.grid(row=0, column=8, sticky='W', padx=5, pady=2)

        self.show_aladin_1 = Button(self.one_label_frame,
                                    text="Show in Aladin",command=self.show_healpix_skymap)
        self.show_aladin_1.grid(row=0, column=15, sticky='E', padx=5, pady=2)

        # 1.1 Enter skymap information: json contour file
        self.entry_json_txt = Label(self.one_label_frame,
                                text="json contour file:")
        self.entry_json_txt.grid(row=1, column=0, padx=5, sticky="E", pady=2)
        
        self.entry_json = Entry(self.one_label_frame, width=40, justify=CENTER)
        self.entry_json.grid(row=1, column=1, columnspan=6, sticky='E', padx=5, pady=2)

##        self.browse_1 = Button(self.one_label_frame,
##                               text="Browse ...", command=self.browsecsv)
##        self.browse_1.grid(row=1, column=8, sticky='W', padx=5, pady=2)

        self.show_aladin_2 = Button(self.one_label_frame,
                                    text="Show in Aladin",command=self.probability_contours)
        self.show_aladin_2.grid(row=1, column=15, sticky='E', padx=5, pady=2)

        # 2. MOC generation
        self.entry_from_txt = Label(self.two_label_frame,
                                    text="from:")     
        self.entry_from_txt.grid(row=2, column=0, sticky='W', padx=5, pady=2)

        from_default = StringVar(self.two_label_frame, value="10") # default

        self.entry_from = Entry(self.two_label_frame, width=5, justify=CENTER,
                                textvariable=from_default)
        self.entry_from.grid(row=2, column=1, columnspan=8, pady=2, sticky='WE')

        self.entry_to_txt = Label(self.two_label_frame,
                                     text="to:")
        self.entry_to_txt.grid(row=2, column=10, sticky='W', padx=5, pady=2)

        to_default = StringVar(self.two_label_frame, value="90")   # default

        self.entry_to = Entry(self.two_label_frame, width=5, justify=CENTER,
                              textvariable=to_default)
        self.entry_to.grid(row=2, column=11, columnspan=8, pady=2, sticky='WE')

        self.entry_step_txt = Label(self.two_label_frame,
                                    text="step:")
        self.entry_step_txt.grid(row=2, column=19, sticky='W', padx=5, pady=2)

        step_default = StringVar(self.two_label_frame, value="10")   # default

        self.entry_step = Entry(self.two_label_frame, width=5, justify=CENTER,
                                textvariable=step_default)
        self.entry_step.grid(row=2, column=20, columnspan=8, pady=2, sticky='WE')
         
        self.show_moc = Button(self.two_label_frame,
                                    text="MOC contour",command=self.moc_contour)
        self.show_moc.grid(row=2, column=50, sticky='E', padx=5, pady=2)

        self.show_intersection = Button(self.two_label_frame,
                                    text="MOC intersection",command=self.moc_intersection)
        self.show_intersection.grid(row=2, column=60, sticky='E', padx=5, pady=2)

        # 3. Site location
        self.entry_latitude_txt = Label(self.three_label_frame,
                                    text="latitude [deg]:")     
        self.entry_latitude_txt.grid(row=3, column=0, sticky='W', padx=5, pady=2)

        latitude_default = StringVar(self.three_label_frame, value="34.2247")   # default
        self.entry_latitude = Entry(self.three_label_frame, width=10, justify=CENTER,
                                    textvariable=latitude_default)
        self.entry_latitude.grid(row=3, column=1, columnspan=8, pady=2, sticky='WE')

        self.entry_longitude_txt = Label(self.three_label_frame,
                                     text="longitude [deg]:")
        self.entry_longitude_txt.grid(row=3, column=10, sticky='W', padx=5, pady=2)

        longitude_default = StringVar(self.three_label_frame, value="-118.0572")   # default
        self.entry_longitude = Entry(self.three_label_frame,width=10, justify=CENTER,
                                     textvariable=longitude_default)
        self.entry_longitude.grid(row=3, column=11, columnspan=8, pady=2, sticky='WE')

        self.entry_altitude_txt = Label(self.three_label_frame,
                                    text="altitude [m]:")
        self.entry_altitude_txt.grid(row=3, column=19, sticky='W', padx=5, pady=2)

        altitude_default = StringVar(self.three_label_frame, value="1742")   # default
        self.entry_altitude = Entry(self.three_label_frame,width=10, justify=CENTER,
                                    textvariable=altitude_default)
        self.entry_altitude.grid(row=3, column=20, columnspan=8, pady=2, sticky='WE')

        # 4. Obs time (UTC)
        self.entry_time_txt = Label(self.four_label_frame,
                                    text=" yyyy-mm-dd hh:mm:ss:")
        
        time_now_utc = Time.now() # getting time now (UTC) as default
        
        time_default = StringVar(self.four_label_frame, value=time_now_utc)   # default

        self.entry_time = Entry(self.four_label_frame, width=25, justify=CENTER,
                                textvariable=time_default)
        
        self.entry_time_txt.grid(row=6, column=2, columnspan=2,
                                 sticky='W', padx=5, pady=2)
        self.entry_time.grid(row=6, column=4, sticky='WE')

        # 5. FoV shape
        var = IntVar()
        
        # box  
        self.radio_box = Radiobutton(self.five_label_frame, text="box =>", state=ACTIVE,
                                     variable=var, value=1, command = self.set_box)
        self.radio_box.grid(row=4, column=1, columnspan=3, pady=2, sticky='WE')

        self.entry_width = Label(self.five_label_frame, text="width [deg]:")
        self.entry_width.grid(row=4, column=10, sticky='W', padx=5, pady=2)

        width_default = StringVar(self.five_label_frame, value="3")   # default 
        self.entry_width = Entry(self.five_label_frame, width=6, justify=CENTER,
                                 textvariable=width_default)
        self.entry_width.grid(row=4, column=15, columnspan=8, pady=2, sticky='WE')

        self.entry_height = Label(self.five_label_frame, text="height [deg]:")
        self.entry_height.grid(row=4, column=30, sticky='W', padx=5, pady=2)

        height_default = StringVar(self.five_label_frame, value="3")   # default 
        self.entry_height = Entry(self.five_label_frame,width=6, justify=CENTER,
                                  textvariable=height_default)                  
        self.entry_height.grid(row=4, column=35, columnspan=8, pady=2, sticky='WE')

        # circle
        self.radio_circle = Radiobutton(self.five_label_frame, text="circle =>",
                                        variable=var, value=2, command = self.set_circle)
        self.radio_circle.grid(row=5, column=1, columnspan=5, pady=2, sticky='WE')

        self.entry_radius = Label(self.five_label_frame, text="radius [deg]:")
        self.entry_radius.grid(row=5, column=10, sticky='W', padx=5, pady=2)

        self.entry_radius = Entry(self.five_label_frame,width=6, justify=CENTER)
        self.entry_radius.grid(row=5, column=15, columnspan=8, pady=2, sticky='WE')

        # 6 query; default GLADE catalog
        self.entry_catalog_txt = Label(self.six_label_frame,
                                    text="Vizier catalog")     
        self.entry_catalog_txt.grid(row=6, column=0, sticky='W', padx=5, pady=2)

        catalog_default = StringVar(self.six_label_frame, value="VII/275/glade1")   # default
        self.entry_catalog = Entry(self.six_label_frame, width=15, justify=CENTER,
                                    textvariable=catalog_default)
        self.entry_catalog.grid(row=6, column=1, columnspan=8, pady=2, sticky='WE')

        self.entry_columns_1_txt = Label(self.six_label_frame,
                                    text="columns & filters")     
        self.entry_columns_1_txt.grid(row=6, column=12, sticky='W', padx=5, pady=2)

        column_1_default = StringVar(self.six_label_frame, value="Dist")   # default
        self.entry_column_1 = Entry(self.six_label_frame, width=8, justify=CENTER,
                                    textvariable=column_1_default)
        self.entry_column_1.grid(row=6, column=18, columnspan=8, pady=2, sticky='WE')

        column_filter_1_default = StringVar(self.six_label_frame, value="<=300")   # default
        self.entry_column_filter_1 = Entry(self.six_label_frame, width=8, justify=CENTER,
                                    textvariable=column_filter_1_default)
        self.entry_column_filter_1.grid(row=6, column=30, columnspan=8, pady=2, sticky='WE')

        column_2_default = StringVar(self.six_label_frame, value="Bmag")   # default
        self.entry_column_2 = Entry(self.six_label_frame, width=8, justify=CENTER,
                                    textvariable=column_2_default)
        self.entry_column_2.grid(row=6, column=50, columnspan=8, pady=2, sticky='WE')

        column_filter_2_default = StringVar(self.six_label_frame, value=">0")   # default
        self.entry_column_filter_2 = Entry(self.six_label_frame, width=8, justify=CENTER,
                                    textvariable=column_filter_2_default)
        self.entry_column_filter_2.grid(row=6, column=70, columnspan=8, pady=2, sticky='WE')

        # 7 Active/Deactive statistic window: Initialize Active __init__
        self.gwsky_basic_active = Radiobutton(self.seven_label_frame, text="Active",variable=self.v, value="A",
                                              command=self.set_active)
        self.gwsky_basic_active.grid(row=7, column=0,sticky='E', padx=0, pady=2)
        
        self.gwsky_basic_deactive = Radiobutton(self.seven_label_frame, text="Deactive",variable=self.v, value="D",
                                                command=self.set_deactive)
        self.gwsky_basic_deactive.grid(row=7, column=1, sticky='E', padx=0, pady=2)
        
        # 8 Launching coverage module. Init. files are created.
        self.init = Button(self.eight_label_frame,
                                    text="Launch tiling coverage",command=self.launch_tiles)
        self.init.grid(row=7, column=75, sticky='E', padx=60, pady=2)

# Action    
    def show_healpix_skymap(self):
        """Inserting a valid LIGO/Virgo probability skymap
                and sending to Aladin plane."""
        try:
            prob = hp.read_map(self.entry_skymap.get(), verbose = False)
            samp.send_file(self.entry_skymap.get()) 
            aladin.rename(self.entry_skymap.get())
        except ValueError as value_error:
            tkMessageBox.showerror ('Error: 1. Enter skymap information',
                                    value_error)
        except IOError as io_error:
            tkMessageBox.showerror ('Error: 1. Enter skymap information',
                                    io_error)
        except SAMPHubError as samphub_error:
            tkMessageBox.showerror ('Error: 1. Enter skymap information',
                                    samphub_error)

##    def browsecsv(self):
##        from tkFileDialog import askopenfilename
##
##        #Tk().withdraw() 
##        filename_from_browser = askopenfilename(filetypes=[("fits files","*.fits.gz")])
##        #fits = filename.read()
##        print filename_from_browser

    def healpix_skymap(self):
        """Inserting a valid LIGO/Virgo probability skymap in GWsky_config."""
        
        try:
            prob = hp.read_map(self.entry_skymap.get(), verbose = False)
        except ValueError as value_error:
            tkMessageBox.showerror ('Error: 1. Enter skymap information',
                                    value_error)
        except IOError as io_error:
            tkMessageBox.showerror ('Error: 1. Enter skymap information',
                                    io_error)
        else:
            npix = len(prob) 
            nside = hp.npix2nside(npix)

            return self.GWsky_config.update({"skymap" : self.entry_skymap.get(),
                                             "nside" : nside})

    def plot_prob_contours(self, data, contour_pieces, percentile_json):
        """Managing the probability-contour file."""       
        
        i = 0
        for i in range(0, contour_pieces):
            contour = data['contours'][i]
            percentile = contour['name']

            if percentile == percentile_json:
                values = contour['coords']
                line = (str(values).replace('[' , '' ).replace(']' , ''))
                aladin.draw_line (line)
        
    def probability_contours(self):
        """Inserting a valid LIGO/Virgo probability contours in json format."""
            
        short_name = self.entry_json.get()
        suffix = short_name[0:5]
        
        aladin.md('JSON' + suffix) # creating a stack folder
        
        try:
            with open(self.entry_json.get()) as data_file:
                data = json.load(data_file)                    
        except IOError as io_error:
            tkMessageBox.showerror ('Error: json contour file', io_error)

        else:
            contour_pieces = len(data['contours'])

            percentile = ('10-percentile','20-percentile','30-percentile','40-percentile',
                          '50-percentile','60-percentile', '70-percentile','80-percentile',
                          '90-percentile')           
            try:
                for percentile_json in percentile:
                    aladin.draw_newtool (percentile_json + suffix) # plotting in different Aladin planes
                    self.plot_prob_contours(data, contour_pieces, percentile_json)
                    aladin.mv(percentile_json + suffix,'JSON' + suffix)
            except SAMPHubError as samphub_error:
                tkMessageBox.showerror ('Error: json contour file', samphub_error)

    def moc_contour(self):
        """Showing MOC plot in a folder."""

        short_name = self.entry_skymap.get()
        suffix = short_name[0:5]
        
        aladin.md('MOC' + suffix) # creating a stack folder
        
        _from = int(self.entry_from.get())
        _to = int(self.entry_to.get())+10
        _step = int(self.entry_step.get())

        # better two "for cicles" in Aladin
        for i in range(_from ,_to, _step):
            
            aladin.cmoc((i/100.0), self.entry_skymap.get(),'moc'+str(i/100.0) + suffix)

        for i in range(int(self.entry_from.get()),(int(self.entry_to.get())+10),
                       int(self.entry_step.get())):
            
            aladin.set_moc('moc'+str(i/100.0)+ suffix)

            aladin.mv('moc'+str(i/100.0)+ suffix,'MOC' + suffix)

    def moc_intersection(self):
        """"MOC intersection with 90% of probability with the VST MOC.
            VST MOC mainly consists of Atlas and Kids survey."""
        
        samp.send_file( 'VST_MOC.fits' )  # sending to Aladin
        aladin.rename('VST_MOC.fits')
        
        aladin.moc_inter('VST_MOC.fits', 'moc0.9', 'INTERSECTION')     

    def telescope_site(self):
        """Inserting the Geodetic coordinates of the Observatory in 
           GWsky_config dict."""

        try:
            float(self.entry_latitude.get())
            float(self.entry_longitude.get())
            float(self.entry_altitude.get())
        except ValueError as value_error:
            message_error= 'insert valid input: latitude [deg], longitude [deg], altitude [deg]'
            tkMessageBox.showerror ('Error: 2. Insert the Geodetic coordinates of the Observatory',
                                    message_error)
            raise
        else:
            return self.GWsky_config.update({"latitude" : float(self.entry_latitude.get()),
                                             "longitude" : float(self.entry_longitude.get()),
                                             "altitude" : float(self.entry_altitude.get())})

    def query_params(self):
        """Inserting the query parameters in GWsky_config dict."""

        return self.GWsky_config.update({"catalog" : self.entry_catalog.get(),
                                         "column_1" : self.entry_column_1.get(),
                                         "column_2" : self.entry_column_2.get(),
                                         "filter_1" : self.entry_column_filter_1.get(),
                                         "filter_2" : self.entry_column_filter_2.get()})

    def starting_time(self):
        """Inserting the observation time yyyy-mm-gg hh:mm:ss."""
        
        try:
            time = Time(self.entry_time.get())
        except ValueError as value_error:
            message_error= 'insert a valid obs time: yyyy-mm-gg hh:mm:ss'
            tkMessageBox.showerror ('Error: 3. Insert the observation time:',
                                    message_error)

        return self.GWsky_config.update({"obs_time" : self.entry_time.get()})

    def set_box(self):
        """Setting shape from Radio Btn - 1: box."""

        self.shape = 1
        return self.shape

    def set_circle(self):
        """Setting shape from Radio Btn - 2: circle."""

        self.shape = 2
        return self.shape

    def set_fov_shape(self):
        if self.shape ==1:
            try:
                float(self.entry_width.get())
                float(self.entry_height.get())
            except ValueError as value_error:
                pass
            else:
                return self.GWsky_config.update({"fov_width" : float(self.entry_width.get()),
                                                 "fov_height" : float(self.entry_height.get()),
                                                 "fov_shape" : 1})
        else:
            try:
                float(self.entry_radius.get())
            except ValueError as value_error:
                pass
            else:
                return self.GWsky_config.update({"fov_width" : float(self.entry_radius.get()),
                                             "fov_height" : float(self.entry_radius.get()),
                                             "fov_radius" : float(self.entry_radius.get()),
                                             "fov_shape" : 2})

    def set_active(self):
        """Setting statistic-window mode from radio button;
           "A": Active"""
      
        self.window_stat = "A"
        return self.window_stat
    
    def set_deactive(self):
        """Setting statistic-window mode from radio button;
           "D": Deactive"""

        self.window_stat = "D"
        return self.window_stat

    def GWsky_basic(self):
        """Setting statistic window: "A": Active; "D": Deactive"""

        if self.window_stat =="A":
            return self.GWsky_config.update({"GWsky_basic" : "A"})
        elif self.window_stat =="D":
            return self.GWsky_config.update({"GWsky_basic" : "D"})
            
    def check_fov_entries(self):
        """Checking the FoV entries: width and height for a box and radius for a circle."""
        
        if self.shape ==1:
            try:
                float(self.entry_width.get())
                float(self.entry_height.get())
            except ValueError as value_error:
                message_error = 'insert valid values: width [deg] and height [deg]'
                tkMessageBox.showerror ('Error: 4. FoV shape',
                                        message_error)
        elif self.shape ==2:
            try:
                float(self.entry_radius.get())
            except ValueError as value_error:
                message_error = 'insert a valid value: radius [deg]'
                tkMessageBox.showerror ('Error: 4. FoV shape',
                                        message_error)
        else:
            message_error = 'select FoV shape: box or circle'
            tkMessageBox.showerror ('Error: 4. FoV shape',
                                    message_error)

    def set_fov_template(self):
        """Downloading footprint template from ggreco77 Git
               - (1) box FoV: (2) circle FoV."""
        
        if self.shape ==1:
            url_id = "https://raw.githubusercontent.com/ggreco77/GWsky/master/footprint_box"
            template_fov_footprint = download_file(url_id, cache=True, timeout=300)

            return self.make_fov_footprint(template_fov_footprint)
            
        else:
            url_id = "https://raw.githubusercontent.com/ggreco77/GWsky/master/footprint_circle"
            template_fov_footprint = download_file(url_id, cache=True, timeout=300)

            return self.make_fov_footprint(template_fov_footprint)

    def make_fov_footprint(self, footprint):
        """Creating the user-defined field-of-view footprint using a template
                  from http://aladin.u-strasbg.fr/footprint_editor/"""

        votable = parse(footprint) # reading footprint template
        table = votable.get_first_table()

        for param in table.params: # retrieving table.params
            param

        # box or circle footprint
        if param.ID == 'radius':
            try:
                param.value = float(self.entry_radius.get())*3600.0
                votable.to_xml('GWsky_fov.vot') # VOTable file output
            except ValueError as value_error:
                raise
        else:
            try:
                data = table.array      
                fov_width_arcsec = float(self.entry_width.get())*3600.0  
                fov_height_arcsec = float(self.entry_height.get())*3600.0
    
                data[0] = - fov_width_arcsec /  2.0,   fov_height_arcsec / 2.0
                data[1] =   fov_width_arcsec /  2.0,   fov_height_arcsec / 2.0
                data[2] =   fov_width_arcsec /  2.0, - fov_height_arcsec / 2.0
                data[3] = - fov_width_arcsec /  2.0, - fov_height_arcsec / 2.0
                votable.to_xml('GWsky_fov.vot') # VOTable file output
            except UnboundLocalError:
                raise
            except ValueError:
                raise
            
        return samp.send_file( 'GWsky_fov.vot' ) 
                

    def _find_highest_pixel(self, infile):
        """Finding the sky position of the highest probability pixel."""

        hpx = hp.read_map(infile, verbose = False)
        npix = len(hpx)
        nside = hp.npix2nside(npix) 

        ipix_max = np.argmax(hpx)
        hpx[ipix_max] 
        theta, phi = hp.pix2ang(nside, ipix_max)
        
        ra_max = np.rad2deg(phi) 
        dec_max = np.rad2deg(0.5 * np.pi - theta)

        return round(ra_max,5), round(dec_max,5)

    def maximum_pixel(self):
        """Return the sky position of the maximum probability pixel."""

        try:
            self._ra_max, self._dec_max = self._find_highest_pixel(self.entry_skymap.get())
                 
            print (' The highest probability pixel is located at ')
            print (' RA =' + str('% .5f' % self._ra_max)+'°' + 'and Dec =' + str('% .5f' % self._dec_max)+'°.')

            return self.GWsky_config.update({"ra_max_pixel" : self._ra_max,
                                             "dec_max_pixel" : self._dec_max})
        except ValueError as value_error:
            raise
        except IOError as io_error:
            raise

    def starting_sky_position(self):
        """Inserting the starting sky coordinates in "GWsky_coords" file;
            default: the maximum probability pixel"""

        try:
            self.GWsky_coords.update({"input_ra" : self._ra_max,
                                      "input_dec" : self._dec_max})

            with open('GWsky_coords', 'wb') as data:
                return pickle.dump(self.GWsky_coords, data)
        except AttributeError as attribute_error:
            raise

    def make_GWsky_config(self):
        """Creating the configuration file "GWsky_config". """

        with open('GWsky_config', 'wb') as data:
            return pickle.dump(self.GWsky_config, data)

    def make_selected_pointing_file(self):
        """Cleaning the file "selected_pointings.txt"."""
        
        open('GWsky_pointings.txt', 'w').close() 
                                                    
    def launch_tiles(self):
        """Command related to the ... Btn"""
        
        self.healpix_skymap()
        self.telescope_site()
        self.query_params()
        self.starting_time()
        
        self.set_fov_shape()
        self.check_fov_entries()
        self.set_fov_template()

        self.maximum_pixel()
        self.starting_sky_position()

        self.GWsky_basic()
        self.make_GWsky_config()
        self.make_selected_pointing_file()
                              
        message= ' The highest probability pixel is located at  RA =' + str('% .5f' % self._ra_max)+'°' + 'and Dec =' + str('% .5f' % self._dec_max)+'°.'
        tkMessageBox.showinfo('User Values has been initialized', message)
        
        import coverage

# running
root = Tk()
root.title('GWsky - User Values ')
uvGUI=UserValuesGUI(root)
root.mainloop()
