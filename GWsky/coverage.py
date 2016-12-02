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

    

import pickle

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import matplotlib.dates as mdates

import dateutil

from Tkinter import *
import tkMessageBox
import tkFont

from math import cos, sin, acos, degrees, radians

import astropy
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import TimeDelta
from astroquery.vizier import Vizier 

import numpy as np
import healpy as hp
import pandas
from scipy.stats import norm

from aladinSAMP import AladinViaSAMP, AladinScriptCommands # class access
samp = AladinViaSAMP()
aladin = AladinScriptCommands()
samp.send_file('user_fov.vot')  # loading user FoV footprint


class Airmass(object):
    """Airmass calculation."""
    
    AIRMASS_MIN = 1
    AIRMASS_MAX = 5.8

    def __init__(self):
        """Initializing with the input values stored in "GWsky_config"
            - observatory site - """

        try:
            with open('GWsky_config', 'rb') as data:
                config_GWsky = pickle.load(data)
            for k, v in config_GWsky.items():
                setattr(self, k, v)
        except IOError as io_error:
            tkMessageBox.showerror('Error', 'Run UserValues module')
                
   
    def airmass(self, ra, dec, lat_input, lon_input, height_input,
                time_input, AIRMASS_MIN, AIRMASS_MAX):
        """ Airmass calculation at a given time in a particular site."""
     
        # Geodetic coordinates of observatory
        observatory = astropy.coordinates.EarthLocation(lat = lat_input*u.deg,
                                                        lon = lon_input*u.deg,
                                                        height = height_input*u.m)
        max_prob_p = SkyCoord(ra = ra*u.deg,
                              dec=dec*u.deg, frame='icrs') # ipix sky coordinates        
        time = Time(time_input) # Time object
        altaz = max_prob_p.transform_to(AltAz(obstime=time,
                                              location=observatory)) # Altitude-Azimuth system
                                                                        #--> Horizontal coordinates     
        airmass_value = altaz.secz # 1/cos(altaz)                      

        if airmass_value <= AIRMASS_MIN or airmass_value >= AIRMASS_MAX:
             airmass_value =  "nan" 
        else:     
             airmass_value = round(airmass_value, 2) 

        return airmass_value

    def airmass_step(self, ra, dec):
        """Airmass calculation at a given time in a particular site in steps of one hour."""

        dt = TimeDelta(3600.0, format = 'sec')
        self.obs_time = Time(self.obs_time)

        airmass_list = []
        time_list = []
        cont = 0
        
        while cont < 10:
            time_input = self.obs_time + cont*dt
            val = self.airmass(ra, dec, self.latitude, self.longitude, self.altitude,
                               time_input, self.AIRMASS_MIN, self.AIRMASS_MAX)           
            airmass_list.append(val)
            time_list.append(str(time_input))
            cont+=1
            
        return airmass_list, time_list


class ShowSkyCoverage(Airmass):
    """Get FoV coverage from a starting pointing; cardinal directions are implemented."""
    
    SHIFT_CORRECTION = 0.00001  # A shift correction of 0.00001 is added
                                 # --> to escape math error during the FoV sequence

    def __init__(self):
        """Initializing with the input values stored in "GWsky_coords"
                - FoV center - """
        
        Airmass.__init__(self)
        
        with open('GWsky_coords', 'rb') as data:
            coords_GWsky = pickle.load(data)

        for k, v in coords_GWsky.items():
            setattr(self, k, v)

    def _ra0ra1_distance(self, A, dec0, dec1):
        """From the angular distance cos(A) =
                        sin(Dec1)sin(Dec2)+cos(Dec1)cos(Dec2)cos(ra1-ra2) and
        thus, A = arccos(A);
                         --------------->
         cos(ra1-ra2) = [cos(A)-sin(dec0)sin(dec1)]/[cos(dec0)cos(dec1)]
         and  thus, acos(cos(ra1-ra2))."""

        dec0, dec1, A = radians(dec0),  radians(dec1), radians(A)
        cos_ra0_ra1 = (cos(A)-sin(dec0)*sin(dec1))/(cos(dec0)*cos(dec1))
        ra0ra1 = degrees(acos(cos_ra0_ra1))

        return  round(ra0ra1, 5)
           
    @classmethod
    def updating_center_coordinates(cls, ra, dec):
        """Updating coordinates of FoV center/getting a new FoV"""

        with open('GWsky_coords', 'rb') as data:
            coords_GWsky = pickle.load(data)
              
        coords_GWsky['input_ra'], coords_GWsky ['input_dec'] = ra, dec

        with open('GWsky_coords', 'wb') as data:
            pickle.dump(coords_GWsky, data)

    def fov_center_separation(self, ra1, dec1, ra2, dec2):
        """Distance between 2 consecutive FoV centers [deg]."""
           
        fov_1 = SkyCoord(ra1, dec1, frame='icrs',unit='deg')
        fov_2 = SkyCoord(ra2, dec2, frame='icrs',unit='deg')
        sep = fov_1.separation(fov_2)

        print ('The Distance between 2 consecutive FoV centers is', sep)

    def are_all_same(self, items):
        """Check if all elements of a list are the same."""
        
        return all(x == items[0] for x in items)

    def probability_inside_box(self, infile, ra_vertices, dec_vertices):
         """Return the probability inside a polygon."""

         prob = hp.read_map(self.skymap, verbose=False)
         
         theta = 0.5 * np.pi - np.deg2rad(dec_vertices)
         phi = np.deg2rad(ra_vertices)
         xyz = hp.ang2vec(theta, phi)

         ipix_poly = hp.query_polygon(self.nside, xyz)
     
         probability_inside_polygon = prob[ipix_poly].sum()
         
         return probability_inside_polygon

    def probability_inside_circle(self, ra, dec, radius):
        """Return the probability inside a circle."""

        prob = hp.read_map(self.skymap, verbose=False)

        theta = 0.5 * np.pi - np.deg2rad(dec)
        phi = np.deg2rad(ra)
        radius = np.deg2rad(radius)

        xyz = hp.ang2vec(theta, phi)
        ipix_disc = hp.query_disc(self.nside, xyz, radius)
        probability_inside_disc = prob[ipix_disc].sum()

        return '%.1e' % probability_inside_disc

    def conditional_distance_fov_center(self, ra, dec):
        """Conditional distance distribution along the line of
            sight of FoV center: see https://arxiv.org/pdf/1605.04242v3.pdf
            section 4.4 for more details."""

        prob, header = hp.read_map(self.skymap, h=True,
                                   verbose=False)
        header = dict(header)

        if header['TFIELDS']==4:
            prob, distmu, distsigma, distnorm = hp.read_map(self.skymap, verbose=False,
                                                            field=[0, 1, 2, 3])
            theta = 0.5 * np.pi - np.deg2rad(dec)
            phi = np.deg2rad(ra)
            ipix = hp.ang2pix(self.nside, theta, phi)

            r = np.linspace(0, 300)
            dp_dr = r**2*distnorm[ipix]*norm(distmu[ipix], distsigma[ipix]).pdf(r)
            
            return r, dp_dr
        else:
            r = "nan"
            dp_dr = "nan"

            return r, dp_dr

    def vertices(self, ra_center, dec_center, fov_base, fov_height):
        """Find the vertices of a FOV given a center position (RA[deg], DEC[deg])
           and the FOV size (FOV_base[deg], FOV_height[deg])

                    ***it will be replaced*** """
        
        A =  fov_base/2.0 # upper vertices
     
        dec0_n = (dec_center) + (fov_height / 2.0)
        dec1_n = dec0_n
         
        offset_n = self._ra0ra1_distance(A, dec0_n, dec1_n)

        dec0_s = (dec_center) - (fov_height / 2.0) # lower vertices
        dec1_s = dec0_s  

        offset_s = self._ra0ra1_distance(A, dec0_s, dec1_s)        

        # FOV vertices 
        vertex_1_ra, vertex_1_dec = (ra_center + offset_n), (dec_center + fov_height / 2.0)   
        vertex_2_ra, vertex_2_dec = (ra_center - offset_n), (dec_center + fov_height / 2.0)
        vertex_3_ra, vertex_3_dec = (ra_center + offset_s), (dec_center - fov_height / 2.0)  
        vertex_4_ra, vertex_4_dec = (ra_center - offset_s), (dec_center - fov_height / 2.0)

        return vertex_1_ra, vertex_1_dec, vertex_2_ra, vertex_2_dec, vertex_3_ra, vertex_3_dec, vertex_4_ra, vertex_4_dec

    def integrated_fov_probability(self, ra, dec, width, height):
        """Return the probability inside a specific FoV."""

        [v1_ra, v1_dec, v2_ra, v2_dec, v3_ra, v3_dec, v4_ra, v4_dec] = self.vertices(ra, dec, self.fov_width, self.fov_height)
        ra_vertices, dec_vertices = [v1_ra, v2_ra,  v4_ra, v3_ra], [v1_dec, v2_dec, v4_dec, v3_dec]
        prob_fov = self.probability_inside_box(self.skymap, ra_vertices, dec_vertices)
        
        return '%.1e' % prob_fov #, v1_ra, v1_dec, v2_ra, v2_dec, v3_ra, v3_dec, v4_ra, v4_dec 

    def query_box(self, ra, dec, base, height, catalog):
        """Vizier.query_region in a square/rectangular FoV."""

        Vizier.ROW_LIMIT = None
        base_str, height_str = str(base) + 'd', str(height) + 'd'
        
        query_result = Vizier.query_region(SkyCoord(ra = ra, dec = dec, unit = (u.deg, u.deg), frame='icrs'),
                                      width = height_str, height = base_str, catalog = [catalog])
        return query_result

    def query_circle(self, ra, dec, radius, catalog):
        """Vizier.query_region in a circle FoV."""
        
        Vizier.ROW_LIMIT = None
        radius_str = str(radius) + 'd'
        
        query_result = Vizier.query_region(SkyCoord(ra = ra, dec = dec, unit = (u.deg, u.deg), frame='icrs'),
                                           radius = radius_str, catalog = [catalog])
        return query_result

    def fov_stats(self, ra, dec, table, integrated_prob, distance_grid, ansatz_distribution):    
        """Managing the values in the window " FoV statistics:
            1) airmass in step of 1h +
            2) query table --> table
            3) prob inside FoV --> prob_fov
            4) FoV coord center --> ra; dec
            5) Conditional Distance Distribution Along a Line of Sight
            """
        
        airmass_values, datestrings = self.airmass_step(ra, dec) # airmass calculation in step of 1h

        same = self.are_all_same(airmass_values) # check if all elements 
                                                  #--> are "nan"
        if same==True:
            tkMessageBox.showerror('Warning',"airmass outside the range of {1 - 5.8}")
            aladin.remove_FoV(ra, dec) # remove FoV center and query
        else:
            time_step = [dateutil.parser.parse(s) for s in datestrings]
            fov_statistics = FoVstatistics()                      # init class for plotting
            fov_statistics.plot_stats(time_step, airmass_values,
                                      ra, dec,
                                      table,
                                      integrated_prob,
                                      distance_grid, ansatz_distribution)

 
    def _main_direction_funs(self, ra, dec):
        """set of functions related to the cardinal directions: north, south, east, west."""

        catalog = "VII/275/glade1" # Vizier catalog--> You can change the default catalog
                                                # read note at the line 700

        if self.GWsky_basic != 'b':
      
            aladin.get_FoV(ra, dec) # get FoV
            
            aladin.draw_string(ra, dec,                     # show FoV center
                                str( ra ) + "/" + str( dec )) #--> in Aladin plane

            if self.fov_shape != 2: # query
                query_result = self.query_box(
                    ra, dec, self.fov_width, self.fov_height, catalog)  # box             
            else:
                query_result = self.query_circle(
                    ra, dec, self.fov_radius, catalog)                  # circle
            
            if self.fov_shape !=2: # integrated FoV probability
                prob_fov = self.integrated_fov_probability(
                    ra, dec, self.fov_width, self.fov_height)           # box
            else:
                prob_fov = self.probability_inside_circle(
                    ra, dec, self.fov_radius)                           # circle
                                                                            
            self.fov_center_separation(self.input_ra, self.input_dec, # Distance between 2
                                        ra, dec)                         #--> consecutive FoV centers

            r, dp_dr = self.conditional_distance_fov_center(ra, dec) # Conditional distance
                                                                    #--> distribution along FoV

            self.fov_stats(ra, dec, query_result, prob_fov, r, dp_dr) #  fov stats: "FoV statistics"

            self.updating_center_coordinates(ra, dec) # cycle variables
            
        else:
            
            aladin.get_FoV(ra, dec) # get FoV
            
            aladin.draw_string(ra, dec,                     # show FoV center
                                str( ra ) + "/" + str( dec )) #--> in Aladin plane
            
            self.fov_center_separation(self.input_ra, self.input_dec, # Distance between 2
                                ra, dec)                         #--> consecutive FoV centers
            
            self.updating_center_coordinates(ra, dec) # cycle variables
            
                           
    def north(self, shift_up_down=0, shift_right_left=0):
        """Showing the FoV coverage in North direction from an input sky position.""" 

        try:
            if shift_right_left > 0:
                shift_east_west = self._ra0ra1_distance((shift_right_left - self.SHIFT_CORRECTION),
                                                        (self.input_dec + self.fov_height + shift_up_down),
                                                        (self.input_dec + self.fov_height + shift_up_down))
                dist = self.input_ra + shift_east_west 
                
            elif shift_right_left < 0 :
                shift_east_west = self._ra0ra1_distance((shift_right_left + self.SHIFT_CORRECTION),
                                                        (self.input_dec + self.fov_height + shift_up_down),
                                                        (self.input_dec + self.fov_height + shift_up_down))                                                   
                dist = self.input_ra - shift_east_west 
                
            else:
                dist = self.input_ra
                                                         
            north_pointing = [(dist),
                              (self.input_dec + self.fov_height + shift_up_down)]
            
            ra, dec = north_pointing[0], north_pointing[1]        
            self._main_direction_funs(ra, dec)
            
        except ValueError as value_error:
            tkMessageBox.showerror ('Error', value_error)   
        
    def south(self, shift_up_down=0, shift_right_left=0):    
        """Showing the FoV coverage in South direction from an input sky position"""

        try:
            if shift_right_left > 0:
                shift_east_west = self._ra0ra1_distance((shift_right_left - self.SHIFT_CORRECTION),
                                                        (self.input_dec - self.fov_height - shift_up_down),
                                                        (self.input_dec - self.fov_height - shift_up_down))
                dist = self.input_ra + shift_east_west 
                
            elif shift_right_left < 0 :
                shift_east_west = self._ra0ra1_distance((shift_right_left + self.SHIFT_CORRECTION),
                                                        (self.input_dec - self.fov_height - shift_up_down),
                                                        (self.input_dec - self.fov_height - shift_up_down))                                                   
                dist = self.input_ra - shift_east_west 
                
            else:
                dist = self.input_ra
           
            south_pointing = [(dist),
                              (self.input_dec - self.fov_height - shift_up_down)]
            
            ra, dec = south_pointing[0], south_pointing[1]
            self._main_direction_funs(ra, dec)
            
        except ValueError as value_error:
            tkMessageBox.showerror ('Error',value_error)
          
    def east(self, shift_up_down=0, shift_right_left=0):
        """Showing the FoV coverage in East direction from an input sky position.
            A shift correction of 0.00001 is added to escape math error."""

        ra_distance = self._ra0ra1_distance((self.fov_width - self.SHIFT_CORRECTION + shift_right_left),
                                            self.input_dec, self.input_dec)      

        east_pointing = [(self.input_ra + ra_distance), (self.input_dec + shift_up_down)]      
        ra, dec = east_pointing[0], east_pointing[1]
        
        self._main_direction_funs(ra, dec)
        
    def west(self, shift_up_down=0, shift_right_left=0):
        """Showing the FoV coverage in West direction from an input sky position.
            A shift correction of 0.00001 is added to escape math error."""
        
        ra_distance = self._ra0ra1_distance((self.fov_width - self.SHIFT_CORRECTION + shift_right_left),
                                            self.input_dec, self.input_dec)
            
        west_pointing = [(self.input_ra - ra_distance), (self.input_dec + shift_up_down)]     
        ra, dec = west_pointing[0], west_pointing[1]
        
        self._main_direction_funs(ra, dec)
        

class ShowSkyCoverageGUI():
    """
        The GUI consists of 9 buttons; 4  directional buttons to move the FoVs
        in cardinal directions (North, South, East, West) 4  buttons to shift the FoVs
        from a consecutive cardinal direction (↕, ↕, ↔, ↔) and 1 button to get a
        new FoV position (Start FoV). ***Input values in deg***.
    """
    
    def __init__(self, tkMainWin):
        frame = Frame(tkMainWin, border=9, bg="dim grey")
        frame.pack()
        bold8 = tkFont.Font(size=8, weight='bold')        

        self.B02 = Button(frame,text="N", font=bold8)   # Btn north
        self.B02.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B02.grid(row=0, column=2)
        
        self.B12 = Button(frame,text="↕↔", fg="grey")   # Btn ↕↔
        self.B12.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B12.grid(row=1, column=2)
  
        self.B30 = Button(frame, text="E", font=bold8)  # Btn East
        self.B30.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B30.grid(row=3,column=0)

        self.B31 = Button(frame,text="↕↔", fg="grey")   # Btn ↕↔
        self.B31.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B31.grid(row=3, column=1)

        self.B32 = Button(frame, text="Start FoV", fg="red4", # Btn Start FoV
                          font=bold8)
        self.B32.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B32.grid(row=3, column=2)

        self.B33 = Button(frame,text="↕↔", fg="grey")   # Btn ↕↔
        self.B33.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B33.grid(row=3, column=3)

        self.B34 = Button(frame, text="W", font=bold8) # Btn West
        self.B34.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B34.grid(row=3, column=4)

        self.B42 = Button(frame,text="↕↔", fg ="grey") # Btn ↕↔
        self.B42.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B42.grid(row=4, column=2)
        
        self.B52 = Button(frame,text="S", font=bold8) # Btn South
        self.B52.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B52.grid(row=5, column=2)

    # Actions    
    def clicked(self, event):
        """Moving in the cardinal and intercardinal directions calling the
           ShowSkyCoverage class."""
        
        run_sequence = ShowSkyCoverage() # ShowSkyCoverage instance
        
        if event.widget == self.B02:
            run_sequence.north()         # north

        if event.widget == self.B12:
            move_fov = ShiftFoV()
            move_fov._north_shift()      # ↕↔
                      
        if event.widget == self.B30:
            run_sequence.east()          # east

        if event.widget == self.B31:
            move_fov = ShiftFoV()
            move_fov._east_shift()       # ↕↔
            
        if event.widget == self.B32: 
            new_starting_fov = StartingFoV() # start FoV
            new_starting_fov

        if event.widget == self.B33:
            move_fov = ShiftFoV()
            move_fov._west_shift()       # ↕↔
            
        if event.widget == self.B34:   
            run_sequence.west()          # west
            
        if event.widget == self.B42:
            move_fov = ShiftFoV()
            move_fov._south_shift()      # ↕↔
            
        if event.widget == self.B52:    
            run_sequence.south()         # south
    

class StartingFoV(Toplevel):
    """Starting a new FoV sequence from an input FoV.
        The window contains 2 keyboard entries: **the sky coordinated in deg**
        and 3 Buttons:
        Show : draw the FoV in Aladin Plane
        No show : not draw the FoV in Aladin Plane
        Close : close the widget"""
    
    def __init__(self):
        Toplevel.__init__(self, border=8, bg="slate grey")

        self.wait_visibility()
        self.wm_attributes('-alpha',0.7) # semi-trasparent windows

        self.title(" Starting FoV")
        self.attributes("-topmost", True)
        
        self.label_1 = Label(self, text="RA (°); DEC (°)", bg="slate grey")
        self.label_1.grid(row=0, sticky=E)
        
        self.entry_1 = Entry(self, width=8)
        self.entry_2 = Entry(self, width=8)

        self.entry_1.grid(row=0, column=1)
        self.entry_2.grid(row=0, column=2)        

        self.checkbox = Button(self, text="Not show",      
                               command=self.not_show_starting_fov)
        self.checkbox.grid(column=1,row=2)                         # Btn Not show

        self.show = Button(self, text='Show',
                           command=self.show_starting_fov)
        self.show.grid(column=2, row=2)                            # Btn Show
        
        self.close = Button(self, text="Close",
                            command=self.close_window)
        self.close.grid(column=3,row=2)                            # Btn Close

    # Actions                                           
    def show_starting_fov(self):
        """Draw the FoV in Aladin plane"""
        
        show_sky_coverage=ShowSkyCoverage() # init class
        try:
            ra_starting, dec_starting = float(self.entry_1.get()), float(self.entry_2.get())

            show_sky_coverage._main_direction_funs(ra_starting, dec_starting)

        except ValueError as value_error:
            tkMessageBox.showerror ('Error',value_error)
        
    def not_show_starting_fov(self):
        """Not Draw the FoV in Aladin plane"""
        
        try:
            ra_starting, dec_starting = float(self.entry_1.get()), float(self.entry_2.get())
            ShowSkyCoverage.updating_center_coordinates(ra_starting, dec_starting)
        except ValueError as value_error:
            tkMessageBox.showerror ('Error',value_error)
            
    def close_window(self):
        self.destroy()

        
class ShiftFoV(Toplevel):
    """Buttons to shift the FoVs from a consecutive cardinal direction (↕, ↕, ↔, ↔);
       one for each direction. The widget contains 2 keyboard entries and 2 Buttons:
       Ok and Close.
    """
    
    def __init__(self):
        Toplevel.__init__(self, border=7, bg="slate grey")
        self.attributes("-topmost", True)
        self.wait_visibility()
        self.wm_attributes('-alpha', 0.8)     
        
        self.label_3 = Label(self, text="↕ (°)",bg="slate grey")  # UpDown Btn
        self.entry_3 = Entry(self, width=6)
        self.label_3.grid(row=0, sticky=E) 
        self.entry_3.grid(row=0, column=1)

        self.label_4 = Label(self, text="↔ (°)",bg="slate grey") # RightLeft Btn
        self.entry_4 = Entry(self, width=6)
        self.label_4.grid(row=0,column=3) 
        self.entry_4.grid(row=0, column=4)

        self.close = Button(self, text="Close",
                            command = self.close_window) # close Btn
        self.close.grid(column=4,row=2)
    
    def _north_shift(self):
        self.title(" Shifting - North")
        self.checkbox = Button(self, text="OK",
                               command=self.shift_calc_north)
        self.checkbox.grid(column=3,row=2)

    def _south_shift(self):
        self.title(" Shifting - South")
        self.checkbox = Button(self, text="OK",
                               command=self.shift_calc_south)
        self.checkbox.grid(column=3,row=2)

    def _east_shift(self):
        self.title(" Shifting - East")
        self.checkbox = Button(self, text="OK",
                               command=self.shift_calc_east)
        self.checkbox.grid(column=3,row=2)

    def _west_shift(self):
        self.title(" Shifting - West")
        self.checkbox = Button(self, text="OK",
                               command=self.shift_calc_west)
        self.checkbox.grid(column=3,row=2)
        
    # Actions
    def shift_calc_north(self):             # north shift
        try:
            shift_up_down = float(self.entry_3.get())
            shift_right_left = float(self.entry_4.get())
            shift = ShowSkyCoverage()
            shift.north(shift_up_down, shift_right_left)
        except ValueError as value_error:
            tkMessageBox.showerror ('Error',value_error)

    def shift_calc_south(self):             # south shift
        try:
            shift_up_down = float(self.entry_3.get())
            shift_right_left = float(self.entry_4.get())
            shift = ShowSkyCoverage()
            shift.south(shift_up_down, shift_right_left)
        except ValueError as value_error:
            tkMessageBox.showerror ('Error',value_error)

    def shift_calc_east(self):              # east shift
        try:
            shift_up_down = float(self.entry_3.get())
            shift_right_left = float(self.entry_4.get())
            shift = ShowSkyCoverage()
            shift.east(shift_up_down, shift_right_left)
        except ValueError as value_error:
            tkMessageBox.showerror ('Error',value_error)

    def shift_calc_west(self):              # west shift
        try:
            shift_up_down = float(self.entry_3.get())
            shift_right_left = float(self.entry_4.get())
            shift = ShowSkyCoverage()
            shift.west(shift_up_down, shift_right_left)
        except ValueError as value_error:
            tkMessageBox.showerror ('Error',value_error)

    def close_window(self):
        self.destroy()


class FoVstatistics(Toplevel):
    """ FoV statistics window: airmass in step of 1h, integrated probability,
        distance between 2 consecutive FoV centers (shown in shell),
        number of galaxies (distance and B mag).
    """
    
    def __init__(self):
        Toplevel.__init__(self)
        self.title("FoV statistics")
        self.attributes("-topmost", True)
              
        self.B00 = Button(self,text="Close and Retain") # Retain Btn
        self.B00.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B00.pack(side=TOP,fill=BOTH)

        self.B01 = Button(self,text="Delate the last drawn FoV") # Delate Btn
        self.B01.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B01.pack(side=TOP,fill=BOTH)
     
    # Actions
    def clicked(self, event):
        """Retain or Delate Button to retain or delate the FoV."""
        
        if event.widget == self.B00:
            """Retain and Close the last FoV."""

            self.destroy()

        if event.widget == self.B01:
            """Delate the last drawn FoV."""

            ssc=ShowSkyCoverage() # init class
            aladin.remove_FoV(ssc.input_ra, ssc.input_dec)
            aladin.remove("Q:"+str(ssc.input_ra)+"/"+str(ssc.input_dec)) # delete FoV and query
            self.destroy()
        
    def plot_stats(self, time_step, airmass_values, ra, dec,
                   table_q, prob_fov, r, dp_dr):
        """Showing the values in the window FoV statistics."""
        
        f = Figure(figsize=(9, 8),
                   facecolor='white')
        
        def airmass_subplot():
            """SubPlot Airmass."""
            
            airmass_plt = f.add_subplot(223)   
        
            ax=plt.gca()
            ax.xaxis.set_major_formatter(
                mdates.DateFormatter('%H:%M:%S'))

            airmass_plt.set_ylabel('airmass')
            airmass_plt.invert_yaxis()
            airmass_plt.grid(True)

            airmass_plt.plot(time_step, airmass_values, "o-", linestyle="--")
            f.autofmt_xdate()

        airmass_subplot = airmass_subplot()
        

        #####################################################################
        #                        NOTE                                       #
        #     The default catalog is GLADE [VII/275/glade1]                 #
        #                                                                   #
        #   You can change the default catalog at the line 312              #
        # e.g from GLADE to GWGC --> "VII/275/glade1" --> "VII/267"         #
        #                                                                   #
        #      You also can change the catalog statistics                   #
        #           query_catalog.hist(table_pandas.Dist.dropna() -->       #
        # query_catalog.hist(table_pandas.**ColumnName**.dropna()           #
        #                                                                   #
        #           query_catalog.hist(table_pandas.Bmag.dropna()) -->      #
        # query_catalog.hist(table_pandas.**ColumnName**.dropna())          #
        #####################################################################


        def subplot_cat_distribution():
            """SubPlot catalogue distributions."""

            try:
                for table_name in table_q.keys():
                    table = table_q[table_name]
                
                table.write("query_items", format = 'votable', overwrite = True)
                samp.send_file("query_items")
                #aladin.rename("Q:"+str(ra)+"/"+str(dec))
                aladin.set_planeID("Q:"+str(ra)+"/"+str(dec))
                aladin.remove('query_items')               # try to do better

                table_pandas = table.to_pandas()
                ra_q = table_pandas._RAJ2000
                dec_q = table_pandas._DEJ2000
        
                query_catalog = f.add_subplot(222)
                query_catalog.hist(table_pandas.Dist.dropna())         # read note above
                query_catalog.set_xlabel('GLADE catalog: Dist [Mpc]')  # read note above
                query_catalog.set_ylabel('Count')                      # read note above
                
                query_catalog = f.add_subplot(224) 
                query_catalog.hist(table_pandas.Bmag.dropna())                   # read note above
                query_catalog.set_xlabel('GLADE catalog: apparent B magnitude')  # read note above
                query_catalog.set_ylabel('Count')                                # read note above
            except UnboundLocalError as unbound_local_error:
                print ('No Galaxies in the selected FoV')              
            finally:
                pass
            
        subplot_cat_distribution = subplot_cat_distribution()                                   

        def subplot_cond_distance():
            """Conditional Distance Distribution Along a Line of Sight."""
            
            conditional_distance_line_sight = f.add_subplot(221) 
            conditional_distance_line_sight.plot(r, dp_dr)
            
            title_string = 'Conditional distance distribution \n along the FoV center'
            conditional_distance_line_sight.set_title(title_string,fontsize=10)
            conditional_distance_line_sight.set_xlabel('distance (Mpc)')
            conditional_distance_line_sight.set_ylabel('prob Mpc$^{-1}$')

        subplot_cond_distance = subplot_cond_distance()

        def draw_area():
            """Drawing Area: FoV Statistics."""

            fov_information_title = "FoV center (ra "+str(ra) + "  " +"dec "+ str(dec)+")" + "; " + "prob: " + str(prob_fov)
            f.suptitle(fov_information_title, fontsize=14)
           
            canvas = FigureCanvasTkAgg(f, self) 
            canvas.show()
            canvas.get_tk_widget().pack(side=BOTTOM, fill=BOTH, expand=True)

            toolbar = NavigationToolbar2TkAgg(canvas, self)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=True)

        draw_area = draw_area()
        
# running
mainWin = Tk()

sscGUI = ShowSkyCoverageGUI(mainWin)
mainWin.title('GWsky')
mainWin.attributes("-topmost", True)

mainWin.wait_visibility(mainWin)
mainWin.wm_attributes('-alpha', 0.8) #semi-trasparent windows

mainWin.mainloop()
