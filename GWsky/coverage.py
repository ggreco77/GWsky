# -*- coding: utf-8 -*-

try:
   import cPickle as pickle
except:
   import pickle

import fileinput
import os
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

from math import cos, sin, acos, asin, atan, degrees, radians

import astropy
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import TimeDelta
from astroquery.vizier import Vizier 

import numpy as np
import healpy as hp
from scipy.stats import norm

from aladinSAMP import AladinViaSAMP, AladinScriptCommands 
samp = AladinViaSAMP()
aladin = AladinScriptCommands()

from utils import Utils
Utils.create_folders(folders=["Queries", "Coords", "FoV"])
Utils.load_user_fov("GWsky_fov.vot")

class UserValues(object):
    """
       Creating a class in which the instance attributes are based on the dictionary
       "GWsky_config" by default. GWsky_config is created and pickled by the "user_values"
       module and will be deleted when the "ShowSkyCoverageGUI" will be closed.
    """
    
    def __init__(self, infile_config="GWsky_config"):      
        """
        "GWsky_config"  contains the keys below:
            self.skymap: a valid LVC skymap in healpix format;
            self.nside:  Resolution parameter of HEALPIX 
            Geodetic coordinates of the Observatory:
                     self.latitude:   latitude[deg];
                     self.longitude:  longitude[deg];
                     self.altitude:   altitude[m];
            self.catalog: catalog name in Vizier code
            self.obs_time: starting time [yyyy-mm-dd hh:mm:ss];
            self.fov_width: field-of-view width [deg];
            self.fov_height: field-of-view height [deg];
            self.fov_radius: field-of-view radius [deg];
            self.ra_max: right ascension of maximum probability pixel [deg];
            self.dec_max: declination of maximum probability pixel [deg];
            self.GWsky_basic: ("A") active statistic window;
                            : ("D") deactive statistic window;
            self.column_1: first catalog column selected by user;
            self.column_2: second catalog column selected by user; 
            self.filter_1: applied filter in colomn 1;
            self.filter_2: applied filter in colomn 2.               
        """

        self.infile_config =  infile_config
        
        with open(infile_config, 'rb') as data:
            config_GWsky = pickle.load(data)

        for k, v in config_GWsky.items():
            setattr(self, k, v)

    def __repr__(self):
        with open(self.infile_config, 'rb') as data:
            config_GWsky = pickle.load(data)
        return str(config_GWsky)

class Airmass(UserValues):
    """Airmass calculation from 1 to 5.8 by default."""
    
    AIRMASS_MIN = 1
    AIRMASS_MAX = 5.8
                  
    def airmass(self, ra, dec, latitude, longitude, height,
                time_input, AIRMASS_MIN, AIRMASS_MAX):
        """ Airmass calculation at a given time in a particular site for an input sky position (ra, dec).
            The airmass is calculated in the range [AIRMASS_MIN, AIRMASS_MAX]."""
     
        observatory = astropy.coordinates.EarthLocation(lat = latitude*u.deg,
                                                        lon = longitude*u.deg,
                                                        height = height*u.m)
        sky_coord = SkyCoord(ra = ra*u.deg,
                              dec=dec*u.deg, frame='icrs')     
        time = Time(time_input) 
        altaz = sky_coord.transform_to(AltAz(obstime=time,
                                              location=observatory))                                                                             
        airmass_value = altaz.secz                  

        if airmass_value < AIRMASS_MIN or airmass_value >= AIRMASS_MAX:
             airmass_value =  "nan" 
        else:     
             airmass_value = round(airmass_value, 2) 

        return airmass_value

    def airmass_step(self, ra, dec):
        """Airmass calculation in step of one hour for an input sky position (ra, dec).
           10 steps are performed: step <10; dt = 1h."""

        dt = TimeDelta(3600.0, format='sec')
        self.obs_time = Time(self.obs_time)

        airmass_list = []
        time_list = []
        step = 0
        
        while step < 10:
            time_input = self.obs_time + step*dt
            val = self.airmass(ra, dec, self.latitude, self.longitude, self.altitude,
                               time_input, self.AIRMASS_MIN, self.AIRMASS_MAX)           
            airmass_list.append(val)
            time_list.append(str(time_input))
            step+=1
            
        return airmass_list, time_list


class ShowSkyCoverage(Airmass, UserValues): 
    """Moving the FoV-footprint coverage by choosing a starting pointing."""
    
    SHIFT_CORRECTION = 0.00001  # A shift correction of 0.00001 is added
                                 # --> to escape math error during the FoV sequence

    def __init__(self, infile_coords='GWsky_coords'):
        """Creating a class in which the instance attributes are based on the dictionary
       "GWsky_coords" by default. GWsky_coords is created and pickled by the "user_values"
       module and will be deleted when the "ShowSkyCoverageGUI" will be closed.
       It contains the keys: "ra", "dec". ra and dec represents the central location of a FoV. 
        
        Starting sky coordinates:
             self.input_ra: right ascension [deg]
             self.input_dec: declination [deg]
         """

        self.infile_coords = infile_coords
        # new entries during the FoV sequence
        self.entries_GWsky_new =[]
        
        UserValues.__init__(self)
        Airmass.__init__(self)
        
        with open(infile_coords, 'rb') as data:  
            coords_GWsky = pickle.load(data)
            
        for k, v in coords_GWsky.items():          
            setattr(self, k, v)

    def __repr__(self):
        with open(self.infile_coords, 'rb') as data:
            coords_GWsky = pickle.load(data)
        return str(coords_GWsky)
            
    def __ra0ra1(self, A, dec0, dec1):
        """From the angular distance:
           cos(A) = sin(Dec1)sin(Dec2)+cos(Dec1)cos(Dec2)cos(ra1-ra2) --> 
           cos(ra1-ra2) = [cos(A)-sin(dec0)sin(dec1)]/[cos(dec0)cos(dec1)]."""

        dec0, dec1, A = radians(dec0),  radians(dec1), radians(A)
        cos_ra0_ra1 = ( cos(A)-sin(dec0)*sin(dec1) )/( cos(dec0)*cos(dec1) )
        ra0ra1 = degrees( acos(cos_ra0_ra1) )

        return  round(ra0ra1, 5)         
    
    def __updating_center_coords(self, ra, dec):
        """Getting/Updating FoV-center (ra, dec) in the dict named by default "GWsky_coords".
.          For each tile across the sky the file is updated."""""

        with open('GWsky_coords', 'rb') as data:
            coords_GWsky = pickle.load(data)
            
        coords_GWsky['input_ra'], coords_GWsky ['input_dec'] = ra, dec

        with open('GWsky_coords', 'wb') as data:
            pickle.dump(coords_GWsky, data)
   
    def __TEST_fov_center_separation(self, ra1, dec1, ra2, dec2):
        """TESTING CODE: Distance between 2 consecutive FoV centers [deg]."""
           
        fov_1 = SkyCoord(ra1, dec1, frame='icrs',unit='deg')
        fov_2 = SkyCoord(ra2, dec2, frame='icrs',unit='deg')
        sep = fov_1.separation(fov_2)

        print ('The distance between 2 consecutive FoV centers is', sep.round(6))

    def __are_all_same(self, items):
        """Check if all elements of a list are the same."""
        
        return all(x == items[0] for x in items)

    def __ipix_sum(self, infile, ra_vertices, dec_vertices):
         """Return the ipix sum inside a polygon."""

         prob = hp.read_map(self.skymap, verbose=False)
         
         theta = 0.5 * np.pi - np.deg2rad(dec_vertices)
         phi = np.deg2rad(ra_vertices)
         xyz = hp.ang2vec(theta, phi)

         ipix_poly = hp.query_polygon(self.nside, xyz)
     
         ipix_sum_polygon = prob[ipix_poly].sum()
         
         return ipix_sum_polygon

    def __vertices(self, ra_center, dec_center, fov_base, fov_height):
        """Finding the vertices of a FoV given the central location (ra[deg], dec[deg])
           and the FoV size (FoV_base [deg], FoV_height [deg])."""
        
        vert_ra, vert_dec = [], []  # ra list,  dec list 
        
        ra_center_rad, dec_center_rad = radians(ra_center), radians(dec_center)

        fov_base_rad, fov_height_rad = radians(fov_base), radians(fov_height)
       
        x = [-fov_base_rad/2, fov_base_rad/2,
             fov_base_rad/2, -fov_base_rad/2]

        y = [fov_height_rad/2, fov_height_rad/2,
             -fov_height_rad/2, -fov_height_rad/2]
        
        for i, j  in zip(x, y):
            arg = -i/(cos(dec_center_rad)-j*sin(dec_center_rad))          
            v_ra = degrees( (ra_center_rad+atan(arg)) )
            
            vert_ra.append(v_ra)
            
            v_dec = degrees( (asin((sin(dec_center_rad)+j*cos(dec_center_rad))/(1+i**2+j**2)**0.5)) )
            
            vert_dec.append(v_dec)

        # test: field-of-view footprint vs vertices function
        #aladin.draw_circle(vert_ra[0], vert_dec[0], size = '2arcmin')
        #aladin.draw_circle(vert_ra[1], vert_dec[1], size = '2arcmin')
        #aladin.draw_circle(vert_ra[2], vert_dec[2], size = '2arcmin')
        #aladin.draw_circle(vert_ra[3], vert_dec[3], size = '2arcmin')
                
        return vert_ra[0], vert_ra[1], vert_ra[3], vert_ra[2], vert_dec[0], vert_dec[1], vert_dec[3], vert_dec[2]

    def __prob_in_box(self, ra, dec, width, height):
        """Return the probability inside a box."""

        v1_ra, v2_ra, v3_ra, v4_ra, v1_dec, v2_dec, v3_dec, v4_dec = self.__vertices(
            ra, dec, self.fov_width, self.fov_height)
        
        ra_vertices, dec_vertices = (
           [v1_ra, v2_ra, v4_ra, v3_ra], [v1_dec, v2_dec, v4_dec, v3_dec])
        
        probability_fov_box = self.__ipix_sum(self.skymap, ra_vertices, dec_vertices)
        
        return '%.1e' % probability_fov_box

    def __prob_in_circle(self, ra, dec, radius):
        """Return the probability inside a circle."""

        prob = hp.read_map(self.skymap, verbose=False)

        theta = 0.5 * np.pi - np.deg2rad(dec)
        phi = np.deg2rad(ra)
        radius = np.deg2rad(radius)

        xyz = hp.ang2vec(theta, phi)
        ipix_disc = hp.query_disc(self.nside, xyz, radius)
        probability_fov_disc = prob[ipix_disc].sum()

        return '%.1e' % probability_fov_disc

    def __conditional_distance_fov_center(self, ra, dec):
        """Conditional distance distribution along the line of sight of a FoV center:
            see https://arxiv.org/pdf/1605.04242v3.pdf - section 4.4 for more details."""

        prob, header = hp.read_map(self.skymap, h=True,
                                   verbose=False)
        header = dict(header)

        if header['TFIELDS']==4:  # 3d skymap
            prob, distmu, distsigma, distnorm = hp.read_map(self.skymap, verbose=False,
                                                            field=range(4))
            theta = 0.5 * np.pi - np.deg2rad(dec)
            phi = np.deg2rad(ra)
            ipix = hp.ang2pix(self.nside, theta, phi)
            
            line_end = header['DISTMEAN'] + (header['DISTSTD']*3)
            r = np.linspace(0, line_end)
            dp_dr = r**2*distnorm[ipix]*norm(distmu[ipix], distsigma[ipix]).pdf(r)
            
            return r, dp_dr
        else:
            r = "nan"
            dp_dr = "nan"

            return r, dp_dr

    def __query_box(self, ra, dec, base, height, catalog):
        """Vizier.query_region in a square/rectangular FoV using astroquery module."""

        base_str, height_str = str(base) + 'd', str(height) + 'd'
        
        catalog_constraints = Vizier(catalog=[self.catalog],
                                     columns=['*', '_RAJ2000', '_DEJ2000', self.column_1, self.column_2],
                                     column_filters={self.column_1 : self.filter_1,
                                                     self.column_2 : self.filter_2})
        
        catalog_constraints.ROW_LIMIT = -1  # completed catalog
        query_result = catalog_constraints.query_region(SkyCoord(
            ra = ra, dec = dec, unit = (u.deg, u.deg), frame='icrs'),
                                                        width = height_str, height = base_str)       
        return query_result

    def __query_circle(self, ra, dec, radius, catalog):
        """Vizier.query_region in a circle FoV  using astroquery module."""
        
        radius_str = str(radius) + 'd'
        
        catalog_constraints = Vizier(catalog=[self.catalog],
                                     columns=['*', '_RAJ2000', '_DEJ2000', self.column_1, self.column_2],
                                     column_filters={self.column_1 : self.filter_1,
                                                     self.column_2 : self.filter_2})
        
        catalog_constraints.ROW_LIMIT = -1  # completed catalog
        query_result = catalog_constraints.query_region(SkyCoord(
            ra = ra, dec = dec, unit = (u.deg, u.deg), frame='icrs'), radius = radius_str)
        
        return query_result

    def __fov_stats(self, ra, dec, table, integrated_prob, distance_grid, ansatz_distribution):    
        """Managing the descriptive statistic window. If the airmass is outside the default range [AIRMASS_MIN, AIRMASS_MAX]
            the window is not opened otherwise the window is shown."""
        
        airmass_values, datestrings = self.airmass_step(ra, dec) 
        same = self.__are_all_same(airmass_values) 
                                                 
        if same == True:
            tkMessageBox.showerror('Warning',"airmass outside the range of {1 - 5.8}")
            aladin.remove_FoV(ra, dec) 
            aladin.remove("C_" + str( ra ) + "/" + str( dec ))  
        else:
            time_step = [dateutil.parser.parse(s) for s in datestrings]
            self.__updating_center_coords(ra,dec)  
            fov_statistics = FoVstatistics()                 
            fov_statistics.plot_stats(time_step, airmass_values,
                                      ra, dec,
                                      table,
                                      integrated_prob,
                                      distance_grid, ansatz_distribution)
    
    def update_pointings_file(self, infile, ra, dec, prob_fov):
         """The central location (ra[deg], dec[deg]) and the integrated probability of
             a selected FoV are saved locally in an external file.
             By default the file is named "GWsky_pointings.txt"."""
           
         with open(infile, 'a') as pointing:
             pointing.write(str(ra) + ' ' + str(dec)+ ' '+ str(prob_fov)+'\n')

    def __query_shape(self, ra, dec, fov_shape):
        """Return the catalog query according the defined-user FoV shape:
                   (1) box and (2) circle. """
        
        if self.fov_shape != 2:  # box FoV
                    query_result = self.__query_box(
                       ra, dec, self.fov_width, self.fov_height, self.catalog)               
        else: # circle FoV
            query_result = self.__query_circle(
               ra, dec, self.fov_radius, self.catalog)
            
        return query_result

    def __prob_shape(self, ra, dec, fov_shape):
        """Return the integrated probability according the defined-user FoV shape:
                   (1) box and (2) circle."""
        
        if self.fov_shape !=2:  # box FoV
            prob_fov = self.__prob_in_box(
               ra, dec, self.fov_width, self.fov_height)                
        else: # circle FoV
            prob_fov = self.__prob_in_circle(
               ra, dec, self.fov_radius)
            
        return  prob_fov

    def __aladin_stack(self, ra, dec):
        """Organizing aladin stack if a FoV footprint is generated."""
        
        aladin.get_FoV(ra, dec)
        aladin.draw_newtool("C_" + str( ra ) + "/" + str( dec ))
        aladin.draw_string(ra, dec,
                           str( ra ) + "/" + str( dec ))
                  
    def pick_coverage(self, ra, dec):
        """Setting GWsky: with statistic window (A); without statistic window (D)."""                                      

        if self.GWsky_basic == "A":  # full version
           
            self.__aladin_stack(ra,dec)

            query_result = self.__query_shape(ra, dec, self.fov_shape)
            prob_fov = self.__prob_shape(ra, dec, self.fov_shape)

            # TEST   
            self.__TEST_fov_center_separation(self.input_ra, self.input_dec, 
                                        ra, dec)                              
            
            r, dp_dr = self.__conditional_distance_fov_center(ra, dec) 
            self.__fov_stats(ra, dec, query_result, prob_fov, r, dp_dr) 
            
        elif self.GWsky_basic == "D":  # basic version-> no Stats win
           
            self.__aladin_stack(ra,dec)

            Utils.move_to_folder(planes=['C_*','P:*'],
                                   folders=['Coords','FoV'])

            #query_result = self.__query_shape(ra, dec, self.fov_shape)
            prob_fov = self.__prob_shape(ra, dec, self.fov_shape)


            # TEST  
            #self.__TEST_fov_center_separation(self.input_ra, self.input_dec,
            #                    ra, dec)                                       

            self.update_pointings_file("GWsky_pointings.txt", ra, dec, prob_fov)
            
    def __intercardinal_distance(self, ra_pos, dec_pos, shift_up_down, shift_right_left):
        """Moving from the fixed cardinal direction using the bi-directional windows;
           shift_up_down ↕ and/or shift_right_left ↔."""

        if shift_right_left > 0:
           shift_east_west = self.__ra0ra1((shift_right_left-self.SHIFT_CORRECTION),
                                                   (dec_pos + self.fov_height + shift_up_down),
                                                   (dec_pos + self.fov_height + shift_up_down))
           dist = ra_pos + shift_east_west 
         
        elif shift_right_left < 0 :
           shift_east_west = self.__ra0ra1((shift_right_left + self.SHIFT_CORRECTION),
                                                   (dec_pos + self.fov_height + shift_up_down),
                                                   (dec_pos + self.fov_height + shift_up_down))
           dist = ra_pos - shift_east_west
         
        else:
           dist = ra_pos

        return dist 

    def __load_entries(self, infile_entries):
        """Opening the file in which the input entries are stored: ra_1 dec_1 ra_2 dec_2 ... ra_n dec_n
           By default the file is named "GWsky_entries". "GWsky_entries" is created from the
           "show_starting_fov" method in "StartingFoV" class.
           A error message invites users to press the "Start FoV" button before carrying out any"""
       
        try:
            with open(infile_entries, 'rb') as data:
                entries_GWsky = pickle.load(data)
                return entries_GWsky
        except IOError as io_error:
            message = "Press the button 'Start FoV' and choose a starting position; \
                        by default the sky position of the max probability pixel is shown"
            tkMessageBox.showerror ('Error', message)
                         
    def north(self, shift_up_down=0, shift_right_left=0):
        """Moving the FoV tiles in North direction from input sky-position(s).
           The related bidirectional button permits small shifts from such pre-defined direction:
           shift_up_down ↕ and/or shift_right_left ↔ """

        entries_GWsky = self.__load_entries("GWsky_entries")
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]  

        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
            dist = self.__intercardinal_distance(float(ra_start), float(dec_start),
                                                 shift_up_down, shift_right_left)
            north_pointing = [(dist),
                               (float(dec_start) + self.fov_height + shift_up_down)]
             
            ra, dec = north_pointing[0], north_pointing[1]
            self.pick_coverage(ra, dec)

            # cycle variables
            new_sky_pos = [ra,dec]
            self.entries_GWsky_new.extend(new_sky_pos)

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)
            
    def south(self, shift_up_down=0, shift_right_left=0):    
        """Moving the FoV tiles in South direction from input sky-position(s).
           The related bidirectional button permits small shifts from such pre-defined direction:
           shift_up_down ↕ and/or shift_right_left ↔"""

        entries_GWsky = self.__load_entries("GWsky_entries")
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]

        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
            dist = self.__intercardinal_distance(float(ra_start), float(dec_start),
                                                 shift_up_down, shift_right_left)
            south_pointing = [(dist), (float(dec_start) - self.fov_height - shift_up_down)]
                    
            ra, dec = south_pointing[0], south_pointing[1]
            self.pick_coverage(ra, dec)
                    
            # cycle variables
            new_sky_pos = [ra,dec]
            self.entries_GWsky_new.extend(new_sky_pos)

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)    
          
    def east(self, shift_up_down=0, shift_right_left=0):
        """Moving the FoV tiles in East direction  from input sky-position(s).
           The related bidirectional button permits small shifts from such pre-defined direction:
           shift_up_down ↕ and/or shift_right_left ↔.
           A shift correction of 0.00001 is added to escape math error."""
        
        entries_GWsky = self.__load_entries("GWsky_entries")
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]

        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):              
            ra_distance = self.__ra0ra1((self.fov_width - self.SHIFT_CORRECTION + shift_right_left),
                                        float(dec_start), float(dec_start))
                
            east_pointing = [(float(ra_start) + ra_distance), (float(dec_start) + shift_up_down)]
            ra, dec = east_pointing[0], east_pointing[1]

            self.pick_coverage(ra, dec)

            # cycle variables
            new_sky_pos = [ra,dec]
            self.entries_GWsky_new.extend(new_sky_pos)

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)
                   
    def west(self, shift_up_down=0, shift_right_left=0):
        """Moving the FoV tiles in West direction  from input sky-position(s).
           The related bidirectional button permits small shifts from such pre-defined direction:
           shift_up_down ↕ and/or shift_right_left ↔.
            A shift correction of 0.00001 is added to escape math error."""
        
        entries_GWsky = self.__load_entries("GWsky_entries")
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]  
      
        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
               
            ra_distance = self.__ra0ra1((self.fov_width - self.SHIFT_CORRECTION + shift_right_left),float(dec_start), float(dec_start))

            west_pointing = [(float(ra_start) - ra_distance), (float(dec_start) + shift_up_down)]
            ra, dec = west_pointing[0], west_pointing[1] 
        
            self.pick_coverage(ra, dec)

            # cycle variables
            new_sky_pos = [ra,dec]
            self.entries_GWsky_new.extend(new_sky_pos)

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)
        

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

        self.B02 = Button(frame,text="N", font=bold8)   
        self.B02.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B02.grid(row=0, column=2)
        
        self.B12 = Button(frame,text="↕↔", fg="grey")  
        self.B12.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B12.grid(row=1, column=2)
  
        self.B30 = Button(frame, text="E", font=bold8)  
        self.B30.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B30.grid(row=3,column=0)

        self.B31 = Button(frame,text="↕↔", fg="grey")   
        self.B31.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B31.grid(row=3, column=1)

        self.B32 = Button(frame, text="Start FoV", fg="red4", 
                          font=bold8)
        self.B32.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B32.grid(row=3, column=2)

        self.B33 = Button(frame,text="↕↔", fg="grey")  
        self.B33.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B33.grid(row=3, column=3)

        self.B34 = Button(frame, text="W", font=bold8) 
        self.B34.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B34.grid(row=3, column=4)

        self.B42 = Button(frame,text="↕↔", fg ="grey") 
        self.B42.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B42.grid(row=4, column=2)
        
        self.B52 = Button(frame,text="S", font=bold8) 
        self.B52.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B52.grid(row=5, column=2)

    # Actions    
    def clicked(self, event):
        """Moving the user-defined FoV footprint."""
        
        run_sequence = ShowSkyCoverage()
        
        if event.widget == self.B02:
            run_sequence.north()         # north

        if event.widget == self.B12:
            move_fov = ShiftFoV()
            move_fov.north_shift()       # ↕↔
                      
        if event.widget == self.B30:
            run_sequence.east()          # east

        if event.widget == self.B31:
            move_fov = ShiftFoV()
            move_fov.east_shift()        # ↕↔
            
        if event.widget == self.B32: 
            new_starting_fov = StartingFoV()  # start FoV
            new_starting_fov

        if event.widget == self.B33:
            move_fov = ShiftFoV()
            move_fov.west_shift()        # ↕↔
            
        if event.widget == self.B34:   
            run_sequence.west()          # west
            
        if event.widget == self.B42:
            move_fov = ShiftFoV()
            move_fov.south_shift()       # ↕↔
            
        if event.widget == self.B52:    
            run_sequence.south()         # south
    

class StartingFoV(Toplevel, UserValues):
    """Starting a sequence from a list of FoV(s). The window contains 1 keyboard entries and 3 Buttons.

        entry:
             ra_1 dec_1 ra_2 dec_2 ra_3 dec_3 ... ra_n dec_n [deg]
             By default: sky coords of maximum probability pixel

        Btns:
           Show : draw user-defined FoV footprint(s) in Aladin Plane(s)
           No show : no draw user-defined FoV footprint(s) in Aladin Plane(s)
           Close : close the widget
        """
    
    def __init__(self):
        Toplevel.__init__(self, border=8, bg="slate grey")
        UserValues.__init__(self)

        # putting the entry value(s) in a list
        self.entries_GWsky=[]  

        self.wait_visibility()
        self.wm_attributes('-alpha',0.7) # transparency

        self.title(" Starting FoV")
        self.attributes("-topmost", True)
        
        self.label_1 = Label(self, text="RA (°) DEC (°)", bg="slate grey")
        self.label_1.grid(row=0, sticky=E, padx=5)

        # default: sky coords of maximum probability pixel
        fov_coords = str(self.ra_max_pixel), str(self.dec_max_pixel) 
        
        max_pixel_default = StringVar(self, value=fov_coords) 
        self.entry_1 = Entry(self, width=22, justify=CENTER,
                             textvariable=max_pixel_default)

        self.entry_1.grid(row=0, column=1)

        self.entryScroll = Scrollbar(self, orient=HORIZONTAL,
                                     command=self.__scrollHandler)
        self.entryScroll.grid(row=1, sticky=E+W)
        self.entry_1['xscrollcommand'] = self.entryScroll.set

        #Btn
        self.checkbox = Button(self, text="Not show",      
                               command=self.no_show_starting_fov)
        self.checkbox.grid(column=2,row=2, sticky=E, padx=2, pady=5)  

        self.show = Button(self, text='Show',
                           command=self.show_starting_fov)
        self.show.grid(column=1, row=2, sticky=E, padx=2, pady=5) 
        
        self.close = Button(self, text="Close",
                            command=self.close_window)  
        self.close.grid(column=3,row=2, sticky=E, padx=2, pady=5)  

    #Actions
    def __scrollHandler(self, *L):
        """Scroll entry in starting FoV window."""
       
        op, howMany = L[0], L[1]

        if op == 'scroll':
            units = L[2]
            self.entry_1.xview_scroll(howMany, units)
        elif op == 'moveto':
            self.entry_1.xview_moveto(howMany)
            
    def __split_entries(self):
        """Splitting the entries in ra and dec; # odd: ra and # even: dec."""
        
        current_fov_coords = self.entry_1.get().replace(';',' ').replace('-',' ').replace(',',' ').split()
        fov_center_ra = current_fov_coords[0::2]
        fov_center_dec = current_fov_coords[1::2]

        return current_fov_coords, fov_center_ra, fov_center_dec
    
    def show_starting_fov(self):
        """Drawing the FoV footprint(s) in the Aladin plane(s).
         By default: sky coords (ra[deg], dec[deg]) of maximum probability pixel."""   
   
        show_sky_coverage = ShowSkyCoverage()
        
        current_fov_coords, fov_center_ra, fov_center_dec = self.__split_entries()
        
        try:
            for ra_starting, dec_starting in zip (fov_center_ra, fov_center_dec):          
                show_sky_coverage.pick_coverage(float(ra_starting), float(dec_starting))
        except ValueError as value_error:
            tkMessageBox.showerror ('Error', value_error)
          
        self.entries_GWsky.extend(current_fov_coords)

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky, data)
            
        self.entries_GWsky=[] # re-init.
        
    def no_show_starting_fov(self):
        """No Draw the FoV footprint(s) in the Aladin plane(s);
               useful to re-initialize the sequence."""

        self.entries_GWsky=[] # re-init.     
        current_fov_coords, fov_center_ra, fov_center_dec = self.__split_entries()

        self.entries_GWsky.extend(current_fov_coords)

        with open('GWsky_entries', 'wb') as data:
            return pickle.dump(self.entries_GWsky, data)
     
    def close_window(self):
        self.destroy()

        
class ShiftFoV(Toplevel):
    """Shifting the FoV footprint(s) from a consecutive cardinal direction (↕, ↕, ↔, ↔);
       The widget contains 2 entries and 2 Buttons."""
    
    def __init__(self):
        Toplevel.__init__(self, border=7, bg="slate grey")
        self.attributes("-topmost", True)
        self.wait_visibility()
        self.wm_attributes('-alpha', 0.8)     
        
        self.label_3 = Label(self, text="↕ (°)",bg="slate grey")  
        self.entry_3 = Entry(self, width=6, justify=CENTER)
        self.label_3.grid(row=0, sticky=E) 
        self.entry_3.grid(row=0, column=1)

        self.label_4 = Label(self, text="↔ (°)",bg="slate grey")  
        self.entry_4 = Entry(self, width=6, justify=CENTER)
        self.label_4.grid(row=0,column=3) 
        self.entry_4.grid(row=0, column=4)

        self.close = Button(self, text="Close",
                            command = self.close_window)  
        self.close.grid(column=4,row=2)
    
    def north_shift(self):
        self.title(" Shifting - North")
        self.checkbox = Button(self, text="OK",
                               command=self.shift_north)
        self.checkbox.grid(column=3,row=2)

    def south_shift(self):
        self.title(" Shifting - South")
        self.checkbox = Button(self, text="OK",
                               command=self.shift_south)
        self.checkbox.grid(column=3,row=2)

    def east_shift(self):
        self.title(" Shifting - East")
        self.checkbox = Button(self, text="OK",
                               command=self.shift_east)
        self.checkbox.grid(column=3,row=2)

    def west_shift(self):
        self.title(" Shifting - West")
        self.checkbox = Button(self, text="OK",
                               command=self.shift_west)
        self.checkbox.grid(column=3,row=2)
        
    # Actions
    def shift_north(self):             
        try:
            shift_up_down = float(self.entry_3.get())
            shift_right_left = float(self.entry_4.get())
            shift = ShowSkyCoverage()
            shift.north(shift_up_down, shift_right_left)
        except ValueError as value_error:
            tkMessageBox.showerror ('Error',value_error)

    def shift_south(self):             
        try:
            shift_up_down = float(self.entry_3.get())
            shift_right_left = float(self.entry_4.get())
            shift = ShowSkyCoverage()
            shift.south(shift_up_down, shift_right_left)
        except ValueError as value_error:
            tkMessageBox.showerror ('Error',value_error)

    def shift_east(self):             
        try:
            shift_up_down = float(self.entry_3.get())
            shift_right_left = float(self.entry_4.get())
            shift = ShowSkyCoverage()
            shift.east(shift_up_down, shift_right_left)
        except ValueError as value_error:
            tkMessageBox.showerror ('Error',value_error)

    def shift_west(self):              
        try:
            shift_up_down = float(self.entry_3.get())
            shift_right_left = float(self.entry_4.get())
            shift = ShowSkyCoverage()
            shift.west(shift_up_down, shift_right_left)
        except ValueError as value_error:
            tkMessageBox.showerror ('Error',value_error)

    def close_window(self):
        self.destroy()


class FoVstatistics(Toplevel, ShowSkyCoverage):
    """ FoV statistics window consists of 3 buttons:
         (1) Confirm the FoV-footprint
         (2) Delete FoV-footprint
         (3) Zoom in FoV. """
    
    def __init__(self):
        Toplevel.__init__(self)
        ShowSkyCoverage.__init__(self)

        self.title("FoV statistics")
        self.attributes("-topmost", True)
              
        self.B00 = Button(self,text="Confirm the pointing in the selected_pointings txt file")  
        self.B00.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B00.pack(side=TOP,fill=BOTH)

        self.B01 = Button(self,text="Delete the FoV")  
        self.B01.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B01.pack(side=TOP,fill=BOTH)

        self.B02 = Button(self,text="Zoom in the FoV")  
        self.B02.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B02.pack(side=TOP,fill=BOTH)

        # default
        fov_center_ra_dec = str(self.input_ra), str(self.input_dec), 'on fly notes, NO DELETE OR CHANGE THE COORDS!'       
        current_fov = StringVar(self, value=fov_center_ra_dec)  
        self.entry_current_fov = Entry(self, width=30, justify=LEFT,
                                       textvariable=current_fov)  
        self.entry_current_fov.pack(side=TOP,fill=BOTH)

    def delete_pointing(self, infile ,ra, dec):
        """Deleting input sky coords from an external file; by default "Pointing.txt" """
        
        for line in fileinput.input(infile, inplace=True):
            if (line.rsplit()[0] != str(ra) or line.rsplit()[1] != str(dec)) :
                print line.rstrip('\n')

    def __rm_from_stack(self, ra, dec):
        """Removing from Aladin stack the associated planes."""
        
        aladin.remove_FoV(ra, dec)       
        aladin.remove("Q:"+ ra +"/"+ dec)           
        aladin.remove("C_" + ra+ "/" + dec)

    def clicked(self, event):
        """Retain or delate the FoV-footprint(s).
            The FoV-center positions are saved in "Pointings txt" """      
        
        if event.widget == self.B00: # Retain and Close the FoV.            
            self.destroy()
            
        if event.widget == self.B01: # Delete FoV
            current_fov_coords = self.entry_current_fov.get().split()
            current_fov_ra, current_fov_dec = current_fov_coords[0], current_fov_coords[1]
            
            self.__rm_from_stack(current_fov_ra,current_fov_dec)
            
            self.delete_pointing(infile="GWsky_pointings.txt",
                                 ra=str(current_fov_ra), dec=str(current_fov_dec))           
            self.destroy()
       
        if event.widget == self.B02:  # Zoom in the  FoV
            aladin.location(str(self.input_ra), str(self.input_dec))             
            aladin.zoom('1x')
            
    # Plots        
    def plot_stats(self, time_step, airmass_values, ra, dec, table_q,
                   prob_fov, r, dp_dr):
        """Showing the plots in the FoV statistic window."""
        
        f = Figure(figsize=(9, 6.2), facecolor='white')
        
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

        def subplot_cat_distribution():
            """SubPlot catalogue distributions."""

            # split entry: ra and dec
            current_fov_coords = self.entry_current_fov.get().split() 
            current_fov_ra, current_fov_dec = current_fov_coords[0], current_fov_coords[1]

            try:
                for table_name in table_q.keys():
                    table = table_q[table_name]
                
                table.write("GWsky_query_items", format = 'votable', overwrite = True)
                samp.send_file("GWsky_query_items")
                #aladin.rename("Q:"+str(ra)+"/"+str(dec))
                aladin.set_planeID("Q:"+str(ra)+"/"+str(dec))    #------------------#
                aladin.remove('GWsky_query_items')               # try to do better #
                                                                 #------------------#              
                aladin.set_planeID("C_" + current_fov_ra + "/" + current_fov_dec)               
                aladin.remove("C_" + current_fov_ra + "/" + current_fov_dec+'~1')

                Utils.move_to_folder(planes=['Q:*','C_*','P:*'],
                                       folders=['Queries','Coords','FoV'])
                
                #table.show_in_browser(jsviewer=True)        
                query_catalog = f.add_subplot(222)

                # column 1
                query_catalog.hist(table[self.column_1].quantity)         
                query_catalog.set_xlabel('cat: ' + self.catalog + ' ' +
                                         'col: ' + self.column_1)
                query_catalog.set_ylabel('Count')                     

                # column 2
                query_catalog = f.add_subplot(224) 
                query_catalog.hist(table[self.column_2].quantity)                   
                query_catalog.set_xlabel('cat: ' + self.catalog + ' ' +
                                         'col: ' + self.column_2) 
                query_catalog.set_ylabel('Count')                                
                
            except UnboundLocalError as unbound_local_error:
                #print ('No Galaxies in the selected FoV')
                tkMessageBox.showerror(' ', unbound_local_error)
            except KeyError as key_error:
                #print ('No catalog column:', key_error)
                tkMessageBox.showerror(' ', key_error)
            finally:
                pass
            
        subplot_cat_distribution = subplot_cat_distribution()                                   

        def subplot_cond_distance():
            """Conditional Distance Distribution Along a Line of Sight (FoV center position)."""
            
            conditional_distance_line_sight = f.add_subplot(221) 
            conditional_distance_line_sight.plot(r, dp_dr)
            
            title_string = 'Conditional distance distribution \n along the FoV center'
            conditional_distance_line_sight.set_title(title_string,fontsize=10)
            conditional_distance_line_sight.set_xlabel('distance (Mpc)')
            conditional_distance_line_sight.set_ylabel('prob Mpc$^{-1}$')

        subplot_cond_distance = subplot_cond_distance()

        def draw_area():
            """Drawing Area of the FoV Statistic window."""

            fov_information_title = "FoV center (ra "+str(ra) + "  " +"dec "+ str(dec)+")" + "; " + "prob: " + str(prob_fov)
            f.suptitle(fov_information_title, fontsize=14)
           
            canvas = FigureCanvasTkAgg(f, self) 
            canvas.show()
            canvas.get_tk_widget().pack(side=BOTTOM, fill=BOTH, expand=True)

            toolbar = NavigationToolbar2TkAgg(canvas, self)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=True)
            
        draw_area = draw_area()

        self.update_pointings_file("GWsky_pointings.txt", ra, dec, prob_fov)

def on_closing():
    """Asking the closure of the coverage window. If "Quit" the files in the list "temp_files" are deleted.
          ***Improving with tempfile module***"""
    
    if tkMessageBox.askokcancel("Quit", "Do you want to quit?"):
        try:
            temp_files=["GWsky_entries", "GWsky_query_items", "GWsky_fov.vot",
                        "GWsky_config", "GWsky_coords"]
            for temp_file in temp_files:
               os.remove(temp_file)
        except OSError:
            pass
        mainWin.destroy()

        
# running
mainWin = Tk()

sscGUI = ShowSkyCoverageGUI(mainWin)
mainWin.title('GWsky')
mainWin.attributes("-topmost", True)

mainWin.wait_visibility(mainWin)
mainWin.wm_attributes('-alpha', 0.8) # transparency

mainWin.protocol("WM_DELETE_WINDOW", on_closing)
mainWin.mainloop()
