# -*- coding: utf-8 -*-

from __future__ import print_function

try:
   import cPickle as pickle
except:
   import pickle

import fileinput
import os
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.dates as mdates

import dateutil

# Python 3 support

from tkinter import *
from tkinter import font, messagebox

from math import cos, sin, acos, asin, atan, degrees, radians, pi

import astropy
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.table import Table

import numpy as np

from mocpy import MOC

from aladinSAMP import AladinScriptCommands 
aladin = AladinScriptCommands()

from localize import LocalizeSources
from load_skymap import LoadSkymap

aladin.setconf_icrsd() # setting ICRSd

from config_values import UserValues

from lvc_skymap import LVCskymap
#lvc = LVCskymap()

from query import Query
query = Query()

from airmass import Airmass
#airmass = Airmass()

from moon import Moon
moon = Moon()

from moc_region import MOC_confidence_region



# creating folders in Aladin planes
from utils import Utils
Utils.create_folders(folders=["Queries", "Coords", "FoV"])
Utils.load_user_fov("GWsky_fov.vot")

# init. votable
table = Table()
table['a'] = [0, 0]
table.write("GWsky_query_items", format =
            'votable', overwrite = True)


# global variable: level of trasparency window
user = UserValues() 
trasparency = user.get_win_trasparency()


class SkyCoverage(object): 
    """Moving the FoV-footprint coverage by choosing a starting pointing."""
    
    SHIFT_CORRECTION = 0.00001  # A shift correction of 0.00001 is added
                                 # --> to escape math error during the FoV sequence

    def __init__(self, infile_coords='GWsky_coords'):
        """Creating a class in which the instance attributes are based on the dictionary
       "GWsky_coords" by default. GWsky_coords is created and pickled by the "config_values"
       module and will be deleted when the "SkyCoverageGUI" will be closed.
       It contains the keys: "ra", "dec". ra and dec represents the central location of a FoV. 
        
        Starting sky coordinates:
             self.input_ra: right ascension [deg]
             self.input_dec: declination [deg]
         """
       
        self.infile_coords = infile_coords
        self.entries_GWsky_new =[] # new entries during the FoV sequence
        
        self.user = UserValues() # composition
        #self.lvc = LVCskymap()
        #self.airmass = Airmass()
        
        with open(infile_coords, 'rb') as data:  
            coords_GWsky = pickle.load(data)
            
        for k, v in coords_GWsky.items():          
            setattr(self, k, v)
            
            
    def ra0ra1(self, A, dec0, dec1):
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

    def __fov_stats(self, ra, dec, moon_illumination, sep_fov_moon, query_result,
                    r, dp_dr, airmass_values, datestrings, prob_fov):    
        """Managing the descriptive statistic window."""
            
        self.__updating_center_coords(ra,dec) # updating file to read value
        
        fov_statistics = FoVstatistics()
        fov_statistics.plot_stats(ra,dec,moon_illumination, sep_fov_moon, query_result,
                                  r, dp_dr,airmass_values, datestrings, prob_fov)
        
        
    def update_pointings_file(self, infile, ra, dec, prob_fov, skymap):
         """The central location (ra[deg], dec[deg]) and the integrated probability of
             a selected FoV are saved locally in an external file.
             By default the file is named "GWsky_pointings.txt"."""
           
         with open(infile, 'a') as pointing:
             pointing.write(str(ra) + ' ' + str(dec)+ ' ' + str(prob_fov) + ' ' + skymap +'\n')

    def __query_shape(self, ra, dec, fov_shape):
        """Return the catalog query according with the defined-user FoV shape:
                   (1) box and (2) circle. """
        
        if self.user.get_fov_shape() != 2:  # box FoV
                    query_result = query.query_box(
                       ra, dec, self.user.get_fov_width(), self.user.get_fov_height(), self.user.get_catalog())              
        else: # circle FoV
            query_result = query.query_circle(
               ra, dec, self.user.get_fov_radius(), self.user.get_catalog())
            
        return query_result

    def __prob_shape(self, ra, dec, fov_shape):
        """Return the integrated probability according with the defined-user FoV shape:
                   (1) box and (2) circle."""
        
        if self.user.get_fov_shape() !=2:  # box FoV
            prob_fov = self.lvc.prob_in_box(
               ra, dec, self.user.get_fov_width(), self.user.get_fov_height())                
        else: # circle FoV
            prob_fov = self.lvc.prob_in_circle(
               ra, dec, self.user.get_fov_radius())
            
        return  prob_fov
                  
    def pick_coverage(self, ra, dec):
        """Setting GWsky: with statistic window (A); without statistic window (D)."""

        self.lvc = LVCskymap()
        self.airmass = Airmass()

        if self.user.get_GWsky_basic() == "A":  # full version 
            
            prob_fov = self.__prob_shape(ra, dec,
                                         self.user.get_fov_shape()) # integrated prob
            
            moon_illumination =  moon.illumination()           # moon_illumination
            
            sep_fov_moon = moon.from_fov(ra, dec)*u.deg        # moon_dist
            sep_fov_moon = sep_fov_moon.round(1)

            query_result = self.__query_shape(ra*u.deg, dec*u.deg,
                                              self.user.get_fov_shape())  # query

            airmass_values, datestrings = self.airmass.airmass_step(ra, dec)
            
            # TEST------------------------------------------------------#    
            #fov_sep = Utils.separation(self.input_ra, self.input_dec,   #
            #                            ra, dec)                        #
            #print ('The distance between 2 consecutive FoV centers is', #
            #       fov_sep.round(6))                                    #
            #-----------------------------------------------------------#
            
            r, dp_dr = self.lvc.conditional_distance_linesight(ra, dec)
            
            self.__fov_stats(ra, dec, moon_illumination, sep_fov_moon, query_result,
                             r, dp_dr, airmass_values, datestrings, prob_fov)          #  Stats win

            del(prob_fov) # memory problem for large nside
            #print (r, dp_dr)
            
        elif self.user.get_GWsky_basic() == "D":  # basic version-> no Stats win
            prob_fov = self.__prob_shape(ra, dec,
                                         self.user.get_fov_shape()) 

            # TEST------------------------------------------------------#  
            #fov_sep = Utils.separation(self.input_ra, self.input_dec,  #
            #                            ra, dec)                       #
            #print ('The distance between 2 consecutive FoV centers is',#
            #       fov_sep.round(6))                                   #
            #-----------------------------------------------------------#
            
            self.update_pointings_file("GWsky_pointings.txt", ra, dec, prob_fov, skymap=self.user.get_skymap())

            del(prob_fov) # memory problem for large nside
            
    def intercardinal_distance(self, ra, dec, shift_up_down, shift_right_left):
        """Moving from the fixed cardinal direction using the bi-directional windows;
           shift_up_down ↕ and/or shift_right_left ↔."""

        if shift_right_left > 0:
           shift_east_west = self.ra0ra1((shift_right_left-self.SHIFT_CORRECTION),
                                                   (dec + self.user.get_fov_height() + shift_up_down),
                                                   (dec + self.user.get_fov_height() + shift_up_down))
           dist = ra + shift_east_west 
         
        elif shift_right_left < 0 :
           shift_east_west = self.ra0ra1((shift_right_left + self.SHIFT_CORRECTION),
                                                   (dec + self.user.get_fov_height() + shift_up_down),
                                                   (dec + self.user.get_fov_height() + shift_up_down))
           dist = ra - shift_east_west
         
        else:
           dist = ra

        return dist 

    def load_entries(self, infile_entries):
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
            messagebox.showerror ('Error', message)
                         
    def north(self, shift_up_down=0, shift_right_left=0):
        """Moving the FoV tiles in North direction."""

        entries_GWsky = self.load_entries("GWsky_entries")
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]  

        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
            dist = self.intercardinal_distance(float(ra_start), float(dec_start),
                                                 shift_up_down, shift_right_left)
            north_pointing = [(dist),
                               (float(dec_start) + self.user.get_fov_height() + shift_up_down)]
             
            ra, dec = round(north_pointing[0], 5), round(north_pointing[1], 5)
            
            assert (dec >= -90 and dec <= 90 and ra >= 0 and ra <= 360), \
                   'Latitude angle(s) must be within -90 deg <= angle <=90 deg, got' + ' ' + str(dec) + ' ' + 'deg'

            aladin.get_FoV(ra, dec)
            self.pick_coverage(ra, dec)
            
            new_sky_pos = [ra,dec] # cycle variables
            self.entries_GWsky_new.extend(new_sky_pos)

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)
            
    def south(self, shift_up_down=0, shift_right_left=0):    
        """Moving the FoV tiles in South direction."""

        entries_GWsky = self.load_entries("GWsky_entries")
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]

        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
            dist = self.intercardinal_distance(float(ra_start), float(dec_start),
                                                 shift_up_down, shift_right_left)
            south_pointing = [(dist), (float(dec_start) - self.user.get_fov_height() - shift_up_down)]
                    
            ra, dec = round(south_pointing[0], 5), round(south_pointing[1], 5)

            assert (dec >= -90 and dec <= 90 and ra >= 0 and ra <= 360), \
                   'Latitude angle(s) must be within -90 deg <= angle <=90 deg, got' + ' ' + str(dec) + ' ' + 'deg'
            
            aladin.get_FoV(ra, dec)
            self.pick_coverage(ra, dec)
                                
            new_sky_pos = [ra,dec] # cycle variables
            self.entries_GWsky_new.extend(new_sky_pos)

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)    
          
    def east(self, shift_up_down=0, shift_right_left=0):
        """Moving the FoV tiles in East direction.
           A shift correction of 0.00001 is added to escape math error."""
        
        entries_GWsky = self.load_entries("GWsky_entries")
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]

        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):              
            ra_distance = self.ra0ra1((self.user.get_fov_width() - self.SHIFT_CORRECTION + shift_right_left),
                                        float(dec_start), float(dec_start))
                
            east_pointing = [(float(ra_start) + ra_distance), (float(dec_start) + shift_up_down)]
            ra, dec = round(east_pointing[0], 5), round(east_pointing[1], 5)

            assert (dec >= -90 and dec <= 90 and ra >= 0 and ra <= 360), \
                   'Longitude angle(s) must be within 0 deg <= angle <=360 deg, got' + ' ' + str(ra) + ' ' + 'deg'

            aladin.get_FoV(ra, dec)
            self.pick_coverage(ra, dec)           

            new_sky_pos = [ra,dec] # cycle variables
            self.entries_GWsky_new.extend(new_sky_pos)

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)
                   
    def west(self, shift_up_down=0, shift_right_left=0):
        """Moving the FoV tiles in West direction.
            A shift correction of 0.00001 is added to escape math error."""
        
        entries_GWsky = self.load_entries("GWsky_entries")
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]  
      
        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
               
            ra_distance = self.ra0ra1((self.user.get_fov_width() - self.SHIFT_CORRECTION + shift_right_left),
                                      float(dec_start), float(dec_start))

            west_pointing = [(float(ra_start) - ra_distance), (float(dec_start) + shift_up_down)]
            ra, dec = round(west_pointing[0], 5), round(west_pointing[1], 5)

            assert (dec >= -90 and dec <= 90 and ra >= 0 and ra <= 360), \
                   'Longitude angle(s) must be within 0 deg <= angle <=360 deg, got' + ' ' + str(ra) + ' ' + 'deg'

            aladin.get_FoV(ra, dec)
            self.pick_coverage(ra, dec)
            
            new_sky_pos = [ra,dec] # cycle variables
            self.entries_GWsky_new.extend(new_sky_pos)

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)
        

class SkyCoverageGUI(Toplevel):
    """Main Windows: Coverage"""
    
    def __init__(self, tkMainWin):
        frame = Frame(tkMainWin, border=9, bg="dim grey")
        frame.pack()
        
        self.B02 = Button(frame,text="N", command=self.north_btn)   
        self.B02.grid(row=0, column=2)
        
        self.B12 = Button(frame,text="↕↔", command=self.north_shift_btn,
                          fg="grey")
        self.B12.grid(row=1, column=2)
  
        self.B30 = Button(frame, text="E", command=self.east_btn)  
        self.B30.grid(row=3,column=0)

        self.B31 = Button(frame,text="↕↔",command=self.east_shift_btn,
                          fg="grey")   
        self.B31.grid(row=3, column=1)

        self.B32 = Button(frame, text="Start FoV", command=self.start_btn,
                          fg="red4",)
        self.B32.grid(row=3, column=2)

        self.B33 = Button(frame,text="↕↔",command=self.west_shift_btn,
                          fg="grey")  
        self.B33.grid(row=3, column=3)

        self.B34 = Button(frame, text="W", command=self.west_btn) 
        self.B34.grid(row=3, column=4)

        self.B42 = Button(frame,text="↕↔",command=self.south_shift_btn,
                          fg ="grey") 
        self.B42.grid(row=4, column=2)
        
        self.B52 = Button(frame,text="S", command=self.south_btn) 
        self.B52.grid(row=5, column=2)

        # Adjustments Btns
        self.B60 = Button(frame,text="↞", command=self.adj_east_btn,
                          fg ="grey", pady=3) 
        self.B60.grid(row=6, column=0)

        self.B61 = Button(frame,text="↠", command=self.adj_west_btn,
                          fg ="grey", pady=3,) 
        self.B61.grid(row=6, column=1)

        self.B62 = Button(frame,text="✓ Accept",fg ="green4", command=self.adj_accept_btn,
                          pady=1, padx=11) 
        self.B62.grid(row=6, column=2)

        self.B63 = Button(frame,text="↟", command=self.adj_north_btn,
                          fg ="grey",pady=3) 
        self.B63.grid(row=6, column=3)

        self.B64 = Button(frame,text="↡", command=self.adj_south_btn,
                          fg ="grey",pady=3) 
        self.B64.grid(row=6, column=4)

        # ▶ Folder
        self.B72 = Button(frame,text="▶ Folder", command=self.folder_btn,
                          fg ="gold4",pady=1,padx=10) 
        self.B72.grid(row=7, column=2)

        # MOC plot default
        self.B82 = Button(frame,text="MOC plot", command=self.moc_plot_btn,
                          fg ="steel blue",pady=1, padx=12) 
        self.B82.grid(row=8, column=2)

        # ⧗ ObsMOC
        self.B80 = Button(frame,text="⧗", command=self.obs_in_MOC_btn,
                          fg ="black",pady=3) 
        self.B80.grid(row=8, column=0)

        # ☰ ObsSkymap
        self.B81 = Button(frame,text="☰", command=self.obs_skymap_btn,
                          fg ="black",pady=3) 
        self.B81.grid(row=8, column=1)

        # ◉ localize sources
        self.B83 = Button(frame,text="◉",
                          fg ="black",pady=3, command=self.localize_btn) 
        self.B83.grid(row=8, column=3)

        # ⬊ load new skymap
        self.B84 = Button(frame,text="⬊", command=self.load_skymap_btn,
                          fg ="black",pady=3) 
        self.B84.grid(row=8, column=4)

           
    # Actions

    # north
    def north_btn(self):
        run_sequence = SkyCoverage() 
        return run_sequence.north()
      
    # ↕↔
    def north_shift_btn(self):
        move_fov = ShiftFoV()
        move_fov.north_shift()              

    # east
    def east_btn(self):
        run_sequence = SkyCoverage()
        run_sequence.east()                 

    # ↕↔
    def east_shift_btn(self):
        move_fov = ShiftFoV()
        move_fov.east_shift()
        
    # start FoV
    def start_btn(self):
        new_starting_fov = StartingFoV() 
        return new_starting_fov           

    # west
    def west_btn(self):
        run_sequence = SkyCoverage()
        run_sequence.west()          

    # ↕↔
    def west_shift_btn(self):
        move_fov = ShiftFoV()
        move_fov.west_shift()        

    # south
    def south_btn(self):
        run_sequence = SkyCoverage()
        run_sequence.south()         

    # ↕↔
    def south_shift_btn(self):
        move_fov = ShiftFoV()
        move_fov.south_shift()       

    # ↞               
    def adj_east_btn(self):           
        adj = Adjustments()
        adj.adj_east()

    # ↠   
    def  adj_west_btn(self):        
        adj = Adjustments()
        adj.adj_west() 

    # ↟                        
    def adj_north_btn(self):       
        adj = Adjustments()
        adj.adj_north()        

    # ↡ 
    def adj_south_btn(self):        
        adj = Adjustments()
        adj.adj_south()

    # ✓ Accept  
    def adj_accept_btn(self):         
        adj = Adjustments()
        adj.adj_accept()      

    # ▶ Folder
    def folder_btn(self):
        Utils.move_to_folder(planes=['Q:*', 'P:*', 'C_*'],
                             folders=['Queries', 'FoV', 'Coords'])
    # ⧗  ObsMOC  
    def obs_in_MOC_btn(self):         
        obs_in_MOC = ObsInMOC()

    # ☰  ObsSkymap
    def obs_skymap_btn(self):         
        obs_skymap = ObsSkymap()

    # MOC plot   
    def moc_plot_btn(self):      
        self.moc = MOC_confidence_region()
        self.user = UserValues()
        self.moc.contour_default(_from=10, _to=100, _step=40,
                                 skymap=self.user.get_skymap())        

    # ◉ localize sources
    def localize_btn(self):        
        localize = LocalizeSources()

    # ⬊ loading a new skymap
    def load_skymap_btn(self):
        self.load = LoadSkymap()            
        
class Adjustments(SkyCoverage):
    """Adjustments FoV position."""

    def __init__ (self):

       SkyCoverage.__init__(self)
       
       self.shift_up = 1     # default adjustments   (up)
       self.shift_down = 1   #      "               (down)
       self.shift_left = 1   #      "               (left)
       self.shift_right = 1  #      "               (right)

       
    def adj_north(self):
        """Adjustments FoV position -> north direction"""
            
        entries_GWsky = self.load_entries("GWsky_entries")        
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]
            
        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
            ra_start, dec_start = float(ra_start), float(dec_start)
               
            aladin.select("P:"+str(ra_start) + ',' + str(dec_start))
                
            dist = self.intercardinal_distance(ra_start, dec_start,
                                               self.shift_up, shift_right_left=0)
            north_adj = [(dist),
                         (dec_start + 0 + self.shift_up)]
             
            ra, dec = north_adj[0], north_adj[1]
                
            aladin.set_target(ra, dec)
            aladin.set_plane_id("P:"+str(ra) + ',' + str(dec))
                
            new_sky_pos = [ra,dec] # cycle variables
            self.entries_GWsky_new.extend(new_sky_pos)
            
            #aladin.remove("Q:"+str(ra_start)+"/"+str(dec_start))
            Utils.delete_pointing(infile="GWsky_pointings.txt",
                                  ra=str(ra_start), dec=str(dec_start))
            
            #aladin.remove("C_" + str(ra_start) + "/" +  str(dec_start))

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)

    def adj_south(self):
        """Adjustments FoV position -> south direction"""
         
        entries_GWsky = self.load_entries("GWsky_entries")        
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]
            
        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
            ra_start, dec_start = float(ra_start), float(dec_start)
               
            aladin.select("P:"+str(ra_start) + ',' + str(dec_start))
                
            dist = self.intercardinal_distance(ra_start, dec_start,
                                               self.shift_down, shift_right_left=0)
            south_adj = [(dist),
                         (dec_start + 0 - self.shift_down)]
             
            ra, dec = south_adj[0], south_adj[1]
                
            aladin.set_target(ra, dec)
            aladin.set_plane_id("P:"+str(ra) + ',' + str(dec))
                
            new_sky_pos = [ra,dec] # cycle variables
            self.entries_GWsky_new.extend(new_sky_pos)
            
            #aladin.remove("Q:"+str(ra_start)+"/"+str(dec_start))
            Utils.delete_pointing(infile="GWsky_pointings.txt",
                                  ra=str(ra_start), dec=str(dec_start))
            #aladin.remove("C_" + str(ra_start) + "/" +  str(dec_start))

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)

    def adj_east(self):
        """Adjustments FoV position -> east direction"""

        entries_GWsky = self.load_entries("GWsky_entries")
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]

        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
            ra_start, dec_start = float(ra_start), float(dec_start)

            aladin.select("P:"+str(ra_start) + ',' + str(dec_start))
            
            ra_distance = self.ra0ra1((0 - self.SHIFT_CORRECTION + self.shift_left),
                                        float(dec_start), float(dec_start))
                          
            east_adj = [(float(ra_start) + ra_distance), (float(dec_start) + 0)]
            ra, dec = round(east_adj[0],5), round(east_adj[1],5)

            aladin.set_target(ra, dec)
            aladin.set_plane_id("P:"+str(ra) + ',' + str(dec))       

            new_sky_pos = [ra,dec] # cycle variables
            self.entries_GWsky_new.extend(new_sky_pos)

            #aladin.remove("Q:"+str(ra_start)+"/"+str(dec_start))
            Utils.delete_pointing(infile="GWsky_pointings.txt",
                                  ra=str(ra_start), dec=str(dec_start))
            
            #aladin.remove("C_" + str(ra_start) + "/" +  str(dec_start))

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)
            
    def adj_west(self):
        """Adjustments FoV position -> west direction"""
         
        entries_GWsky = self.load_entries("GWsky_entries")
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]

        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
            ra_start, dec_start = float(ra_start), float(dec_start)

            aladin.select("P:"+str(ra_start) + ',' + str(dec_start))
             
            ra_distance = self.ra0ra1((0 - self.SHIFT_CORRECTION + self.shift_right),
                                        float(dec_start), float(dec_start))
            
            aladin.select("P:"+str(ra_start) + ',' + str(dec_start))
                
            west_adj = [(float(ra_start) - ra_distance), (float(dec_start) + 0)]
            ra, dec = west_adj[0], west_adj[1]

            aladin.set_target(ra, dec)
            aladin.set_plane_id("P:"+str(ra) + ',' + str(dec))       

            new_sky_pos = [ra,dec] # cycle variables
            self.entries_GWsky_new.extend(new_sky_pos)

            #aladin.remove("Q:"+str(ra_start)+"/"+str(dec_start))
            Utils.delete_pointing(infile="GWsky_pointings.txt",
                                  ra=str(ra_start), dec=str(dec_start))

            #aladin.remove("C_" + str(ra_start) + "/" +  str(dec_start))

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)

    def adj_accept(self):
        """Confirming the adjustments FoV position -> open statistic win."""
         
        entries_GWsky = self.load_entries("GWsky_entries")        
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]
                   
        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
            ra_start, dec_start = float(ra_start), float(dec_start)
                    
            self.pick_coverage(float(ra_start), float(dec_start))

            new_sky_pos = [ra_start,dec_start]
            self.entries_GWsky_new.extend(new_sky_pos)

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)

    
class StartingFoV(Toplevel):
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
        
        self.user = UserValues()

        # putting the entry value(s) in a list
        self.entries_GWsky=[]  

        self.wait_visibility()

        self.user = UserValues()
        # get trasparency windows
        #trasparency = self.user.get_win_trasparency()
        self.wm_attributes('-alpha', trasparency)   

        self.title(" Starting FoV")
        self.attributes("-topmost", True)
        
        self.label_1 = Label(self, text="RA (°) DEC (°)", bg="slate grey")
        self.label_1.grid(row=0, column=0, sticky=E, pady=0)

        # default: sky coords of maximum probability pixel
        fov_coords = str(self.user.get_ra_max_pixel()), str(self.user.get_dec_max_pixel()) 
        
        max_pixel_default = StringVar(self, value=fov_coords) 
        self.entry_1 = Entry(self, width=30, justify=CENTER,
                             textvariable=max_pixel_default)

        self.entry_1.grid(row=0, padx=15, column=1)

        self.entryScroll = Scrollbar(self, orient=HORIZONTAL,
                                     command=self.__scrollHandler)
        self.entryScroll.grid(row=1, column=1, sticky=E+W)
        self.entry_1['xscrollcommand'] = self.entryScroll.set

        #Btns
        self.show = Button(self, text='Show',
                           command=self.show_starting_fov)
        self.show.grid(column=2, row=0, sticky=W, padx=2, pady=5)
        
        self.checkbox = Button(self, text="Not show",      
                               command=self.no_show_starting_fov)
        self.checkbox.grid(column=3,row=0, sticky=E, padx=2, pady=5)

        self.browse_1 = Button(self, text="Browse ...",
                               command=self.browsecsv)
        self.browse_1.grid(row=0, column=4, padx=2, pady=2)

        #self.close = Button(self, text="Obs",  fg='dark green',
        #                    command=self.obs)  
        #self.close.grid(column=4,row=0, sticky=W, padx=2, pady=5)

        self.close = Button(self, text="Close",
                            command=self.close_window)  
        self.close.grid(column=6,row=0, sticky=E, padx=2, pady=5)


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
        
        current_fov_coords = self.entry_1.get().replace(';',' ').replace(',',' ').split()
        fov_center_ra = current_fov_coords[0::2]
        fov_center_dec = current_fov_coords[1::2]

        return current_fov_coords, fov_center_ra, fov_center_dec
    
    def show_starting_fov(self):
        """Drawing the FoV footprint(s) in the Aladin plane(s).
         By default: sky coords (ra[deg], dec[deg]) of maximum probability pixel."""   
   
        show_sky_coverage = SkyCoverage()
        
        current_fov_coords, fov_center_ra, fov_center_dec = self.__split_entries()
        
        try:
            for ra_starting, dec_starting in zip (fov_center_ra, fov_center_dec):
                aladin.get_FoV(float(ra_starting), float(dec_starting))
                show_sky_coverage.pick_coverage(float(ra_starting), float(dec_starting))
                              
        except ValueError as value_error:
            messagebox.showerror ('Error', value_error)
          
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

    def browsecsv(self):
        from tkFileDialog import askopenfilename

        #Tk().withdraw() 
        filename_from_browser = askopenfilename(filetypes=[("fits files","*.fits.gz")])
        #fits = filename.read()
        print (filename_from_browser)

     
    def close_window(self):
        self.destroy()

class ObsSkymap(Toplevel):
    """Initializi"""

    def __init__(self):
        Toplevel.__init__(self, border=8, bg="slate grey")
        
        self.user = UserValues()
        self.moc = MOC_confidence_region()

        self.observatory = astropy.coordinates.EarthLocation(
           lat=self.user.get_latitude()*u.deg, lon=self.user.get_longitude()*u.deg,
           height=self.user.get_altitude()*u.m)

        self.obs_time = Time(self.user.get_obs_time())
        self.nside = self.user.get_nside()

        self.wait_visibility()
        
        self.user = UserValues()
        # get trasparency windows
        #trasparency = self.user.get_win_trasparency()
        self.wm_attributes('-alpha', trasparency)   

        self.title("ObsSkymap --" + "Observability" + " " + "starting from" + \
                   " " + self.user.get_obs_time() + " " + "UT")
        self.attributes("-topmost", True)

        self.bttn_clicks = 0 # counter ">>" Button
        
        # first label
        self.label_1 = Label(self, text="Show the airmass regions over the skymap from 1 to 4 in step of 1",
                             bg="slate grey")
        self.label_1.grid(row=0, column=0, sticky=E, padx=0)

        #Btn
        self.show = Button(self, text='Show',
                           command=self.snapshot_airmass)
        self.show.grid(column=4, row=0, sticky=W, padx=2, pady=5)
        
        self.forward = Button(self, text=">>",      
                               command=self.snapshot_airmass_up)
        self.forward.grid(column=5,row=0, sticky=E, padx=2, pady=5)

        self.close = Button(self, text="Close",
                            command=self.close_window)  
        self.close.grid(column=7,row=0, sticky=E, padx=2, pady=5)

    #Actions
    def update_count(self):
        """Return the time in step of 1h when the button ">>" is clicked."""
   
        self.bttn_clicks += 1
              
        dt = TimeDelta(3600.0, format='sec')
        update_time = int(self.bttn_clicks) * dt
        
        obs_time = Time(self.user.get_obs_time())
        time_update = obs_time + update_time
        
        return time_update

    def from_ipixs_to_moc(self, time_input, percentage = 0.9):
        """Return ipix table with the associated airmass"""

        prob = self.moc.read_prob(self.user.get_skymap())
        print (self.user.get_skymap())
        #percentage = 
        #float(self.entry_percentage.get())/100.0
        
        ipixs = self.moc.ipixs_in_percentage(prob, percentage )
        nside = int(self.user.get_nside())
        
        ra, dec = self.moc.sky_coords(ipixs, nside)
        
        sky_coord = SkyCoord(ra = ra*u.deg,dec=dec*u.deg, frame='icrs')
        altaz = sky_coord.transform_to(AltAz(obstime=time_input, location=self.observatory))
        airmass_values = altaz.secz
        
        contour_ipix = Table([ ra, dec, airmass_values, ipixs ],
                             names = ('RA[deg]', 'DEC[deg]', 'airmass', 'ipix'),
                             meta = {'ipix': 'ipix table'})             # astropy table       

        return contour_ipix


    def snapshot_airmass(self):
        """Showing an airmass snapshot at a given time."""
        
        #obs_time = Time(self.user.get_obs_time())
        nside = self.user.get_nside()
        
        contour_ipix = self.from_ipixs_to_moc(self.obs_time)
        
        airmass_start = [1, 2, 3, 4]
        airmass_end = [2, 3, 4, 5.8]

        aladin.md( "airmass@"+str(self.obs_time))
        
        for i,j in zip(airmass_start, airmass_end):
         
            snap_1 = contour_ipix[(contour_ipix['airmass'] >= i) & (contour_ipix['airmass'] < j) ]
            moc_order = self.moc.moc_order(nside)


            snap_1['RA[deg]'].unit = 'deg'
            snap_1['DEC[deg]'].unit = 'deg'

            moc = MOC.from_lonlat( snap_1['RA[deg]'], snap_1['DEC[deg]'],
                                  moc_order )                # moc creation
            moc.write( 'snap_airmass_'+'initial', format = 'fits', write_to_file = True)     # fits file

            if len(snap_1)!=0:
                aladin.send_file('snap_airmass_'+'initial')
                aladin.rename(str(i)+"=<airmass<"+str(j))
                #aladin.set_planeID(str(i)+"=<airmass<"+str(j))
                aladin.mv(str(i)+'=<airmass<'+str(j), '"'+'airmass@'+str(self.obs_time)+'"' )
                

    def snapshot_airmass_up(self):
        """Showing an airmass snapshot increasing the time in step of 1 hour."""

        time_update = self.update_count()
        
        nside = self.user.get_nside()
        
        contour_ipix = self.from_ipixs_to_moc(time_update)
        #print (time_update)

        airmass_start = [1, 2, 3, 4]
        airmass_end = [2, 3, 4, 5.8]
        snaps = ['s1', 's2', 's3','s4']

        aladin.md( "airmass@"+str(time_update))
        
        for i,j, snap in zip(airmass_start, airmass_end, snaps):
   
            snap = contour_ipix[(contour_ipix['airmass'] >= i) & (contour_ipix['airmass'] < j) ]
            moc_order = self.moc.moc_order(nside)


            snap['RA[deg]'].unit = 'deg'
            snap['DEC[deg]'].unit = 'deg'

            moc = MOC.from_lonlat( snap['RA[deg]'], snap['DEC[deg]'],
                                  moc_order )                # moc creation
            moc.write( 'snap', format = 'fits', write_to_file = True )     # fits file

            if len(snap)!=0:
                aladin.send_file('snap')
                aladin.set_planeID(' '+str(i)+"<airmass<"+str(j)+'@'+str(time_update))
                aladin.mv('"'+str(i)+'<airmass<'+str(j)+'@'+str(time_update)+'"', '"'+'airmass@'+str(time_update)+'"' )

    def close_window(self):
        self.destroy()



class ObsInMOC(Toplevel):
    """The class is designed to define the sky area"""

    def __init__(self):
        Toplevel.__init__(self, border=8, bg="slate grey")
        
        self.user = UserValues()
        
        self.moc = MOC_confidence_region() #  non ha senso

        self.observatory = astropy.coordinates.EarthLocation(
           lat=self.user.get_latitude()*u.deg, lon=self.user.get_longitude()*u.deg,
           height=self.user.get_altitude()*u.m)

        self.wait_visibility()
        
        self.user = UserValues()
        # get trasparency windows
        #trasparency = self.user.get_win_trasparency()
        self.wm_attributes('-alpha', trasparency)   

        self.title("Observability" + " " + "starting from" + " " + self.user.get_obs_time() + " " + "UT")
        self.attributes("-topmost", True)

        self.bttn_clicks = 0 # counter ">>" Button
        
        # first label
        self.label_1 = Label(self, text="Show the region in the",
                             bg="slate grey")
        self.label_1.grid(row=0, column=0, sticky=E, padx=0)

        moc_value = 90  # default     
        moc_default = StringVar(self, value=moc_value)
        
        self.entry_percentage = Entry(self, width=5, justify=CENTER,
                             textvariable=moc_default)
        self.entry_percentage.grid(row=0, padx=2, column=1)

        # second label
        self.label_2 = Label(self, text="% MOC in which the airmass is ≤",
                             bg="slate grey")
        self.label_2.grid(row=0, column=2, sticky=E, pady=0)

        airmass_value = "2.5" # default
        airmass_default = StringVar(self, value=airmass_value)
        
        self.entry_airmass = Entry(self, width=5, justify=CENTER,
                             textvariable=airmass_default)
        self.entry_airmass.grid(row=0, padx=2, column=3)

        #Btn
        self.show = Button(self, text='Show',
                           command=self.moc_obs)
        self.show.grid(column=4, row=0, sticky=W, padx=2, pady=5)
        
        self.moon = Button(self, text="Moon",
                            command=self.get_moon_position)  
        self.moon.grid(column=6,row=0, sticky=W, padx=2, pady=5) 
        
        self.forward = Button(self, text=">>",      
                               command=self.moc_obs_update)
        self.forward.grid(column=5,row=0, sticky=E, padx=2, pady=5)

        self.close = Button(self, text="Close",
                            command=self.close_window)  
        self.close.grid(column=7,row=0, sticky=E, padx=2, pady=5)


    #Actions
    def update_count(self):
        """Return the time in step of 1h when the button ">>" is clicked."""
   
        self.bttn_clicks += 1
              
        dt = TimeDelta(3600.0, format='sec')
        update_time = int(self.bttn_clicks) * dt
        
        obs_time = Time(self.user.get_obs_time())
        time_update = obs_time + update_time
        
        return time_update

    def from_ipixs_to_moc(self, time_input):
        """Return ipix table with the associated airmass"""

        prob = self.moc.read_prob(self.user.get_skymap())
        percentage = float(self.entry_percentage.get())/100.0
        
        ipixs = self.moc.ipixs_in_percentage(prob, percentage )
        nside = int(self.user.get_nside())
        
        ra, dec = self.moc.sky_coords(ipixs, nside)
        
        sky_coord = SkyCoord(ra = ra*u.deg,dec=dec*u.deg, frame='icrs')
        altaz = sky_coord.transform_to(AltAz(obstime=time_input, location=self.observatory))
        airmass_values = altaz.secz
        
        contour_ipix = Table([ ra, dec, airmass_values, ipixs ],
                             names = ('RA[deg]', 'DEC[deg]', 'airmass', 'ipix'),
                             meta = {'ipix': 'ipix table'})             # astropy table       

        mask = (contour_ipix['airmass']) >= 1 # clearing
        obs1 = contour_ipix[mask]

        mask2 = (obs1['airmass']) <= float(self.entry_airmass.get())  # airmass user values
        obs = obs1[mask2]

        obs['RA[deg]'].unit = 'deg'
        obs['DEC[deg]'].unit = 'deg'        

        # TEST
        #print obs, "sono qui"
        
        nside = self.user.get_nside()

        #if len(obs)!=0:
        moc_order = self.moc.moc_order(nside)
        moc = MOC.from_lonlat( obs['RA[deg]'], obs['DEC[deg]'],
                               moc_order )                # moc creation
            
        moc.write( 'obs_airmass_', format = 'fits',write_to_file = True )     # fits file

        if len(obs)!=0:
           aladin.send_file('obs_airmass_')
         
        return obs, moc
      
    def moc_obs(self):
        """Return the MOC region in which the airmass is <= the airmass value defined by the user."""

        time_start = Time(self.user.get_obs_time())
        ipixs_start, moc_start = self.from_ipixs_to_moc(time_start)

        percentage = float(self.entry_percentage.get())/100.0

        if len(ipixs_start) !=0:
            aladin.rename('obs_airmass_'+self.entry_airmass.get()\
                      +'MOC_'+str(percentage)+ '@' + str(time_start.isot))

            # printing area
            square_degrees_sphere = (360.0**2)/pi           
            area_sq2 = round( ( moc_start.sky_fraction * square_degrees_sphere ), 1 )
            
            messagebox.showinfo('MOC visibility'+ '@' + str(time_start.isot),
                                 '   sky coverage: ' + str(area_sq2)+ ' ' + 'sq. deg')          
        else:
            messagebox.showinfo('MOC visibility'+ '@' + str(time_start.isot),
                                  'No region for the selected  airmass')

        #print time_start
        #print type (time_start)

        # TEST
        #print (time_start)

    def moc_obs_update(self):
        """Return the MOC region in which the airmass is <= the airmass value defined by the user."""

        time_update = self.update_count()
        ipixs_update, moc_update = self.from_ipixs_to_moc(time_update)

        percentage = float(self.entry_percentage.get())/100.0

        if len(ipixs_update) != 0:
            aladin.rename('obs_airmass_'+self.entry_airmass.get()+\
                          'MOC_'+str(percentage)+ '@' + str(time_update.isot))

            # printing area
            square_degrees_sphere = (360.0**2)/pi
            area_sq2 = round( ( moc_update.sky_fraction * square_degrees_sphere ), 1 )
            
            messagebox.showinfo('MOC visibility'+ '@' + str(time_update.isot),
                                 '   sky coverage: ' + str(area_sq2)+ ' ' + 'sq. deg')
            
        else:
            messagebox.showinfo('MOC visibility'+ '@' + str(time_update.isot),
                                  'No region for the selected  airmass')

        # TEST
        #print (time_update)

    def get_moon_position(self):
        moon = Moon()
        moon.sky_position()

        #moon.steps()

    def close_window(self):
        self.destroy()

        
class ShiftFoV(Toplevel):
    """Shifting the FoV footprint(s) from a consecutive cardinal direction (↕, ↕, ↔, ↔);
       The widget contains 2 entries and 2 Buttons."""
    
    def __init__(self):
        Toplevel.__init__(self, border=7, bg="slate grey")
        self.attributes("-topmost", True)
        self.wait_visibility()

        self.user = UserValues()
        # get trasparency windows
        #trasparency = self.user.get_win_trasparency()
        self.wm_attributes('-alpha', trasparency)     
        
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
            shift = SkyCoverage()
            shift.north(shift_up_down, shift_right_left)
        except ValueError as value_error:
            messagebox.showerror ('Error',value_error)

    def shift_south(self):             
        try:
            shift_up_down = float(self.entry_3.get())
            shift_right_left = float(self.entry_4.get())
            shift = SkyCoverage()
            shift.south(shift_up_down, shift_right_left)
        except ValueError as value_error:
            messagebox.showerror ('Error',value_error)

    def shift_east(self):             
        try:
            shift_up_down = float(self.entry_3.get())
            shift_right_left = float(self.entry_4.get())
            shift = SkyCoverage()
            shift.east(shift_up_down, shift_right_left)
        except ValueError as value_error:
            messagebox.showerror ('Error',value_error)

    def shift_west(self):              
        try:
            shift_up_down = float(self.entry_3.get())
            shift_right_left = float(self.entry_4.get())
            shift = SkyCoverage()
            shift.west(shift_up_down, shift_right_left)
        except ValueError as value_error:
            messagebox.showerror ('Error',value_error)

    def close_window(self):
        self.destroy()

class FoVstatistics(Toplevel, SkyCoverage):
    """ FoV statistics window consists of 3 buttons:

         (1) Confirm the FoV-footprint
         (2) Delete FoV-footprint
         (3) Zoom in FoV. """
    
    def __init__(self):
        Toplevel.__init__(self)
        SkyCoverage.__init__(self)

        self.title("FoV statistics")
        self.attributes("-topmost", True)
              
        self.B00 = Button(self,text="Confirm the pointing in the GWsky_pointings txt file")  
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
        
        self.entry_current_fov.bind("<Key>", lambda e: "break") # saving entries
        self.entry_current_fov.pack(side=TOP,fill=BOTH)

    def __rm_from_stack(self, ra, dec):
        """Removing from Aladin stack the associated planes."""
        
        aladin.remove_FoV(ra, dec)       
        aladin.remove("Q:"+ ra +"/"+ dec)           
        #aladin.remove("C_" + ra+ "/" + dec)

    def clicked(self, event):
        """Retain/Delete/Zoom a FoV-footprint.
           The FoV-center positions are saved in "GWsky_pointings.txt" """      
        
        if event.widget == self.B00: # Retain and Close the FoV. 
            self.destroy()
            
        if event.widget == self.B01: # Delete FoV
            current_fov_coords = self.entry_current_fov.get().split() # getting entries
            current_fov_ra, current_fov_dec = current_fov_coords[0], current_fov_coords[1]
            
            self.__rm_from_stack(current_fov_ra, current_fov_dec)
            
            Utils.delete_pointing(infile="GWsky_pointings.txt",
                                 ra=str(current_fov_ra), dec=str(current_fov_dec))           
            self.destroy()
       
        if event.widget == self.B02:  # Zoom in the  FoV
            aladin.location(str(self.input_ra), str(self.input_dec))             
            aladin.zoom('1x')

    def __are_all_same(self, items):
        """Check if all elements of a list are the same."""
        
        return all(x == items[0] for x in items)
      
    # Plots        
    def plot_stats(self, ra,dec, moon_illumination, sep_fov_moon, query_result,
                   r, dp_dr, airmass_values, datestrings, prob_fov):
        """Showing the plots in the FoV statistic window."""
        
        f = Figure(figsize=(9, 5.2), facecolor='white')
        f.subplots_adjust(left=.13, bottom=.16, right=.93, top=.84, wspace=.26, hspace=.3)
        

        def airmass_subplot():
            """SubPlot Airmass."""
            
            airmass_plt = f.add_subplot(223)
            airmass_plt.set_ylabel('airmass')
            
            same = self.__are_all_same(airmass_values)

            if same !=True:
                time_step = [dateutil.parser.parse(s) for s in datestrings]
            
                #ax=plt.gca()
                #ax.xaxis.set_major_formatter(
                #mdates.DateFormatter('%H:%M:%S'))

                airmass_plt.set_ylabel('airmass')
                airmass_plt.set_xlabel('Universal Time')
                airmass_plt.invert_yaxis()
                airmass_plt.grid(True)
           
                airmass_plt.plot(time_step, airmass_values, "o-", linestyle="--")
                f.autofmt_xdate()
            else:
                messagebox.showinfo('Warning',"airmass outside the range of 1 - 5.8")
                pass

        airmass_subplot = airmass_subplot()

        def subplot_cat_distrib():
            """SubPlot catalog distributions."""

            # setting column 1
            query_catalog_1 = f.add_subplot(222)
            query_catalog_1.set_xlabel('cat: ' + self.user.get_catalog() + ' ' +
                                       'col: ' + self.user.get_column_1())
            query_catalog_1.set_ylabel('Count')    

            # setting column 2
            query_catalog_2 = f.add_subplot(224)
            query_catalog_2.set_xlabel('cat: ' + self.user.get_catalog() + ' ' +
                                       'col: ' + self.user.get_column_2())
            query_catalog_2.set_ylabel('Count')

            try:
                # filtering and hist
                mask_1 = query_result[self.user.get_column_1()] > float(self.user.get_filter_1())
                query_result_filtered_1 = query_result[mask_1]

                mask_1_1 = query_result_filtered_1[self.user.get_column_1()] < float(self.user.get_filter_2())
                query_result_filtered_2 = query_result_filtered_1[mask_1_1]

                query_catalog_1.hist(query_result_filtered_2[self.user.get_column_1()])
                
                query_catalog_2.hist(query_result_filtered_2[self.user.get_column_2()])

                # creating vo table for each footprint
                query_result_filtered_2.write("GWsky_query_items", format =
                            'votable', overwrite = True)

                #aladin.send_file("GWsky_query_items")
                #aladin.rename("Q:"+str(ra)+"/"+str(dec))
         
            except KeyError as key_error:
                messagebox.showerror(' Error: no key found', key_error)
            except ValueError:
                c1, c2  = np.array(query_result_filtered_2[self.user.get_column_1()]), \
                              np.array(query_result_filtered_2[self.user.get_column_2()])
                newc1, newc2 = c1[~np.isnan(c1)], \
                                   c2[~np.isnan(c2)]
                query_catalog_1.hist(newc1)
                query_catalog_2.hist(newc2)
                         
        subplot_cat_distrib = subplot_cat_distrib()       

        def subplot_cond_dist():
            """Conditional Distance Distribution Along a Line of Sight (FoV center position)."""
            
            conditional_distance_line_sight = f.add_subplot(221) 
            conditional_distance_line_sight.plot(r, dp_dr)
            
            title_string = 'Conditional distance distribution \n along the FoV center'
            conditional_distance_line_sight.set_title(title_string,fontsize=10)
            conditional_distance_line_sight.set_xlabel('distance (Mpc)')
            conditional_distance_line_sight.set_ylabel('prob Mpc$^{-1}$')

        subplot_cond_dist = subplot_cond_dist()

        def draw_area():
            """Drawing Area of the FoV Statistic window."""

            fov_information_title = "FoV center (ra "+str(ra) + "  " +"dec "+ str(dec)+")" + "; " + "prob: " + str(prob_fov)+ \
                                    ";" + " " + "Moon" + " " + "(illumi.:" + " " + str(moon_illumination) + " " + \
                                    "dist.:" + " " + str(sep_fov_moon) + ")"
            
            f.suptitle(fov_information_title, fontsize=10)
           
            canvas = FigureCanvasTkAgg(f, self) 
            canvas.draw()
            canvas.get_tk_widget().pack(side=BOTTOM, fill=BOTH, expand=True)

            toolbar = NavigationToolbar2Tk(canvas, self)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=True)
            
        draw_area = draw_area()

        self.update_pointings_file("GWsky_pointings.txt", ra, dec, prob_fov, skymap=self.user.get_skymap())

        aladin.send_file("GWsky_query_items")
        aladin.rename("Q:"+str(ra)+"/"+str(dec))
        
        #aladin.draw_newtool("C_" + str( ra ) + "/" + str( dec ))
        #aladin.draw_string(ra, dec, str( ra ) + "/" + str( dec ))
        print (str(sep_fov_moon))
        print (sep_fov_moon.round(1))

def on_closing():
    """Asking the closure of the coverage window. If "Quit" the files in the list "temp_files" are deleted.
          ***Improving with tempfile module***"""
    
    if messagebox.askokcancel("Quit", "Do you want to quit?"):
        try:
            temp_files=["GWsky_entries", "GWsky_query_items", "GWsky_fov.vot",
                        "GWsky_config", "GWsky_coords", "obs_airmass_", "snap_airmass_initial", "snap",]
            for temp_file in temp_files:
               os.remove(temp_file)
        except OSError:
            pass
        mainWin.destroy()




# running
# ***TO BE COMPLETED***

mainWin = Tk()
sscGUI = SkyCoverageGUI(mainWin)

mainWin.title('GWsky')
mainWin.attributes("-topmost", True)

mainWin.wait_visibility(mainWin)
mainWin.wm_attributes('-alpha', trasparency) 

mainWin.protocol("WM_DELETE_WINDOW", on_closing)
mainWin.mainloop()
