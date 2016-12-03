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
import pickle

from astropy.vo.samp import SAMPHubError # no work in py3
from astropy.time import Time
from astropy.io.votable import parse
from astropy.io.votable.tree import VOTableFile, Resource, Field, Param
from astropy.utils.data import download_file

from aladinSAMP import AladinViaSAMP, AladinScriptCommands # class access
samp = AladinViaSAMP()
aladin = AladinScriptCommands()
        
    
class UserInput(object):
    """Insert user values at the command prompt to create the configuration files."""

    def __init__(self):
        """
            skymap : a valid LVC skymap in healpix format;
            prob_contours : a valid probability contours file in json format;
            fov_width : field-of-view width in degrees;
            fov_height : field-of-view height in degrees;
            fov_radius : field-of-view radius in degrees;
            latitude : Geodetic coordinates of the Observatory [deg];
            longitude :              "                         [deg];
            altitude :               "                          [m];
            obs_time : starting time yyyy-mm-dd hh:mm:ss;
            input_ra : right ascension in degrees;             
            input_dec : declination in degrees;                
            GWsky_coords : dict to store the sky coords;
            GWsky_config : dict to store the input values;
            GWsky_base : basic version of GWsky by typing "b";
            ra_max : right ascension of maximum probability pixel;   [deg]
            dec_max : declination of maximum probability pixel.      [deg]
        """
        
        self._skymap = str()
        self._prob_contours = str()       
        self._fov_width = float()
        self._fov_height = float()
        self._fov_radius = float()
        self._latitude = float()
        self._longitude = float()
        self._altitude = float()
        self._obs_time = str()
        self._input_ra = float()
        self._input_dec = float()
        self.GWsky_coords = dict()
        self.GWsky_config = dict()
        self.GWsky_base = str()
        self._ra_max = float()
        self._dec_max = float()

    def display_welcome(self):
        print('                 *************************           ')
        print('                 Welcome to GWsky -v3.5.1            ')
        print('                 *************************           ')
        print('                                                     ')
        print('')
        print (' Launch Aladin sky Atlas to run the script! ')
        print ('             Aladin java applet ')
        print (' http://aladin.u-strasbg.fr/java/nph-aladin.pl ')
        print ('             Aladin Desktop ')
        print (' http://aladin.u-strasbg.fr/java/nph-aladin.pl?frame=downloading ')
        print ('            and check the SAMP connection:')
        print ( ' Aladin bar --> Interop --> Connect with SAMP. ')
        print('')             
        print('  The *UserValues* module creates 1 configuration file:    ')
        print('  GWsky_config                                             ')
        print (' The configuration file is read by the module *coverage*. ')
        print('')

    def display_SampError(self):
        print('')
        print (' Launch Aladin sky Atlas to run the script! ')
        print ('             Aladin java applet ')
        print (' http://aladin.u-strasbg.fr/java/nph-aladin.pl ')
        print ('             Aladin Desktop ')
        print (' http://aladin.u-strasbg.fr/java/nph-aladin.pl?frame=downloading ')
        print ('            and check the SAMP connection:')
        print ( ' Aladin bar --> Interop --> Connect with SAMP. ')
        print('')

    def healpix_skymap(self):
        """Inserting a valid LIGO/Virgo probability skymap."""
        
        while True:
            try:
                print ('--> Load a probability skymap')
                self._skymap = raw_input( ' from a local directory : ' ).strip( )
                print('')
                prob = hp.read_map(self._skymap, verbose = False)
            except IOError as io_error:
                print (io_error)
            except ValueError as value_error:
                print (value_error)
            else:
                break

        npix = len(prob) # get healpix resolution
        nside = hp.npix2nside(npix)

        return self.GWsky_config.update({"skymap" : self._skymap,
                                       "nside" : nside})

    def update_GWsky_config(self):
        """Add input skymap in GWsky_config file."""
        
        return self.GWsky_config.update({"skymap" : self._skymap})
        
    def display_healpix_skymap(self):
        """Displaying the healpix skymap in Aladin plane."""

        while True:
            try:
                return samp.send_file(self._skymap)
            except SAMPHubError as samphub_error:
                print(samphub_error, self.display_SampError())
                time.sleep(5)
            else:
                break

    def plot_prob_contours(self, data, contour_pieces, percentile_json):
        """Managing the probability-contour file."""
        
        i = 0
        for i in range(0, contour_pieces):
            contour = data['contours'][i]
            percentile = contour['name']

            if percentile == percentile_json:
                values = contour['coords']
                line = (str(values).replace('[' , '' ).replace(']' , ''))
                aladin.draw_line (line) # sending to Aladin plane
            
    def probability_contours(self):
        """Inserting a valid LIGO/Virgo probability contours in json format."""
        
        while True:
            try:
                print ('--> Load a valid contour-line file in json format')
                self._prob_contours = raw_input( ' : ' ).strip( )
                print('')
                with open(self._prob_contours) as data_file:
                    data = json.load(data_file)                    
            except IOError as io_error:
                print (io_error)
            else:
                break

        contour_pieces = len(data['contours'])

        percentile = ('10-percentile','20-percentile','30-percentile','40-percentile',
                      '50-percentile','60-percentile', '70-percentile','80-percentile',
                      '90-percentile')

        for percentile_json in percentile:
            aladin.draw_newtool (percentile_json) # plotting in different Aladin planes
            self.plot_prob_contours(data, contour_pieces, percentile_json)

    def fov_shape(self):
        """Choosing the FoV shape: box (1) or circle (2)."""

        while True:
            try:
                print ('--> FoV shape: type (1) for a box and (2) for a circle')
                shape = input( ' : ' )
            
                # loading footprint template from Git--> box or circle
                if shape==1: 
                    url_id = "https://raw.githubusercontent.com/ggreco77/GWsky/master/footprint_box"
                    template_fov_footprint = download_file(url_id, cache=True, timeout=300)
                    self._fov_box()
                    self._make_fov_footprint(template_fov_footprint)
                elif shape==2:
                    url_id = "https://raw.githubusercontent.com/ggreco77/GWsky/master/footprint_circle"
                    template_fov_footprint = download_file(url_id, cache=True, timeout=300)
                    self._fov_circle()
                    self._make_fov_footprint(template_fov_footprint)                 
                else:
                    continue
            except SyntaxError as io_error:
                    print (io_error)
            except NameError as value_error:
                    print (value_error)
            else:
                break
                 
    def _fov_box(self):
        """Inserting the instrument FoV size (width and height) in degrees."""
        
        while True:
            try:
                print(' ')
                print ('--> Insert the size of your Field of View (FoV)')
                self._fov_width, self._fov_height = input( ' width [deg], height [deg] : ' )
                print('')
            except SyntaxError as syntax_error:
               print(syntax_error)
            except TypeError as type_error:
               print(type_error)
            except NameError as name_error:
               print(name_error) 
            except ValueError as value_error:
               print(value_error) 
            else:
               break

        return self.GWsky_config.update({"fov_width" : self._fov_width,
                                       "fov_height" : self._fov_height,
                                       "fov_shape" : 1})

    def _fov_circle(self):
        """Inserting the instrument FoV radius in degrees."""
        
        while True:
            try:
                print(' ')
                print ('--> Insert the radius of your Field of View (FoV)')
                self._fov_radius = input( ' radius [deg] : ' )
                print('')
            except SyntaxError as syntax_error:
               print(syntax_error)
            except TypeError as type_error:
               print(type_error)
            except NameError as name_error:
               print(name_error) 
            except ValueError as value_error:
               print(value_error) 
            else:
               break
    
        return self.GWsky_config.update({"fov_width" : self._fov_radius,
                                       "fov_height" : self._fov_radius,
                                       "fov_radius" : self._fov_radius,
                                       "fov_shape" : 2})
    
    def _make_fov_footprint(self, footprint):
        """Making the instrument field of view footprint using a template
                  from http://aladin.u-strasbg.fr/footprint_editor/"""

        votable = parse(footprint) # reading footprint template
        table = votable.get_first_table()

        for param in table.params: # retrieving table.params
            param

        # box or circle footprint
        if param.ID == 'radius':
            param.value = self._fov_radius*3600.0
            votable.to_xml('user_fov.vot') # VOTable file output
        else:
            data = table.array      
            fov_width_arcsec = self._fov_width*3600.0  
            fov_height_arcsec = self._fov_height*3600.0
    
            data[0] = - fov_width_arcsec /  2.0,   fov_height_arcsec / 2.0
            data[1] =   fov_width_arcsec /  2.0,   fov_height_arcsec / 2.0
            data[2] =   fov_width_arcsec /  2.0, - fov_height_arcsec / 2.0
            data[3] = - fov_width_arcsec /  2.0, - fov_height_arcsec / 2.0
            votable.to_xml('user_fov.vot') # VOTable file output       

        return samp.send_file( 'user_fov.vot' ) # sending to Aladin
          
    def telescope_site(self):
        """Inserting the Geodetic coordinates of the Observatory."""
        
        while True:
            try:
                self._latitude, self._longitude, self._altitude =input(
                    '--> Insert the Geodetic coordinates of the Observatory \n latitude[deg], longitude[deg], altitude[m]: ' )
                print('')
            except SyntaxError as syntax_error:
                print (syntax_error)
            except TypeError as type_error:
                print (type_error)
            except NameError as name_error:
                print (name_error)
            except ValueError as value_error:
                print (value_error) 
            else:
                break

        return self.GWsky_config.update({"latitude" : self._latitude,
                                           "longitude" : self._longitude,
                                           "altitude" : self._altitude })

    def starting_time(self):
        """Inserting the observation time yyyy-mm-gg hh:mm:ss."""
        
        while True:
            try:
                print ('--> Insert the observation time ') 
                self._obs_time = raw_input( '  yyyy-mm-dd hh:mm:ss : ' ).strip()
                print('')
                time = Time(self._obs_time)
            except ValueError as value_error:
                print (value_error)
            else:
                break
            
        return self.GWsky_config.update({"obs_time" : self._obs_time})

    def GWsky_basic(self):
        """GWsky basic version: no statistic window"""
        
        while True:
            try:
                print ('--> Type "b" for GWsky basic (no statistic window) ')
                print ('--> otherwise any key ')
                self.GWsky_basic = raw_input( ' : ' ).strip()
                print('')
            except ValueError as value_error:
                print (value_error)
            else:
                break
            
        return self.GWsky_config.update({"GWsky_basic" : self.GWsky_basic})

    def _find_highest_pixel(self, infile):
        """Finding the sky position of the highest probability pixel."""

        hpx = hp.read_map(infile, verbose = False)
        npix = len(hpx)
        nside = hp.npix2nside(npix) # resolution for the HEALPix map

        ipix_max = np.argmax(hpx)# pixel position
        hpx[ipix_max] 
        theta, phi = hp.pix2ang(nside, ipix_max)
        
        ra_max = np.rad2deg(phi) # sky coordinates
        dec_max = np.rad2deg(0.5 * np.pi - theta)

        return round(ra_max,5), round(dec_max,5)

    def maximum_pixel(self):
        """Return the sky position of the maximum probability pixel"""
        
        self._ra_max, self._dec_max = self._find_highest_pixel(self._skymap)
        
        print (' The highest probability pixel is located at ')
        print (' RA =' + str('% .5f' % self._ra_max)+'°' + 'and Dec =' + str('% .5f' % self._dec_max)+'°.')

        return self.GWsky_config.update({"ra_max_pixel" : self._ra_max,    
                                           "dec_min_pixel" : self._dec_max})

    def sky_position(self):
        """Init. sky coordinates in "GWsky_coords" file."""

        self.GWsky_coords.update({"input_ra" : self._ra_max,
                                  "input_dec" : self._dec_max})

        with open('GWsky_coords', 'wb') as data:
            return pickle.dump(self.GWsky_coords, data)

    def site(self):
        """Creating the configuration file "GWsky_config". """

        with open('GWsky_config', 'wb') as data:
            return pickle.dump(self.GWsky_config, data)                        

    def display_note(self):
        print(' NOTE: the default catalog is GLADE [VII/275/glade1] ')
        print('')
        print(' You can change the default catalog at the line 312 in **coverage** module ')
        print (' e.g from GLADE to GWGC --> "VII/275/glade1" --> "VII/267" ')
        print('')
        print (' You can also change the catalog statistics ')
        print (' LINE 792: query_catalog.hist(table_pandas.Dist.dropna() --> ')
        print (' query_catalog.hist(table_pandas.**ColumnName**.dropna() ')
        print ('')
        print(' LINE 797: query_catalog.hist(table_pandas.Bmag.dropna()) --> ')
        print (' query_catalog.hist(table_pandas.**ColumnName**.dropna()) ')
        print(' ')
        print (' Launch the GUI by typing: ')
        print (' from GWsky import coverage ')

        
# running module
user_input = UserInput()

user_input.display_welcome()
user_input.healpix_skymap()
user_input.update_GWsky_config()
user_input.display_healpix_skymap()
user_input.probability_contours()
user_input.fov_shape()
user_input.telescope_site()
user_input.starting_time()
user_input.GWsky_basic()
user_input.maximum_pixel()
user_input.sky_position()
user_input.site()
aladin.draw_newtool('FoVcenters') # new plane: FoV sky coords
user_input.display_note()
