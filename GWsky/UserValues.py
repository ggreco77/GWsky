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
from mocpy import MOC

import json
import time
import pickle

from astropy.vo.samp import SAMPHubError # no work in py3
from  astropy.time import Time
from astropy.io.votable import parse

from aladinSAMP import AladinViaSAMP, AladinScriptCommands # class access
samp = AladinViaSAMP()
aladin = AladinScriptCommands()


class SkymapCollectionMethods(object):
    """Some utilities for working with a probability skymap.

     Methods :
              MOC_confidence_region : to extract a multi-order coverage map
              (MOC) at a selected probability level.

              query_catalogs : query a list of Vizier catalog(s) in
              a MOC region (default: 90%).

              find_highest_pixel : find the highest probability pixel.
              """


    def MOC_confidence_region(self, infile, percentage, short_name):
        """
        Multi-Order coverage map (MOC) of sky area enclosed within a contour plot
        at a given confidence level.
    
        Input:
             infile: healpix format
                     LVC probability sky map
             percentage: float
                      probability percentage of the enclosed area  
             short_name: str
                     output file name
     
        Output: fits format
                     MOC map named "short_name"_"percentage"      
        """
        
        hpx = hp.read_map(infile, verbose = False) # reading skymap
        npix = len(hpx)
        nside = hp.npix2nside(npix)
 
        sort = sorted(hpx, reverse = True)
        cumsum = np.cumsum(sort)
        index, value = min(enumerate(cumsum),
                           key = lambda x: abs(x[1] - percentage))
 
        index_hpx = range(0, len(hpx)) # finding ipix indices confined 
        hpx_index = np.c_[hpx, index_hpx]   #--> in a given percentage

        sort_2array = sorted(hpx_index,
                             key = lambda x: x[0], reverse = True)
        value_contour = sort_2array[0:index]

        j = 1 
        table_ipix_contour = [ ]

        for i in range (0, len(value_contour)):
            ipix_contour = int(value_contour[i][j])
            table_ipix_contour.append(ipix_contour)
          
        theta, phi = hp.pix2ang(nside, table_ipix_contour) # from index to polar coordinates
        ra = np.rad2deg(phi) 
        dec = np.rad2deg(0.5 * np.pi - theta)
 
        from astropy.table import Table # creating an astropy.table:
                                            #--> contour_ipix      
        contour_ipix = Table([ra, dec], names = ('RA[deg]', 'DEC[deg]'), 
                         meta = {'ipix': 'ipix table'})     
    
        from math import log
        moc_order = int(log(nside, 2)) # setting MOC order

        # creating a MOC map from the contour_ipix table
        moc = MOC.from_table(contour_ipix, 'RA[deg]',
                             'DEC[deg]', moc_order)

        moc_file = short_name + '_MOC_' + str(percentage) 
        moc.write(moc_file, format = 'fits') # writing MOC file in fits

        return moc_file
        
    def query_catalogs(self, moc_file):
        """Query a list of selected Vizier catalog(s) in 90% MOC
                    (Multi-Order Coverage) region."""

        catalogs = ['VII/267/gwgc'] # selected catalog(s)

        moc = MOC.from_file(moc_file)  # loading MOC coverage
        
        import warnings # querying from MOC ignoring 
        with warnings.catch_warnings():     #--> astropy.io.votable.exceptions
            warnings.simplefilter("ignore")
            for catalog in catalogs:
                print ("")
                table = moc.query_vizier_table(catalog,
                                               max_rows = 2000000) # set rows limits
                #print (table)
                catalog_renamed = catalog.replace('/', '_')
                table.write(catalog_renamed + 'MOC_query',
                            format = 'votable', overwrite = True)
                samp.send_file(catalog_renamed + 'MOC_query')

    def find_highest_pixel(self, infile):
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
        
    
class UserInput(SkymapCollectionMethods):
    """Insert at the command prompt user values to create the configuration
       file **GWskyConfigure**.
       The coverage** module is initialized by the GWskyConfigure."""

    def __init__(self):
        """
            skymap : a valid LVC skymap in healpix format;
            prob_contours : a valid probability contours file in json format;
            fov_width : field-of-view width in degrees;
            fov_height : field-of-view height in degrees;
            latitude : Geodetic coordinates of the Observatory [deg];
            longitude :              "                         [deg];
            altitude :               "                          [m];
            obs_time : starting time yyyy-mm-dd hh:mm:ss;
            input_ra : right ascension in degrees;
            input_dec : declination in degrees;
            GWskyConfigure : dict to store the user values;
            percentage : probability percentage of the enclosed area.
        """
        
        self._skymap = str()
        self._prob_contours = str()       
        self._fov_width = float()
        self._fov_height = float()
        self._latitude = float()
        self._longitude = float()
        self._altitude = float()
        self._obs_time = str()
        self._input_ra = float()
        self._input_dec = float()
        self.GWskyConfigure = dict()
        self.percentage = 0.9

    def healpix_skymap(self):
        """Inserting a valid LIGO/Virgo probability skymap."""
        
        while True:
            try:
                print ('--> Load a probability skymap')
                self._skymap = raw_input( ' : ' ).strip( )
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
        
        return self.GWskyConfigure.update({"skymap" : self._skymap,
                                           "prob" : prob,
                                           "nside" : nside})

    def is_3d_skymap(self):
        """Check if the healpix skymap is a 3d-skymap and read all columns."""

        prob, header = hp.read_map(self._skymap, h=True, verbose=False)
        header = dict(header)
        
        if header['TFIELDS']==4:
            distmu, distsigma, distnorm = hp.read_map(self._skymap, field=[1, 2, 3])
            
            return self.GWskyConfigure.update({"distmu" : distmu,
                                    "distsigma" : distsigma,
                                    "distnorm" : distnorm})
        else:
            return self.GWskyConfigure.update({"distmu" : [ ],
                                               "distsigma" : [ ],
                                               "distnorm" : [ ]})
        
    def display_healpix_skymap(self):
        """Displaying the healpix skymap in Aladin plane."""

        while True:
            try:
                return samp.send_file(self._skymap)
            except SAMPHubError as samphub_error:
                print(samphub_error, 
                      ' launch Aladin sky Atlas to run the script ' 
                      ' http://aladin.u-strasbg.fr/java/nph-aladin.pl '
                      ' and check the SAMP connection:'
                      ' Aladin bar --> Interop --> Connect with SAMP. ')
                print('')
                time.sleep(5)
            else:
                break

    def catalogs(self):
        """Loading a list of Vizier catalogs in 90% confidence region in MOC
           (Multi-Order Coverage) representation."""

        print(" Loading selected catalogs in the MOC region (90%): GWGC ")

        self.MOC_confidence_region(self._skymap, self.percentage,
                                   short_name = self._skymap)
        moc_region = self._skymap + '_MOC_' + str(self.percentage)
        
        return self.query_catalogs(moc_region)

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
           
    def fov_size(self):
        """Inserting the instrument FoV size (width and height) in degrees."""
        
        while True:
            try:
                print ('--> Insert the size of your Field of View (FoV)')
                self._fov_width, self._fov_height = input( ' width[deg], height[deg] : ' )
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

        return self.GWskyConfigure.update({"fov_width" : self._fov_width,
                                           "fov_height" : self._fov_height})
    
    def make_FoV_footprint(self):
        """Making the instrument field of view footprint using a generic template
                  from http://aladin.u-strasbg.fr/footprint_editor/"""
                  
        # loading footprint template from Git
        from astropy.utils.data import download_file
        url_id = "https://raw.githubusercontent.com/ggreco77/GWsky/master/template_fov_footprint"
        template_fov_footprint = download_file(url_id, cache=True, timeout=300)

        votable = parse(template_fov_footprint) # reading the Instrument 
        table = votable.get_first_table()               #--> Footprint file template
        data = table.array      

        fov_width_arcsec = self._fov_width*3600.0  # from degree to arcsecond
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

        return self.GWskyConfigure.update({"latitude" : self._latitude,
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
            
        return self.GWskyConfigure.update({"obs_time" : self._obs_time})

    def maximum_pixel(self):
        """Return the sky position of the maximum probability pixel"""
        
        ra_max, dec_max = self.find_highest_pixel(self._skymap)
        
        print (' The highest probability pixel is located  at RA =' 
               + str('% .5f' % ra_max)+'°' + 'and Dec =' + str('% .5f' % dec_max)+'°.')

        return self.GWskyConfigure.update({"input_ra" : ra_max,    
                                           "input_dec" : dec_max})

    def sky_position(self):
        """Inserting a sky position in degrees: RA[deg], DEC[deg]."""
        
        while True:
            try:
                print ('--> Specify a sky position (in degrees) to generate a FOV sequence.')
                self._input_ra, self._input_dec = input(' RA[deg], DEC[deg]: ')
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
    
        return self.GWskyConfigure.update({"input_ra" : self._input_ra,    
                                           "input_dec" : self._input_dec})

    def display_fov_user(self):
        """Displaying the FoV footprint at the user sky position."""
        
        return aladin.get_FoV(self._input_ra, self._input_dec)   
   
    def make_configure_file(self):
        """Creating the configuration file "GWskyConfigure". """
        
        with open('GWskyConfigure', 'wb') as data:
            return pickle.dump(self.GWskyConfigure, data)

    def show_configure_file(self):
        """Showing the values in the configuration file GWskyConfigure."""
        
        with open('GWskyConfigure', 'rb') as data:
            GWskyConfigure = pickle.load(data)
            
        for k, v in sorted(GWskyConfigure.items()):
            print ('GWskyConfigure:')
            print str(k) + (' : ') + str(v)
            print('')
            

# running module
user_input = UserInput()
user_input.healpix_skymap()
user_input.display_healpix_skymap()
user_input.is_3d_skymap()
#user_input.catalogs()
user_input.probability_contours()
user_input.fov_size()
user_input.make_FoV_footprint()
user_input.telescope_site()
user_input.starting_time()
user_input.maximum_pixel()
#user_input.sky_position() 
#user_input.display_fov_user()
user_input.make_configure_file()
user_input.show_configure_file()
aladin.draw_newtool('FoVcenters') # new plane: FoV sky coords
