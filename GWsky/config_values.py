#!/usr/bin/env python
# -*- coding: utf-8 -*-

try:
   import cPickle as pickle
except:
   import pickle


class UserValues(object):
    """
       Creating a class in which the instance attributes are based on the dictionary
       'GWsky_config' by default. GWsky_config is created and pickled by the 'user_values' *** ***
       module and will be deleted when the program will be closed.
    """
    
    def __init__(self, infile_config='GWsky_config'):      
        """
        'GWsky_config'  contains the keys below:
        
            self.skymap: a valid LVC skymap in healpix format;
            self.nside:  Resolution parameter of HEALPIX 
            Geodetic coordinates of the Observatory:
                     self.latitude:   latitude[deg];
                     self.longitude:  longitude[deg];
                     self.altitude:   altitude[m];
            self.catalog: catalog name in Vizier code
            self.obs_time: starting time [yyyy-mm-dd hh:mm:ss];
            self.fov_width: field-of-view width [deg] -->> if selected;
            self.fov_height: field-of-view height [deg] -->> if selected;
            self.fov_radius: field-of-view radius [deg] -->> if selected;
            self.fov_shape: (1) box or (2) circle;
            self.ra_max_pixel: right ascension of maximum probability pixel [deg];
            self.dec_max_pixel: declination of maximum probability pixel [deg];
            self.GWsky_basic: ('A') active statistic window;
                            : ('D') deactive statistic window;
            self.trasparency: window trasparency
            self.column_1: first catalog column selected by user;
            self.column_2: second catalog column selected by user; 
            self.filter_1: applied filter in colomn 1;
            self.filter_2: applied filter in colomn 2.               
        """
      
        self.infile_config = infile_config
        
        with open(self.infile_config, 'rb') as data:
            config_GWsky = pickle.load(data)

        for k, v in config_GWsky.items():
            setattr(self, k, v)

    def set_infile_config(self, new_infile_config):
        self.infile_config = new_infile_config

    def get_infile_config(self):
         return self.infile_config

    def set_skymap(self, new_skymap):
         self.skymap = new_skymap
    
    def get_skymap(self):
        return self.skymap

    def set_nside(self, new_nside):
        self.nside = new_nside
        
    def get_nside(self):
        return self.nside

    def set_latitude(self, new_latitude):
        self.latitude = new_latitude
        
    def get_latitude(self):
        return self.latitude

    def set_longitude(self, new_longitude):
        self.longitude = new_longitude
        
    def get_longitude(self):
        return self.longitude

    def set_altitude(self, new_altitude):
        self.altitude = new_altitude
    
    def get_altitude(self):
        return self.altitude

    def set_catalog(self, new_catalog):
        self.catalog = new_catalog
        
    def get_catalog(self):
        return self.catalog

    def set_obs_time(self, new_obs_time):
        self.obs_time = new_obs_time
        
    def get_obs_time(self):
        return self.obs_time

    def set_fov_width(self, new_fov_width):
        self.fov_width = new_fov_width
        
    def get_fov_width(self):
        return self.fov_width

    def set_fov_height(self, new_fov_height):
        self.fov_height = new_fov_height
        
    def get_fov_height(self):
        return self.fov_height

    def set_fov_radius(self, new_fov_radius):
         self.fov_radius = new_fov_radius
         
    def get_fov_radius(self):
        return self.fov_radius

    def set_ra_max_pixel(self, new_ra_max_pixel):
        self.ra_max_pixel = new_ra_max_pixel
        
    def get_ra_max_pixel(self):
        return self.ra_max_pixel

    def set_dec_max_pixel(self, new_dec_max_pixel):
        self.dec_max_pixel = new_dec_max_pixel
        
    def get_dec_max_pixel(self):
        return self.dec_max_pixel

    def set_GWsky_basic(self, new_GWsky_basic):
        self.GWsky_basic = new_GWsky_basic
        
    def get_GWsky_basic(self):
        return self.GWsky_basic

    def set_win_trasparency(self):
        return self.trasparency

    def get_win_trasparency(self):
        return self.trasparency

    def set_column_1(self, new_column_1):
        self.column_1 = new_column_1
        
    def get_column_1(self):
        return self.column_1  

    def set_column_2(self, new_column_2):
        self.column_2 = new_column_2
        
    def get_column_2(self):
        return self.column_2

    def set_get_filter_1(self, new_filter_1):
        self.filter_1 = new_filter_1
        
    def get_filter_1(self):
        return self.filter_1
    
    def set__filter_2(self, new_filter_2):
        self.filter_2 =  new_filter_2
        
    def get_filter_2(self):
        return self.filter_2

    def set_fov_shape(self, new_fov_shape):
        self.fov_shape = self.new_fov_shape

    def get_fov_shape(self):
       return self.fov_shape
           
    def __repr__(self): 
        with open(self.infile_config, 'rb') as data:
            config_GWsky = pickle.load(data)
        return str(config_GWsky)
