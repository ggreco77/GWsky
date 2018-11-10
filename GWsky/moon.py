from __future__ import print_function

import astropy
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_moon, get_sun
from astropy.time import TimeDelta
from astropy.time import Time

import numpy as np

from .config_values import UserValues
from .aladinSAMP import AladinViaSAMP, AladinScriptCommands

from .utils import Utils

class Moon(object):

    def __init__(self):
        """Getting user-values from config_values module."""
        
        self.user = UserValues()
        self.aladin = AladinScriptCommands()
        
        self.latitude = self.user.get_latitude()
        self.longitude = self.user.get_longitude()
        self.altitude = self.user.get_altitude()
        self.obs_time = Time(self.user.get_obs_time())

        self.dt = TimeDelta(7200.0, format='sec')

        self.step = 0
        self.end_step = 11

    def get_location(self):

        observatory = astropy.coordinates.EarthLocation(lat = self.latitude*u.deg,
                                                        lon = self.longitude*u.deg,
                                                        height = self.altitude*u.m)
        return observatory

    def get_time(self):

        self.time = Time(self.obs_time)

        return self.time


    def moon_on_sky(self):

        self.get_location()
        #self.get.time()
        self.sky_position()

    def steps(self):
        """Moon position in step of one hour for an input sky position (ra, dec).
           10 steps are performed: step <10; dt = 1h."""

        #obs_time = Time(self.obs_time)
        self.aladin.draw_newtool("Moon")
        self.time = Time(self.obs_time)
        self.observatory = self.get_location()
        
        while self.step < self.end_step:
            time_update = self.time + self.step*self.dt
            
            position_moon = get_moon(time_update, self.observatory)
            
            #val = self.airmass(ra, dec, self.altitude, self.longitude, self.altitude,
            #                   time_input)           
            #self.airmass_list.append(val)
            #self.time_list.append(str(time_input))
            self.step+=1
            
            self.aladin.draw_string(position_moon.ra, position_moon.dec, "MOON"+ "-->" + str(time_update.isot))
            #print str(time_update.isot) 
            #print  position_moon.ra, position_moon.ra
            
    def illumination(self):
        """Return the fraction of the moon illumination.
           Modified version of astroplan project."""
        
        sun = get_sun(self.obs_time)
        observatory = self.get_location()
        moon = get_moon(self.obs_time, observatory)
        #print moon
        
        elongation = sun.separation(moon)
        i = np.arctan2(sun.distance*np.sin(elongation), moon.distance - sun.distance*np.cos(elongation))
        k = (1 + np.cos(i))/2.0
        
        return round(k.value, 2)

    def from_fov(self, ra_fov_center, dec_fov_center):
        """ Return the Moon position over the sky."""

        observatory = self.get_location()
        moon = get_moon(self.obs_time, observatory)
        
        distance = Utils.separation(ra_fov_center, dec_fov_center, moon.ra, moon.dec)
        #print moon.ra*u.deg, moon.dec*u.deg                
        return distance.deg

        

    def sky_position(self):
        """Plot the Moon position on the Aladin plane."""

        time = Time(self.obs_time)
        observatory = self.get_location()
        position_moon = get_moon(time, observatory)

        illumination = self.illumination()

        #self.aladin.draw_string(position_moon.ra, position_moon.dec, "MOON position")
        self.aladin.draw_moon(position_moon.ra, position_moon.dec, illumination)       



    
