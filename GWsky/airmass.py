import astropy
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import TimeDelta
from astropy.time import Time

import dateutil

from config_values import UserValues

class Airmass(object):
    """Airmass calculation from 1 to 5.8 by default."""

    def __init__(self):
        """Getting user-values from config_values module."""
        
        self.user = UserValues()
        
        self.latitude = self.user.get_latitude()
        self.longitude = self.user.get_longitude()
        self.altitude = self.user.get_altitude()
        self.obs_time = self.user.get_obs_time()

        self.dt = TimeDelta(3600.0, format='sec') 
        
        self.airmass_list = []
        self.time_list = []
        self.step = 0
        self.end_step = 11
                     
    def airmass(self, ra, dec, lat, lon, height,
                time_input, airmass_min=1, airmass_max=5.8):
        """ Airmass calculation at a given time in a particular site for an input sky position (ra, dec) in degrees.
            The airmass is calculated in the range [airmass_min, airmass_max]."""
     
        observatory = astropy.coordinates.EarthLocation(lat = self.latitude*u.deg,
                                                        lon = self.longitude*u.deg,
                                                        height = self.altitude*u.m)
        sky_coord = SkyCoord(ra = ra*u.deg,
                              dec=dec*u.deg, frame='icrs')     
        time = Time(time_input) 
        altaz = sky_coord.transform_to(AltAz(obstime=time,
                                              location=observatory))                                                                             
        airmass_value = altaz.secz                  

        if airmass_value < airmass_min or airmass_value > airmass_max:
            airmass_value =  "nan"
            return airmass_value
        else:     
            airmass_value = round(float(airmass_value), 2)
            return airmass_value

    def airmass_step(self, ra, dec):
        """Airmass calculation in step of one hour for an input sky position (ra, dec) in degrees.
           10 steps are performed: step <10; dt = 1h."""

        obs_time = Time(self.obs_time)
        
        while self.step < self.end_step:
            time_input = obs_time + self.step*self.dt
            val = self.airmass(ra, dec, self.altitude, self.longitude, self.altitude,
                               time_input)           
            self.airmass_list.append(val)
            self.time_list.append(str(time_input))
            self.step+=1
            
        return self.airmass_list, self.time_list
