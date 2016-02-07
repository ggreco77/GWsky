
def airmass( ra, dec, lat_input, lon_input, height_input, time_input = '', airmass_min = 1, airmass_max = 6 ):
     

     """

          Airmass calculation at a given time in a particular site.

     """

     

     import astropy
     from   astropy import units as u
     from  astropy.time import Time
     from  astropy.coordinates import SkyCoord, EarthLocation, AltAz

     # Geodetic coordinates of observatory
     observatory = astropy.coordinates.EarthLocation( lat = lat_input*u.deg,lon = lon_input*u.deg,height = height_input*u.m )


     # Sky coordinates of an ipix
          # (the maximum ipix is suggested during the script)
     max_prob_p = SkyCoord( ra = ra*u.deg, dec=dec*u.deg, frame='icrs' )

     # Time object
     time = Time( time_input )
     
     # coordinate or frame in the Altitude-Azimuth system (Horizontal coordinates).
     altaz = max_prob_p.transform_to( AltAz( obstime=time, location=observatory ) )

     # 1/cos(altaz)
     airmass_value = altaz.secz
     
     
     if airmass_value <= airmass_min or airmass_value >= airmass_max:
           print '---'
          #print "---> The airmass of the FOV center is outside the range 1 < airmass < 6."

     else:     
          print "---> The airmass of the FOV center is ", str( round( airmass_value,2 ) )+'.'


def airmass_step(ra, dec, lat_input, lon_input, height_input):

     """

      Airmass calculation at a given time in a particular site
      in steps of one hour

     """
     
     print '-----------------------------------------------------------------------'
     print '----------------------Airmass Table------------------------------------'
     print '-----------------------------------------------------------------------'

     from astropy.time import TimeDelta
     from  astropy.time import Time
     import pickle
     from airmass import airmass

     # read parameters from config_GWsky.ini
     with open ( 'config_GWsky.ini', 'rb' ) as data:
          config_GWsky = pickle.load ( data )
    
     # reading time from the input
     time_user = config_GWsky[ 'input_time' ]

     dt = TimeDelta( 3600.0, format = 'sec' )
     time_user = Time( time_user )

     cont = 0
     while cont < 10:
          time_input = time_user + cont*dt
          print time_input,
          airmass( ra, dec, lat_input, lon_input, height_input, time_input, airmass_min = 1, airmass_max = 5.8 )
          cont+=1
     print '-------------------------------------------------------------------------'
