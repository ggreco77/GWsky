
def airmass(ra,dec,lat_input,lon_input,height_input,time_input='',airmass_min=1,airmass_max=5.8):
     

     """

          Airmass calculation at a given time in a particular site.

     """

     

     import astropy
     from   astropy import units as u
     from  astropy.time import Time
     from  astropy.coordinates import SkyCoord, EarthLocation, AltAz

     # Geodetic coordinates of observatory
     observatory = astropy.coordinates.EarthLocation(lat=lat_input*u.deg,lon=lon_input*u.deg,height=height_input*u.m)


     # Sky coordinate of an ipix
          # (the maximum ipix in the script is suggested)
     max_prob_p = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')

     # Time object
     time = Time(time_input)
     
     # coordinate or frame in the Altitude-Azimuth system (Horizontal coordinates).
     altaz = max_prob_p.transform_to(AltAz(obstime=time,location=observatory))

     #1/cos(altaz)
     airmass_value = altaz.secz
     
     
     if airmass_value <= airmass_min or airmass_value >= airmass_max:
          print "---> The airmass of the FOV center is outside the range 1 < airmass < 5.8. <--- "

     else:     
          print "---> The airmass of the FOV center is ", str(round(airmass_value,2))+'. <---'

     


     
      
      
