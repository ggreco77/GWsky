# -*- coding: utf-8 -*-

def find_highest_probability_pixel( sky_map ):

     """

     finding the highest probability pixel RA[deg], DEC[deg].

     """
 
     
     import healpy as hp
     import numpy as np

     # read probability skymap
     hpx = hp.read_map( sky_map, verbose = False )

     # number of pixels
     npix = len( hpx )

     # nside: resolution for the HEALPix map
     nside = hp.npix2nside( npix )

     # pixel position
     ipix_max = np.argmax( hpx )
     hpx[ ipix_max ]

     # sky coordinate
     theta, phi = hp.pix2ang( nside, ipix_max )

     ra_max = np.rad2deg( phi )

     dec_max = np.rad2deg( 0.5 * np.pi - theta )

     print ' or insert the highest probability pixel located \n at RA =', str(( '% .5f' % ra_max))+'°', 'and Dec =',str(( '% .5f' % dec_max))+'°.' 
