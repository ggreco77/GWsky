def table_ipix_percentage( infile, percentage ):
          
     """

     Saving in ascii format the Healpix pixel (ipix) confined in a given probability percentage.

     """

     import healpy as hp
     import numpy as np

     #reading skymap
     hpx = hp.read_map( infile, verbose = False )

     # number of pixels
     npix = len( hpx )

     # nside: resolution for the HEALPix map
     nside = hp.npix2nside( npix )
 
     # sorting probability array
     sort = sorted( hpx, reverse = True )

     # cumulative sum 
     cumsum = np.cumsum( sort )

     # ipix index in percentage
     index, value = min( enumerate( cumsum ), key = lambda x: abs( x[1] - percentage ) )

     # ----- find ipix index confined in  percentage ------------------------

     #ipix index in hpx array
     index_hpx = range( 0, len( hpx ) )
     
     # hpx and ipix index 2darray
     hpx_index = np.c_[ hpx, index_hpx ]

     # sorting 2d array hpx and ipix index
     sort_2array = sorted( hpx_index, key = lambda x: x[0], reverse = True )

     # finding ipix confined in  percentage
     value_contour = sort_2array[ 0:index ]

     # extracting ipix index and put in table_ipix_contour
     j = 1 
     table_ipix_contour = []

     for i in range ( 0, len( value_contour ) ):
          ipix_contour = int( value_contour[i][j] )
          table_ipix_contour.append( ipix_contour )

     # -----------------------------------------------------------------

     
     # from index to polar coordinates
     theta, phi = hp.pix2ang( nside, table_ipix_contour )

     # converting these to right ascension and declination in degrees
     ra = np.rad2deg( phi )
     dec = np.rad2deg( 0.5 * np.pi - theta )

     # saving table, ra and dec of ipix
     np.savetxt( 'contour_ipix.out', np.c_[ra,dec],delimiter=" ", fmt='%1.5f,%1.5f', header='RA[deg],DEC[deg]', comments=' ' )
