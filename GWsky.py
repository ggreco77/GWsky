#! /usr/bin/env python
# -*- coding: utf-8 -*-  

#-------------------------GWsky (v2)----------------------------------------
#                                                                          |
#  This interactive script is free software: you can redistribute it       |
#  and/or modify it under the terms of the BSD license.                    |
#                                                                          |
#     This script is distributed in the hope that it will be useful,       |
#     but WITHOUT ANY WARRANTY; without even the implied warranty of       |
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                 |
#            See BSD license for more details.                             |
#                                                                          |
#---------------------------------------------------------------------------

#-------------------Launching the Script------------------------------------
#                        GWsky (v2)                                        |
#---------------------------------------------------------------------------
#                                                                          |
# from idle:                                                               |
#    >>> import GWsky                                                      |
#    >>> GWsky.main()                                                      |
#                                                                          |
# from terminal: ./GWsky                                                   |
# if ./GWsky: Permission denied; type: chmod u+x GWsky                     |
#---------------------------------------------------------------------------

__author__ = "Giuseppe Greco <giuseppe.greco@uniurb.it>"
 

def main():

	   
     """

     The interactive script GWsky (v2) defines a sequence of Fields of View (FoV)
     from a fixed position over the sky.

     North/South/East/West directions are allowed.

     The results are displayed in Aladin Sky Atlas using the |SAMPIntegratedClient| class
     (http://aladin.u-strasbg.fr/).

     The airmass at the FoV center and the integrated probability (%)
     are provided during the FoV sequence.

     ***Running it***
     from idle:
          >>> import GWsky
          >>> GWsky.main()
          
     from terminal: ./GWsky
     if ./GWsky: Permission denied; type chmod u+x GWsky  


     SUMMARY OF DEPENDENCIES:
     -----------------------
     
     from astropy.vo.samp import SAMPIntegratedClient
     from astropy.vo.samp import SAMPHubError
     from  astropy.time import Time
     from astroquery.vizier import Vizier
     import astropy.coordinates as coord
     import astropy.units as u
     import healpy 
     import urlparse
     import os.path
     from astropy.io.votable import parse
     from math import sin, cos, acos, degrees, radians
     import healpy
     import numpy
     import time                                                                   
     import sys
     import pickle
     import time

     Input Parameters :
     ------------------

     sky_map : str
              LVC probability skymap in healpix format.

     percentage : float
             probability percentage to determine the area (in square degrees) confined in it.

     FOV_base, FOV_height : float
              size of Field of View . 

     time_input : str
             the time in the format "2012-7-12 23:00:00"

     lat_input : float
             Geodetic coordinates of the Observatory: latitude (deg).

     lon_input : float
             Geodetic coordinates of the Observatory: longitude (deg).

     height_input : float
             Geodetic coordinates of the Observatory: altitude (m).

     ra : float
            right ascention of FoV center (deg).

     dec : float
           declination of FoV center (deg).


     Return :
     --------

     area_probability : float
          the area (in square degrees) confined in a given probability percentage.

     percentage_poly : float
          the integrated probability in a specific FOV (%).

     airmass : float
          the airmass of the FOV center.

     contour_ipix.out: list
          the table that contained the pixels confined
               in a given probability percentage.

     ra_max : float
          right ascention of the highest probability pixel.

     dec_max : float
          declination of the highest probability pixel.

     instrument_FOV.vot : VOTABLE
          Instrument Footprint Editor from
               http://aladin.u-strasbg.fr/footprint_editor/
            
     N/S/E/W/R/Q : str
          N/S/E/W: a set of command line to add a contiguous FOVs in North/South/East/West  (N/S/E/W) directions;
          R: to insert a new FOV center RA[deg], DEC[deg] to begin a new sequence in N/S/E/W directions;
          Q: quit
     
     """
     
     import time

     import pickle
     
     # External Dependencies
     from astropy.vo.samp import SAMPHubError

     from  astropy.time import Time
     
     import healpy as hp

     # GWsky functions
     
     import aladinSAMP
     
     from print_area_prob import print_area_prob

     from table_ipix_percentage import table_ipix_percentage

     from highest_probability_pixel import highest_probability_pixel

     from instrument_FOV import instrument_FOV

     from gw_sequence import add_FOV

     from progress_bar import progress_bar
     
     
     # load a probability Healpix sky map
     while True: 
         try:
              print ' Load a probability skymap; '
              sky_map = raw_input( ' HEALPix format: ' ).strip( )
              hpx_test = hp.read_map( sky_map,verbose=False )
         except IOError as io_error:
               print ''
               print '    -------------------------------------------------------------'
               print '               ***Oops!  No sky map has been found.***'
               print '    -------------------------------------------------------------'
               print ''
               print '', io_error
               print ''
         except ValueError as value_error:  
               print ''
               print '    -------------------------------------------------------------'
               print '               ***Oops!  No sky map has been found.***'
               print '    -------------------------------------------------------------'
               print ''
               print '', value_error
         else:
              break
          
     progress_bar()
     
     # SAMPIntegratedClient
     while True:     
          try:
               aladinSAMP.send_file( sky_map )
          except SAMPHubError as samphub_error:
               print ''
               print '    -------------------------------------------------------------'
               print '       ***Open Aladin Sky Atlas before running the script***'
               print '                    http://aladin.u-strasbg.fr/'
               print '    -------------------------------------------------------------'
               print ''
               print'', samphub_error
               print ''
               time.sleep ( 5 )
          else:
               break

     # probability percentage 
     while True:
          try:
               print(' Calculate the area confined in a given probability percentage; ')
               percentage = input(' (0 < p > 1) : ')
               while True:
                    if percentage > 1 or percentage < 0:
                         print ''
                         percentage = input( 'Please, insert a number from 0 to 1: ' )
                         print ''
                    else:
                         break
          except SyntaxError as syntax_error:
               print ''
               print ''
               print ' Oops: try again! ', syntax_error
               print ''
               print ''
          except NameError as name_error:
               print ''
               print ''
               print ' Oops: try again! ', name_error
               print ''
               print ''
          else:
               break

     # print the area confined in a given percentage of probability.
     print_area_prob( sky_map, percentage )

     print ' The table that contained those pixels is displayed in Aladin  \n < contour_ipix.out >'

     # build the table that contained those pixels 
     table_ipix_percentage( sky_map, percentage )

     # send to Aladin the contour_ipix.out table
     aladinSAMP.send_file( 'contour_ipix.out' )


     progress_bar()

     # insert the FOV size
     while True:     
          try:          
               print ' Insert the size of your Field of View (FoV); ' 
               FOV_base, FOV_height = input( ' x[deg], y[deg] : ' )
          except SyntaxError as syntax_error:
               print ''
               print '    -------------------------------------------------------------'
               print '    ***Please insert the size of your Field Of View in degrees*** '
               print '                   e.g. for a 3° x 2° FOV: 3,2    '
               print '    -------------------------------------------------------------'
               print '', syntax_error
          except TypeError as type_error:
               print ''
               print '    -------------------------------------------------------------'
               print '    ***Please insert the size of your Field Of View in degrees*** '
               print '                   e.g. for a 3° x 2° FOV: 3,2    '
               print '    -------------------------------------------------------------'
               print '', type_error
          except NameError as name_error:
               print ''
               print '    -------------------------------------------------------------'
               print '    ***Please insert the size of your Field Of View in degrees*** '
               print '                   e.g. for a 3° x 2° FOV: 3,2    '
               print '    -------------------------------------------------------------'
               print '', name_error
          else:
               break


     progress_bar()

     #catalogue VizieR Queries
     print 'Specify the ID of a catalog for VizieR Query; '
     print ' e.g. in the case of Gravitational Wave Galaxy Catalog type: VII/267/gwgc '
     catalog = raw_input (' : ')
     
     progress_bar()

     # modify the template VOTable of Instrument Footprint Editor
     instrument_FOV( FOV_base, FOV_height )

     # INPUT values for airmass calculation
     while True:
          try:
               print ' For airmass calculation insert the observation time ' 
               time_input = raw_input( '  2012-7-12 23:00:00 : ' ).strip()
               time = Time( time_input )
          except ValueError as value_error:
               print ''
               print ''
               print 'Please insert the time in the format  "2012-7-12 23:00:00" ', value_error
               print ''
               print ''
          else:
               break
     
     while True:
          try:    
               lat_input, lon_input, height_input = input( ' and the Geodetic coordinates of the Observatory \n latitude[deg], longitude[deg], altitude[m]: ' )
          except SyntaxError:
               print ''
               print '    -------------------------------------------------------------'
               print '    ***Please insert the Geodetic coordinates of the Observatory*** '
               print '        e.g. for Cerro Paranal Observatory: -24.63, -70.40, 2635     '
               print '    -------------------------------------------------------------'
               print ''
          except TypeError:
               print ''
               print '    -------------------------------------------------------------'
               print '    ***Please insert the Geodetic coordinates of the Observatory*** '
               print '        e.g. for Cerro Paranal Observatory: -24.63, -70.40, 2635     '
               print '    -------------------------------------------------------------'
               print ''
          except NameError as name_error:
               print ''
               print '    -------------------------------------------------------------'
               print '    ***Please insert the Geodetic coordinates of the Observatory*** '
               print '        e.g. for Cerro Paranal Observatory: -24.63, -70.40, 2635     '
               print '    -------------------------------------------------------------'
               print ''
          else:
               break

     progress_bar()

# display a possible EM candidate
     while True:
          print ' Do you want to display a possible EM candidate? '
          print ' Type Y to add it or any key to exit '
          option = raw_input( ' : ' ).strip()

          if option.strip().capitalize() == 'Y':
               while True:
                    try:
                         print ' Insert an ID '
                         ID = raw_input ( ' : ' )
          
                         print ' and the sky position in degrees; RA[deg], DEC[deg] '
                         ra_candidate, dec_candidate = input( ' : ' )
                    except SyntaxError as syntax_error:
                         print ''
                         print 'Oops, try again', syntax_error
                         print
                    except TypeError as type_error:
                         print ''
                         print 'Oops, try again', type_error
                         print
                    except ValueError as value_error:
                         print ''
                         print 'Oops, try again', value_error
                         print
                    except NameError as name_error:
                         print ''
                         print 'Oops, try again', name_error
                         print
                    else:
                         break

               # building command script for Aladin
               position_candidate = [ ra_candidate, dec_candidate ]
               position_candidate = ' , '.join( map ( str, position_candidate ) )
               EM_candidate = 'draw string' + '( ' + position_candidate + ', ' + ID + ')'
               aladinSAMP.send_script( EM_candidate )
          
               size = '10arcmin'
               EM_candidate_tag = 'draw red circle' + '( ' + position_candidate + ', ' + size + ')'
               aladinSAMP.send_script( EM_candidate_tag )

               progress_bar()
               print''        
          else:
              break
     

     progress_bar()


     print ' Specify a sky position to generate a FOV sequence,'

     # find the highest probability pixel RA[deg], DEC[deg]
     highest_probability_pixel( sky_map )

     print ''

     while True:
          try:
               ra, dec = input( ' RA[deg], DEC[deg]: ' )
          except SyntaxError:
               print ''
               print '    -------------------------------------------------------------'
               print '    ***Please specify a sky position to generate a FOV sequence*** '
               print '                      e.g. -24.63, -70.40   '
               print '    -------------------------------------------------------------'
               print ''
          except TypeError:
               print ''
               print '    -------------------------------------------------------------'
               print '    ***Please specify a sky position to generate a FOV sequence*** '
               print '                      e.g. -24.63, -70.40   '
               print '    -------------------------------------------------------------'
               print ''
          except NameError:
               print ''
               print '    -------------------------------------------------------------'
               print '    ***Please specify a sky position to generate a FOV sequence*** '
               print '                      e.g. -24.63, -70.40   '
               print '    -------------------------------------------------------------'
               print ''
          except ValueError:
               print ''
               print '    -------------------------------------------------------------'
               print '    ***Please specify a sky position to generate a FOV sequence*** '
               print '                      e.g. -24.63, -70.40   '
               print '    -------------------------------------------------------------'
               print ''
          else:
               break

     print ''


     #-------------------------------------------------------------
    
     # write the input parameters in a config file "config_GWsky"
     config_GWsky = {
                     'sky_map' : 'sky_map',
                     'percentage' : 'percentage', 
                     'FOV_base' : 'FOV_base',
                     'FOV_height' : 'FOV_height',
                     'catalog' : 'catalog',
                     'time_input' : 'time_input',
                     'lat_input' : 'lat_input',
                     'lon_input' : 'lon_input',
                     'height_input' : 'height_input',
                     'ra' : 'ra',
                     'dec' : 'dec'
                     }
     
     config_GWsky [ 'sky_map' ] = sky_map
     config_GWsky [ 'percentage' ] = percentage
     config_GWsky [ 'FOV_base' ] = FOV_base
     config_GWsky [ 'FOV_height' ] = FOV_height
     config_GWsky [ 'catalog' ] = catalog
     config_GWsky [ 'time_input' ] = time_input
     config_GWsky [ 'lat_input' ] = lat_input
     config_GWsky [ 'lon_input' ] = lon_input
     config_GWsky [ 'height_input' ] = height_input
     config_GWsky [ 'ra' ] = ra
     config_GWsky [ 'dec' ] = dec
     
     with open ( 'config_GWsky.ini', 'wb' ) as data:
               pickle.dump ( config_GWsky, data )


     # sequence of FOVs from a predetermined sky position
     add_FOV ( sky_map, FOV_base, FOV_height, ra, dec )
