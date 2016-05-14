#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-  
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
     North/NorthWest/West/SouthWest/South/SouthEast/East/NorthEast directions are allowed.

     """
     
     import time
     import pickle
     
     # External Dependencies

     from astropy.vo.samp import SAMPHubError
     from  astropy.time import Time
     import healpy as hp

     # GWsky functions
     
     import aladinSAMP
     import aladin_console
     from print_area_prob import print_area_prob
     from table_ipix_percentage import table_ipix_percentage
     from find_highest_probability_pixel import find_highest_probability_pixel
     from instrument_FOV import instrument_FOV
     from FOV_sequence import add_FOV
     from progress_bars import progress_bar
     
     # load a probability healpix sky map
     while True: 
         try:
              print ' Load a probability skymap '
              input_skymap = raw_input( ' : ' ).strip( )
              hpx_test = hp.read_map( input_skymap, verbose = False )
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
     
     # sending the healpix skymap to Aladin plane
     while True:     
          try:
               aladinSAMP.send_file( input_skymap )
          except SAMPHubError as samphub_error:
               print ''
               print '    -------------------------------------------------------------'
               print '       ***Launch Aladin Sky Atlas for running the script***'
               print '                    http://aladin.u-strasbg.fr/'
               print '    -------------------------------------------------------------'
               print ''
               print'', samphub_error
               print ''
               time.sleep ( 5 )
          else:
               break


     # inserting a value of probability percentage 
     while True:
          try:
               print( ' Provide the probability value (< 1) to define the confidence region; ' )
               print ( '(e.g. 1 for all the sky, 0.5 to manage only the 50% confidence region)' )
               percentage = input(' : ')
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

     # printing the area confined in a given percentage of probability.
     print_area_prob( input_skymap, percentage )
     print ''
     print ' The table containing the pixels is displayed in Aladin plane  \n < contour_ipix.out >'

     # building the table that contained those pixels 
     table_ipix_percentage( input_skymap, percentage )

     # sending to Aladin plane the contour_ipix.out table
     aladinSAMP.send_file( 'contour_ipix.out' )


     progress_bar()

     # inserting the FOV size
     while True:     
          try:          
               print ' Insert the size of your Field of View (FoV); ' 
               fov_width, fov_height = input( ' x[deg], y[deg] : ' )
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
          except ValueError as value_error:
               print ''
               print '    -------------------------------------------------------------'
               print '    ***Please insert the size of your Field Of View in degrees*** '
               print '                   e.g. for a 3° x 2° FOV: 3,2    '
               print '    -------------------------------------------------------------'
               print '', value_error
          else:
               break


     progress_bar()

     catalog = ''
     while catalog == '':
          print 'Specify the ID of a catalog for VizieR Query; '
          print ' e.g. in the case of Gravitational Wave Galaxy Catalog type: VII/267/gwgc '
          print ''
          print ' CDS Catalogues: http://cdsarc.u-strasbg.fr/cats/U.htx '
          catalog = raw_input(' : ').strip()
     
     
     # sending the selected catalog to Aladin plane
     aladin_console.get_VizieR( catalog )


     progress_bar()

     # modifing the template VOTable of Instrument Footprint Editor
     instrument_FOV( fov_width, fov_height )

     # inserting input values for airmass calculation
     while True:
          try:
               print ' For airmass calculation insert the observation time ' 
               input_time = raw_input( '  2012-7-12 23:00:00 : ' ).strip()
               time = Time( input_time )
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
               input_latitude, input_longitude, input_altitude = input( ' and the Geodetic coordinates of the Observatory \n latitude[deg], longitude[deg], altitude[m]: ' )
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
          except NameError:
               print ''
               print '    -------------------------------------------------------------'
               print '    ***Please insert the Geodetic coordinates of the Observatory*** '
               print '        e.g. for Cerro Paranal Observatory: -24.63, -70.40, 2635     '
               print '    -------------------------------------------------------------'
               print ''
          except ValueError:
               print ''
               print '    -------------------------------------------------------------'
               print '    ***Please insert the Geodetic coordinates of the Observatory*** '
               print '        e.g. for Cerro Paranal Observatory: -24.63, -70.40, 2635     '
               print '    -------------------------------------------------------------'
               print ''
          else:
               break

     progress_bar()

     # displaying a possible EM candidate
     while True:
          print ' Do you want to display a possible EM candidate? '
          print ' Type Y to add it or any key to exit '
          option = raw_input( ' : ' ).strip()

          if option.strip().capitalize() == 'Y':
               while True:
                    try:
                         print ' Insert an ID '
                         ID = raw_input ( ' : ' ).strip()
          
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

     	       # sending EM-candidate position to Aladin plane         
               aladin_console.draw_string( ra_candidate, dec_candidate, ID )

               # sending circle to Aladin plane - tag the EM-candidate position
               aladin_console.draw_circle( ra_candidate, dec_candidate, size = '10arcmin' )
     
               progress_bar()
               print''        
          else:
              break
          

     progress_bar()
     
     # specifying a sky position to generate a FOV sequence
     print ' Specify a sky position to generate a FOV sequence,'

     # the highest probability pixel RA[deg], DEC[deg] is suggested
     find_highest_probability_pixel( input_skymap )

     print ''

     while True:
          try:
               input_ra, input_dec = input( ' RA[deg], DEC[deg]: ' )
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

    
     # config file "config_GWsky"

     ids = [ 'input_skymap', 'percentage', 'fov_width', 'fov_height', 'catalog', 'input_time', 'input_latitude',
            'input_longitude', 'input_altitude', 'input_ra', 'input_dec', 'ra_n', 'dec_n', 'ra_s', 'dec_s', 'ra_e', 'dec_e',
            'ra_w', 'dec_w', 'ra_nw', 'dec_nw', 'ra_sw', 'dec_sw', 'ra_ne', 'dec_ne', 'ra_se', 'dec_se', 'ra_last', 'dec_last',
             'ra_last2', 'dec_last2' ]

     input_values =[ input_skymap, percentage, fov_width, fov_height, catalog, input_time, input_latitude, input_longitude,
                     input_altitude,input_ra, input_dec,input_ra, input_dec,input_ra, input_dec,input_ra, input_dec,input_ra,
                     input_dec, input_ra, input_dec, input_ra, input_dec, input_ra, input_dec, input_ra, input_dec,
                     input_ra, input_dec, input_ra, input_dec ]
     

     config_GWsky = dict ( zip( ids, input_values ) )
     
     with open( 'config_GWsky.ini', 'wb' ) as data:
               pickle.dump( config_GWsky, data )


     # building a sequence of FoVs from a predetermined sky position
               # e.g. from the highest probability pixel
     add_FOV( input_skymap, fov_width, fov_height, input_ra, input_dec )
