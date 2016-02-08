# -*- coding: utf-8 -*-

def towards_north( input_skymap, fov_width, fov_height, catalog,
                  input_latitude, input_longitude, input_altitude, init_ra, init_dec):

     import pickle
     import aladin_console
     import query_FOV
     from airmass import airmass_step
     from spherical_distance import shift_direction_nord_sud
     from progress_bars import print_separation_line, print_aladin_plane_x

     with open ( 'config_GWsky.ini', 'rb' ) as data:
          config_GWsky = pickle.load ( data )

     ra_n, dec_n = float( config_GWsky[ init_ra ] ), float( config_GWsky[ init_dec ] )

     # Undo "U"
     config_GWsky [ 'ra_last2' ] = ra_n
     config_GWsky [ 'dec_last2' ] =  dec_n
     
     shift_n = shift_direction_nord_sud( dec_n )
               
     pointing_n = ( dec_n + fov_height ) + shift_n
               
     # variable cycle
     dec_n = pointing_n

     [ v1_ra, v1_dec, v2_ra, v2_dec, v3_ra, v3_dec,
      v4_ra, v4_dec ] = vertices_FOV( ra_n, pointing_n, fov_width, fov_height )
               
     print '---> The new FOV is centered at RA = ', str(('% .5f' % ra_n))+'°', 'and DEC =', str(('%.5f' % pointing_n))+'°. <---'
     print ' '
               
     # vectors of ra and dec FOV 
     ra_vertices, dec_vertices = [ v1_ra, v2_ra,  v4_ra, v3_ra ], [ v1_dec, v2_dec, v4_dec, v3_dec ]
               
     prob_fov = probability_inside_FOV( input_skymap, ra_vertices, dec_vertices )

     print '---> The integrated probability in the FOV is', str(('% .4e' % prob_fov))+'%. <---'
     print ''

     # airmass in steps of one hour
     airmass_step( ra_n, pointing_n, input_latitude, input_longitude, input_altitude )
                      
     print_aladin_plane_x( )
               
     # sending FoV
     aladin_console.get_FoV( ra_n, pointing_n )

     # sending the string of probability percentage
     aladin_console.draw_string_float( ra_n, pointing_n, prob_fov )
     
     # VizieR Queries (astroquery.vizier)
     query_FOV.selected_catalog( ra_n, pointing_n, fov_width, fov_height, catalog )

     print_separation_line( )

     # re-initializing variables
     config_GWsky [ init_ra ], config_GWsky [ init_dec ] = ra_n, dec_n

     # calling last FOV "L"
     config_GWsky [ 'ra_last' ] = ra_n
     config_GWsky [ 'dec_last' ] =  dec_n

      
     with open( 'config_GWsky.ini', 'wb' ) as data:
                    pickle.dump( config_GWsky, data )
                    


def towards_south( input_skymap, fov_width, fov_height, catalog,
                  input_latitude, input_longitude, input_altitude, init_ra, init_dec ):

     import pickle
     import aladin_console
     import query_FOV
     from airmass import airmass_step
     from spherical_distance import shift_direction_nord_sud
     from progress_bars import print_separation_line, print_aladin_plane_x

     with open ( 'config_GWsky.ini', 'rb' ) as data:
          config_GWsky = pickle.load ( data )

     ra_s, dec_s = float( config_GWsky[ init_ra ] ), float( config_GWsky[ init_dec ] )

     # Undo "U"
     config_GWsky [ 'ra_last2' ] = ra_s
     config_GWsky [ 'dec_last2' ] =  dec_s
     
     shift_s = shift_direction_nord_sud( dec_s )
               
     pointing_s = ( dec_s - fov_height) - shift_s

     # variable cycle
     dec_s = pointing_s

     [ v1_ra, v1_dec, v2_ra, v2_dec, v3_ra, v3_dec,
      v4_ra, v4_dec ] = vertices_FOV( ra_s, pointing_s, fov_width, fov_height )
               
     print '---> The new FOV is centered at RA = ', str(('% .5f' % ra_s))+'°', 'and DEC =', str(('%.5f' % pointing_s))+'°. <---'
     print ''

     # vectors of ra and dec FOV corners
     ra_vertices, dec_vertices = [ v1_ra, v2_ra,  v4_ra, v3_ra ], [ v1_dec, v2_dec, v4_dec, v3_dec ]
               
     prob_fov = probability_inside_FOV( input_skymap, ra_vertices, dec_vertices )

     print '---> The integrated probability in the FOV is', str(('% .4e' % prob_fov ))+'%. <---'
     print ''

     # airmass in steps of one hour
     airmass_step( ra_s, pointing_s, input_latitude, input_longitude, input_altitude )
               
     print_aladin_plane_x( )
               
     # sending FoV
     aladin_console.get_FoV( ra_s, pointing_s )

     # sending the string of probability percentage
     aladin_console.draw_string_float( ra_s, pointing_s, prob_fov )
     
     # VizieR Queries (astroquery.vizier)
     query_FOV.selected_catalog( ra_s, pointing_s, fov_width, fov_height, catalog )

     print_separation_line( )

     # re-initializing variables
     config_GWsky [ init_ra ], config_GWsky [ init_dec ] = ra_s, dec_s

     # calling last FOV "L"
     config_GWsky [ 'ra_last' ] = ra_s
     config_GWsky [ 'dec_last' ] =  dec_s
     

     with open( 'config_GWsky.ini', 'wb' ) as data:
                    pickle.dump( config_GWsky, data )

   
                    
def towards_east( input_skymap, fov_width, fov_height, catalog,
                  input_latitude, input_longitude, input_altitude, init_ra, init_dec ):

     import pickle
     import aladin_console
     import query_FOV
     from airmass import airmass_step
     from spherical_distance import shift_direction_east_west
     from spherical_distance import ra0ra1_distance
     from progress_bars import print_separation_line, print_aladin_plane_x

     with open ( 'config_GWsky.ini', 'rb' ) as data:
          config_GWsky = pickle.load ( data )

     input_dec = float( config_GWsky[ 'input_dec' ] )

     ra_e, dec_e = float( config_GWsky[ init_ra ] ), float( config_GWsky[ init_dec ] )


     # Undo "U"
     config_GWsky [ 'ra_last2' ] = ra_e
     config_GWsky [ 'dec_last2' ] = dec_e


     shift_east = shift_direction_east_west( dec_e )

     offset_east_west = ra0ra1_distance( fov_width , input_dec, input_dec )
     
     pointing_e = ( ra_e + offset_east_west ) + shift_east

     # variable cycle
     ra_e = pointing_e

     [ v1_ra, v1_dec, v2_ra, v2_dec, v3_ra, v3_dec,
      v4_ra, v4_dec ] = vertices_FOV( pointing_e, dec_e, fov_width, fov_height )

     print '---> The new FOV is centered at RA = ', str(('% .5f' % pointing_e))+'°', 'and DEC =', str(('%.5f' % dec_e))+'°. <---'
     print ''

     # vectors of ra and dec FOV corners
     ra_vertices, dec_vertices = [ v1_ra, v2_ra,  v4_ra, v3_ra ], [ v1_dec, v2_dec, v4_dec, v3_dec ]
               
     prob_fov = probability_inside_FOV( input_skymap, ra_vertices, dec_vertices )

     print '---> The integrated probability in the FOV is', str(('% .4e' % prob_fov ))+'%. <---'
     print ''

     # airmass in steps of one hour
     airmass_step( pointing_e, dec_e, input_latitude, input_longitude, input_altitude )

     print_aladin_plane_x( )
               
     # sending FoV
     aladin_console.get_FoV( pointing_e, dec_e )

     # sending the string of probability percentage
     aladin_console.draw_string_float( pointing_e, dec_e, prob_fov )

     # VizieR Queries (astroquery.vizier)
     query_FOV.selected_catalog( pointing_e, dec_e, fov_width, fov_height, catalog )

     print_separation_line( )

     # re-initializing variables
     config_GWsky [ init_ra ], config_GWsky [ init_dec ] = ra_e, dec_e

     # calling last FOV "L"
     config_GWsky [ 'ra_last' ] = ra_e
     config_GWsky [ 'dec_last' ] = dec_e
     

     with open( 'config_GWsky.ini', 'wb' ) as data:
                    pickle.dump( config_GWsky, data )


                    
def towards_west( input_skymap, fov_width, fov_height, catalog,
                  input_latitude, input_longitude, input_altitude, init_ra, init_dec ):

     import pickle
     import aladin_console
     import query_FOV
     from airmass import airmass_step
     from spherical_distance import shift_direction_east_west
     from spherical_distance import ra0ra1_distance
     from progress_bars import print_separation_line, print_aladin_plane_x

     with open ( 'config_GWsky.ini', 'rb' ) as data:
          config_GWsky = pickle.load ( data )

     input_dec = float( config_GWsky[ 'input_dec' ] )

     ra_w, dec_w = float( config_GWsky[ init_ra ] ), float( config_GWsky[ init_dec ] )
     

     # Undo "U"
     config_GWsky [ 'ra_last2' ] = ra_w
     config_GWsky [ 'dec_last2' ] = dec_w


     shift_west = shift_direction_east_west( dec_w )
     
     offset_east_west = ra0ra1_distance( fov_width , input_dec, input_dec )
               
     pointing_w = ( ra_w - offset_east_west ) - shift_west

     # variable cycle
     ra_w = pointing_w

     [ v1_ra, v1_dec, v2_ra, v2_dec, v3_ra, v3_dec,
      v4_ra, v4_dec ] = vertices_FOV( pointing_w, dec_w, fov_width, fov_height )

     print '---> The new FOV is centered at RA = ', str(('% .5f' % pointing_w))+'°', 'and DEC =', str(('%.5f' % dec_w))+'°. <---'
     print ''

     # vectors of ra and dec FOV corners
     ra_vertices, dec_vertices = [ v1_ra, v2_ra,  v4_ra, v3_ra ], [ v1_dec, v2_dec, v4_dec, v3_dec ]
               
     prob_fov = probability_inside_FOV( input_skymap, ra_vertices, dec_vertices )

     print '---> The integrated probability in the FOV is', str(('% .4e' % prob_fov ))+'%. <---'
     print ''

     # airmass in steps of one hour
     airmass_step( pointing_w, dec_w, input_latitude, input_longitude, input_altitude )

     print_aladin_plane_x( )
               
     # sending FoV
     aladin_console.get_FoV( pointing_w, dec_w )

     # sending the string of probability percentage
     aladin_console.draw_string_float( pointing_w, dec_w, prob_fov )

     # VizieR Queries (astroquery.vizier)
     query_FOV.selected_catalog( pointing_w, dec_w, fov_width, fov_height, catalog )

     print_separation_line( )
     
     # re-initializing variables
     config_GWsky [ init_ra ], config_GWsky [ init_dec ] = ra_w, dec_w

     # calling last FOV "L"
     config_GWsky [ 'ra_last' ] = ra_w
     config_GWsky [ 'dec_last' ] = dec_w
    

     with open( 'config_GWsky.ini', 'wb' ) as data:
                    pickle.dump( config_GWsky, data )



def towards_north_east( input_skymap, fov_width, fov_height, catalog,
                  input_latitude, input_longitude, input_altitude, init_ra, init_dec ):

     import pickle
     import aladin_console
     import query_FOV
     from airmass import airmass_step
     from spherical_distance import ra0ra1_distance
     from spherical_distance import shift_direction_east_west
     from spherical_distance import shift_direction_nord_sud
     from spherical_distance import shape_FOV
     from progress_bars import print_separation_line, print_aladin_plane_x

     with open ( 'config_GWsky.ini', 'rb' ) as data:
          config_GWsky = pickle.load ( data )

     input_dec = float( config_GWsky[ 'input_dec' ] )
     
     ra_ne, dec_ne = float( config_GWsky [ 'ra_ne' ] ), float( config_GWsky [ 'dec_ne' ] )   # reading variables            

     # Undo "U"
     config_GWsky [ 'ra_last2' ] = ra_ne
     config_GWsky [ 'dec_last2' ] = dec_ne


     dec1 =  dec_ne
     dec2 = dec1 + fov_height

     shift_x_ne, shift_y_ne = shift_direction_east_west( input_dec ), shift_direction_nord_sud( input_dec )
               
               
     # square or rectangular FoV
     A = shape_FOV( fov_width, fov_height )

     offset_nord_est = ra0ra1_distance( A, dec1, dec2 )

     pointing_ne_ra, pointing_ne_dec = ( ra_ne + offset_nord_est ) + shift_x_ne , ( dec_ne + fov_height ) + shift_y_ne
               
                           
     # variable cycle; not saved
     ra_ne, dec_ne = pointing_ne_ra, pointing_ne_dec
               
     [ v1_ra, v1_dec, v2_ra, v2_dec, v3_ra, v3_dec,
      v4_ra, v4_dec ] = vertices_FOV ( pointing_ne_ra, pointing_ne_dec, fov_width, fov_height )

     print '---> The new FOV is centered at RA = ', str(('% .5f' % pointing_ne_ra))+'°', 'and DEC =', str(('%.5f' % pointing_ne_dec))+'°. <---'
     print ''

     # vectors of ra and dec FoV corners
     ra_vertices, dec_vertices = [ v1_ra, v2_ra,  v4_ra, v3_ra ], [ v1_dec, v2_dec, v4_dec, v3_dec ]
               
     prob_fov = probability_inside_FOV( input_skymap, ra_vertices, dec_vertices )

     print '---> The integrated probability in the FOV is', str(('% .4e' % prob_fov))+'%. <---'
     print ''

     # airmass in steps of one hour
     airmass_step( pointing_ne_ra, pointing_ne_dec, input_latitude, input_longitude, input_altitude )

     print_aladin_plane_x( )
          
     # sending FoV
     aladin_console.get_FoV( pointing_ne_ra, pointing_ne_dec )

     # sending the string of probability percentage
     aladin_console.draw_string_float( pointing_ne_ra, pointing_ne_dec, prob_fov )

     # VizieR Queries (astroquery.vizier)
     query_FOV.selected_catalog( pointing_ne_ra, pointing_ne_dec, fov_width, fov_height, catalog )

     print_separation_line( )

     # re-initializing variables
     config_GWsky [ init_ra ], config_GWsky [ init_dec ] = ra_ne, dec_ne    

     # calling last FOV "L"
     config_GWsky [ 'ra_last' ] = ra_ne
     config_GWsky [ 'dec_last' ] = dec_ne
     

     with open( 'config_GWsky.ini', 'wb' ) as data:
                    pickle.dump( config_GWsky, data )


                    
def towards_north_west( input_skymap, fov_width, fov_height, catalog,
                  input_latitude, input_longitude, input_altitude, init_ra, init_dec ):

     
     import pickle
     import aladin_console
     import query_FOV
     from airmass import airmass_step
     from spherical_distance import ra0ra1_distance
     from spherical_distance import shift_direction_east_west
     from spherical_distance import shift_direction_nord_sud
     from spherical_distance import shape_FOV
     from progress_bars import print_separation_line, print_aladin_plane_x

     with open ( 'config_GWsky.ini', 'rb' ) as data:
          config_GWsky = pickle.load ( data )

     input_dec = float( config_GWsky[ 'input_dec' ] )

     ra_nw, dec_nw = float( config_GWsky [ 'ra_nw' ] ), float( config_GWsky [ 'dec_nw' ] )   # reading variables           


     # Undo "U"
     config_GWsky [ 'ra_last2' ] = ra_nw
     config_GWsky [ 'dec_last2' ] = dec_nw

     
     # square or rectangular FoV
     A = shape_FOV( fov_width, fov_height )
                    
     dec1 = dec_nw
     dec2 = dec1 + fov_height

     shift_x_nw,  shift_y_nw = shift_direction_east_west( input_dec ), shift_direction_nord_sud( input_dec )
                              
     offset_nord_west = ra0ra1_distance( A, dec1, dec2 )

     pointing_nw_ra, pointing_nw_dec = ( ra_nw - offset_nord_west ) - shift_x_nw, ( dec_nw +  fov_height ) + shift_y_nw
                           
     # variable cycle
     ra_nw, dec_nw = pointing_nw_ra, pointing_nw_dec
               
               
     [ v1_ra, v1_dec, v2_ra, v2_dec, v3_ra, v3_dec,
      v4_ra, v4_dec ] = vertices_FOV( pointing_nw_ra, pointing_nw_dec, fov_width, fov_height )


     print '---> The new FOV is centered at RA = ', str(('% .5f' % pointing_nw_ra))+'°', 'and DEC =', str(('%.5f' % pointing_nw_dec))+'°. <---'
     print ''

     # vectors of ra and dec FoV corners
     ra_vertices, dec_vertices = [ v1_ra, v2_ra,  v4_ra, v3_ra ], [ v1_dec, v2_dec, v4_dec, v3_dec ]
               
     prob_fov = probability_inside_FOV( input_skymap, ra_vertices, dec_vertices )

     print '---> The integrated probability in the FOV is', str(('% .4e' % prob_fov ))+'%. <---'
     print ''

     # airmass in steps of one hour
     airmass_step( pointing_nw_ra, pointing_nw_dec, input_latitude, input_longitude, input_altitude )

     print_aladin_plane_x( )
          
     # sending FoV
     aladin_console.get_FoV( pointing_nw_ra, pointing_nw_dec )

     # sending the string of probability percentage
     aladin_console.draw_string_float( pointing_nw_ra, pointing_nw_dec, prob_fov )

     # VizieR Queries (astroquery.vizier)
     query_FOV.selected_catalog( pointing_nw_ra, pointing_nw_dec, fov_width, fov_height, catalog )

     print_separation_line( )

     # re-initializing variables
     config_GWsky [ init_ra ], config_GWsky [ init_dec ] = ra_nw, dec_nw

     # calling last FOV "L"
     config_GWsky [ 'ra_last' ] = ra_nw
     config_GWsky [ 'dec_last' ] = dec_nw

     
     with open( 'config_GWsky.ini', 'wb' ) as data:
                    pickle.dump( config_GWsky, data )

                        
                    
def towards_south_east( input_skymap, fov_width, fov_height, catalog,
             input_latitude, input_longitude, input_altitude, init_ra, init_dec ):

     
     import pickle
     import aladin_console
     import query_FOV
     from airmass import airmass_step
     from spherical_distance import ra0ra1_distance
     from spherical_distance import shift_direction_east_west
     from spherical_distance import shift_direction_nord_sud
     from spherical_distance import shape_FOV
     from progress_bars import print_separation_line, print_aladin_plane_x

     with open ( 'config_GWsky.ini', 'rb' ) as data:
          config_GWsky = pickle.load ( data )

     input_dec = float( config_GWsky[ 'input_dec' ] )

     ra_se, dec_se = float( config_GWsky [ 'ra_se' ] ), float( config_GWsky [ 'dec_se' ] )   # reading variables           


     # Undo "U"
     config_GWsky [ 'ra_last2' ] = ra_se
     config_GWsky [ 'dec_last2' ] = dec_se

     
     dec1 = dec_se
     dec2 = dec1 - fov_height

     shift_x_se, shift_y_se = shift_direction_east_west( input_dec ), shift_direction_nord_sud( input_dec )
               
               
     # square or rectangular FoV
     A = shape_FOV( fov_width, fov_height )

     offset_sud_est = ra0ra1_distance( A, dec1, dec2 )

     pointing_se_ra, pointing_se_dec = ( ra_se + offset_sud_est ) + shift_x_se , ( dec_se - fov_height ) - shift_y_se  

                          
     # variable cycle; not saved
     ra_se, dec_se = pointing_se_ra, pointing_se_dec

     [ v1_ra, v1_dec, v2_ra, v2_dec, v3_ra, v3_dec,
      v4_ra, v4_dec ] = vertices_FOV ( pointing_se_ra, pointing_se_dec, fov_width, fov_height )

     print '---> The new FOV is centered at RA = ', str(('% .5f' % pointing_se_ra))+'°', 'and DEC =', str(('%.5f' % pointing_se_dec))+'°. <---'
     print ''

     # vectors of ra and dec FoV corners
     ra_vertices, dec_vertices = [ v1_ra, v2_ra,  v4_ra, v3_ra ], [ v1_dec, v2_dec, v4_dec, v3_dec ]
               
     prob_fov = probability_inside_FOV( input_skymap, ra_vertices, dec_vertices )

     print '---> The integrated probability in the FOV is', str(('% .4e' % prob_fov))+'%. <---'
     print ''

     # airmass in steps of one hour
     airmass_step( pointing_se_ra, pointing_se_dec, input_latitude, input_longitude, input_altitude )

     print_aladin_plane_x( )
          
     # sending FoV
     aladin_console.get_FoV( pointing_se_ra, pointing_se_dec )

     # sending the string of probability percentage
     aladin_console.draw_string_float( pointing_se_ra, pointing_se_dec, prob_fov )

     # VizieR Queries (astroquery.vizier)
     query_FOV.selected_catalog( pointing_se_ra, pointing_se_dec, fov_width, fov_height, catalog )

     print_separation_line( )

     # re-initializing variables
     config_GWsky [ init_ra ], config_GWsky [ init_dec ] = ra_se, dec_se   

     # calling last FOV "L"
     config_GWsky [ 'ra_last' ] = ra_se
     config_GWsky [ 'dec_last' ] = dec_se


     with open( 'config_GWsky.ini', 'wb' ) as data:
                    pickle.dump( config_GWsky, data )


               
def towards_south_west( input_skymap, fov_width, fov_height, catalog,
             input_latitude, input_longitude, input_altitude, init_ra, init_dec ):


     import pickle
     import aladin_console
     import query_FOV
     from airmass import airmass_step
     from spherical_distance import ra0ra1_distance
     from spherical_distance import shift_direction_east_west
     from spherical_distance import shift_direction_nord_sud
     from spherical_distance import shape_FOV
     from progress_bars import print_separation_line, print_aladin_plane_x

     with open ( 'config_GWsky.ini', 'rb' ) as data:
          config_GWsky = pickle.load ( data )

     input_dec = float( config_GWsky[ 'input_dec' ] )

     ra_sw, dec_sw = float( config_GWsky [ 'ra_sw' ] ), float( config_GWsky [ 'dec_sw' ] )   # reading variables
                       

     # Undo "U"
     config_GWsky [ 'ra_last2' ] = ra_sw
     config_GWsky [ 'dec_last2' ] = dec_sw

     dec1 =  dec_sw
     dec2 = dec1 - fov_height

     shift_x_sw, shift_y_sw = shift_direction_east_west( input_dec ), shift_direction_nord_sud( input_dec )
               
               
     # square or rectangular FoV
     A = shape_FOV( fov_width, fov_height )

     offset_sud_west = ra0ra1_distance( A, dec1, dec2 )

     pointing_sw_ra, pointing_sw_dec = ( ra_sw - offset_sud_west ) - shift_x_sw , ( dec_sw - fov_height ) - shift_y_sw   
               
                          
     # variable cycle; not saved
     ra_sw, dec_sw = pointing_sw_ra, pointing_sw_dec

     [ v1_ra, v1_dec, v2_ra, v2_dec, v3_ra, v3_dec,
      v4_ra, v4_dec ] = vertices_FOV( pointing_sw_ra, pointing_sw_dec, fov_width, fov_height )

     print '---> The new FOV is centered at RA = ', str(('% .5f' % pointing_sw_ra))+'°', 'and DEC =', str(('%.5f' % pointing_sw_dec))+'°. <---'
     print ''

     # vectors of ra and dec FoV corners
     ra_vertices, dec_vertices = [ v1_ra, v2_ra,  v4_ra, v3_ra ], [ v1_dec, v2_dec, v4_dec, v3_dec ]
               
     prob_fov = probability_inside_FOV( input_skymap, ra_vertices, dec_vertices )

     print '---> The integrated probability in the FOV is', str(('% .4e' % prob_fov ))+'%. <---'
     print ''

     # airmass in steps of one hour
     airmass_step( pointing_sw_ra, pointing_sw_dec, input_latitude, input_longitude, input_altitude )

     print_aladin_plane_x( )
          
     # sending FoV
     aladin_console.get_FoV( pointing_sw_ra, pointing_sw_dec )

     # sending the string of probability percentage
     aladin_console.draw_string_float( pointing_sw_ra, pointing_sw_dec, prob_fov )

     # VizieR Queries (astroquery.vizier)
     query_FOV.selected_catalog( pointing_sw_ra, pointing_sw_dec, fov_width, fov_height, catalog )

     print_separation_line( )

     # re-initializing variables
     config_GWsky [ init_ra ], config_GWsky [ init_dec ] = ra_sw, dec_sw  


     # calling last FOV "L"
     config_GWsky [ 'ra_last' ] = ra_sw
     config_GWsky [ 'dec_last' ] = dec_sw
     

     with open( 'config_GWsky.ini', 'wb' ) as data:
                    pickle.dump( config_GWsky, data )



def user_FOV( input_skymap, fov_width, fov_height, catalog,
             input_latitude, input_longitude, input_altitude, init_ra, init_dec ):

     import numpy as np
     import pickle
     import aladinSAMP
     import aladin_console
     import query_FOV
     from airmass import airmass_step
     from progress_bars import print_separation_line, print_aladin_plane


     with open ( 'config_GWsky.ini', 'rb' ) as data:
          config_GWsky = pickle.load ( data )

     input_ra, input_dec = float( config_GWsky[ init_ra ] ), float( config_GWsky[ init_dec ] )

     # Undo "U"
     config_GWsky [ 'ra_last2' ] = input_ra
     config_GWsky [ 'dec_last2' ] = input_dec

     [ v1_ra, v1_dec, v2_ra, v2_dec, v3_ra, v3_dec,
      v4_ra, v4_dec ] = vertices_FOV( input_ra, input_dec, fov_width, fov_height )

     # vectors of ra and dec FoV corners
     ra_vertices, dec_vertices = [ v1_ra, v2_ra,  v4_ra, v3_ra ], [ v1_dec, v2_dec, v4_dec, v3_dec ]

     prob_fov = probability_inside_FOV( input_skymap, ra_vertices, dec_vertices )

     
     print ' The integrated probability in the selected FOV \n at RA = ', str(input_ra)+'°', 'and DEC =',  str(input_dec)+'°', 'is', str(('% .4e' % prob_fov)) +'%.'
     print ''

     # airmass in steps of one hour
     airmass_step( input_ra, input_dec, input_latitude, input_longitude, input_altitude )
          
     print_aladin_plane( )
     
     # sending FoV
     aladin_console.get_FoV( input_ra, input_dec )

     # sending the string of probability percentage
     aladin_console.draw_string_float( input_ra, input_dec, prob_fov )
     
     # VizieR Queries (astroquery.vizier)
     query_FOV.selected_catalog( input_ra, input_dec, fov_width, fov_height, catalog )

     print_separation_line( )

     # calling last FOV "L"
     config_GWsky [ 'ra_last' ] = input_ra
     config_GWsky [ 'dec_last' ] = input_dec


     with open( 'config_GWsky.ini', 'wb' ) as data:
                    pickle.dump( config_GWsky, data )



def probability_inside_FOV( input_skymap, ra_vertices, dec_vertices ):

     """

        Give the probability contained in a specific FOV
     
     """

     
     import healpy as hp
     import numpy as np

     hpx = hp.read_map( input_skymap, verbose = False )

     npix = len( hpx )

     nside = hp.npix2nside( npix )

     theta = 0.5 * np.pi - np.deg2rad( dec_vertices )

     phi = np.deg2rad( ra_vertices )
     
     xyz = hp.ang2vec( theta, phi )

     ipix_poly = hp.query_polygon( nside, xyz )

     # probability contains in a specific FOV
     probability_inside_polygon = hpx[ipix_poly].sum()

     return probability_inside_polygon
     

def vertices_FOV( ra_center, dec_center, FOV_base, FOV_height ):

     """

     Find the vertices of a FOV given a center position (RA[deg], DEC[deg])
     and the FOV size (FOV_base[deg], FOV_height[deg])

                ***it will be replaced***
                
     """

     from math import sin, cos, acos, degrees, radians
     from spherical_distance import ra0ra1_distance


     # upper vertices
     A =  FOV_base/2.0
     
     dec0_N = ( dec_center ) + ( FOV_height / 2.0 )
     dec1_N = dec0_N
         
     offset_N = ra0ra1_distance( A, dec0_N, dec1_N )

     # lower vertices
     dec0_S = ( dec_center ) - ( FOV_height / 2.0  ) 
     dec1_S = dec0_S  

     offset_S = ra0ra1_distance( A, dec0_S, dec1_S )
     
     
     # FOV vertices 
     vertex_1_ra, vertex_1_dec = ( ra_center + offset_N ), ( dec_center + FOV_height / 2.0 )   
     vertex_2_ra, vertex_2_dec = ( ra_center - offset_N ), ( dec_center + FOV_height / 2.0 )
     vertex_3_ra, vertex_3_dec = ( ra_center + offset_S ), ( dec_center - FOV_height / 2.0 )  
     vertex_4_ra, vertex_4_dec = ( ra_center - offset_S ), ( dec_center - FOV_height / 2.0 )


     return vertex_1_ra, vertex_1_dec, vertex_2_ra, vertex_2_dec,  vertex_3_ra, vertex_3_dec, vertex_4_ra, vertex_4_dec    

def replaced( ):

     '''

     the new right ascension and declination values are iteratively replaced for a new cycle
     
     '''
     
     import pickle

     while True:
          try:
               ra_new, dec_new = input( ' Insert a new FOV center RA[deg], DEC[deg]: ' )
          except SyntaxError:
               print ''
               print '    --------------------------------------------------------'
               print '        ***Please specify a new sky position*** '
               print '                  ex.: -24.63, -70.40   '
               print '    --------------------------------------------------------'
               print ''
          except TypeError:
               print ''
               print '    --------------------------------------------------------'
               print '        ***Please specify a new sky position*** '
               print '                  ex.: -24.63, -70.40   '
               print '    --------------------------------------------------------'
               print ''
          except NameError:
               print ''
               print '    --------------------------------------------------------'
               print '        ***Please specify a new sky position *** '
               print '                  ex.: -24.63, -70.40   '
               print '    --------------------------------------------------------'
               print ''
          except ValueError:
               print ''
               print '    --------------------------------------------------------'
               print '        ***Please specify a new sky position *** '
               print '                  ex.: -24.63, -70.40   '
               print '    --------------------------------------------------------'
               print ''
          else:
               break
     
     with open ( 'config_GWsky.ini', 'rb' ) as data:
          config_GWsky = pickle.load( data )
 
     config_GWsky[ 'input_ra' ] = ra_new
     config_GWsky[ 'input_dec' ] = dec_new

     config_GWsky[ 'ra_n' ] = ra_new
     config_GWsky[ 'dec_n' ] = dec_new

     config_GWsky[ 'ra_s' ] = ra_new
     config_GWsky[ 'dec_s' ] = dec_new

     config_GWsky[ 'ra_e' ] = ra_new
     config_GWsky[ 'dec_e' ] = dec_new

     config_GWsky[ 'ra_w' ] = ra_new
     config_GWsky[ 'dec_w' ] = dec_new

     config_GWsky[ 'ra_nw' ] = ra_new
     config_GWsky[ 'dec_nw' ] = dec_new

     config_GWsky[ 'ra_sw' ] = ra_new
     config_GWsky[ 'dec_sw' ] = dec_new

     config_GWsky[ 'ra_ne' ] = ra_new
     config_GWsky[ 'dec_ne' ] = dec_new

     config_GWsky[ 'ra_se' ] = ra_new
     config_GWsky[ 'dec_se' ] = dec_new

     with open( 'config_GWsky.ini', 'wb' ) as data:
                    pickle.dump( config_GWsky, data )


def add_FOV( skymap, FOV_base, FOV_height, ra, dec ):

     """
          Define a sequence of instrument FOVs from a predetermined sky position:
          North/NorthWest/West/SouthWest/South/SouthEast/East/NorthEast directions.
     
     """

     import pickle
     
     # reading parameters from config_GWsky.ini
     
     with open ( 'config_GWsky.ini', 'rb' ) as data:
          config_GWsky = pickle.load ( data )
          

     input_latitude = float ( config_GWsky [ 'input_latitude' ] )
     input_longitude = float ( config_GWsky [ 'input_longitude' ] )
     input_altitude = float ( config_GWsky [ 'input_altitude' ] )
     fov_width = float ( config_GWsky [ 'fov_width' ] )
     fov_height = float ( config_GWsky [ 'fov_height' ] )
     catalog = config_GWsky[ 'catalog' ]
     input_skymap = config_GWsky[ 'input_skymap' ]


     #---------------- FoV at the sky coordinates insert from users; (ra dec)---------------------#
     user_FOV( input_skymap, fov_width, fov_height, catalog,
               input_latitude, input_longitude, input_altitude, init_ra = 'input_ra', init_dec = 'input_dec' )
    
     
     # -----------sequence of instrument FOVs from a predetermined sky position---------------------#
      #----North/NorthWest/West/SouthWest/South/SouthEast/East/NorthEast directions are permitted----#

     while True:

          print ' '
          print ' '
          print ' Press (N/NW/W/SW/S/SE/E/NE) to add FOVs in '
          print ' North/NorthWest/West/SouthWest/South/SouthEst/Est/NorthEst; '
          print ' Press R to add a new FOV center; '
          print ' Press Q to exit. '      
          direction = raw_input( ' : ' )

          while True:
               if direction == '':
                    print ' '
                    print ' '
                    print ' Press (N/NW/W/SW/S/SE/E/NE) to add FOVs in '
                    print ' North/NorthWest/West/SouthWest/South/SouthEast/East/NorthEast; '
                    print ' Press R to add a new FOV center; '
                    print ' Press Q to exit. '      
                    direction = raw_input( ' : ' )
               else:
                    break
         
           
          if direction.strip().capitalize() == 'N':
               
               towards_north( input_skymap, fov_width, fov_height, catalog,
                             input_latitude, input_longitude, input_altitude, init_ra = 'ra_n', init_dec = 'dec_n' )


          elif direction.strip().capitalize() == 'S':

               towards_south( input_skymap, fov_width, fov_height, catalog,
                             input_latitude, input_longitude, input_altitude, init_ra = 'ra_s', init_dec = 'dec_s' )

   
          elif direction.strip().capitalize() == 'E':

               towards_east( input_skymap, fov_width, fov_height, catalog,
                             input_latitude, input_longitude, input_altitude, init_ra = 'ra_e', init_dec = 'dec_e' )              
               
                                                                      
          elif direction.strip().capitalize() == 'W':

               towards_west( input_skymap, fov_width, fov_height, catalog,
                             input_latitude, input_longitude, input_altitude, init_ra = 'ra_w', init_dec = 'dec_w' )


          elif direction[0].strip().capitalize() == 'N' and direction[1].strip().capitalize() == 'W':

               towards_north_west( input_skymap, fov_width, fov_height, catalog,
                             input_latitude, input_longitude, input_altitude, init_ra = 'ra_nw', init_dec = 'dec_nw' )

               
          elif direction[0].strip().capitalize() == 'S' and direction[1].strip().capitalize() == 'W':

               towards_south_west( input_skymap, fov_width, fov_height, catalog,
                             input_latitude, input_longitude, input_altitude, init_ra = 'ra_sw', init_dec = 'dec_sw' )

               
               
          elif direction[0].strip().capitalize() == 'N' and direction[1].strip().capitalize() == 'E':

               towards_north_east( input_skymap, fov_width, fov_height, catalog,
                             input_latitude, input_longitude, input_altitude, init_ra = 'ra_ne', init_dec = 'dec_ne' )

               
               
          elif direction[0].strip().capitalize() == 'S' and direction[1].strip().capitalize() == 'E':

               towards_south_east( input_skymap, fov_width, fov_height, catalog,
                             input_latitude, input_longitude, input_altitude, init_ra = 'ra_se', init_dec = 'dec_se' )

               
 
          # running a new cycle from a  new position
          elif direction.strip().capitalize() == 'R':
     
               replaced( ) # the ra and dec are iteratively replaced for a new cycle

               import pickle
     
               # reading parameters from config_GWsky.ini
     
               with open ( 'config_GWsky.ini', 'rb' ) as data:
                    config_GWsky = pickle.load ( data )
          

               input_latitude = float ( config_GWsky [ 'input_latitude' ] )
               input_longitude = float ( config_GWsky [ 'input_longitude' ] )
               input_altitude = float ( config_GWsky [ 'input_altitude' ] )
               fov_width = float ( config_GWsky [ 'fov_width' ] )
               fov_height = float ( config_GWsky [ 'fov_height' ] )
               catalog = config_GWsky[ 'catalog' ]
               input_skymap = config_GWsky[ 'input_skymap' ]
    

               user_FOV( input_skymap, fov_width, fov_height, catalog,
               input_latitude, input_longitude, input_altitude, init_ra = 'input_ra', init_dec = 'input_dec' )

               print ' '
               print ' '
               print ' Press (N/NW/W/SW/S/SE/E/NE) to add FOVs in '
               print ' North/NorthWest/West/SouthWest/South/SouthEast/East/NorthEast; '
               #print ' Press R to add a new FOV center; '
               print ' Press Q to exit. '      
               direction = raw_input( ' : ' )

               while True:
                    if direction == '':
                         print ' '
                         print ' '
                         print ' Press (N/NW/W/SW/S/SE/E/NE) to add FOVs in '
                         print ' North/NorthWest/West/SouthWest/South/SouthEast/East/NorthEast; '
                         print ' Press R to add a new FOV center; '
                         print ' Press Q to exit. '      
                         direction = raw_input( ' : ' )
                    else:
                         break
         
           
               if direction.strip().capitalize() == 'N':
                    towards_north( input_skymap, fov_width, fov_height, catalog,
                                  input_latitude, input_longitude, input_altitude, init_ra = 'ra_n', init_dec = 'dec_n' )


               elif direction.strip().capitalize() == 'S':
                    towards_south( input_skymap, fov_width, fov_height, catalog,
                                   input_latitude, input_longitude, input_altitude, init_ra = 'ra_s', init_dec = 'dec_s' )

   
               elif direction.strip().capitalize() == 'E':
                    towards_east( input_skymap, fov_width, fov_height, catalog,
                                  input_latitude, input_longitude, input_altitude, init_ra = 'ra_e', init_dec = 'dec_e' )              
               
                                                                      
               elif direction.strip().capitalize() == 'W':
                    towards_west( input_skymap, fov_width, fov_height, catalog,
                                  input_latitude, input_longitude, input_altitude, init_ra = 'ra_w', init_dec = 'dec_w' )


               elif direction[0].strip().capitalize() == 'N' and direction[1].strip().capitalize() == 'W':
                    towards_north_west( input_skymap, fov_width, fov_height, catalog,
                                       input_latitude, input_longitude, input_altitude, init_ra = 'ra_nw', init_dec = 'dec_nw' )


               elif direction[0].strip().capitalize() == 'S' and direction[1].strip().capitalize() == 'W':
                    towards_south_west( input_skymap, fov_width, fov_height, catalog,
                                       input_latitude, input_longitude, input_altitude, init_ra = 'ra_sw', init_dec = 'dec_sw' )

               
               elif direction[0].strip().capitalize() == 'N' and direction[1].strip().capitalize() == 'E':
                    towards_north_east( input_skymap, fov_width, fov_height, catalog,
                                       input_latitude, input_longitude, input_altitude, init_ra = 'ra_ne', init_dec = 'dec_ne' )

                
               elif direction[0].strip().capitalize() == 'S' and direction[1].strip().capitalize() == 'E':
                   towards_south_east( input_skymap, fov_width, fov_height, catalog,
                                      input_latitude, input_longitude, input_altitude, init_ra = 'ra_se', init_dec = 'dec_se' )


          # running a new cycle from an existing FOV
          elif direction.strip().capitalize() == 'F':
               
               replaced( ) # the ra and dec are iteratively replaced for a new cycle

               import pickle
     
               # reading parameters from config_GWsky.ini
     
               with open ( 'config_GWsky.ini', 'rb' ) as data:
                    config_GWsky = pickle.load ( data )
          

               input_latitude = float ( config_GWsky [ 'input_latitude' ] )
               input_longitude = float ( config_GWsky [ 'input_longitude' ] )
               input_altitude = float ( config_GWsky [ 'input_altitude' ] )
               fov_width = float ( config_GWsky [ 'fov_width' ] )
               fov_height = float ( config_GWsky [ 'fov_height' ] )
               catalog = config_GWsky[ 'catalog' ]
               input_skymap = config_GWsky[ 'input_skymap' ]
    

               print ' '
               print ' '
               print ' Press (N/NW/W/SW/S/SE/E/NE) to add FOVs in '
               print ' North/NorthWest/West/SouthWest/South/SouthEast/East/NorthEast; '
               print ' Press R to add a new FOV center; '
               print ' Press Q to exit. '      
               direction = raw_input( ' : ' )

               while True:
                    if direction == '':
                         print ' '
                         print ' '
                         print ' Press (N/NW/W/SW/S/SE/E/NE) to add FOVs in '
                         print ' North/NorthWest/West/SouthWest/South/SouthEast/East/NorthEast; '
                         print ' Press R to add a new FOV center; '
                         print ' Press Q to exit. '      
                         direction = raw_input( ' : ' )
                    else:
                         break
         
           
               if direction.strip().capitalize() == 'N':
                    towards_north( input_skymap, fov_width, fov_height, catalog,
                                  input_latitude, input_longitude, input_altitude, init_ra = 'ra_n', init_dec = 'dec_n' )


               elif direction.strip().capitalize() == 'S':
                    towards_south( input_skymap, fov_width, fov_height, catalog,
                                   input_latitude, input_longitude, input_altitude, init_ra = 'ra_s', init_dec = 'dec_s' )

   
               elif direction.strip().capitalize() == 'E':
                    towards_east( input_skymap, fov_width, fov_height, catalog,
                                  input_latitude, input_longitude, input_altitude, init_ra = 'ra_e', init_dec = 'dec_e' )              
               
                                                                      
               elif direction.strip().capitalize() == 'W':
                    towards_west( input_skymap, fov_width, fov_height, catalog,
                                  input_latitude, input_longitude, input_altitude, init_ra = 'ra_w', init_dec = 'dec_w' )


               elif direction[0].strip().capitalize() == 'N' and direction[1].strip().capitalize() == 'W':
                    towards_north_west( input_skymap, fov_width, fov_height, catalog,
                                       input_latitude, input_longitude, input_altitude, init_ra = 'ra_nw', init_dec = 'dec_nw' )


               elif direction[0].strip().capitalize() == 'S' and direction[1].strip().capitalize() == 'W':
                    towards_south_west( input_skymap, fov_width, fov_height, catalog,
                                       input_latitude, input_longitude, input_altitude, init_ra = 'ra_sw', init_dec = 'dec_sw' )

               
               elif direction[0].strip().capitalize() == 'N' and direction[1].strip().capitalize() == 'E':
                    towards_north_east( input_skymap, fov_width, fov_height, catalog,
                                       input_latitude, input_longitude, input_altitude, init_ra = 'ra_ne', init_dec = 'dec_ne' )

                
               elif direction[0].strip().capitalize() == 'S' and direction[1].strip().capitalize() == 'E':
                   towards_south_east( input_skymap, fov_width, fov_height, catalog,
                                      input_latitude, input_longitude, input_altitude, init_ra = 'ra_se', init_dec = 'dec_se' )




             # running a new cycle from the last FOV
          elif direction.strip().capitalize() == 'L':
               
               # ra and dec are iteratively replaced for a new cycle

               import pickle

               with open ( 'config_GWsky.ini', 'rb' ) as data:
                    config_GWsky = pickle.load( data )

               ra_last, dec_last = float ( config_GWsky [ 'ra_last' ] ), float ( config_GWsky [ 'dec_last' ] )
 
               config_GWsky[ 'input_ra' ] = ra_last
               config_GWsky[ 'input_dec' ] = dec_last

               config_GWsky[ 'ra_n' ] = ra_last
               config_GWsky[ 'dec_n' ] = dec_last

               config_GWsky[ 'ra_s' ] = ra_last
               config_GWsky[ 'dec_s' ] = dec_last

               config_GWsky[ 'ra_e' ] = ra_last
               config_GWsky[ 'dec_e' ] = dec_last

               config_GWsky[ 'ra_w' ] = ra_last
               config_GWsky[ 'dec_w' ] = dec_last

               config_GWsky[ 'ra_nw' ] = ra_last
               config_GWsky[ 'dec_nw' ] = dec_last

               config_GWsky[ 'ra_sw' ] = ra_last
               config_GWsky[ 'dec_sw' ] = dec_last

               config_GWsky[ 'ra_ne' ] = ra_last
               config_GWsky[ 'dec_ne' ] = dec_last

               config_GWsky[ 'ra_se' ] = ra_last
               config_GWsky[ 'dec_se' ] = dec_last

               with open( 'config_GWsky.ini', 'wb' ) as data:
                    pickle.dump( config_GWsky, data )

               
     
               # reading parameters from config_GWsky.ini
     
               with open ( 'config_GWsky.ini', 'rb' ) as data:
                    config_GWsky = pickle.load ( data )
          

               input_latitude = float ( config_GWsky [ 'input_latitude' ] )
               input_longitude = float ( config_GWsky [ 'input_longitude' ] )
               input_altitude = float ( config_GWsky [ 'input_altitude' ] )
               fov_width = float ( config_GWsky [ 'fov_width' ] )
               fov_height = float ( config_GWsky [ 'fov_height' ] )
               catalog = config_GWsky[ 'catalog' ]
               input_skymap = config_GWsky[ 'input_skymap' ]
    

               print ' '
               print ' '
               print ' Press (N/NW/W/SW/S/SE/E/NE) to add FOVs in '
               print ' North/NorthWest/West/SouthWest/South/SouthEast/East/NorthEast; '
               print ' Press R to add a new FOV center; '
               print ' Press Q to exit. '      
               direction = raw_input( ' : ' )

               while True:
                    if direction == '':
                         print ' '
                         print ' '
                         print ' Press (N/NW/W/SW/S/SE/E/NE) to add FOVs in '
                         print ' North/NorthWest/West/SouthWest/South/SouthEast/East/NorthEast; '
                         print ' Press R to add a new FOV center; '
                         print ' Press Q to exit. '      
                         direction = raw_input( ' : ' )
                    else:
                         break
         
           
               if direction.strip().capitalize() == 'N':
                    towards_north( input_skymap, fov_width, fov_height, catalog,
                                  input_latitude, input_longitude, input_altitude, init_ra = 'ra_n', init_dec = 'dec_n' )

                  

               elif direction.strip().capitalize() == 'S':
                    towards_south( input_skymap, fov_width, fov_height, catalog,
                                   input_latitude, input_longitude, input_altitude, init_ra = 'ra_s', init_dec = 'dec_s' )

                   
   
               elif direction.strip().capitalize() == 'E':
                    towards_east( input_skymap, fov_width, fov_height, catalog,
                                  input_latitude, input_longitude, input_altitude, init_ra = 'ra_e', init_dec = 'dec_e' )              


                   
               elif direction.strip().capitalize() == 'W':
                    towards_west( input_skymap, fov_width, fov_height, catalog,
                                  input_latitude, input_longitude, input_altitude, init_ra = 'ra_w', init_dec = 'dec_w' )

                         

               elif direction[0].strip().capitalize() == 'N' and direction[1].strip().capitalize() == 'W':
                    towards_north_west( input_skymap, fov_width, fov_height, catalog,
                                       input_latitude, input_longitude, input_altitude, init_ra = 'ra_nw', init_dec = 'dec_nw' )

                    
                         

               elif direction[0].strip().capitalize() == 'S' and direction[1].strip().capitalize() == 'W':
                    towards_south_west( input_skymap, fov_width, fov_height, catalog,
                                       input_latitude, input_longitude, input_altitude, init_ra = 'ra_sw', init_dec = 'dec_sw' )
                    
                     
               
               elif direction[0].strip().capitalize() == 'N' and direction[1].strip().capitalize() == 'E':
                    towards_north_east( input_skymap, fov_width, fov_height, catalog,
                                       input_latitude, input_longitude, input_altitude, init_ra = 'ra_ne', init_dec = 'dec_ne' )

                    
                
               elif direction[0].strip().capitalize() == 'S' and direction[1].strip().capitalize() == 'E':
                   towards_south_east( input_skymap, fov_width, fov_height, catalog,
                                      input_latitude, input_longitude, input_altitude, init_ra = 'ra_se', init_dec = 'dec_se' )



              # running a new cycle from the second-last FOV
          elif direction.strip().capitalize() == 'U':
               
               # ra and dec are iteratively replaced for a new cycle

               import pickle

               with open ( 'config_GWsky.ini', 'rb' ) as data:
                    config_GWsky = pickle.load( data )

               ra_last2, dec_last2 = float ( config_GWsky [ 'ra_last2' ] ), float ( config_GWsky [ 'dec_last2' ] )
 
               config_GWsky[ 'input_ra' ] = ra_last2
               config_GWsky[ 'input_dec' ] = dec_last2

               config_GWsky[ 'ra_n' ] = ra_last2
               config_GWsky[ 'dec_n' ] = dec_last2

               config_GWsky[ 'ra_s' ] = ra_last2
               config_GWsky[ 'dec_s' ] = dec_last2

               config_GWsky[ 'ra_e' ] = ra_last2
               config_GWsky[ 'dec_e' ] = dec_last2

               config_GWsky[ 'ra_w' ] = ra_last2
               config_GWsky[ 'dec_w' ] = dec_last2

               config_GWsky[ 'ra_nw' ] = ra_last2
               config_GWsky[ 'dec_nw' ] = dec_last2

               config_GWsky[ 'ra_sw' ] = ra_last2
               config_GWsky[ 'dec_sw' ] = dec_last2

               config_GWsky[ 'ra_ne' ] = ra_last2
               config_GWsky[ 'dec_ne' ] = dec_last2

               config_GWsky[ 'ra_se' ] = ra_last2
               config_GWsky[ 'dec_se' ] = dec_last2

               with open( 'config_GWsky.ini', 'wb' ) as data:
                    pickle.dump( config_GWsky, data )

               
     
               # reading parameters from config_GWsky.ini
     
               with open ( 'config_GWsky.ini', 'rb' ) as data:
                    config_GWsky = pickle.load ( data )
          

               input_latitude = float ( config_GWsky [ 'input_latitude' ] )
               input_longitude = float ( config_GWsky [ 'input_longitude' ] )
               input_altitude = float ( config_GWsky [ 'input_altitude' ] )
               fov_width = float ( config_GWsky [ 'fov_width' ] )
               fov_height = float ( config_GWsky [ 'fov_height' ] )
               catalog = config_GWsky[ 'catalog' ]
               input_skymap = config_GWsky[ 'input_skymap' ]
    

               print ' '
               print ' '
               print ' Press (N/NW/W/SW/S/SE/E/NE) to add FOVs in '
               print ' North/NorthWest/West/SouthWest/South/SouthEast/East/NorthEast; '
               print ' Press R to add a new FOV center; '
               print ' Press Q to exit. '      
               direction = raw_input( ' : ' )

               while True:
                    if direction == '':
                         print ' '
                         print ' '
                         print ' Press (N/NW/W/SW/S/SE/E/NE) to add FOVs in '
                         print ' North/NorthWest/West/SouthWest/South/SouthEast/East/NorthEast; '
                         print ' Press R to add a new FOV center; '
                         print ' Press Q to exit. '      
                         direction = raw_input( ' : ' )
                    else:
                         break
         
           
               if direction.strip().capitalize() == 'N':
                    towards_north( input_skymap, fov_width, fov_height, catalog,
                                  input_latitude, input_longitude, input_altitude, init_ra = 'ra_n', init_dec = 'dec_n' )

                  

               elif direction.strip().capitalize() == 'S':
                    towards_south( input_skymap, fov_width, fov_height, catalog,
                                   input_latitude, input_longitude, input_altitude, init_ra = 'ra_s', init_dec = 'dec_s' )

                   
   
               elif direction.strip().capitalize() == 'E':
                    towards_east( input_skymap, fov_width, fov_height, catalog,
                                  input_latitude, input_longitude, input_altitude, init_ra = 'ra_e', init_dec = 'dec_e' )              


                   
               elif direction.strip().capitalize() == 'W':
                    towards_west( input_skymap, fov_width, fov_height, catalog,
                                  input_latitude, input_longitude, input_altitude, init_ra = 'ra_w', init_dec = 'dec_w' )

                         

               elif direction[0].strip().capitalize() == 'N' and direction[1].strip().capitalize() == 'W':
                    towards_north_west( input_skymap, fov_width, fov_height, catalog,
                                       input_latitude, input_longitude, input_altitude, init_ra = 'ra_nw', init_dec = 'dec_nw' )

                    
                         

               elif direction[0].strip().capitalize() == 'S' and direction[1].strip().capitalize() == 'W':
                    towards_south_west( input_skymap, fov_width, fov_height, catalog,
                                       input_latitude, input_longitude, input_altitude, init_ra = 'ra_sw', init_dec = 'dec_sw' )
                    
                     
               
               elif direction[0].strip().capitalize() == 'N' and direction[1].strip().capitalize() == 'E':
                    towards_north_east( input_skymap, fov_width, fov_height, catalog,
                                       input_latitude, input_longitude, input_altitude, init_ra = 'ra_ne', init_dec = 'dec_ne' )

                    
                
               elif direction[0].strip().capitalize() == 'S' and direction[1].strip().capitalize() == 'E':
                   towards_south_east( input_skymap, fov_width, fov_height, catalog,
                                      input_latitude, input_longitude, input_altitude, init_ra = 'ra_se', init_dec = 'dec_se' )



                         
          # closing cycle          
          elif direction.strip().capitalize() == 'Q':
               break
