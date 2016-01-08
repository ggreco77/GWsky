# -*- coding: utf-8 -*-

def percentage_FOV( sky_map,ra_vertices,dec_vertices ):

     """

     Give the probability contains in a specific FOV
     
     """

     global  percentage_poly
     
     import healpy as hp
     import numpy as np

     hpx = hp.read_map( sky_map,verbose = False )

     npix = len( hpx )

     nside = hp.npix2nside( npix )

     theta = 0.5 * np.pi - np.deg2rad( dec_vertices )

     phi = np.deg2rad( ra_vertices )
     
     xyz = hp.ang2vec( theta, phi )

     ipix_poly = hp.query_polygon( nside, xyz )

     # probability contains in a specific FOV
     percentage_poly = hpx[ipix_poly].sum()


def vertices_FOV( ra, dec, FOV_base, FOV_height ):

     """

        Find the vertices of a FOV given a sky center position (RA[deg], DEC[deg])
        and the FOV size (deg).
        
     """

     from math import sin, cos, acos, degrees, radians

     # FOV vertices global variables
     global vertex_1_ra
     global vertex_1_dec
     global vertex_2_ra
     global vertex_2_dec
     global vertex_3_ra
     global vertex_3_dec
     global vertex_4_ra
     global vertex_4_dec
     
     delta0_rads = radians( dec ) + radians( FOV_height / 2.0 ) 
     delta1_rads = delta0_rads 

     x_vertices_FOV = ( cos(radians(FOV_base/2.0))-sin(delta0_rads)*sin(delta1_rads))/(cos(delta0_rads)*cos(delta1_rads ))

     offset_N = degrees( acos(x_vertices_FOV ))

     delta0_rads = radians( dec ) - radians( FOV_height / 2.0  ) 
     delta1_rads = delta0_rads 

     x_vertices_FOV = ( cos(radians(FOV_base/2.0))-sin(delta0_rads)*sin(delta1_rads))/(cos(delta0_rads)*cos(delta1_rads ))

     offset_S = degrees( acos(x_vertices_FOV ))

     # FOV vertices 
     vertex_1_ra, vertex_1_dec = ra + offset_N, dec + FOV_height / 2.0   

     vertex_2_ra, vertex_2_dec = ra - offset_N, dec + FOV_height / 2.0

     vertex_3_ra, vertex_3_dec = ra + offset_S, dec - FOV_height / 2.0   

     vertex_4_ra, vertex_4_dec = ra - offset_S, dec - FOV_height / 2.0



def add_FOV( sky_map, FOV_base, FOV_height, ra, dec ):

     """
          Define a sequence of instrument FOVs from a predetermined sky position.
          North/South/East/West directions.
     
     """
     
     import aladinSAMP
     
     from airmass import airmass
     
     import healpy as hp
     
     import numpy as np
     
     from math import sin, cos, acos, degrees, radians
     
     import pickle

     # read parameters from config_GWsky.ini
     with open ( 'config_GWsky.ini', 'rb' ) as data:
          config_GWsky = pickle.load ( data )

     ra = float(config_GWsky[ 'ra' ])

     dec = float(config_GWsky[ 'dec' ])

     lat_input = float ( config_GWsky [ 'lat_input' ] )

     lon_input = float ( config_GWsky [ 'lon_input' ] )

     height_input = float ( config_GWsky [ 'height_input' ] )

     time_input = config_GWsky[ 'time_input' ]

     catalog = config_GWsky[ 'catalog' ]

     
     delta0_rads = radians(dec)

     delta1_rads = delta0_rads 

     x_est_west = (cos(radians(FOV_base))-sin(delta0_rads)*sin(delta1_rads))/(cos(delta0_rads)*cos(delta1_rads))

     offset_est_west = degrees(acos(x_est_west))

     # initializing variables
     p0_N = [ ra, dec ]
     p0_S = [ ra, dec ]
     p0_O = [ ra, dec ]
     p0_E = [ ra, dec ]


     #---------------- FOV at the sky coordinates insert from users; (ra dec)---------------------#


     # Find the vertices of FOV 
     vertices_FOV ( ra, dec, FOV_base, FOV_height )

     # sky coord FOV center
     coord = [ ra, dec ]

     # FOV corners
     vtx_1 =  [ vertex_1_ra, vertex_1_dec ]
     vtx_2 =  [ vertex_2_ra, vertex_2_dec ]
     vtx_3 =  [ vertex_3_ra, vertex_3_dec ]
     vtx_4 =  [ vertex_4_ra, vertex_4_dec ]

     # vectors of ra and dec FOV corners
     ra_vertices = [ vertex_1_ra, vertex_2_ra,  vertex_4_ra, vertex_3_ra ]
     dec_vertices = [ vertex_1_dec, vertex_2_dec, vertex_4_dec, vertex_3_dec ]

     percentage_FOV ( sky_map,ra_vertices,dec_vertices )
                                                                                                                               
     print ' The integrated probability in the selected FOV \n at RA = ', str(ra)+'°', 'and DEC =',  str(dec)+'°', 'is', str(('% .4e' % percentage_poly))+'%.'
     print ''

     airmass ( ra, dec, lat_input, lon_input, height_input, time_input, airmass_min = 1, airmass_max = 5.8 )
          
     print ''
     print '---[The FOV is displayed in Aladin plane as **pointing**]---'

     # building command script for Aladin: get FoV(pointing)
     FOV_center = [ ra, dec ] 
     FOV_center = '  '.join(map(str, FOV_center))
     file_script_FOV = 'get FoV(pointing)' + ' ' + FOV_center
     
     aladinSAMP.send_script ( file_script_FOV )
     print ''
     
     # building command script for Aladin: draw string (ra, dec, "percentage_poly")
     FOV_center = [ ra + (FOV_height/4.0), dec + (FOV_height/4.0) ]
     FOV_center = ' , '.join(map(str, FOV_center))
     draw_string_percentage_poly = 'draw string' + '( ' + FOV_center +','+str(('% .1e' % percentage_poly))+'%)'
     

     aladinSAMP.send_script( draw_string_percentage_poly )


     # VizieR Queries (astroquery.vizier)

     from astroquery.vizier import Vizier

     import astropy.coordinates as coord

     import astropy.units as u

     FOV_base_str = str(FOV_base)+'d'
     FOV_height_str = str(FOV_height)+'d'

     result = Vizier.query_region(coord.SkyCoord(ra = ra, dec = dec,unit = (u.deg, u.deg),frame='icrs'),
                                  width = FOV_height_str,height = FOV_base_str,catalog = [ catalog ])

     result.pprint()
     
     try:
          print(result[ catalog ])
     except TypeError as type_error:
          print 'No Galaxy find in the selected FoV using the Catalog ID:', catalog
     finally:
          pass


     print '     ===================================================================== '

    
     
     # -----------sequence of instrument FOVs from a predetermined sky position------------------#
     #-----------------North/South/East/West directions are permitted-------------------#

     while True:
     
          print ' '
          print ' '
          direction=raw_input(' Press (N/S/E/W) to add a contiguous FOVs in North/South/East/West directions. \n Press R to add a new FOV center; \n Press Q to exit. ') 
          print ' '
          print ' '
           
          if direction.strip().capitalize() == 'N':
               
               add_pointing_N_dec = p0_N[1] + FOV_height 

               # variable cycle; not saved
               p0_N[1] = add_pointing_N_dec

               vertices_FOV ( ra, add_pointing_N_dec, FOV_base, FOV_height )

               # sky coord FOV center
               coord = [ ra, add_pointing_N_dec ]

               # FOV corners
               vtx_1 =  [ vertex_1_ra, vertex_1_dec ]
               vtx_2 =  [ vertex_2_ra, vertex_2_dec ]
               vtx_3 =  [ vertex_3_ra, vertex_3_dec ]
               vtx_4 =  [ vertex_4_ra, vertex_4_dec ]
               
               print '---> The new FOV is centered at RA = ', str(('% .5f' % ra))+'°', 'and DEC =', str(('%.5f' % add_pointing_N_dec))+'°. <---'
               print ' '
               
               #vectors of ra and dec FOV corners
               ra_vertices = [ vertex_1_ra, vertex_2_ra,  vertex_4_ra, vertex_3_ra ]
               dec_vertices= [ vertex_1_dec, vertex_2_dec, vertex_4_dec, vertex_3_dec ]
               
               percentage_FOV( sky_map, ra_vertices, dec_vertices )

               print '---> The integrated probability in the FOV is', str(('% .4e' % percentage_poly))+'%. <---'
               print ''
               
               airmass( ra, add_pointing_N_dec, lat_input, lon_input, height_input, time_input, airmass_min=1, airmass_max=5.8 )

               print ''
               print '---[The FOV is displayed in Aladin plane as **pointing~x**]---'

               print ''

               # building command script for Aladin: get FoV(pointing)
               FOV_center = [ ra, add_pointing_N_dec ] 
               FOV_center = '  '.join(map(str, FOV_center))
               file_script_FOV = 'get FoV(pointing)' + ' ' + FOV_center
               
               aladinSAMP.send_script( file_script_FOV )

               # building command script for Aladin: draw string (ra, dec, "percentage_poly")
               FOV_center = [ra+(FOV_height/4.0), add_pointing_N_dec+(FOV_height/4.0)]
               FOV_center = ' , '.join(map(str, FOV_center))
               draw_string_percentage_poly = 'draw string' + '( ' + FOV_center +','+str(('% .1e' % percentage_poly))+'%)'
               
               aladinSAMP.send_script( draw_string_percentage_poly )

               # VizieR Queries (astroquery.vizier)
               import astropy.coordinates as coord

               FOV_base_str = str(FOV_base)+'d'
               FOV_height_str = str(FOV_height)+'d'
     
               result = Vizier.query_region(coord.SkyCoord(ra=ra, dec=add_pointing_N_dec,unit=(u.deg, u.deg),frame='icrs'),
                                  width=FOV_height_str,height=FOV_base_str,catalog=[ catalog ])

               result.pprint()
     
               try:
                    print(result[ catalog ])
               except TypeError as type_error:
                    print 'No Galaxy find in the selected FoV using the Catalog ID:', catalog
               finally:
                    pass
               
               print ''
               print '     ===================================================================== '
               
               
          elif direction.strip().capitalize() == 'S':

               add_pointing_S_dec = p0_S[1] - FOV_height 

               # variable cycle; not saved
               p0_S[1] = add_pointing_S_dec

               vertices_FOV(ra, add_pointing_S_dec, FOV_base, FOV_height)

               # sky coord FOV center
               coord = [ ra, add_pointing_S_dec ]

               # FOV corners
               vtx_1 =  [ vertex_1_ra, vertex_1_dec ]
               vtx_2 =  [ vertex_2_ra, vertex_2_dec ]
               vtx_3 =  [ vertex_3_ra, vertex_3_dec ]
               vtx_4 =  [ vertex_4_ra, vertex_4_dec ]
               
               print '---> The new FOV is centered at RA = ', str(('% .5f' % ra))+'°', 'and DEC =', str(('%.5f' % add_pointing_S_dec))+'°. <---'
               print ''

               #vectors of ra and dec FOV corners
               ra_vertices = [ vertex_1_ra, vertex_2_ra,  vertex_4_ra, vertex_3_ra ]
               dec_vertices= [ vertex_1_dec, vertex_2_dec, vertex_4_dec, vertex_3_dec ]
               
               percentage_FOV ( sky_map, ra_vertices, dec_vertices )

               print '---> The integrated probability in the FOV is', str(('% .4e' % percentage_poly))+'%. <---'
               print ''

               airmass ( ra, add_pointing_S_dec, lat_input, lon_input, height_input, time_input, airmass_min = 1, airmass_max = 5.8 )

               print ''
               print '---[The FOV is displayed in Aladin plane as **pointing~x**]---'

               print ''

               # building command script for Aladin: get FoV(pointing)
               FOV_center = [ ra, add_pointing_S_dec ] 
               FOV_center = '  '.join(map(str, FOV_center))
               file_script_FOV = 'get FoV(pointing)' + ' ' + FOV_center
     
               aladinSAMP.send_script( file_script_FOV )

               # building command script for Aladin: draw string (ra, dec, "percentage_poly")
               FOV_center = [ra+(FOV_height/4.0), add_pointing_S_dec+(FOV_height/4.0)]
               FOV_center = ' , '.join(map(str, FOV_center))
               draw_string_percentage_poly = 'draw string' + '( ' + FOV_center +','+str(('% .1e' % percentage_poly))+'%)'

               aladinSAMP.send_script( draw_string_percentage_poly )

               # VizieR Queries (astroquery.vizier)
               import astropy.coordinates as coord
               
               FOV_base_str = str(FOV_base)+'d'
               FOV_height_str = str(FOV_height)+'d'
     
               result = Vizier.query_region(coord.SkyCoord(ra=ra, dec=add_pointing_S_dec,unit=(u.deg, u.deg),frame='icrs'),
                                  width=FOV_height_str,height=FOV_base_str,catalog=[ catalog ])

               result.pprint()
     
               try:
                    print(result[ catalog ])
               except TypeError as type_error:
                    print 'No Galaxy find in the selected FoV using the Catalog ID:', catalog
               finally:
                    pass
               
               print ''
               print '     ===================================================================== '
               print ''

          
          elif direction.strip().capitalize() == 'E':
               
               add_pointing_E_ra = p0_E[0] + offset_est_west

               # variable cycle; not saved
               p0_E[0] = add_pointing_E_ra

               vertices_FOV( add_pointing_E_ra, dec, FOV_base, FOV_height )

               # sky coord FOV center
               coord = [ add_pointing_E_ra, dec ]

               # FOV corners
               vtx_1 =  [ vertex_1_ra, vertex_1_dec ]
               vtx_2 =  [ vertex_2_ra, vertex_2_dec ]
               vtx_3 =  [ vertex_3_ra, vertex_3_dec ]
               vtx_4 =  [ vertex_4_ra, vertex_4_dec ]

               print '---> The new FOV is centered at RA = ', str(('% .5f' % add_pointing_E_ra))+'°', 'and DEC =', str(('%.5f' % dec))+'°. <---'
               print ''

               #vectors of ra and dec FOV corners
               ra_vertices = [ vertex_1_ra, vertex_2_ra,  vertex_4_ra, vertex_3_ra ]
               dec_vertices= [ vertex_1_dec, vertex_2_dec, vertex_4_dec, vertex_3_dec ]
               
               percentage_FOV( sky_map, ra_vertices, dec_vertices )

               print '---> The integrated probability in the FOV is', str(('% .4e' % percentage_poly))+'%. <---'
               print ''

               airmass( add_pointing_E_ra, dec, lat_input, lon_input, height_input, time_input, airmass_min = 1, airmass_max = 5.8 )

               print ''
               print '---[The FOV is displayed in Aladin plane as **pointing~x**]---'


               # build command script for Aladin: get FoV(pointing)
               FOV_center = [ add_pointing_E_ra, dec ] 
               FOV_center = '  '.join(map(str, FOV_center))
               file_script_FOV = 'get FoV(pointing)' + ' ' + FOV_center
     
               aladinSAMP.send_script( file_script_FOV )
               print ''

               # building command script for Aladin: draw string (ra, dec, "percentage_poly")
               FOV_center = [add_pointing_E_ra + (FOV_height/4.0), dec + (FOV_height/4.0)]
               FOV_center = ' , '.join(map(str, FOV_center))
               draw_string_percentage_poly = 'draw string' + '( ' + FOV_center +','+str(('% .1e' % percentage_poly))+'%)'

               aladinSAMP.send_script( draw_string_percentage_poly )

               # VizieR Queries (astroquery.vizier)
               import astropy.coordinates as coord

               FOV_base_str = str(FOV_base)+'d'
               FOV_height_str = str(FOV_height)+'d'
     
               result = Vizier.query_region(coord.SkyCoord(ra=add_pointing_E_ra, dec=dec,unit=(u.deg, u.deg),frame='icrs'),
                                  width=FOV_height_str,height=FOV_base_str,catalog=[ catalog ])

               result.pprint()
     
               try:
                    print( result[ catalog ])
               except TypeError as type_error:
                    print 'No Galaxy find in the selected FoV using the Catalog ID:', catalog
               finally:
                    pass

               
               print '     ===================================================================== '
               print ''

                                                                      
          elif direction.strip().capitalize() == 'W':

               add_pointing_O_ra = p0_O[0] - offset_est_west

               # variable cycle; not saved
               p0_O[0] = add_pointing_O_ra

               vertices_FOV ( add_pointing_O_ra, dec, FOV_base,FOV_height )

               # sky coord FOV center
               coord = [ add_pointing_O_ra, dec ]

               # FOV corners
               vtx_1 =  [ vertex_1_ra, vertex_1_dec ]
               vtx_2 =  [ vertex_2_ra, vertex_2_dec ]
               vtx_3 =  [ vertex_3_ra, vertex_3_dec ]
               vtx_4 =  [ vertex_4_ra, vertex_4_dec ]

               print '---> The new FOV is centered at RA = ', str(('% .5f' % add_pointing_O_ra))+'°', 'and DEC =', str(('%.5f' % dec))+'°. <---'
               print ''

               #vectors of ra and dec FOV corners
               ra_vertices = [ vertex_1_ra, vertex_2_ra,  vertex_4_ra, vertex_3_ra ]
               dec_vertices= [ vertex_1_dec, vertex_2_dec, vertex_4_dec, vertex_3_dec ]
               
               percentage_FOV( sky_map, ra_vertices, dec_vertices )

               print '---> The integrated probability in the FOV is', str(('% .4e' % percentage_poly))+'%. <---'
               print ''

               airmass( add_pointing_O_ra, dec, lat_input, lon_input, height_input, time_input, airmass_min = 1, airmass_max = 5.8 )

               print ''
               print '---[The FOV is displayed in Aladin plane as **pointing~x**]---'

               # build command script for Aladin: get FoV(pointing)
               FOV_center = [ add_pointing_O_ra, dec ] 
               FOV_center = '  '.join(map(str, FOV_center))
               file_script_FOV = 'get FoV(pointing)' + ' ' + FOV_center
     
               aladinSAMP.send_script( file_script_FOV )
               print ''

               # building command script for Aladin: draw string (ra, dec, "percentage_poly")
               FOV_center = [ add_pointing_O_ra+(FOV_height/4.0), dec+(FOV_height/4.0) ]
               FOV_center = ' , '.join(map(str, FOV_center))
               draw_string_percentage_poly = 'draw string' + '( ' + FOV_center +','+str(('% .1e' % percentage_poly))+'%)'

               aladinSAMP.send_script( draw_string_percentage_poly )

               # VizieR Queries (astroquery.vizier)
               import astropy.coordinates as coord

               FOV_base_str = str(FOV_base)+'d'
               FOV_height_str = str(FOV_height)+'d'
     
               result = Vizier.query_region(coord.SkyCoord(ra=add_pointing_O_ra, dec=dec,unit=(u.deg, u.deg),frame='icrs'),
                                  width=FOV_height_str,height=FOV_base_str,catalog=[ catalog ])


               result.pprint()
     
               try:
                    print( result[ catalog ])
               except TypeError as type_error:
                    print 'No Galaxy find in the selected FoV using the Catalog ID:', catalog
               finally:
                    pass
               
               print '     ===================================================================== '
               print ''
          
 
          # run a new cycle
          elif direction.strip().capitalize() == 'R':
               while True:
                    try:
                         ra_new, dec_new = input(' Insert a new FOV center RA[deg], DEC[deg]: ')
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


               # the new ra and dec are iteratively replaced
               config_GWsky[ 'ra' ] = ra_new
               config_GWsky[ 'dec' ] = dec_new
     
               with open( 'config_GWsky.ini', 'wb' ) as data:
                    pickle.dump ( config_GWsky, data )
          
     
               add_FOV( sky_map,FOV_base,FOV_height,ra,dec )

          # closes all cycles one after another
          elif direction.strip().capitalize() == 'Q':
               break
     
