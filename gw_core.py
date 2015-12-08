# -*- coding: utf-8 -*-



def percentage_FOV(infile,ra_vertices,dec_vertices):

     """

     Give the probability contains in a specific FOV
     
     """

     global  percentage_poly
     
     import healpy as hp
     import numpy as np

     # read probability skymap
     hpx = hp.read_map(infile,verbose=False)

     # number of pixels
     npix = len(hpx)

     # nside: resolution for the HEALPix map
     nside = hp.npix2nside(npix)

     # Spherical polar coordinates of vertices in radians
     theta = 0.5 * np.pi - np.deg2rad(dec_vertices)

     phi = np.deg2rad(ra_vertices)
     
     xyz = hp.ang2vec(theta, phi)

     # Array of indices of pixels inside polygon
     ipix_poly = hp.query_polygon(nside, xyz)

     # probability contains in a specific FOV
     percentage_poly=hpx[ipix_poly].sum()


def vertices_FOV(ra,dec,FOV_base,FOV_height):

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

     # new FOV_base; FOV_height -----> FOV_height =  FOV_height
     
     #-------------------------------------------------------------------------
     # same declination is assumed; from degree to rad (+ FOV_height / 2.0)
     delta0_rads = radians(dec) + radians(FOV_height / 2.0) 
     delta1_rads = delta0_rads 

     # x = cos(alfa0_rads-alfa1_rads) in _vertices_FOV function
     x_vertices_FOV=(cos(radians(FOV_base/2.0))-sin(delta0_rads)*sin(delta1_rads))/(cos(delta0_rads)*cos(delta1_rads))

     # offset to obtain a specific distance in degrees on the sky
     offset_N = degrees(acos(x_vertices_FOV))
     #--------------------------------------------------------------------------

     

     #-------------------------------------------------------------------------
     # same declination is assumed; from degree to rad (- FOV_height / 2.0)
     delta0_rads = radians(dec) - radians(FOV_height / 2.0  ) 
     delta1_rads = delta0_rads 

     # x = cos(alfa0_rads-alfa1_rads) in _vertices_FOV function
     x_vertices_FOV=(cos(radians(FOV_base/2.0))-sin(delta0_rads)*sin(delta1_rads))/(cos(delta0_rads)*cos(delta1_rads))

     # offset to obtain a specific distance in degrees on the sky
     offset_S = degrees(acos(x_vertices_FOV))
     #--------------------------------------------------------------------------

     # FOV vertices 
     vertex_1_ra, vertex_1_dec = ra + offset_N, dec + FOV_height / 2.0   

     vertex_2_ra, vertex_2_dec = ra - offset_N, dec + FOV_height / 2.0

     vertex_3_ra, vertex_3_dec = ra + offset_S, dec - FOV_height / 2.0   

     vertex_4_ra, vertex_4_dec = ra - offset_S, dec - FOV_height / 2.0



def add_FOV(infile,FOV_base,FOV_height,ra,dec):

     """
          Define a sequence of instrument FOVs from a predetermined sky position.
          North/South/East/West directions.
          The FOV centers and the corners are displayed in Aladin Sky Atlas.
     
     """
     
     import send_Aladin_image as sAi
     import send_Aladin_script as sAs
     from airmass import airmass
     
     import healpy as hp
     import numpy as np
     from math import sin, cos, acos, degrees, radians
     
     
     # read value from config file "config_GWsky"
     from configobj import ConfigObj
     config = ConfigObj('config_GWsky')

     infile = config['sky_map']
     
     time_input = config['time_input']
     
     lat_input =  config['lat_input']
     lat_input=float(lat_input)
     
     lon_input =  config['lon_input']
     lon_input = float(lon_input)
     
     height_input =  config['height_input']
     height_input=float(height_input)

     ra =  config['ra']
     ra=float(ra)
     
     dec =  config['dec']
     dec=float(dec)
     
     # from degree to rad     
     delta0_rads = radians(dec)

     delta1_rads = delta0_rads 

     # x=cos(alfa0_rads-alfa1_rads)
     x=(cos(radians(FOV_height))-sin(delta0_rads)*sin(delta1_rads))/(cos(delta0_rads)*cos(delta1_rads))

     x_est_west=(cos(radians(FOV_base))-sin(delta0_rads)*sin(delta1_rads))/(cos(delta0_rads)*cos(delta1_rads))

     # offset to obtain a specific distance in degrees on the sky
     offset = degrees(acos(x))

     offset_est_west = degrees(acos(x_est_west))

     

     # initializing variables
     p0_N = [ra, dec]
     p0_S = [ra, dec]
     p0_O = [ra, dec]
     p0_E = [ra, dec]


     #---------------- FOV at the sky coordinates insert from users; (ra dec)---------------------#


     # Find the vertices of FOV 
     vertices_FOV(ra,dec,FOV_base,FOV_height)

     # sky coord FOV center
     coord = [ra, dec]

     # FOV corners
     vtx_1 =  [vertex_1_ra, vertex_1_dec]
     vtx_2 =  [vertex_2_ra, vertex_2_dec]
     vtx_3 =  [vertex_3_ra, vertex_3_dec]
     vtx_4 =  [vertex_4_ra, vertex_4_dec]

     # save the center and the corners of initial FOV in "Initial_FOV.out" and send the sky positions to Aladin
     np.savetxt('starting_FOV.out',[coord,vtx_1,vtx_2,vtx_3,vtx_4],fmt='%1.5f,%1.5f',header='RA[deg],DEC[deg]',comments=' ')
     sAi.send_Aladin_image('starting_FOV.out')

     # vectors of ra and dec FOV corners
     ra_vertices = [vertex_1_ra, vertex_2_ra,  vertex_4_ra, vertex_3_ra]
     dec_vertices= [vertex_1_dec, vertex_2_dec, vertex_4_dec, vertex_3_dec]

     # probability contains in a specific FOV
     percentage_FOV(infile,ra_vertices,dec_vertices)

     
                                                                                                                               
     print ' The integrated probability in the selected FOV \n at RA = ', str(ra)+'°', 'and DEC =',  str(dec)+'°', 'is', str(('% .4e' % percentage_poly))+'%.'
     print ''

     # airmass calculation at a given time in a particular site
     airmass(ra, dec, lat_input, lon_input, height_input, time_input, airmass_min=1, airmass_max=5.8)
     #print 'dimmi qualcosa', airmass_value
     
     
     print ''
     print '---[The FOV is displayed in Aladin plane as **starting_FOV**]---'

     # build command script for Aladin: get FoV(pointing)
     FOV_center = [ra, dec] 
     FOV_center = '  '.join(map(str, FOV_center))
     file_script_FOV = 'get FoV(pointing)' + ' ' + FOV_center
     
     # send script to Aladin
     sAs.send_Aladin_script(file_script_FOV)
     print ''

     
     # build command script for Aladin: draw string (ra, dec, "percentage_poly")
     FOV_center = [ra+(FOV_height/4.0), dec+(FOV_height/4.0)]
     FOV_center = ' , '.join(map(str, FOV_center))
     draw_string_percentage_poly = 'draw string' + '( ' + FOV_center +','+str(('% .1e' % percentage_poly))+'%)'

     # send script to Aladin
     #send_Aladin_script(draw_string_percentage_poly)
     sAs.send_Aladin_script(draw_string_percentage_poly)


     # VizieR Queries (astroquery.vizier)
     #import FOV_query as Fq

     #FOV_base_str=str(FOV_base)+'d'
     #FOV_height_str=str(FOV_height)+'d'
     
     #Fq.FOV_query(ra=ra,dec=dec,width=FOV_base_str,height=FOV_height_str,catalog_id='GWGC')



     print '     ===================================================================== '

    
     
     # -----------sequence of instrument FOVs from a predetermined sky position------------------#
     #-----------------North/South/East/West directions are permitted-------------------#
    
     while True:
     
          print ' '
          print ' '
          direction=raw_input(' Press (N/S/E/W) to add a contiguous FOVs in North/South/East/West directions. \n Press R to add a new FOV center; \n Press Q to exit (closes all cycles one after another). ') 
          print ' '
          print ' '
           
          if direction == 'N':

               # set Dec coordinate for North direction
               add_pointing_N_dec = p0_N[1] + FOV_height 

               # variable cycle; not saved
               p0_N[1] = add_pointing_N_dec


               # Find the vertices of FOV 
               vertices_FOV(ra, add_pointing_N_dec, FOV_base, FOV_height)


               # sky coord FOV center
               coord = [ra, add_pointing_N_dec]


               # FOV corners
               vtx_1 =  [vertex_1_ra, vertex_1_dec]
               vtx_2 =  [vertex_2_ra, vertex_2_dec]
               vtx_3 =  [vertex_3_ra, vertex_3_dec]
               vtx_4 =  [vertex_4_ra, vertex_4_dec]
               
          
               # save the center and the corners of selected FOV in "pointing.out" and send the sky positions to Aladin
               np.savetxt('pointing',[coord,vtx_1,vtx_2,vtx_3,vtx_4],fmt='%1.5f,%1.5f',header='RA[deg],DEC[deg]',comments=' ')
               sAi.send_Aladin_image('pointing')
               
               print '---> The new FOV is centered at RA = ', str(('% .5f' % ra))+'°', 'and DEC =', str(('%.5f' % add_pointing_N_dec))+'°. <---'
               print ' '
               

               #vectors of ra and dec FOV corners
               ra_vertices = [vertex_1_ra, vertex_2_ra,  vertex_4_ra, vertex_3_ra]
               dec_vertices= [vertex_1_dec, vertex_2_dec, vertex_4_dec, vertex_3_dec]
               
               # probability contains in a specific FOV
               percentage_FOV(infile, ra_vertices, dec_vertices)

               print '---> The integrated probability in the FOV is', str(('% .4e' % percentage_poly))+'%. <---'
               print ''
               
               # airmass calculation at a given time in a particular site
               airmass(ra, add_pointing_N_dec, lat_input, lon_input, height_input, time_input, airmass_min=1, airmass_max=5.8)

               print ''
               print '---[The FOV is displayed in Aladin plane as **pointing~x**]---'

               print ''

               # build command script for Aladin: get FoV(pointing)
               FOV_center = [ra, add_pointing_N_dec] 
               FOV_center = '  '.join(map(str, FOV_center))
               file_script_FOV = 'get FoV(pointing)' + ' ' + FOV_center
               
     
               # send script to Aladin
               sAs.send_Aladin_script(file_script_FOV)



               # build command script for Aladin: draw string (ra, dec, "percentage_poly")
               FOV_center = [ra+(FOV_height/4.0), add_pointing_N_dec+(FOV_height/4.0)]
               FOV_center = ' , '.join(map(str, FOV_center))
               draw_string_percentage_poly = 'draw string' + '( ' + FOV_center +','+str(('% .1e' % percentage_poly))+'%)'

               # send script to Aladin
               #send_Aladin_script(draw_string_percentage_poly)
               sAs.send_Aladin_script(draw_string_percentage_poly)

               
               print ''
               print '     ===================================================================== '
               
               
          elif direction == 'S':

               # set Dec coordinate for South direction
               add_pointing_S_dec = p0_S[1] - FOV_height 

               # variable cycle; not saved
               p0_S[1] = add_pointing_S_dec

               # Find the vertices of FOV
               vertices_FOV(ra, add_pointing_S_dec, FOV_base, FOV_height)

               # sky coord FOV center
               coord = [ra, add_pointing_S_dec]

               # FOV corners
               vtx_1 =  [vertex_1_ra, vertex_1_dec]
               vtx_2 =  [vertex_2_ra, vertex_2_dec]
               vtx_3 =  [vertex_3_ra, vertex_3_dec]
               vtx_4 =  [vertex_4_ra, vertex_4_dec]

               # save the center and the corners of selected FOV in "pointing.out" and send the sky positions to Aladin
               np.savetxt('pointing',[coord,vtx_1,vtx_2,vtx_3,vtx_4],fmt='%1.5f,%1.5f',header='RA[deg],DEC[deg]',comments=' ')
               sAi.send_Aladin_image('pointing')
               
               print '---> The new FOV is centered at RA = ', str(('% .5f' % ra))+'°', 'and DEC =', str(('%.5f' % add_pointing_S_dec))+'°. <---'
               print ''

               #vectors of ra and dec FOV corners
               ra_vertices = [vertex_1_ra, vertex_2_ra,  vertex_4_ra, vertex_3_ra]
               dec_vertices= [vertex_1_dec, vertex_2_dec, vertex_4_dec, vertex_3_dec]
               
               # probability contains in a specific FOV
               percentage_FOV(infile, ra_vertices, dec_vertices)

               print '---> The integrated probability in the FOV is', str(('% .4e' % percentage_poly))+'%. <---'
               print ''

               # airmass calculation at a given time in a particular site
               airmass(ra, add_pointing_S_dec, lat_input, lon_input, height_input, time_input, airmass_min=1, airmass_max=5.8)

               print ''
               print '---[The FOV is displayed in Aladin plane as **pointing~x**]---'

               print ''

               # build command script for Aladin: get FoV(pointing)
               FOV_center = [ra, add_pointing_S_dec] 
               FOV_center = '  '.join(map(str, FOV_center))
               file_script_FOV = 'get FoV(pointing)' + ' ' + FOV_center
     
               # send script to Aladin
               sAs.send_Aladin_script(file_script_FOV)


               # build command script for Aladin: draw string (ra, dec, "percentage_poly")
               FOV_center = [ra+(FOV_height/4.0), add_pointing_S_dec+(FOV_height/4.0)]
               FOV_center = ' , '.join(map(str, FOV_center))
               draw_string_percentage_poly = 'draw string' + '( ' + FOV_center +','+str(('% .1e' % percentage_poly))+'%)'

               # send script to Aladin
               #send_Aladin_script(draw_string_percentage_poly)
               sAs.send_Aladin_script(draw_string_percentage_poly)

               
               print ''
               print '     ===================================================================== '
               print ''

          
          elif direction == 'E':

               # set Dec coordinate for east direction
               add_pointing_E_ra = p0_E[0] + offset_est_west

               # variable cycle; not saved
               p0_E[0] = add_pointing_E_ra

               # Find the vertices of FOV
               vertices_FOV(add_pointing_E_ra, dec, FOV_base, FOV_height)

               # sky coord FOV center
               coord = [add_pointing_E_ra, dec]

               # FOV corners
               vtx_1 =  [vertex_1_ra, vertex_1_dec]
               vtx_2 =  [vertex_2_ra, vertex_2_dec]
               vtx_3 =  [vertex_3_ra, vertex_3_dec]
               vtx_4 =  [vertex_4_ra, vertex_4_dec]


               # save the center and the corners of selected FOV in "pointing.out" and send the sky positions to Aladin
               np.savetxt('pointing',[coord,vtx_1,vtx_2,vtx_3,vtx_4],fmt='%1.5f,%1.5f',header='RA[deg],DEC[deg]',comments=' ')
               sAi.send_Aladin_image('pointing')

               print '---> The new FOV is centered at RA = ', str(('% .5f' % add_pointing_E_ra))+'°', 'and DEC =', str(('%.5f' % dec))+'°. <---'
               print ''

               #vectors of ra and dec FOV corners
               ra_vertices = [vertex_1_ra, vertex_2_ra,  vertex_4_ra, vertex_3_ra]
               dec_vertices= [vertex_1_dec, vertex_2_dec, vertex_4_dec, vertex_3_dec]
               
               # probability contains in a specific FOV
               percentage_FOV(infile, ra_vertices, dec_vertices)

               print '---> The integrated probability in the FOV is', str(('% .4e' % percentage_poly))+'%. <---'
               print ''

               # airmass calculation at a given time in a particular site
               airmass(add_pointing_E_ra, dec, lat_input, lon_input, height_input, time_input, airmass_min=1, airmass_max=5.8)

               print ''
               print '---[The FOV is displayed in Aladin plane as **pointing~x**]---'


               # build command script for Aladin: get FoV(pointing)
               FOV_center = [add_pointing_E_ra, dec] 
               FOV_center = '  '.join(map(str, FOV_center))
               file_script_FOV = 'get FoV(pointing)' + ' ' + FOV_center
     
               # send script to Aladin
               sAs.send_Aladin_script(file_script_FOV)
               print ''

               # build command script for Aladin: draw string (ra, dec, "percentage_poly")
               FOV_center = [add_pointing_E_ra+(FOV_height/4.0), dec+(FOV_height/4.0)]
               FOV_center = ' , '.join(map(str, FOV_center))
               draw_string_percentage_poly = 'draw string' + '( ' + FOV_center +','+str(('% .1e' % percentage_poly))+'%)'

               # send script to Aladin
               #send_Aladin_script(draw_string_percentage_poly)
               sAs.send_Aladin_script(draw_string_percentage_poly)



               
               print '     ===================================================================== '
               print ''

                                                                      
          elif direction == 'O':

               # set Dec coordinate for West direction
               add_pointing_O_ra = p0_O[0] - offset_est_west

               # variable cycle; not saved
               p0_O[0] = add_pointing_O_ra

               # Find the vertices of FOV
               vertices_FOV(add_pointing_O_ra, dec, FOV_base,FOV_height)

               # sky coord FOV center
               coord = [add_pointing_O_ra, dec]

               # FOV corners
               vtx_1 =  [vertex_1_ra, vertex_1_dec]
               vtx_2 =  [vertex_2_ra, vertex_2_dec]
               vtx_3 =  [vertex_3_ra, vertex_3_dec]
               vtx_4 =  [vertex_4_ra, vertex_4_dec]


               # save the center and the corners of selected FOV in "pointing.out" and send the sky positions to Aladin
               np.savetxt('pointing',[coord,vtx_1,vtx_2,vtx_3,vtx_4],fmt='%1.5f,%1.5f',header='RA[deg],DEC[deg]',comments=' ')
               sAi.send_Aladin_image('pointing')

               print '---> The new FOV is centered at RA = ', str(('% .5f' % add_pointing_O_ra))+'°', 'and DEC =', str(('%.5f' % dec))+'°. <---'
               print ''

               #vectors of ra and dec FOV corners
               ra_vertices = [vertex_1_ra, vertex_2_ra,  vertex_4_ra, vertex_3_ra]
               dec_vertices= [vertex_1_dec, vertex_2_dec, vertex_4_dec, vertex_3_dec]
               
               # probability contains in a specific FOV
               percentage_FOV(infile, ra_vertices, dec_vertices)

               print '---> The integrated probability in the FOV is', str(('% .4e' % percentage_poly))+'%. <---'
               print ''

               # airmass calculation at a given time in a particular site
               airmass(add_pointing_O_ra, dec, lat_input, lon_input, height_input, time_input, airmass_min=1, airmass_max=5.8)

               print ''
               print '---[The FOV is displayed in Aladin plane as **pointing~x**]---'


               # build command script for Aladin: get FoV(pointing)
               FOV_center = [add_pointing_O_ra, dec] 
               FOV_center = '  '.join(map(str, FOV_center))
               file_script_FOV = 'get FoV(pointing)' + ' ' + FOV_center
     
               # send script to Aladin
               sAs.send_Aladin_script(file_script_FOV)
               print ''

               # build command script for Aladin: draw string (ra, dec, "percentage_poly")
               FOV_center = [add_pointing_O_ra+(FOV_height/4.0), dec+(FOV_height/4.0)]
               FOV_center = ' , '.join(map(str, FOV_center))
               draw_string_percentage_poly = 'draw string' + '( ' + FOV_center +','+str(('% .1e' % percentage_poly))+'%)'

               # send script to Aladin
               #send_Aladin_script(draw_string_percentage_poly)
               sAs.send_Aladin_script(draw_string_percentage_poly)


               
               print '     ===================================================================== '
               print ''
 
           # run a new cycle
          elif direction == 'R':

               ra_new, dec_new = input(' Insert a new FOV center RA[deg], DEC[deg]: ')

               # the new ra and dec are iteratively replaced 
               config['ra'] = ra_new
               config.write()

               config['dec'] = dec_new
               config.write()

               
               add_FOV(infile, FOV_height,ra,dec)

          # closes all cycles one after another
          elif direction == 'Q':
               break
               
          else:
               print ' '
               print '     *** Please, type < N: North; S: South; E: East; W: West; R: new FOV center; Q: exit > ***'

     
                        




