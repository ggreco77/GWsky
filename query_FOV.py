# -*- coding: utf-8 -*-
def selected_catalog( ra, dec, FOV_base, FOV_height, catalog ):

     '''

     Vizier.query_region in selected FoV. The hist of the apparent blue magnitude
     and the sky coordinates are provided if GWGC catalog is selected

     '''

     from astroquery.vizier import Vizier 

     import astropy.coordinates as coord
     import astropy.units as u

     import aladinSAMP

     # setting rows limit
     Vizier.ROW_LIMIT = None

     # no reduced query size (rq)
     reduce_query_size = 1
     FOV_base_reduced, FOV_height_reduced = FOV_base * reduce_query_size, FOV_height * reduce_query_size

     FOV_base_str, FOV_height_str = str( FOV_base_reduced ) + 'd', str( FOV_height_reduced ) + 'd'

     result = Vizier.query_region(coord.SkyCoord(ra = ra, dec = dec, unit = (u.deg, u.deg), frame='icrs'),
                                  width = FOV_height_str, height = FOV_base_str, catalog = [ catalog ])

     result.pprint()
     print result.values()
     

# --------------------------------------------------------------------------------
#                           work in progress
#---------------------------------------------------------------------------------

     #if result != [] and catalog == 'VII/267/gwgc':
          

          #import pandas
          #import numpy as np
          #import matplotlib.pyplot as plt

          # converting to pandas
          #for table_name in result.keys():
          #     table = result[ table_name ]

          #table_pandas = table.to_pandas()

          # sky coordinates
          #ra_query_FOV = table_pandas.RAJ2000
     
          #dec_query_FOV = table_pandas.DEJ2000

          # hists of Dist, BMAG and Bmag
          #try:
               #plt.hist(table_pandas.Dist.dropna())
               #plt.title('Histogram of Distance (Mpc)')
               #plt.xlabel('Dist')
               #plt.ylabel('Count')
     
               #plt.show()
          
     
               #plt.hist(table_pandas.BMAG.dropna())
               #plt.title('Histogram of absolute blue magnitude')
               #plt.xlabel('BMAG')
               #plt.ylabel('Count')
     
               #plt.show()

               #plt.hist(table_pandas.Bmag.dropna())
               #plt.title('Histogram of the apparent blue magnitude')
               #plt.xlabel('Bmag')
               #plt.ylabel('Count')

               #plt.show()

          #except KeyError as key_error:
          #     print '',  key_error
          #finally:
          #     pass
     
          # sending to Aladin plane
          #np.savetxt('FOV_query.out',np.c_[ra_query_FOV,dec_query_FOV],delimiter=" ",fmt='%1.5f,%1.5f',header='RA[deg],DEC[deg]',comments=' ')

          #aladinSAMP.send_file( 'FOV_query.out' )
