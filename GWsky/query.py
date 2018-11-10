
#from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
from astropy import units as u

#from astroquery.vizier import Vizier
#import astropy.units as u
import astropy.coordinates as coord

class Query(object):

   
    """Querying galaxy within an user-defined FoV."""

    def query_box(self,ra, dec, base, height, catalog): #constrain
        """Vizier.query_region in a square/rectangular FoV using astroquery module."""

        base_str, height_str = str(base) + 'deg', str(height) + 'deg'
        catalog = str(catalog)
        
        #catalog_constraints = Vizier(catalog=[catalog],
        #                             columns=['*', '_RAJ2000', '_DEJ2000', column_1],
        #                             column_filters={column_1 : filter_1})
        
        #catalog_constraints.ROW_LIMIT = -1  # completed catalog
        query_result = Vizier.query_region(coord.SkyCoord(ra=ra, dec=dec,
                                                          frame='icrs'), width=base_str, height=height_str,
                                           catalog=[catalog])[0]
        #print(query_result)
        return query_result

#query_box(20,20,20,20,"VII/281/glade2")

##
    def query_circle(self, ra, dec, radius, catalog):
        """Vizier.query_region in a circle FoV  using astroquery module."""

        
##        
        radius_str = str(radius) + 'd'
        catalog = str(catalog)
##        
##        #catalog_constraints = Vizier(catalog=[catalog],
##        #                             columns=['*', 'RAJ2000', 'DEJ2000', column_1],
##        #                             column_filters={column_1 : filter_1})
##        
##        #catalog_constraints.ROW_LIMIT = -1  # completed catalog
        query_result = Vizier.query_region(coord.SkyCoord(ra = ra, dec = dec, unit = (u.deg, u.deg),
                                                          frame='icrs'), radius = radius_str,
                                           catalog=[catalog])[0]
                                           
        
        return query_result
##
##
#query=Query()
#query.query_box(20,20,20,20,"VII/275/glade1")


