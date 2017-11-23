
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
from astropy import units as u

class Query(object):
    """Galaxy query inside a user-defined FoV."""

    def query_box(self, ra, dec, base, height, catalog, column_1, filter_1):
        """Vizier.query_region in a square/rectangular FoV using astroquery module."""

        base_str, height_str = str(base) + 'd', str(height) + 'd'
        
        catalog_constraints = Vizier(catalog=[catalog],
                                     columns=['*', '_RAJ2000', '_DEJ2000', column_1],
                                     column_filters={column_1 : filter_1})
        
        catalog_constraints.ROW_LIMIT = -1  # completed catalog
        query_result = catalog_constraints.query_region(SkyCoord(
            ra = ra, dec = dec, unit = (u.deg, u.deg), frame='icrs'),
                                                        width = height_str, height = base_str)       
        return query_result

    def query_circle(self, ra, dec, radius, catalog, column_1, filter_1):
        """Vizier.query_region in a circle FoV  using astroquery module."""
        
        radius_str = str(radius) + 'd'
        
        catalog_constraints = Vizier(catalog=[catalog],
                                     columns=['*', '_RAJ2000', '_DEJ2000', column_1],
                                     column_filters={column_1 : filter_1})
        
        catalog_constraints.ROW_LIMIT = -1  # completed catalog
        query_result = catalog_constraints.query_region(SkyCoord(
            ra = ra, dec = dec, unit = (u.deg, u.deg), frame='icrs'), radius = radius_str)
        
        return query_result
