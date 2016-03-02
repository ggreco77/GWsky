# GWsky: tiling the skymap in Fields of View
                             
GWsky is a Python script to generate iteratively pointings given a specific Field of View (FoV).
The script aims to split up into several independent areas the large error boxes of gravitational wave triggers.

GWsky defines a sequence of FoVs from a fixed position over the sky, *e.g*., starting from the highest probability pixel. The intercardinal and cardinal directions are allowed. The results are displayed in Aladin Sky Atlas (http://aladin.u-strasbg.fr/) using the SAMPIntegratedClient class. The airmass and the integrated probability are provided in real time. Moreover, specifying the ID of a catalog, a query request to the Vizier database is sent; the relative items are listed in each FoV. 
    
The FoVs are evenly spaced assuming that the shortest angular distance between two points on the celestial sphere is measured along a great circle that passes through both of them:

                            cosθ=sinδ1sinδ2+cosδ1cosδ2cos(α1−α2), 
where (α1,δ1) and (α2,δ2) are the right ascension and declination of the two points on the sky.


**Running it**

    from idle:  
              >>> import GWsky 
              >>> GWsky.main() 
    
    from terminal: ./GWsky
    if ./GWsky: Permission denied; type: chmod u+x GWsky


**Input Parameters:**


    sky_map : str 
    LVC probability skymap in healpix format

    percentage : float
    probability percentage to determine the area (in square degrees) confined in it

    FOV_base   : float
    FOV_height : float
    size of Field of View in degrees
    
    catalog : str
    ID of a catalog for VizieR Query
    
    time_input : str
    the time in the format "2012-7-12 23:00:00"

    lat_input : float
    Geodetic coordinates of the Observatory: latitude (deg)

    lon_input : float
    Geodetic coordinates of the Observatory: longitude (deg)

    height_input : float
    Geodetic coordinates of the Observatory: altitude (m)

    ra : float
    right ascention of FoV center (deg)

    dec : float
    declination of FoV center (deg)
    
    position_candidate : float
    sky position in degrees; RA[deg], DEC[deg];
    a specific object to display in the Aladin plan
    
    ID : str
    object ID to display in the Aladin plan


**Return:**

    area_probability : float
    the area (in square degrees) confined in a given probability percentage

    percentage_poly : float
    the integrated probability in a specific FOV (%)

    airmass : float
    the airmass of the FOV center

    contour_ipix.out: list
    the table that contained the pixels confined in a given probability percentage

    ra_max : float
    right ascention of the highest probability pixel
    
    dec_max : float
    declination of the highest probability pixel

    instrument_FOV.vot : VOTABLE
    Instrument Footprint Editor from http://aladin.u-strasbg.fr/footprint_editor/
    
    result : <class 'astroquery.utils.commons.TableList'>
    result of the query

    N/S/E/W/R/Q : str
    N/S/E/W: a set of command line to add contiguous FOVs in North/South/East/West  (N/S/E/W) directions;
    R: to insert a new FOV center RA[deg], DEC[deg] for a new sequence in N/S/E/W directions;
    Q: quit
    
    option : str
    Y or any key: a set of command line to add interactively any specific objects
    
**SUMMARY OF DEPENDENCIES**

    from astroquery.vizier import Vizier
    import astropy.coordinates 
    import astropy.units 
    from astropy.vo.samp import SAMPIntegratedClient
    import urlparse
    import os.path
    from astropy.io.votable import parse
    import healpy
    import numpy

