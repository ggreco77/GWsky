# GWsky
The interactive script GWsky (v2) defines a sequence of Fields of View (FoV)
from a fixed position over the sky. North/South/East/West directions are allowed.
The results are displayed in Aladin Sky Atlas using the |SAMPIntegratedClient| class
(http://aladin.u-strasbg.fr/).
The airmass at the FoV center and the integrated probability (%)
are provided during the FoV sequence.

**Running it**

from idle:  execfile('GWsky.py')

from terminal: python GWsky.py


**Input Parameters:**


sky_map : str 

*LVC probability skymap in healpix format*

percentage : float

*probability percentage to determine the area (in square degrees) confined in it*

FOV_size : float

*size of Field-of-View Instrument*

time_input : str

*the time in the format "2012-7-12 23:00:00"* 

lat_input : float

*Geodetic coordinates of the Observatory: latitude (deg)*

lon_input : float

*Geodetic coordinates of the Observatory: longitude (deg)*

height_input : float

*Geodetic coordinates of the Observatory: altitude (m)*

ra : float

*right ascention of FoV center (deg)*

dec : float

*declination of FoV center (deg)*


**Return:**

area_probability : float

*the area (in square degrees) confined in a given probability percentage*

percentage_poly : float

*the integrated probability in a specific FOV (%)*

airmass : float

*the airmass of the FOV center*

contour_ipix.out: list

*the table that contained the pixels confined in a given probability percentage*

ra_max : float

*right ascention of the highest probability pixel*

dec_max : float

*declination of the highest probability pixel*

instrument_FOV.vot : VOTABLE

*Instrument Footprint Editor from http://aladin.u-strasbg.fr/footprint_editor/*

N/S/E/W/R/Q : str

    N/S/E/W: a set of command line to add a contiguous FOVs in North/South/East/West  (N/S/E/W) directions;
  
    R: to insert a new FOV center RA[deg], DEC[deg] to begin a new sequence in N/S/E/W directions;
  
    Q: quit
     



