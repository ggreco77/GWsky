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


__author__ = "Giuseppe Greco <giuseppe.greco@uniurb.it>"

"""

The interactive script GWsky (v2) defines a sequence of Fields of View (FoV)
from a fixed position over the sky.

North/South/East/West directions are allowed.

The results are displayed in Aladin Sky Atlas using the |SAMPIntegratedClient| class
(http://aladin.u-strasbg.fr/).

The airmass at the FoV center and the integrated probability (%)
are provided during the FoV sequence.

***Running it***
from idle:  execfile('GWsky.py')
from terminal: python GWsky.py


SUMMARY OF DEPENDENCIES:
-----------------------
from astropy.vo.samp import SAMPIntegratedClient
import urlparse
import os.path
import config
----> install "sudo apt-get install python-configobj"
from astropy.io.votable import parse
from math import sin, cos, acos, degrees, radians
import healpy
import numpy
import time                                                                   
import sys 

Input Parameters:
-----------

sky_map : str
	 LVC probability skymap in healpix format.

percentage : float
        probability percentage to determine the area (in square degrees) confined in it.

FOV_size : float
	 size of Field-of-View Instrument (square). 

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

# Script functions
import send_Aladin_image as sAi
import send_Aladin_script as sAs

import print_area_prob as pap
import table_ipix_percentage as tip

import highest_probability_pixel as hpp
import instrument_FOV as iF

import airmass as airmass
from gw_core import *

import progress_bar as pb


     # write and read parameters from a config file
from configobj import ConfigObj

# create an empty config file "config_GWsky"
config = ConfigObj(infile='config_GWsky', create_empty=True)

# add values to config file "config_GWsky"
config.filename = 'config_GWsky'

# read value from config file "config_GWsky"
config = ConfigObj('config_GWsky')




#-----------------Launching the Script-------------------------------
#                     GWsky (v2)                                    |
#--------------------------------------------------------------------
# from idle:  execfile('GWsky_v2.py')
# from terminal: python GWsky_v2.py


# load a probability sky map
sky_map = ""
while sky_map=="":
         sky_map = raw_input(' Load a probability skymap: \n ')
         

# write and read parameters from a config file: sky_map
config['sky_map'] = sky_map
infile = config['sky_map']
config.write()

# Send to Aladin Sky Atlas the probability skymap using SAMP client
try:
     sAi.send_Aladin_image(sky_map)
except:
     print ''
     print '    -------------------------------------------------------------'
     print '    ***Open Aladin Sky Atlas before running the script***'
     print '    -------------------------------------------------------------'
     print ''
     raise

pb.progress_bar()

# probability percentage
percentage=input(' Determine the area (in square degrees) \n confined in a given probability percentage. \n Insert a number from 0 to 1: ')

# write and read parameter "percentage"
     # from the config file "config_GWsky"
config['percentage'] = percentage
percentage = config['percentage'] 
config.write()

pb.progress_bar()

# print the area (in square degrees) confined
     # in a given percentage of probability.
pap.print_area_prob(sky_map,percentage)


pb.progress_bar()

print ' The table that contained those pixels is displayed in Aladin \n < contour_ipix.out >'


#build the table that contained those pixels 
tip.table_ipix_percentage(sky_map,percentage)

# send to Aladin the contour_ipix.out table
sAi.send_Aladin_image('contour_ipix.out')

pb.progress_bar()


try:
     #rectangle FOV
     FOV_ba se, FOV_height = input(' Insert the size of your FOV instrument [deg]: ')

# write and read in config file "config_GWsky" the parameter "FOV_base"
     config['FOV_base'] = FOV_base
     FOV_base = config['FOV_base']
     config.write()
     
# write and read in config file "config_GWsky" the parameter "FOV_height"
     config['FOV_height'] = FOV_height
     FOV_height = config['FOV_height']
     config.write()

except:
     print ''
     print '    -------------------------------------------------------------'
     print '    ***Please insert the size of your Field Of View in degrees*** '
     print '                   ex. for a 3° x 2° FOV: 3,2    '
     print '    -------------------------------------------------------------'
     print ''
     raise

pb.progress_bar()

# Modify the file output of Instrument Footprint Editor
iF.instrument_FOV(FOV_base, FOV_height)

#progress_bar()

# INPUT values for airmass calculation
time_input=''                                            
while time_input=='':                                    
     time_input=raw_input(' For airmass calculation insert the time in the format \n "2012-7-12 23:00:00": ')

# write in config file "config_GWsky" the parameter time_input
config['time_input'] = time_input
time_input=config['time_input']
config.write()

try:
    
     lat_input, lon_input, height_input = input(' and the Geodetic coordinates of the Observatory \n latitude[deg], longitude[deg], altitude[m]: ')
except:
     print ''
     print '    -------------------------------------------------------------'
     print '    ***Please insert the Geodetic coordinates of the Observatory*** '
     print '        ex. for Cerro Paranal Observatory: -24.63, -70.40, 2635     '
     print '    -------------------------------------------------------------'
     print ''
     raise

# write in config file "config_GWsky" the parameters: time_input,
      #lat_input, lon_input,  height_input
config['time_input'] = time_input
time_input = config['time_input']
config.write()

config['lat_input'] = lat_input
lat_input = config['lat_input']
config.write()

config['lon_input'] = lon_input
lon_input = config['lon_input']
config.write()

config['height_input'] = height_input
height_input = config['height_input']
config.write()


pb.progress_bar()


# sequence of instrument FOVs from a predetermined sky position
     # North/South/East/West directions are permitted


# input FOV center to generate a sequence of FOV
print ' Specify a sky position to generate a FOV sequence,'

# print ' or insert the highest probability pixel located
     # \n at RA =', str(( '% .5f' % ra_max))+'°', 'and Dec =',str(( '% .5f' % dec_max))+'°.' 

#find the highest probability pixel RA[deg], DEC[deg]
hpp.highest_probability_pixel(sky_map)

print ''

ra,dec = input(' RA[deg], DEC[deg]: ')

print ''

# write in config file "config_GWsky" the parameters ra and dec
config['ra'] = ra
ra =  config['ra']
config.write()

config['dec'] = dec
dec =  config['dec']
config.write()


add_FOV(infile,FOV_base,FOV_height,ra,dec)

pb.progress_bar



#---------------------------end-------------------------------------------------
