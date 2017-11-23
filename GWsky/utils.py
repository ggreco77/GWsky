# -*- coding: utf-8 -*-

from __future__ import print_function
import sys

from aladinSAMP import AladinViaSAMP, AladinScriptCommands
samp = AladinViaSAMP()
aladin = AladinScriptCommands()


class Utils(object):

    @classmethod
    def create_folders(cls, folders=[]):
        """Creating aladin stack folders."""
         
        for folder in folders:
            aladin.md(folder)

    @classmethod
    def load_user_fov(cls, footprint):
        """Loading user-defined FoV footprint."""
         
        samp.send_file(footprint)

    @classmethod
    def check_version(cls):
        """Checking the python version."""

        if sys.version[0] == "3":                                            
            print('  ===================================================')   
            print('  ***Please, run the script in  python 2.7.x***')         
            print('  ===================================================')   
            sys.exit()

    @classmethod
    def move_to_folder(cls, planes=[], folders=[]):
        """Moving to folders specific Aladin planes"""

        for plane, folder in zip(planes, folders):
            aladin.mv(plane, folder)

    @classmethod
    def delete_pointing(cls, infile ,ra, dec):
        """Deleting input sky coords from an external file; by default "Pointing.txt" """
        import fileinput
        
        for line in fileinput.input(infile, inplace=True):
            if (line.rsplit()[0] != str(ra) or line.rsplit()[1] != str(dec)) :
                print (line.rstrip('\n'))

    @classmethod
    def separation(self, ra1, dec1, ra2, dec2):
        """Return the distance between 2 sky positions [deg]."""
        
        from astropy.coordinates import SkyCoord
        
        pos_1 = SkyCoord(ra1, dec1, frame='icrs',unit='deg')
        pos_2 = SkyCoord(ra2, dec2, frame='icrs',unit='deg')
        
        sep = pos_1.separation(pos_2)
        
        return sep

        #print ('The distance between 2 consecutive FoV centers is', sep.round(6))


            
