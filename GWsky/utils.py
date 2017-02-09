# -*- coding: utf-8 -*-

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


            