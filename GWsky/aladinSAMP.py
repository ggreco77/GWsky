#!/usr/bin/env python
#-*- coding: iso-8859-15 -*-

import os.path
import sys

if sys.version > '3':
    import urllib.parse
else:
    import urlparse

from astropy.vo.samp import SAMPIntegratedClient


class AladinViaSAMP(object):

    def __init__(self):
        self._client = SAMPIntegratedClient()
        
    def send_file(self, infile=str()):
        """Sending a file (image or table) to Aladin Sky Atlas using the SAMPIntegratedClient class.
             http://docs.astropy.org/en/stable/vo/samp/example_table_image.html
        """   
     
        self._client.connect()
        params = {}

        if sys.version > '3':
            params["url"] = urllib.parse.urljoin('file:', os.path.abspath(infile))
        else:
            params["url"] = urlparse.urljoin('file:', os.path.abspath(infile))
            
        message = {}
        message[ "samp.mtype" ] = "image.load.fits"
        message[ "samp.params" ] = params
     
        self._client.notify_all(message)
        self._client.disconnect()

    def move(self, script=str()):
        """Sending a script to Aladin Sky Atlas using the SAMPIntegratedClient class.
           http://docs.astropy.org/en/stable/vo/samp/example_table_image.html
         """

        self._client.connect()

        params = {}
        message = {} 
        message[ "samp.mtype" ] = "coord.pointAt.sky"
        message[ "samp.params" ] = { "script" : script }  

        self._client.notify_all(message)
        self._client.disconnect()

    def send_script_command(self, script=str()):
        """Sending a script to Aladin Sky Atlas using the SAMPIntegratedClient class.
           http://docs.astropy.org/en/stable/vo/samp/example_table_image.html
         """

        self._client.connect()

        params = {}
        message = {} 
        message[ "samp.mtype" ] = "script.aladin.send"
        message[ "samp.params" ] = { "script" : script }  

        self._client.notify_all(message)
        self._client.disconnect()


class AladinScriptCommands(AladinViaSAMP):
    """A set of the main script commands for Aladin console.
        http://aladin.u-strasbg.fr/java/AladinScriptManual.gml"""

    def cview(self, url):
        """Creation of view: url."""
    
        cview_str = 'cview' + ' ' + url     
        return self.send_script_command(cview_str)

    def draw_circle(self, x, y, size = '10arcmin'):
        """Draw circle (x, y, size)."""
     
        position = [ x, y ]
        position = ' , '.join(map(str, position))  
        draw_circle_str = 'draw red circle' + '( ' + position + ', ' + size + ')'

        return self.send_script_command(draw_circle_str)

    def draw_line(self, line_values):
        """Draw line: .line(x1,y1,x2,y2,...[,text])."""
    
        draw_line_str = 'draw' + ' ' + 'line' + ' ' + line_values   
        return self.send_script_command(draw_line_str)

    def draw_newtool(self, name):
        """Create manually a new plane: draw newtool(name)."""
    
        draw_newtool_str = 'draw' + ' ' + 'newtool' + ' ' + name    
        return self.send_script_command(draw_newtool_str)

    def draw_string(self, x, y, text):
        """Draw string (x, y, text)."""

        position = [x, y]
        position = ' , '.join(map(str, position))
        draw_string_str = 'draw string' + '( ' + position + ', ' + text + ')'

        return self.send_script_command(draw_string_str)

    def draw_string_float(self, x, y, number):
        """Draw string (x, y, number)."""
     
        position = [ x, y ]
        position = ' , '.join(map(str, position))  
        draw_string_number = 'draw string' + '( ' + position +','+str(('% .1e' % number))+'%)'
     
        return send_script_command(self, draw_string_number)

    def get_FoV(self, x, y):
        """P_ra_dec = get FoV(pointing)."""
     
        position = [x, y] 
        position = '  '.join(map(str, position))

        plane_name = 'P:'+ str(x) + ',' + str(y) 
        get_fov_str =  plane_name + '= get FoV(pointing)' + ' ' + position

        return self.send_script_command (get_fov_str)

    def get_VizieR(self, catalog):
        """get VizieR(catalog,allsky)."""  
     
        get_vizier_str = 'get VizieR(' + catalog + ',' + 'allsky' + ')' 
        return self.send_script_command(get_vizier_str)

    def rename(self, plane):
        """Rename plane."""
        
        rename_plane_str = 'rename' + ' ' + plane      
        return self.send_script_command(rename_plane_str)

    def remove_FoV(self, x, y):
        """Remove Field of View"""

        position = [x, y]
        position = '  '.join(map(str, position))
        
        plane_name = 'P:'+ str(x) + ',' + str(y)
        remove_fov_str= 'rm' + ' ' + plane_name

        return self.send_script_command(remove_fov_str)

    def set_planeID(self, name):
        """Set planeID = todo"""

        set_plane_str = 'set' + ' ' + 'planeID=' + name
        return self.send_script_command(set_plane_str)

    def remove(self, plane):
        """Remove plane"""

        rm_str = 'rm' + ' ' + plane
        return self.send_script_command(rm_str)

    def cmoc(self, threshold, skymap, plane_name):

        cmoc_str = plane_name+'='+'cmoc' + ' ' + '-threshold=' + str(threshold) + ' ' + skymap
          
        return self.send_script_command(cmoc_str)

    def set_moc(self, plane):
        """"""
        set_srt = 'set'+' '+plane +' '+'drawing=+perimeter,-border,-fill'

        return self.send_script_command(set_srt)

    def md(self, name):
        """"""

        md_srt = 'md'+ ' '+name

        return self.send_script_command(md_srt)

    def mv(self, plane, folder):
        """"""

        mv_srt = 'mv'+ ' ' + plane + ' ' + folder

        return self.send_script_command(mv_srt)

    def moc_inter(self, moc_1, moc_2, name):
        """"""

        cmoc_inter_str = name + '=' + '' +'cmoc' +' ' + '-inter' + ' ' + moc_1 + ' ' + moc_2

        return self.send_script_command(cmoc_inter_str)

    def zoom(self, factor):

        zoom_str = 'zoom' + ' ' + factor

        return self.send_script_command(zoom_str)

    def location(self, ra, dec):

        location_str = ra + ' ' + dec

        return self.send_script_command(location_str)


    def set_target(self, ra, dec):

        set_target_str = 'set' + ' ' + 'Target' + ' ' + '=' + ' ' + str(ra) + ' ' + str(dec)

        return self.send_script_command(set_target_str)

    def select (self, plane):

        select_str = 'select' + ' ' + str(plane)
        return self.send_script_command(select_str)

    def set_plane_id(self, name):

        set_plane_id_str = 'set PlaneID='+name
        return self.send_script_command(set_plane_id_str)

        
    def setconf_icrsd(self):
        
        setconf_icrsd_str = "setconf" + " " + "frame=ICRSd"
        return self.send_script_command(setconf_icrsd_str)

    def collapse(self, folder):

        collapse_str = "collapse" + " " + folder
        return self.send_script_command(collapse_str)


    def expand(self, folder):
        """"""
        expand_str = "expand" + " " + folder
        return self.send_script_command(expand_str)


    def set_color(self, plane, color):
        """"""
        set_color = "set" + " " + plane + " " + "color =" + color
        return self.send_script_command(set_color)

    def set_roll(self, roll_value):
        """"""
        
        set_roll = "set Roll =" + str(roll_value)
        return self.send_script_command(set_roll)


#    def draw_moon(self, ra, dec, illumination):
#        """"""
        
#        draw_moon_str = 'draw rgb(224,224,224) tag(' + str(ra) + str(dec) + ',' + 'MOON ,20 ,60 , bigreticle ,15)'


    def draw_moon(self, x, y, illumination):
        """Draw string for Moon illumination."""

        position = [x, y]
        position = ' , '.join(map(str, position))
        
        #draw_string_moon = 'draw string' + '( ' + position + ', ' + text + ')'
        
        draw_moon_str = 'draw rgb(224,224,224) tag(' + position + ',' + 'MOON' + '-->' + str(illumination) +',20 ,60 , bigreticle ,15)'
        
        return self.send_script_command(draw_moon_str)  
        


    def draw_airmass(self, x, y, airmass_step):
        """Draw string for Moon illumination."""

        position = [x, y]
        position = ' , '.join(map(str, position))
        
        #draw_string_moon = 'draw string' + '( ' + position + ', ' + text + ')'
        
        draw_airmass_str = 'draw rgb(224,224,224) tag(' + position + ',' + 'Airmass' + str(airmass_step) +',20 ,60 , bigreticle ,15)'
        
        return self.send_script_command(draw_airmass_str)


    def draw_source(self, x, y, label):
        """Draw string for Moon illumination."""

        position = [x, y]
        position = ' , '.join(map(str, position))
        
        #draw_string_moon = 'draw string' + '( ' + position + ', ' + text + ')'
        
        draw_source_str = 'draw' + ' ' + 'tag(' + position + ','  + label +',20 ,60 , bigreticle ,15)'
        
        return self.send_script_command(draw_source_str)

    def load(self, image, name):
         """"""
         load_str = name + '=' + ' '+ 'load'+ ' ' + image

         return self.send_script_command(load_str)
       
        




    
   
