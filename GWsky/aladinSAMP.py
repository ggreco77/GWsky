#-*- coding: utf-8 -*-


from astropy.vo.samp import SAMPIntegratedClient
import urlparse
import os.path

class AladinViaSAMP(object):

    def __init__(self):
        self._client = SAMPIntegratedClient()
        
    def send_file(self, infile=str()):
        """Sending a file (image or table) to Aladin Sky Atlas using the SAMPIntegratedClient class.
             http://docs.astropy.org/en/stable/vo/samp/example_table_image.html
        """   
     
        self._client.connect()

        params = {}       
        params[ "url" ] = urlparse.urljoin( 'file:',
				 os.path.abspath( infile ) )
        message = {}
        message[ "samp.mtype" ] = "image.load.fits"
        message[ "samp.params" ] = params
     
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

        position = [x, y]
        position = '  '.join(map(str, position))
        
        plane_name = 'P:'+ str(x) + ',' + str(y)
        remove_fov_str= 'rm' + ' ' + plane_name

        return self.send_script_command(remove_fov_str)

        
        




    
   
