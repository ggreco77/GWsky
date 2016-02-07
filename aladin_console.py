# -*- coding: utf-8 -*-

def get_FoV( x, y ):

     '''

     building command script for Aladin console: "get FoV(pointing)".
     It is sent via SAMP to Aladin console.

     '''
     
     import aladinSAMP
     
     position = [ x, y ] 
     position = '  '.join(map(str, position))
     
     FoV_pointing = 'get FoV(pointing)' + ' ' + position

     aladinSAMP.send_script ( FoV_pointing )


def  draw_string( x, y, text ):

     '''

     building command script for Aladin console: "draw string (x, y, text )".
     It is sent via SAMP to Aladin console.

     '''

     import aladinSAMP

     position = [ x, y]
     position = ' , '.join(map(str, position))

     draw_str = 'draw string' + '( ' + position + ', ' + text + ')' # string

     aladinSAMP.send_script( draw_str )


def  draw_string_float( x, y, number ):

     '''

     building command script for Aladin console: "draw string ( x, y, number )".
     It is sent via SAMP to Aladin console; the parameter num is a float.

     '''
     
     import aladinSAMP

     position = [ x, y ]
     position = ' , '.join(map(str, position))
     
     draw_string_number = 'draw string' + '( ' + position +','+str(('% .1e' % number))+'%)' # float
     
     aladinSAMP.send_script( draw_string_number )


def draw_circle( x, y, size = '10arcmin' ):

     '''

     building command script for Aladin console: "draw circle ( x, y, size )".
     It is sent via SAMP to Aladin console.
       
     '''
     
     import aladinSAMP

     position = [ x, y ]
     position = ' , '.join(map(str, position))
     
     draw_circ = 'draw red circle' + '( ' + position + ', ' + size + ')'

     aladinSAMP.send_script( draw_circ )
     

def get_VizieR(catalog):
     
     '''

      building command script for Aladin console: "get VizieR(catalog,allsky)".
      It is sent via SAMP to Aladin console.
      
     '''

     import aladinSAMP
     
     get_vizier = 'get VizieR(' + catalog + ',' + 'allsky' + ')'
     
     aladinSAMP.send_script( get_vizier )
