def ra0ra1_distance( A, dec0, dec1 ):

     """

     from the angular distance
     cos(A) = sin(Dec1)sin(Dec2)+cos(Dec1)cos(Dec2)cos(ra1-ra2) and thus, A = arccos(A);

                                 --------------->

     cos(ra1-ra2) = [cos(A)-sin(dec0)sin(dec1)]/[cos(dec0)cos(dec1)]
     and  thus, acos(cos(ra1-ra2))
      
     """
     
     from math import cos, sin, acos, degrees, radians
     
     dec0, dec1, A = radians( dec0 ),  radians ( dec1 ), radians( A )
     
     cos_ra0_ra1 = ( cos( A ) - sin( dec0 ) * sin( dec1 ) ) / ( cos( dec0 ) * cos( dec1 ) )

     ra0ra1 = degrees( acos( cos_ra0_ra1 ))

     return  ra0ra1

#ValueError


def shift_direction_east_west( dec ):

     '''

     shifting the est and the west direction
     by an input value
     
     '''
     
     while True:
          try:
               print " Shift the RA position (+/-) "
               print " otherwise type '0' "
               shift = float( input ( " [deg] : " ) )
          except SyntaxError as syntax_error :
               print '', syntax_error
          except NameError as name_error:
               print '', name_error
          except TypeError as type_error:
               print '', type_error
          except ValueError as value_error:
               print '', value_error
          else:
               break
                   
     if shift < 0:
          shift = ( ra0ra1_distance( shift, dec, dec )) * (-1)
     elif shift == 0:
          shift = 0
     else:
          shift = ra0ra1_distance( shift, dec, dec )
                    
     return shift


def shift_direction_nord_sud( dec ):

     '''

     shifting the nord and the sud direction
     by an input value
     
     '''

     while True:
          try:
               print " Shift the DEC position (+/-) "
               print " otherwise type '0' "
               shift = float( input ( " [deg] : " ) )
          except SyntaxError as syntax_error :
                print '', syntax_error
          except NameError as name_error:
               print '', name_error
          except TypeError as type_error:
               print '', type_error
          except ValueError as value_error:
               print '', value_error
          else:
               break
          
     return shift


def shape_FOV( base, width ):

     '''

     checking the FOV shape: square or rectangular FOV
     to determine the diagonal A
      
     '''  
     
     from math import sqrt
     
     if width != base:
          diagonal = sqrt( width**2 + base**2 )
     else:
          diagonal = base * sqrt( 2.0 )

     return diagonal
