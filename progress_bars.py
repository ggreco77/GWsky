def progress_bar():

     import time                                                                         
     import sys                                                                          
                                                                                 
     toolbar_width = 15                                                               
                                                                                 
                                                                    
     sys.stdout.flush()                                
     sys.stdout.write("\b" * (toolbar_width+1))        
                                                  
     for i in xrange(toolbar_width):                   
         time.sleep(0.005)                               

                                      
         sys.stdout.write("\b")                         
         sys.stdout.flush()                            
                                                       
     sys.stdout.write(" ")
     print ''

def print_separation_line():
     print ''
     print '     ===================================================================== '
     print ''


def print_aladin_plane():
     print '***********************----------------------***********************'
     print '   [The FOV is displayed in Aladin plane as **pointing**]           '
     print '***********************----------------------***********************'


def print_aladin_plane_x():
     print '***********************----------------------***********************'
     print '   [The FOV is displayed in Aladin plane as **pointing~x**]         '
     print '***********************----------------------***********************'


