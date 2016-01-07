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
