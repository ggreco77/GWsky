def progress_bar():

     import time                                                                         
     import sys                                                                          
                                                                                 
     toolbar_width = 35                                                                
                                                                                 
     # setup toolbar                                                                   
     sys.stdout.write("[%s]" % (" " * toolbar_width))   
     sys.stdout.flush()                                
     sys.stdout.write("\b" * (toolbar_width+1))        
                                                  
     for i in xrange(toolbar_width):                   
         time.sleep(0.1)                               

         # update the bar                              
         sys.stdout.write("*")                         
         sys.stdout.flush()                            
                                                       
     sys.stdout.write("\n")    
