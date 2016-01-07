def send_file(infile):
     
     """

         Sending a file to Aladin Sky Atlas using the SAMPIntegratedClient class.
             http://docs.astropy.org/en/stable/vo/samp/example_table_image.html

     """
     
     from astropy.vo.samp import SAMPIntegratedClient
     
     client = SAMPIntegratedClient()
     client.connect()

     params = {}
     import urlparse
     import os.path
     params["url"] = urlparse.urljoin('file:',
				 os.path.abspath(infile))

     message = {}
     message["samp.mtype"] = "image.load.fits"
     message["samp.params"] = params
     
     client.notify_all(message)

     client.disconnect()


def send_script(file_script):

     """

         Sending a script to Aladin Sky Atlas using the SAMPIntegratedClient class.
               http://docs.astropy.org/en/stable/vo/samp/example_table_image.html

     """

     from astropy.vo.samp import SAMPIntegratedClient
     
     client = SAMPIntegratedClient()
     client.connect()

     params = {}
     message = {} 
     message["samp.mtype"] = "script.aladin.send"
     message["samp.params"] = {"script":file_script}  

     client.notify_all(message)

     client.disconnect()



