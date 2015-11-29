def send_Aladin_image(infile):
     
     """

         Send to Aladin Sky Atlas tables or fits images using SAMP client.

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
