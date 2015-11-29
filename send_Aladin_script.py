def send_Aladin_script(file_script):

     """

         Send to Aladin Sky Atlas an Aladin script using SAMP client.

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

