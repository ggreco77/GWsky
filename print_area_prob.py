def print_area_prob(infile,percentage):
     
     """
          Find the area (in square degrees) confined in a given percentage of probability.
          
     """


     import healpy as hp
     import numpy as np
     
     # read healpix sky map
     hpx = hp.read_map(infile,verbose=False)            

     # sort probability array
     sort = sorted(hpx, reverse=True)

     cumsum = np.cumsum(sort)

     # ipix index in percentage
     index,value=min(enumerate(cumsum), key=lambda x: abs(x[1]-percentage)) 

     # Give pixel area (deg) given nside 
     #nside = hp.pixelfunc.get_nside(hpx)
     #area_pixel = hp.nside2pixarea(nside, degrees = True)

     # area pixel
     npix = len(hpx)
     sky_area = 4 * 180**2 / np.pi
     area_pixel=sky_area / npix

     # area confined in a given percentage of probability (deg)
     area_probability = area_pixel * index

     # 2 digit
     area_probability = round(area_probability,2)

     # value in % 
     percentage_frac = percentage * 100

     # from float to string
     percentage_frac = str(percentage_frac)

     print "The area confined in", percentage_frac+'% of probability is', area_probability, "deg^2."


                      
