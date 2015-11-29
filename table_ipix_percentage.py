def table_ipix_percentage(infile,percentage):
          
     """

     Save in ascii format the table that contained the pixels confined in a given probability percentage


     """


     import healpy as hp
     import numpy as np

     #read skymap
     hpx = hp.read_map(infile,verbose=False)

     # number of pixels
     npix = len(hpx)

     # nside: resolution for the HEALPix map
     nside = hp.npix2nside(npix)

     # sort probability array
     sort = sorted(hpx, reverse=True)

     # cumulative sum 
     cumsum = np.cumsum(sort)

     # ipix index in percentage
     index,value = min(enumerate(cumsum), key=lambda x: abs(x[1] - percentage))

     # find ipix confined in  percentage
     
     value_contour = sort[0:index]

     hpx_list = list(hpx)

     i = 0
     table_ipix_contour = []

     for i in value_contour:
             ipix_contour = hpx_list.index(i)
             table_ipix_contour.append(ipix_contour)


     # from index to polar coordinates
     theta,phi=hp.pix2ang(nside, table_ipix_contour)

     # convert these to right ascension and declination in degrees
     ra = np.rad2deg(phi)
     dec = np.rad2deg(0.5 * np.pi - theta)

     # save table, ra and dec of ipix
     np.savetxt('contour_ipix.out',np.c_[ra,dec],delimiter=" ",fmt='%1.5f,%1.5f',header='RA[deg],DEC[deg]',comments=' ')

