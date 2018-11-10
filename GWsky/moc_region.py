# -*- coding: utf-8 -*-

import healpy as hp
import numpy as np
from mocpy import MOC

from math import log

from lvc_skymap import LVCskymap

class MOC_confidence_region(object):
    """Multi-Order coverage map (MOC) of sky areas enclosed within a contour plot
    at a given confidence level."""


    def read_prob(self, infile):
        """Reading healpix skymap.
        
        Input parameters
        ----------------
        infile : string
              LVC probability sky localization in healpix format
              
        Return
        -------
        hpx : list
            1D array of values (probability stored in each pixel)
        """      
        
        hpx = hp.read_map(infile, verbose = False)
        
        return hpx
 
    def ipixs_in_percentage(self, hpx, percentage):
        """Finding the ipix indices confined in a given percentage.
        
        Input parameters
        ----------------
        hpx : numpy array
            1D array of values (probability stored in each pixel)
            
        percentage : float
                 fractional percentage from 0 to 1  
        
        Return
        ------- 
        ipixs : numpy array
              indices of pixels
        """

        # ranked the healpix pixels from most probable to least,  and finally counted how many
        # pixels  summed  to  a  given  total  probability.
        # see https://arxiv.org/pdf/1404.5623.pdf

        cumsum = np.sort(hpx)[::-1].cumsum() 
        how_many_ipixs, cut_percentage = min(enumerate(cumsum),
                                             key = lambda x: abs(x[1] - percentage))

        del(cumsum)
        
        index = np.arange(0, len(hpx))        
        hpx_index = np.c_[hpx, index]
        
        sort = hpx_index[hpx_index[:, 0].argsort()[::-1]]
        ipixs = sort[0:how_many_ipixs, [1]].astype(int)

        return ipixs

##    def __ipix_box(self, ra_vertices, dec_vertices):
##        """Return the ipix inside a polygon."""
##
##        ## TO BE COMPLETED
##        
##        NSIDE=512   # fixed nside resolution
##        theta = 0.5 * np.pi - np.deg2rad(dec_vertices)
##        phi = np.deg2rad(ra_vertices)
##        xyz = hp.ang2vec(theta, phi)
##        ipix_poly = hp.query_polygon(NSIDE, xyz)
##
##        return ipix_poly, NSIDE

    def ipix_in_box(self, ra, dec, width, height):
        """Return the probability inside a box."""
        
        self.lvc_skymap = LVCskymap()
        
        v1_ra, v2_ra, v3_ra, v4_ra, v1_dec, v2_dec, v3_dec, v4_dec = self.lvc_skymap.vertices(
            ra, dec, width, height)
        
        ra_vertices, dec_vertices = (
           [v1_ra, v2_ra, v4_ra, v3_ra], [v1_dec, v2_dec, v4_dec, v3_dec])
        
        ipix_fov_box, NSIDE = self.__ipix_box(ra_vertices, dec_vertices)
        
        return ipix_fov_box, NSIDE
        
    def ipix_within_circle(self, ra_vertices, dec_vertices):
        pass
         
    def sky_coords(self, ipixs, nside):
        """Converting the ipix into right ascension and declination in degrees
        
        Return
        ------- 
        contour_ipix : list
                    sky coords in degrees
        """
       
        # from index to polar coordinates
        theta, phi = hp.pix2ang(nside, ipixs)

        # converting these to right ascension and declination in degrees
        ra = np.rad2deg(phi)
        dec = np.rad2deg(0.5 * np.pi - theta)

        
        # creating an astropy.table with RA[deg] and DEC[deg]
        #contour_ipix = Table([ra, dec], names = ('RA[deg]', 'DEC[deg]'), 
        #                     meta = {'ipix': 'ipix table'})
     
        return ra, dec
   
    def moc_order(self, nside):
        """Setting MOC order.
        
        Return
        ------- 
        moc_order : int
              
        """       
        
        order = int(log( nside, 2))
     
        return order

    def create_moc(self):
        """Creating a MOC map from the contour_ipix table."""
        
        self.moc = MOC.from_table(self.contour_ipix, 'RA[deg]', 'DEC[deg]',
                                  self.moc_order)

        return self.moc

    def write_moc(self, percentage, short_name):
        """Writing MOC file in fits format.
        
        Input parameters
        ----------------
        percentage : float
                 fractional percentage from 0 to 1 converted into a string
        short_name : str
                 file output
        """
        
        return self.moc.write(short_name + '_MOC_' + str(percentage), format = 'fits')
    
    def contour_plot(self, infile, percentage, short_name=''):
        """Creating/Writing a MOC contour region at a fixed level of probability.
        
        Input parameters
        ---------------
        infile : string
              LVC probability sky localization in healpix format
        percentage : float
                 fractional percentage from 0 to 1
        """
        
        self.read_skymap(infile)
        self.ipixs_in_percentage(percentage)
        self.sky_coords()
        self.moc_order()
        self.create_moc()

        return self.write_moc(percentage, short_name)

    def contour_default(self, _from, _to, _step, skymap=""):
        """Creating & Showing MOC plots (from 10% to 90% in step of 10%) in a folder."""

        import time
        
        from aladinSAMP import AladinScriptCommands 
        aladin = AladinScriptCommands()

        colors=["#ff0000","#ffaa00 ","#aaff00","#00ff00","#00ffa9",
                "#00a9ff","#0000ff","#aa00ff","#ff00aa"] # skymap viewer color
        
        #short_name = skymap
        suffix = skymap[0:]
        
        aladin.md('MOC' + suffix) # creating a stack folder
        aladin.remove('MOC' + suffix + '~1') # removing multiple copy of the folder

        for i, color in zip(np.arange(_from, _to, _step),colors):
            aladin.cmoc((i/100.0), skymap,'moc'+str(i/100.0) + suffix)

            time.sleep(1) # random break to organize the aladin planes
            
            plane = 'moc'+str(i/100.0) + suffix        
            aladin.set_color(plane, color)
            aladin.set_moc('moc'+str(i/100.0)+ suffix)
            aladin.mv('moc'+str(i/100.0)+ suffix,'MOC' + suffix)
                

        
        
        
       



