
import healpy as hp
import numpy as np
from math import cos, sin, acos, asin, atan, degrees, radians
from scipy.stats import norm

from .config_values import UserValues

from .aladinSAMP import AladinViaSAMP, AladinScriptCommands 
samp = AladinViaSAMP() # potresti toglierlo, ereditato da aladin cone jupyter
aladin = AladinScriptCommands()


class Is3d(object):
    """Determining if it is a 3d skymap by reading the TFIELDS field in the header."""

    def __init__(self, infile):
        self.infile = infile # LVC skymap
    
    def get_header(self):
        """Reading TFIELDS in the header."""
           
        prob, header = hp.read_map(self.infile, h=True, verbose=False)
        header = dict(header)
        tfield = header['TFIELDS']

        return tfield, header

    def get_values(self, tfield):
        """Getting the values stored in the healpix plane(s)."""
        
        if tfield == 4:
            prob, distmu, distsigma, distnorm = hp.read_map(
                self.infile, verbose=False, h=False, field=range(4))

            return prob, distmu, distsigma, distnorm
        else:
            prob = hp.read_map(self.infile, verbose=False, h=False)
            distmu, distsigma, distnorm = [], [], []
            
            return prob, distmu, distsigma, distnorm
        
class LVCskymap(object):
    """A set of methods for working with LVC healpix (3d) skymaps."""

    def __init__ (self):
        
        self.user = UserValues() #comp.
        self.skymap = self.user.get_skymap()
        print (self.skymap) # TEST
        self.nside = self.user.get_nside()

        self.is_3d = Is3d(self.skymap) #comp --> eredita
        self.tfield, self.header = self.is_3d.get_header()
        self.prob, self.distmu, self.distsigma, self.distnorm = self.is_3d.get_values(
            self.tfield)

    
    def vertices(self, ra_center, dec_center, fov_base, fov_height):
        """Finding the vertices of a FoV given the central location (ra[deg], dec[deg])
           and the FoV size (FoV_base [deg], FoV_height [deg])."""
        
        vert_ra, vert_dec = [], []  # ra list,  dec list 
        
        ra_center_rad, dec_center_rad = radians(ra_center), radians(dec_center)

        fov_base_rad, fov_height_rad = radians(fov_base), radians(fov_height)
       
        x = [-fov_base_rad/2, fov_base_rad/2,
             fov_base_rad/2, -fov_base_rad/2]

        y = [fov_height_rad/2, fov_height_rad/2,
             -fov_height_rad/2, -fov_height_rad/2]
        
        for i, j  in zip(x, y):
            arg = -i/(cos(dec_center_rad)-j*sin(dec_center_rad))          
            v_ra = degrees( (ra_center_rad+atan(arg)) )
            
            vert_ra.append(v_ra)
            
            v_dec = degrees( (asin((sin(dec_center_rad)+j*cos(dec_center_rad))/(1+i**2+j**2)**0.5)) )
            
            vert_dec.append(v_dec)

        # test: field-of-view footprint vs vertices function
        #aladin.draw_circle(vert_ra[0], vert_dec[0], size = '15arsec')
        #aladin.draw_circle(vert_ra[1], vert_dec[1], size = '15arsec')
        #aladin.draw_circle(vert_ra[2], vert_dec[2], size = '15arsec')
        #aladin.draw_circle(vert_ra[3], vert_dec[3], size = '15arsec')
        aladin.draw_polygon(vert_ra[0], vert_dec[0],vert_ra[1], vert_dec[1],
                            vert_ra[2], vert_dec[2], vert_ra[3], vert_dec[3])

        #print round(vert_ra[0],6), round(vert_dec[0],6), round(vert_ra[1],6), round(vert_dec[1],6),round(vert_ra[2],6), round(vert_dec[2],6), round(vert_ra[3],6), round(vert_dec[3],6)
        #print vert_ra[0], vert_ra[1], vert_ra[3], vert_ra[2], vert_dec[0], vert_dec[1], vert_dec[3], vert_dec[2]
        return vert_ra[0], vert_ra[1], vert_ra[3], vert_ra[2], vert_dec[0], vert_dec[1], vert_dec[3], vert_dec[2]
        #print vert_ra[0], vert_ra[1], vert_ra[3], vert_ra[2], vert_dec[0], vert_dec[1], vert_dec[3], vert_dec[2]

##    def vertices(self, ra_center, dec_center, fov_base, fov_height):
##        """Finding the vertices of a FoV given the central location (ra[deg], dec[deg])
##           and the FoV size (FoV_base [deg], FoV_height [deg])."""
##        
##        vert_ra, vert_dec = [], []  # ra list,  dec list 
##        
##        ra_center_rad, dec_center_rad = radians(ra_center), radians(dec_center)
##
##        fov_base_rad, fov_height_rad = radians((fov_base+0.00069)/2.0), radians((fov_height+0.00069)/2.0) #0.0017119592190950605
##       
##        x = [-fov_base_rad, fov_base_rad,
##             fov_base_rad, -fov_base_rad]
##
##        y = [fov_height_rad, fov_height_rad,
##             -fov_height_rad, -fov_height_rad]
##
##
##
####        Dato un FOV quadrato di lato L e orientato come dicevo nella mia mail precedente, es:
####       deg2rad=pi/180;
####
####      L=3*deg2rad;
####
####       In coord. ortogonali i vertici del quadrato in senso orario partendo da quello in alto a sx, assumendo
##        #che l'origine sia nel punto centrale del FOV, sono:
####
####
####   X1=(-L/2,L/2,L/2,-L/2);
####
####    Y1=(L/2,L/2,-L/2,-L/2);
####
####
####In coord. angolari diventano:
####
####
####arg1=-X1/(cos(d0)-Y1*sin(d0));
####
####a1=(a0+atan(arg1))/deg2rad
####
####d1=(asin( (sin(d0)+Y1*cos(d0))/(1+X1.^2+Y1.^2).^0.5))/deg2rad
####
####dove (a0,d0) sono le coordinate del centro del FOV
####
####Ho fatto delle prove con Aladin e a me torna ma fai un check...
##
##        
##        for i, j  in zip(x, y):
##            arg = -i/(cos(dec_center_rad)-j*sin(dec_center_rad))          
##            v_ra = degrees( (ra_center_rad+atan(arg)) )
##            
##            vert_ra.append(v_ra)
##            
##            v_dec = degrees( (asin((sin(dec_center_rad)+j*cos(dec_center_rad))/(1+i**2+j**2)**0.5)) )
##            
##            vert_dec.append(v_dec)
##
##        # test: field-of-view footprint vs vertices function
##        aladin.draw_circle(vert_ra[0], vert_dec[0], size = '5arcmin')
##        aladin.draw_circle(vert_ra[1], vert_dec[1], size = '5arcmin')
##        aladin.draw_circle(vert_ra[2], vert_dec[2], size = '5arcmin')
##        aladin.draw_circle(vert_ra[3], vert_dec[3], size = '5arcmin')
##
##        #print round(vert_ra[0],6), round(vert_dec[0],6), round(vert_ra[1],6), round(vert_dec[1],6),round(vert_ra[2],6), round(vert_dec[2],6), round(vert_ra[3],6), round(vert_dec[3],6)
##        #print vert_ra[0], vert_ra[1], vert_ra[3], vert_ra[2], vert_dec[0], vert_dec[1], vert_dec[3], vert_dec[2]
##        return vert_ra[0], vert_ra[1], vert_ra[3], vert_ra[2], vert_dec[0], vert_dec[1], vert_dec[3], vert_dec[2]
##        #print vert_ra[0], vert_ra[1], vert_ra[3], vert_ra[2], vert_dec[0], vert_dec[1], vert_dec[3], vert_dec[2]



    def __ipix_sum(self, ra_vertices, dec_vertices):
        """Return the ipix sum inside a polygon."""
        
        theta = 0.5 * np.pi - np.deg2rad(dec_vertices)
        phi = np.deg2rad(ra_vertices)
        xyz = hp.ang2vec(theta, phi)

        ipix_poly = hp.query_polygon(self.nside, xyz)
     
        ipix_sum_polygon = self.prob[ipix_poly].sum()
         
        return ipix_sum_polygon

    def prob_in_box(self, ra, dec, width, height):
        """Return the probability inside a box."""

        v1_ra, v2_ra, v3_ra, v4_ra, v1_dec, v2_dec, v3_dec, v4_dec = self.vertices(
            ra, dec, width, height)
        
        ra_vertices, dec_vertices = (
           [v1_ra, v2_ra, v4_ra, v3_ra], [v1_dec, v2_dec, v4_dec, v3_dec])
        
        probability_fov_box = self.__ipix_sum(ra_vertices, dec_vertices)
        
        return '%.1e' % probability_fov_box

    def prob_in_circle(self, ra, dec, radius):
        """Return the probability inside a circle."""

        theta = 0.5 * np.pi - np.deg2rad(dec)
        phi = np.deg2rad(ra)
        radius = np.deg2rad(radius)

        xyz = hp.ang2vec(theta, phi)
        ipix_disc = hp.query_disc(self.nside, xyz, radius)
        probability_fov_disc = self.prob[ipix_disc].sum()

        aladin.draw_circle(ra, dec, size = str(self.user.get_fov_radius())+' deg')

        return '%.1e' % probability_fov_disc

    def conditional_distance_linesight(self, ra, dec):
        """Conditional distance distribution along the line of sight of a FoV center:
            see https://arxiv.org/pdf/1605.04242v3.pdf - section 4.4 for more details."""

        if self.tfield==4:  # 3d skymap
            theta = 0.5 * np.pi - np.deg2rad(dec)
            phi = np.deg2rad(ra)
            ipix = hp.ang2pix(self.nside, theta, phi)

            line_end = self.header['DISTMEAN'] + (self.header['DISTSTD']*4)
            r = np.linspace(0, line_end)
            dp_dr = r**2*self.distnorm[ipix]*norm(self.distmu[ipix], self.distsigma[ipix]).pdf(r)
            
            return r, dp_dr
        else: 
            r = "nan"
            dp_dr = "nan"

            return r, dp_dr

class healpixIpix(object):

    @classmethod            
    def find_ipix(cls, ra, dec, nside):
        """Finding ipix."""
        
        theta = 0.5 * np.pi - np.deg2rad(float(dec))
        phi = np.deg2rad(float(ra))
        ipix = hp.ang2pix(nside, theta, phi)
        
        return ipix
