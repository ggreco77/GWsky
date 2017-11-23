# -*- coding: utf-8 -*-

#from __future__ import print_function

# py2 and py3 compatibility
try:
    from Tkinter import *
except ImportError:
    from tkinter import *
    
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp

from moc_region import MOC_confidence_region
moc = MOC_confidence_region()

from msg_pyversion import MSG

from lvc_skymap import LVCskymap
from lvc_skymap import healpixIpix

from aladinSAMP import AladinScriptCommands 
aladin = AladinScriptCommands()

from config_values import UserValues

# global variable: transparency level of windows
user = UserValues() 
trasparency = user.get_win_trasparency()

class LocalizeSources(Toplevel):
    """The class is designed to answer if an astrophysical source falls within a
        specific level of probability (MOC contour plot)."""
    
    def __init__(self):
        Toplevel.__init__(self, border=8, bg="slate grey")

        # get trasparency windows from global variable
        self.wait_visibility()
        self.wm_attributes('-alpha', trasparency)   

        self.title(" Localize Sources in probability skymap")
        self.attributes("-topmost", True)

        # label 1
        self.label_1 = Label(self, text="Is/Are the [ID Source(s) RA (°) DEC (°)]",
                             bg="slate grey")
        self.label_1.grid(row=0, column=0, sticky=E, pady=0)

        # showing example entries
        help_string = ("source_1 079.91935 +43.34139; source_2 063.87703 +33.36612 ")
        
        help_string_default = StringVar(self, value=help_string) 
        self.entry_sources = Entry(self, width=30, justify=CENTER,
                             textvariable=help_string_default)
        self.entry_sources.grid(row=0, padx=15, column=1)

        # label 2
        self.label_2 = Label(self, text="whithin the", bg="slate grey")
        self.label_2.grid(row=0, column=3, sticky=E, pady=0)

        moc_value = 90  # default     
        moc_default = StringVar(self, value=moc_value)
        
        self.entry_percentage = Entry(self, width=5, justify=CENTER,
                             textvariable=moc_default)
        self.entry_percentage.grid(row=0, padx=2, column=5)

        # label 3
        self.label_2 = Label(self, text="% MOC?  ",
                             bg="slate grey")
        self.label_2.grid(row=0, column=6, sticky=E, pady=0)

        # label 4
        folder = "transients"  # default     
        folder_default = StringVar(self, value=folder)
        
        self.entry_folder = Entry(self, width=15, justify=CENTER,
                             textvariable=folder_default)
        self.entry_folder.grid(row=1, padx=2, column=0)

        # label 3.1
        self.label_3 = Label(self, text="  Folder: ",justify=LEFT,
                             bg="slate grey")
        self.label_3.grid(row=1, column=0, sticky=W, pady=0)

        self.entryScroll = Scrollbar(self, orient=HORIZONTAL,
                                     command=self.__scrollHandler)
        self.entryScroll.grid(row=1, column=1, sticky=E+W)
        self.entry_sources['xscrollcommand'] = self.entryScroll.set

        # -------------Btns--------------#

        # Ask
        self.show = Button(self, text='Ask',
                           command=self.in_skymap)
        self.show.grid(column=7, row=0, sticky=W, padx=2, pady=5)

        # Pinpoint
        self.checkbox = Button(self, text="Pinpoint", fg='black',     
                               command=self.pinpoint)
        self.checkbox.grid(column=8,row=0, sticky=E, padx=2, pady=5)

        # Dist
        self.checkbox = Button(self, text="Dist", fg='black',     
                               command=self.cond_distance_source_for)
        self.checkbox.grid(column=15,row=0, sticky=E, padx=2, pady=5)

        # Close
        self.close = Button(self, text="Close",
                            command=self.close_window)  
        self.close.grid(column=25,row=0, sticky=W, padx=2, pady=5)

    #Actions
    def __scrollHandler(self, *L):
        """Scroll entry."""
       
        op, howMany = L[0], L[1]

        if op == 'scroll':
            units = L[2]
            self.entry_sources.xview_scroll(howMany, units)
        elif op == 'moveto':
            self.entry_sources.xview_moveto(howMany)
           
    def __split_entries_3(self):
        """Splitting the entries in 'id source', 'ra' and 'dec'."""
        
        entry_sources = self.entry_sources.get().replace(';',' ').replace(',',' ').split()

        # defined lists
        label = entry_sources[::3]
        source_ra = entry_sources[1::3]
        source_dec = entry_sources[2::3]

        # Check if the lists are of the same length
        if len(label)==len(source_ra)==len(source_dec):

            return label[0::], source_ra[0::], source_dec[0::]
        else:
            msg_err = 'ENTER: id source ra[deg] dec[deg]'
            MSG.split_entries_3(msg_err)
                              
    def in_skymap(self):
        """Checking if an object falls in a given probability level defined by an user.
           List of sources are inserted in 'self.entry_sources'."""

        aladin.md(self.entry_folder.get()) # creating folder defined by user
        aladin.remove(self.entry_folder.get() + '~1') # removing multiple copy of the folder
                                                      # TO DO BETTER
        # (re)initialization: skymap/nside
        self.user = UserValues()
        skymap = self.user.get_skymap()
        nside = int(self.user.get_nside())

        # getting probability array
        prob = moc.read_prob(skymap)

        # user input: MOC confidence level in percentage
        percentage = float(self.entry_percentage.get())/100.0

        # splitting the entries              
        labels, ra_transients, dec_transients = self.__split_entries_3()       
      
        for ra_transient, dec_transient, label in zip(
           ra_transients, dec_transients, labels):

            # aladin stack organization: draw source position/move in folder
            aladin.draw_newtool(label)
            aladin.draw_source(ra_transient, dec_transient, label)
            aladin.mv(label, self.entry_folder.get())

            # from sky coords to ipix
            ipixs = moc.ipixs_in_percentage(prob, percentage)
            
            try:
                ipix = healpixIpix.find_ipix(ra_transient, dec_transient,
                                             nside)
            
                is_there = ipix in ipixs  # is the ipix within the MOC contour plot defined by user?

                if is_there is True: # ipix found
                    res_true = ("The sky coord" + " " + "ra="+str(ra_transient)+"°," + " " + "dec="+str(dec_transient)+"°" + " " + \
                                "(labels: " + label+")" + " " +  "lies within the" + " " + str(percentage*100)+'%' + " " + "c.l.\n" +"["+skymap+"]")
                    MSG.in_skymap_true(label, res_true)
                else:
                    res_false = ("The sky coord" + " " + "ra="+str(ra_transient)+"°," + " " +"dec="+str(dec_transient)+"°" + " " + \
                                 "(labels: " + label+")" + " " + "is outside the" + " " + str(percentage*100)+'%' + " " + "c.l.\n" + "["+skymap+"]")
                    MSG.in_skymap_false(label, res_false)

            except ValueError as value_error:
                MSG.value_error(value_error)

    def pinpoint(self):
        """Initialize Pinpoint class."""
        
        pinpoint_localize = Pinpoint()
        return pinpoint_localize

    def cond_distance_source_for(self):
        """Plot of the conditional distance distribution  along the line of sight for a list of sources."""
        
        self.lvc = LVCskymap() # basic module for handling LVC skymaps
        
        labels, ra_transients, dec_transients = self.__split_entries_3()
        
        plt.ion()
        for label, ra_transient, dec_transient in zip(
            labels, ra_transients, dec_transients):

           try:
               r, dp_dr = self.lvc.conditional_distance_linesight(
                   float(ra_transient), float(dec_transient))
           
               self.__cond_distance_source(label, r, dp_dr)
               
           except ValueError as value_error:
               MSG.value_error(value_error)          

    def __cond_distance_source(self, label, r, dp_dr):
        """Plot of the conditional distance distribution  along the line of sight."""
        
        # (re)initialization: skymap
        self.user = UserValues() 
        skymap=self.user.get_skymap()
        
        fig, ax = plt.subplots()

        ax.plot(r, dp_dr)
        title_string = label + ':' + ' '+ ' \n conditional distance distribution along the line of sight \n' + '['+skymap+']'
        ax.set_title(title_string,fontsize=10)
        ax.set_xlabel('distance (Mpc)')
        ax.set_ylabel('prob Mpc$^{-1}$')

        plt.show()

    def close_window(self):
        """Closing window"""
        
        return self.destroy()

class Pinpoint(Toplevel):
    """The class is designed to determine in which level of probability a source is localized."""

    def __init__(self):
        Toplevel.__init__(self, border=8, bg="slate grey")
        
        # get trasparency windows from global variable
        self.wait_visibility()
        self.wm_attributes('-alpha', trasparency)   

        self.title(" Pinpoint Localization")
        self.attributes("-topmost", True)

        # label 1
        self.label_1 = Label(self,
                             text=" In which level of probability the source(s) falls/fall", bg="slate grey")
        self.label_1.grid(row=0, column=0, sticky=E, pady=0)

        # showing example entries
        help_string = ("source_3 085.91935 +33.34139; source_4 073.87703 +23.36612 ")
        
        help_string_default = StringVar(self, value=help_string) 
        self.entry_pin = Entry(self, width=30, justify=CENTER,
                             textvariable=help_string_default)
        self.entry_pin.grid(row=0, padx=15, column=1)

        # label 2
        self.label_2 = Label(self, text=" from  ", bg="slate grey")
        self.label_2.grid(row=0, column=3, sticky=E, pady=0) 

        from_default = StringVar(self, value="10") # default

        self.entry_from = Entry(self, width=5, justify=CENTER,
                                textvariable=from_default)
        self.entry_from.grid(row=0, column=4, columnspan=8, pady=2, sticky='WE')

        # label 3
        self.label_3 = Label(self, text=" to  ", bg="slate grey")
        self.label_3.grid(row=0, column=20, sticky=E, pady=0) 

        to_default = StringVar(self, value="90") # default

        self.entry_to = Entry(self, width=5, justify=CENTER,
                                textvariable=to_default)
        self.entry_to.grid(row=0, column=25, columnspan=8, pady=2, sticky='WE')

        # label 4
        self.label_4 = Label(self, text=" grid  ", bg="slate grey")
        self.label_4.grid(row=0, column=40, sticky=E, pady=0) 

        grid_default = StringVar(self, value="10") # default

        self.entry_grid = Entry(self, width=5, justify=CENTER,
                                textvariable=grid_default)
        self.entry_grid.grid(row=0, column=45, columnspan=8, pady=2, sticky='WE')

        # label 4
        folder = "pinpoint"  # default     
        folder_default = StringVar(self, value=folder)
        
        self.entry_folder = Entry(self, width=20, justify=CENTER,
                             textvariable=folder_default)
        self.entry_folder.grid(row=1, padx=2, column=0)

        # label 3.1
        self.label_3 = Label(self, text="  Folder:",justify=LEFT,
                             bg="slate grey")
        self.label_3.grid(row=1, column=0, sticky=W, pady=0)
        
        self.entryScroll = Scrollbar(self, orient=HORIZONTAL,
                                     command=self.__scrollHandler)
        self.entryScroll.grid(row=1, column=1, sticky=E+W)
        self.entry_pin['xscrollcommand'] = self.entryScroll.set

        # -------------Btns--------------#

        # Do
        self.show = Button(self, text='Do',
                           command=self.pinpoint_for)
        self.show.grid(column=55, row=0, sticky=W, padx=2, pady=5)

        # Close
        self.close = Button(self, text="Close",
                            command=self.close_window)  
        self.close.grid(column=70,row=0, sticky=W, padx=2, pady=5)
        
    #Actions
    def __scrollHandler(self, *L):
        """Scroll entry."""
       
        op, howMany = L[0], L[1]

        if op == 'scroll':
            units = L[2]
            self.entry_pin.xview_scroll(howMany, units)
        elif op == 'moveto':
            self.entry_pin.xview_moveto(howMany)
          
    def __split_entries_3(self):
        """Splitting the entries in 'id source', 'ra' and 'dec'."""
        
        entry_sources = self.entry_pin.get().replace(';',' ').replace(',',' ').split()
                                                          # TRY TO DO BETTER: numpy!!!
        # defined lists
        label = entry_sources[::3]
        source_ra = entry_sources[1::3]
        source_dec = entry_sources[2::3]

        # Check if the lists are of the same length
        if len(label)==len(source_ra)==len(source_dec):
            return label[0::], source_ra[0::], source_dec[0::]
        else:
            msg_err = 'ENTER: id source ra[deg] dec[deg]'
            MSG.split_entries_3(msg_err)
    
    def pinpoint_for(self):
        """Finding in which confidence level the sources fall.
           List of sources are inserted in 'self.entry_pin'."""

        aladin.md(self.entry_folder.get()) # creating folder defined by user
        aladin.remove(self.entry_folder.get() + '~1') # removing multiple copy of the folder
                                                       ## TRY TO DO BETTER!!!!
        # splitting the entries
        labels, ra_transients, dec_transients = self.__split_entries_3()
        
        for ra_transient, dec_transient, label in zip(
            ra_transients, dec_transients, labels):

            # aladin stack organization: draw source position/move in folder                       
            aladin.draw_newtool(label)
            aladin.draw_source(ra_transient, dec_transient, label)
            aladin.mv(label, self.entry_folder.get()) 

            try:
                self.pinpoint(ra_transient, dec_transient, label)
            except ValueError as value_error:
                MSG.value_error(value_error)
            
    def pinpoint(self, ra_transient, dec_transient, label):        
        """Finding in which confidence level the source falls.
        
        Input parameters
        ---------------
        ra_transient, dec_transient : float
              sky coordinates in degrees
        label : string
              id source transient
        """
        
        # (re)initialization: skymap/nside
        self.user = UserValues()
        skymap = self.user.get_skymap() 
        nside = int(self.user.get_nside())

        # getting probability array
        prob = moc.read_prob(skymap)

        # user input values: from/to/resolution
        from_percentage = float(self.entry_from.get())/100.0
        to_percentage = float(self.entry_to.get())/100.0
        resolution_percentage = float(self.entry_grid.get())/100.0

        # from sky coords to ipix
        ipix = healpixIpix.find_ipix(ra_transient, dec_transient,
                                     nside)            

        find = "n"         
        while from_percentage <= to_percentage or find =="y":
            ipixs = moc.ipixs_in_percentage(prob, from_percentage)        
            is_there = ipix in ipixs  # is the ipix within the MOC contour plot defined by user?
            
            if is_there != True: # ipix not found                            
                from_percentage = from_percentage + resolution_percentage                
            else:
                find = "y"  # ipix found    
                res_yes = ("The sky coord" + " " + "ra="+str(ra_transient)+"°," + ' ' + "dec="+str(dec_transient)+"°"+" " + "(label:" + label+")" \
                           "lies within the" + " " + str(from_percentage*100)+'%' + " " + "c.l.\n" +"["+skymap+"]")
                MSG.pinpoint_find(label, res_yes)

                return find
            
        # ipix not found [from_percentage -- to_percentage]
        from_percentage = to_percentage
        
        res_no = ("The sky coord" + " " + "ra="+str(ra_transient)+"°," + ' '
                 + "dec="+str(dec_transient)+"°"+" " + "(label:" + label+")" \
                 + " " + "is not localized within the" + " " + str(from_percentage*100)+'%' + " " + "c.l.\n" +"["+skymap+"]")

        MSG.pinpoint_nofind(label, res_no)
                           
    def close_window(self):
        """Closing window"""
        
        return self.destroy()
