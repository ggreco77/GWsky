# -*- coding: utf-8 -*-

from __future__ import print_function

try:
   import cPickle as pickle
except:
   import pickle
   
# Python 3 support
try:
    from Tkinter import *
    import tkMessageBox
    import tkFont
except ImportError:
    from tkinter import *
    from tkinter import font, messagebox

import healpy as hp
import numpy as np


from aladinSAMP import AladinScriptCommands 
aladin = AladinScriptCommands()

from config_values import UserValues

# global variable: level of trasparency window
user = UserValues() 
trasparency = user.get_win_trasparency()


class LoadSkymap(Toplevel):
    """Loading a new skymap."""
    
    def __init__(self):
        Toplevel.__init__(self, border=8, bg="slate grey")
        
        self.user = UserValues()

        # get trasparency windows
        self.wait_visibility()
        self.wm_attributes('-alpha', trasparency)   

        self.title("Load a new skymap")
        self.attributes("-topmost", True)
        
        self.label_1 = Label(self, text="LVC skymap", bg="slate grey")
        self.label_1.grid(row=0, column=0, sticky=E, pady=0)

        # default: input skymap
        skymap_input = StringVar(self, value=self.user.get_skymap()) 
        self.entry_new_skymap = Entry(self, width=30, justify=CENTER,
                             textvariable=skymap_input)
        self.entry_new_skymap.grid(row=0, padx=15, column=1)

        #Btns
        self.show = Button(self, text='Load',
                           command=self.new_skymap)
        self.show.grid(column=2, row=0, sticky=W, padx=2, pady=5)
        
        self.close = Button(self, text="Close",
                            command=self.close_window)  
        self.close.grid(column=5,row=0, sticky=E, padx=2, pady=5)

    
    def new_skymap(self):
        """Loading a new LVC skymap."""
        
        try:
            aladin.send_file(self.entry_new_skymap.get()) 
            aladin.rename(self.entry_new_skymap.get())
        except ValueError as value_error:
            tkMessageBox.showerror ('Load a new skymap',
                                    value_error)
        except IOError as io_error:
            tkMessageBox.showerror ('Load a new skymap',
                                    io_error)

        def update_GWsky_config():
            """Updating GWsky_config file: coords max probability pixel and nside."""
        
            prob = hp.read_map(self.entry_new_skymap.get(), verbose = False)      

            # update nside
            npix = len(prob)
            nside = hp.npix2nside(npix)

            # update coord. maximum prob pixel
            ipix_max = np.argmax(prob)
            prob[ipix_max] 
            theta, phi = hp.pix2ang(nside, ipix_max)
        
            ra_max = round(np.rad2deg(phi), 5) 
            dec_max = round(np.rad2deg(0.5 * np.pi - theta), 5)
        
            with open('GWsky_config', 'rb') as data:
               config_GWsky = pickle.load(data)
            
            config_GWsky['skymap'], config_GWsky['nside'],config_GWsky['ra_max_pixel'],config_GWsky['dec_max_pixel']=\
                                    self.entry_new_skymap.get(), nside, ra_max, dec_max

            with open('GWsky_config', 'wb') as data:
               pickle.dump(config_GWsky, data)

        update_GWsky_config()
        
    def close_window(self):
        return self.destroy()
