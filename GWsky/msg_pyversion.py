# -*- coding: utf-8 -*-

import os
import sys

from tkinter import *
from tkinter import font, messagebox

class MSG(object):
    """Message according with the python version: tkMessageBox/messagebox in localize module."""

    @classmethod
    def split_entries_3(cls, msg_err):
        """Error message for __split_entries_3 method."""

        if sys.version_info[0] !=3:
            return tkMessageBox.showerror('check your list', msg_err)
        else:
            return messagebox.showerror('check your list', msg_err)
        
    @classmethod
    def in_skymap_true(cls, label, res_true):
        """Showing info message for 'in_skymap' method: true."""
        
        if sys.version_info[0] !=3:
            tkMessageBox.showinfo("Localize Result: " + label, res_true)
        else:
            messagebox.showinfo("Localize Result: " + label, res_true)
            
    @classmethod    
    def in_skymap_false(cls, label, res_false):
        """Showing info message for 'in_skymap' method: false."""

        if sys.version_info[0] !=3:
            tkMessageBox.showinfo("Localize Results: " + label, res_false)
        else:
            messagebox.showinfo("Localize Results: " + label, res_false)
            
    @classmethod
    def value_error(cls, value_error):
        """Show error message for 'in_skymap' method: ValueError."""
        
        if sys.version_info[0] !=3:
            tkMessageBox.showerror ('check your list', value_error)
        else:
            messagebox.showerror ('check your list', value_error)

    @classmethod
    def pinpoint_find(cls, label, res_yes):
        """Showing info message for 'pinpoint_find' method: find/yes."""
        if sys.version_info[0] !=3:
            tkMessageBox.showinfo("Localize Result: " + label, res_yes)
        else:
            messagebox.showinfo("Localize Result: " + label, res_yes)

    @classmethod
    def pinpoint_nofind(cls, label, res_no):
        """Showing info message for 'pinpoint_find' method: no find/no."""
        if sys.version_info[0] !=3:       
            tkMessageBox.showinfo("Localize Result: " + label, res_no)
        else:
            messagebox.showinfo("Localize Result: " + label, res_no)


        
