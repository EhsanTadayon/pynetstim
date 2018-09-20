"""
Main stimulation project functions and classes
"""

from .utils import Coords



class StimProject(object):
    
    """ main stimulation project class"""
    def __init__(self, subject, bids_dir, fs_dir):
        
        self.subject = subject
        self.bids_dir = bids_dir
        self.fs_dir = fs_dir
        
    def add_stim_coord(self,stim_coord_file):
        
        self.stim_coord_file = stim_coord_file
        
        