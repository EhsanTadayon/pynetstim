"""
classes and functions that implement Target and Sample :
Author: Ehsan Tadayon, M.D. [sunny.tadayon@gmail.com / stadayon@bidmc.harvard.edu]
"""



from coordinates import FreesurferCoords, _FreesurferCoord
from surface import Surf, FreesurferSurf
import numpy as np
import os


### TO DO:
### make coords and targets iterable 

class Targets(FreesurferCoords):
    
    def __init__(self, coords, subject, subjects_dir, names=None, colors=None, directions=None):
        
        FreesurferCoords.__init__(self,coords, subject, subjects_dir, guess_hemi=True, names=names, colors=colors)
        
        if directions is None:
            self.directions = np.empty((self.npoints,9),dtype='object')
        else:
            self.directions = directions
            
              
    def __iter__(self):
        return self
        
    def next(self):
        """ make the Targets iterable, calls another class _Target that will keep all the necessary information for that point. Note that Targets bear multiple points and _Target only keep one point. This is to have the ability to speed up the calculations using linear algebra. """
        if self._count>=self.npoints:
            self._count = 0
            raise StopIteration
        else:
            t = _Target(ras_coord = self.coordinates['ras_coords'][self._count,:],
                           voxel_coord = self.coordinates['voxel_coords'][self._count],
                           ras_tkr_coord = self.coordinates['ras_tkr_coords'][self._count,:],
                           fsvoxel_coord = self.coordinates['fsvoxel_coords'][self._count],
                           hemi = self.hemis[self._count],
                           roi = self.rois[self._count],
                           name = self.names[self._count],
                           color = self.colors[self._count],
                           direction= self.directions[self._count,:])
                       
            self._count +=1
            return t
    
    
    def __getitem__(self,idx):
        t = _Target(ras_coord = self.coordinates['ras_coords'][idx,:],
                        voxel_coord = self.coordinates['voxel_coords'][idx],
                        ras_tkr_coord = self.coordinates['ras_tkr_coords'][idx,:],
                        fsvoxel_coord = self.coordinates['fsvoxel_coords'][idx],
                        hemi = self.hemis[idx],
                        roi = self.rois[idx],
                        name = self.names[idx],
                        color = self.colors[idx],
                        direction = self.directions[idx,:])
        return t
    
    
    

class _Target(_FreesurferCoord):
    def __init__(self,ras_coord,voxel_coord, ras_tkr_coord, fsvoxel_coord, hemi, roi, name, color, direction):
        _FreesurferCoord.__init__(self, ras_coord, voxel_coord, ras_tkr_coord, fsvoxel_coord, hemi, roi, name, color)
        self.direction = direction
        
        
    

    
        

        
    
        
        
        