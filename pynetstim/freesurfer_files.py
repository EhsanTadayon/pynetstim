"""
classes that implement surface and freesurfer surface 
Author: Ehsan Tadayon [sunny.tadayon@gmail.com/stadayon@bidmc.harvard.edu]
"""

import nibabel as nib
import numpy as np
from scipy.spatial.distance import cdist
import os

class Surf():
    
    def __init__(self,surf_file):
        
        self.surf_file = os.path.abspath(surf_file)
        self.vertices, self.faces = self.read_geometry()

    def read_geometry(self):
        vertices, faces = nib.freesurfer.read_geometry(self.surf_file)
        return vertices, faces
    
    def project_coords(self,coords):
        coords = np.atleast_2d(coords)
        indices = np.argmin(cdist(self.vertices, coords), axis=0)
        return indices, self.vertices[indices,:]
        
    def get_coords(self,vertices_num):
        
        return self.vertices[vertices_num]
    
class FreesurferSurf(Surf):
    
    def __init__(self,hemi,surf, subject,subjects_dir):
        
        self.subject = subject
        self.subjects_dir = os.path.abspath(subjects_dir)
        self.surf = surf
        self.hemi = hemi
        
        surf_file = '{subjects_dir}/{subject}/surf/{hemi}.{surf}'.format(subjects_dir=subjects_dir, subject = subject, surf=surf, hemi=hemi)
        Surf.__init__(self,surf_file)