"""
classes that implement freesurfer files
Author: Ehsan Tadayon [sunny.tadayon@gmail.com/stadayon@bidmc.harvard.edu]
"""

import nibabel as nib
import numpy as np
import os

####################################################################
######                  Surface classes
####################################################################

class Surf():
    
    def __init__(self,surf_file):
        
        self.surf_file = os.path.abspath(surf_file)
        self.vertices, self.faces = self.read_geometry()

    def read_geometry(self):
        vertices, faces = nib.freesurfer.read_geometry(self.surf_file)
        return vertices, faces
        
    def get_coords(self,vertices_num):
        return self.vertices[vertices_num,:]
    


class FreesurferSurf(Surf):
    
    def __init__(self,hemi,surf, subject,subjects_dir):
        
        self.subject = subject
        self.subjects_dir = os.path.abspath(subjects_dir)
        self.surf = surf
        self.hemi = hemi
        
        surf_file = '{subjects_dir}/{subject}/surf/{hemi}.{surf}'.format(subjects_dir=subjects_dir, subject = subject, surf=surf, hemi=hemi)
        Surf.__init__(self,surf_file)
        

####################################################################
######                  Annotation classes
####################################################################

class Annot(object):
    
    def __init__(self,hemi, annot, subject, subjects_dir):
        
        self.hemi = hemi
        self.annot = annot
        self.subject = subject
        self.subjects_dir = subjects_dir
        self._labels, self._ctab, self._structures = nib.freesurfer.read_annot('{subjects_dir}/{subject}/label/{hemi}.{annot}.annot'.format(subjects_dir=subjects_dir, subject=subject, annot=annot, hemi=self.hemi))
       
        
    def get_ctab(self):
        return self._ctab
        
    def get_labels(self):
        return self._labels
        
    def get_structures(self):
        return self._structures
        
    def get_vertices_colors(self,vertices):
        
        labels = self._labels[vertices]
        colors = [self._ctab[x][:3]/255.0 for x in labels]
        return colors
        
    def get_vertices_names(self,vertices):
        
        labels = self._labels[vertices]
        structures = [self._structures[x] for x in labels]
        return structures