"""
Implements Freesurfer annotation file
"""
import nibabel as nib
import numpy as np
from surface import FreesurferSurf



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
        
    def get_vertices(self,vertices):
        
        labels = self._labels[vertices]
        colors = [self._ctab[x][:3]/255.0 for x in labels]
        structures = [self._structures[x] for x in labels]
        return colors,structures
        
    def map_coords(self, coords, map_surface='white'):
        
        # first map your coords to surface
        if map_surface in ['white','pial']:
            map_surface = FreesurferSurf(self.hemi, map_surface, self.subject, self.subjects_dir)
            
        elif isinstance(map_surface,FreesurferSurf):
            pass
        else:
            raise ValueError('map_surface should be either string(white, pial) or FreesurferSurf instance')
        
        mapped_vertices_indices, mapped_vertices_coords = map_surface.project_coords(coords)
        return self.get_vertices(mapped_vertices_indices)