#### Surf module which will take care of freesurfer surface,annot and label files
### Author: Ehsan Tadayon, M.D. [sunny.tadayon@gmail.com / stadayon@bidmc.harvard.edu]
### Start date: Sep 16, 2018

import nibabel as nib
import numpy as np
from scipy.spatial.distance import cdist
import os
## to do list:




############################## Surf class

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
        return self.vertices[indices,:]
    
class FreesurferSurf(Surf):
    
    def __init__(self,hemi,surf, subject,subjects_dir):
        
        self.subject = subject
        self.subjects_dir = os.path.abspath(subjects_dir)
        self.surf = surf
        self.hemi = hemi
        
        surf_file = '{subjects_dir}/{subject}/surf/{hemi}.{surf}'.format(subjects_dir=subjects_dir, subject = subject, surf=surf, hemi=hemi)
        Surf.__init__(self,surf_file)
        
        
############################# Working with coordinates######################
### two classes have been implemented to work with coordinates. 1) Coords class which is more general class that could take care of coordinate in any volume 2) FreesurferCooords which inherits from Coords but requires a subject that has been processed by recon-all.


class Coords():
    
    def __init__(self,coords, img_file, type='ras'):
        """ coordinate class"""
        
        self.img_file = img_file
        self.type = type
        self.img = nib.load(img_file)
        self.vox2ras = self.img.affine
        self.ras2vox = np.linalg.inv(self.vox2ras)
        self.npoints = coords.shape[0]
        self.coordinates = {}
        self._affineM = np.hstack((coords, np.ones((coords.shape[0],1)))).T
        if type=='ras':
            self.coordinates['ras_coords'] = coords
            self.coordinates['voxel_coords'] = np.round(np.dot(self.ras2vox,self._affineM).T[:,:3])
        else:
            self.coordinates['voxel_coords'] = coords
            self.coordinates['ras_coords'] = np.dot(self.vox2ras,self._affineM).T[:,:3]
    
    def img2imgcoord(self,img):
        pass
    
    
class FreesurferCoords(Coords):
    
    def __init__(self,coords,subject,subjects_dir,guess_hemi=True):
        """ Freesurfer Coordinate class"""
        
        self.subjects_dir = subjects_dir
        self.subject = subject
        rawavg_file = '{subjects_dir}/{subject}/mri/rawavg.mgz'.format(subjects_dir=subjects_dir,subject=subject)

        Coords.__init__(self, coords, rawavg_file,type='ras')
        
        orig_file = '{subjects_dir}/{subject}/mri/orig.mgz'.format(subjects_dir=subjects_dir,subject=subject)
        self.orig_img = nib.freesurfer.load(orig_file)
        self.ras2fsvox = self.orig_img.header.get_ras2vox()
        self.fsvox2ras_tkr = self.orig_img.header.get_vox2ras_tkr()
        self.coordinates['ras_tkr_coords'] = np.dot(self.fsvox2ras_tkr,np.dot(self.ras2fsvox,self._affineM)).T[:,:3]
        self.coordinates['fsvoxel_coords'] = np.round(np.dot(self.ras2fsvox,self._affineM).T[:,:3])
        self.coordinates['talairach'] = self._get_talairach_coords()
        
        if guess_hemi:
            self._guess_hemi()
            
    def _guess_hemi(self):
        self.hemis = []
        for s in np.arange(self.npoints):
            if self.coordinates['fsvoxel_coords'][s,0]> 128: 
                 self.hemis.append('lh')   
            elif self.coordinates['fsvoxel_coords'][s,0] < 128:
                 self.hemis.append('rh')   
            else:
                 raise 'Could not determine hemisphere'
                    
        
    def _read_talaraich_transformation(self):
        """ read talairach transformation from freesurfer talairach.xfm output"""
        
        fname = '{subjects_dir}/{subject}/mri/transforms/talairach.xfm'.format(subjects_dir=self.subjects_dir,
                                                                               subject=self.subject)
        f = open(fname,'r').read().split('\n')
        
        ### cleaning rows of file
        def clean_number(x):
            if ';' in x:
                return float(x[0:-1])
            else:
                return float(x)
            
        rows = []
        for i in [5,6,7]:
            row = f[i].split(' ')
            row = map(clean_number,row)
            rows.append(row)
            
        return np.array(rows)
    
    def _get_talairach_coords(self):
        """ transforms the coordinates by talairach transform matrix from freesurfer talairach.xfm"""
        talairach_tr = self._read_talaraich_transformation()
        return np.dot(talairach_tr,self._affineM).T[:,:3]
        


        
