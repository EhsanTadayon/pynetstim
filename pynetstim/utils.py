#### Surf module which will take care of freesurfer surface,annot and label files
### Author: Ehsan Tadayon, M.D. [sunny.tadayon@gmail.com / stadayon@bidmc.harvard.edu]
### Start date: Sep 16, 2018

import nibabel as nib
import numpy as np
from scipy.spatial.distance import cdist
import os
from nipype.interfaces.fsl import WarpPoints
import warnings

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
        
        
################################ Freesurfer annotation       
    
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
    

        
        
############################# Working with coordinates######################
### 


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
            
        self._count=0
     
    
    def img2imgcoord(self, dest_img, warp_file):
       
        np.savetxt('./temp_coords.txt',self.coordinates['ras_coords'])
        warppoints = WarpPoints()
        warppoints.inputs.in_coords = './temp_coords.txt'
        warppoints.inputs.src_file = self.img_file
        warppoints.inputs.dest_file = dest_img
        warppoints.inputs.warp_file = warp_file
        warppoints.inputs.coord_mm = True
        res = warppoints.run()
        res = np.loadtxt('./temp_coords_warped.txt')
        
        ## removing the files
        os.remove('./temp_coords.txt')
        os.remove('./temp_coords_warped.txt')
        
        return res
        
    def __iter__(self):
        return self
    
    def next(self):
        
        if self._count>=self.npoints:
            self._count = 0
            raise StopIteration
            
        else:
            c = _Coord(ras_coord,voxel_coord)
            self._count+=1
            return c
        
    
    
class MNICoords(Coords):
    
    def __init__(self,coords,mni_template='MNI152_T1_2mm.nii.gz',mni_directory=os.environ['FSLDIR']):
        
        mni_file = os.path.join(mni_directory,mni_tmeplate)
        Coords.__init__(self,coords,mni_file)
        
            
        
    
    
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
        self.hemis_not_determined = []
        for s in np.arange(self.npoints):
            if self.coordinates['fsvoxel_coords'][s,0]> 128: 
                self.hemis.append('lh')   
            elif self.coordinates['fsvoxel_coords'][s,0] < 128:
                self.hemis.append('rh')   
            else:
                warnings.warn('Could not determine hemisphere for point {x},{y},{z}. Right hemisphere has been chosen arbitrarily for this point. Manually set the hemisphere for this point by calling set_hemi_manually!'.format(x=self.coordinates['ras_coords'][s,0], y=self.coordinates['ras_coords'][s,1], z=self.coordinates['ras_coords'][s,2]))
                
                self.hemis_not_determined.append(s)
                self.hemis.append('rh')
                 
        self.hemis = np.array(self.hemis)
        
    def set_hemi_manually(self, n, hemi):
        self.hemis[n] = hemi
                    
        
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
            row = [clean_number(x) for x in row]
            rows.append(row)
            
        return np.array(rows)
    
    def _get_talairach_coords(self):
        """ transforms the coordinates by talairach transform matrix from freesurfer talairach.xfm"""
        talairach_tr = self._read_talaraich_transformation()
        return np.dot(talairach_tr,self._affineM).T[:,:3]
        
        
    def __iter__(self):
        return self
        
    def next(self):
        
        if self._count>=self.npoints:
            self._count = 0
            raise StopIteration
        else:
            t = _FreesurferCoord(ras_coord = self.coordinates['ras_coords'][self._count,:],
                           voxel_coord = self.coordinates['voxel_coords'][self._count],
                           ras_tkr_coord = self.coordinates['ras_tkr_coords'][self._count,:],
                           fsvoxel_coord = self.coordinates['fsvoxel_coords'][self._count],
                           hemi = self.hemis[self._count])
                       
            self._count +=1
            return t
            


#### class _Coord and _FreesurferCoord

class _Coord(object):
    def __init__(self, ras_coord, voxel_coord):
        self.ras_coord = ras_coord
        self.voxel_coord = voxel_coord
        


class _FreesurferCoord(object):
    def __init__(self,ras_coord, voxel_coord, ras_tkr_coord, fsvoxel_coord, hemi):
        self.ras_coord = ras_coord
        self.voxel_coord = voxel_coord
        self.ras_tkr_coord = ras_tkr_coord
        self.fsvoxel_coord = fsvoxel_coord
        self.hemi = hemi

    
        
        
        





        
    
    
        
    


        
