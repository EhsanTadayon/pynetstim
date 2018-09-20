#### Surf module which will take care of freesurfer surface,annot and label files
### Author: Ehsan Tadayon, M.D. [sunny.tadayon@gmail.com / stadayon@bidmc.harvard.edu]
### Start date: Sep 16, 2018

import nibabel as nib
import numpy as np
from scipy.spatial.distance import cdist

############################## Freesurfer files


class FreesurferFile(object):
    def __init__(self, subject, fs_subjects_dir):
        self.subject = subject
        self.fs_subjects_dir = fs_subjects_dir
        

class Surf(FreesurferFile):
    
    def __init__(self,hemi, surf, subject, fs_subjects_dir):
        
        FreesurferFile.__init__(self, subject,fs_subjects_dir)
        self.surf = surf
        self.hemi = hemi
        
    def read_geometry(self):
        
        vertices, faces = nib.freesurfer.read_geometry('{fs_subjects_dir}/{subject}/surf/{hemi}.{surf}'.format(fs_subjects_dir = self.fs_subjects_dir, 
                                                                                                    subject = self.subject, surf = self.surf, hemi=self.hemi))
        return vertices, faces
        
    def get_closest_vertex(self,tkr_ras_coord):
        
        vertices, faces = self.read_geometry()
        d = []
        for i in range(vertices.shape[0]):
            d.append(np.linalg.norm(vertices[i,:] - tkr_ras_coord))
            
        vtx_id = np.argmin(d)
        vtx_coord = vertices[vtx_id,:]
        return vtx_id, vtx_coord
        
        
############################# Working with coordinates


class Coords(object):
    
    """ A class that gets RAS coordinate as input and enables converting between RAS, tkrRAS and voxel coordinates in native and freesurfer space
    
    """
    def __init__(self,ras_coords, subject, fs_subjects_dir=None):
        
        self.subject = subject
        self.fs_subjects_dir = fs_subjects_dir
        
        ras_coords = np.atleast_2d(ras_coords)
        if ras_coords.shape[1]!=3:
            raise ValueError('RAS coordinate should be 3-Dimensional.')
         
         
        ## number of points provided  
        self.npoints = ras_coords.shape[0]
        
        ## get the voxel coordinates in the native space
        self.coords = {}
        self.coords['ras_coords'] = ras_coords
        self._pointsM = np.hstack((self.coords['ras_coords'], np.ones((self.npoints, 1)))).T
        self.coords['voxels'] = self._get_native_voxels()
        
            
        ### get the freesurfer corresponding voxels and tkrRAS coordinate    
        if fs_subjects_dir:
            
            self.coords['talairachs'] = self._get_talairach_coords()
            self.coords['freesurfer_voxels'] = self._get_freesurfer_voxels()
            self.coords['tkr-ras_coords'] = self._get_freesurfer_tkr_ras()
            
            
            self.hemis = []
            for s in np.arange(self.npoints):
                
                if self.coords['freesurfer_voxels'][s,:][0] > 128: 
                    self.hemis.append('lh')
                
                elif self.coords['freesurfer_voxels'][s,:][0] < 128:
                    self.hemis.append('rh')
                
                else:
                    raise 'Could not determine hemisphere'
                    
                    
        
    def _get_native_voxels(self):
        rawObj = nib.load('{fs_subjects_dir}/{subject}/mri/rawavg.mgz'.format(fs_subjects_dir = self.fs_subjects_dir,
                                                                          subject = self.subject))
        ras2vox = rawObj.header.get_ras2vox()                                                        
        native_voxels = np.dot(ras2vox,self._pointsM).T[:,:3].astype(np.double)
        native_voxels = np.around(native_voxels)
        return native_voxels
    
        
    def _get_freesurfer_tkr_ras(self):

        origObj = nib.load('{fs_subjects_dir}/{subject}/mri/orig.mgz'.format(fs_subjects_dir = self.fs_subjects_dir,
                                                                          subject =self.subject))
        ras2vox = origObj.header.get_ras2vox()
        vox2ras_tkr = origObj.header.get_vox2ras_tkr()
        freesurfer_ras_tkr = np.dot(vox2ras_tkr,np.dot(ras2vox,self._pointsM)).T[:,:3]
        return freesurfer_ras_tkr
        
        
    
    def _get_freesurfer_voxels(self):
        
        origObj = nib.load('{fs_subjects_dir}/{subject}/mri/orig.mgz'.format(fs_subjects_dir = self.fs_subjects_dir,
                                                                         subject = self.subject))
        ras2vox = origObj.header.get_ras2vox()

        fs_voxel = np.dot(ras2vox,self._pointsM).T[:,:3]
        fs_voxel = np.round(fs_voxel)
        return fs_voxel
        
    
    def _read_talaraich_transformation(self):
        
        fname = '{fs_subjects_dir}/{subject}/mri/transforms/talairach.xfm'.format(fs_subjects_dir=self.fs_subjects_dir,subject=self.subject)
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
        
        talairach_tr = self._read_talaraich_transformation()
        return np.dot(talairach_tr,self._pointsM).T[:,:3]
        
        
    def map_to_surface(self, surf):
        
        if not self.fs_subjects_dir:
            raise ValueError('Freesurfer fs_subjects_dir has not been provided!')
        
        mapped_coords = np.zeros((self.npoints,3))
        lh_idx = np.array(self.hemis)=='lh'
        rh_idx = np.logical_not(lh_idx)
        point_coords = np.atleast_2d(self.coords['tkr-ras_coords'])
        
        if surf in ['pial','white']:
            lh_surf_coords = Surf('lh', surf, self.subject, self.fs_subjects_dir).read_geometry()[0]
            rh_surf_coords = Surf('rh',surf,self.subject, self.fs_subjects_dir).read_geometry()[0]
            
            lh_vtx_id = np.argmin(cdist(lh_surf_coords, point_coords[lh_idx,:]), axis=0)
            rh_vtx_id = np.argmin(cdist(rh_surf_coords, point_coords[rh_idx,:]), axis=0)

            mapped_coords[lh_idx,:] = lh_surf_coords[lh_vtx_id,:]
            mapped_coords[rh_idx,:] = rh_surf_coords[rh_vtx_id,:]
            
            
        elif surf=='inflated':
            
            lh_white_coords = Surf('lh', 'white', self.subject, self.fs_subjects_dir).read_geometry()[0]
            rh_white_coords = Surf('rh', 'white',self.subject, self.fs_subjects_dir).read_geometry()[0]
            
            lh_vtx_id = np.argmin(cdist(lh_white_coords, point_coords[lh_idx,:]), axis=0)
            rh_vtx_id = np.argmin(cdist(rh_white_coords, point_coords[rh_idx,:]), axis=0)
            
            lh_inflated_coords = Surf('lh', 'inflated', self.subject, self.fs_subjects_dir).read_geometry()[0]
            rh_inflated_coords = Surf('rh', 'inflated', self.subject, self.fs_subjects_dir).read_geometry()[0]
            
            mapped_coords[lh_idx,:] = lh_inflated_coords[lh_vtx_id,:]
            mapped_coords[rh_idx,:] = lh_inflated_coords[rh_vtx_id,:]
            
        else:
            
            raise ValueError('Surface should be only white, pial or inflated')
        
        return mapped_coords


            