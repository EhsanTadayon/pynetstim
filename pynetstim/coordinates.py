"""
classes to work with coordinates in native, freesurfer and MNI space
Author: Ehsan Tadayon, M.D. [sunny.tadayon@gmail.com / stadayon@bidcm.harvard.edu]
"""

import nibabel as nib
import numpy as np
import os
import nipype.pipeline.engine as pe
from nipype.interfaces.fsl import WarpPoints,  Reorient2Std
from nipype.interfaces.freesurfer import MRIConvert,Label2Label
import warnings
from mne.label import grow_labels
from nipype import Node, Workflow
from .freesurfer_files import Surf, FreesurferSurf, Annot
from nipype.interfaces.fsl import FLIRT,FNIRT
from .image_manipulation import img2img_register, mri_label2vol, img2img_coord_register
from scipy.spatial.distance import cdist
import pandas as pd
import shutil
from nibabel.freesurfer import read_label
from mne import Label
from collections import defaultdict


#############################################################################################################################
#                                                          Coords class 
#############################################################################################################################
class Coords(object):
    
    def __init__(self,coords, img_file, subject=None, coord_type='ras', working_dir=None, **traits):
        """ coordinate class
        
        Parameters
        ===============
        
        coords: numpy array of shape: npoints x 3; 
        img_file: volume file where the coordinates come from
        subject: name of subject
        coord_type: either "ras" or "voxel"
        working_dir: the path to working_dir, where the results will be written to
        **traits: other traits of the coordinates; the traits size should be equal to the number of npoints
    
        """
        
        coords = np.atleast_2d(coords)
        self.img_file = img_file
        self.subject=subject
        self.coord_type = coord_type
        self.img = nib.load(img_file)
        self.working_dir = working_dir
        self.vox2ras = self.img.affine
        self.ras2vox = np.linalg.inv(self.vox2ras)
        self.npoints = coords.shape[0]
        self.coordinates = {}
        affineM= self._to_affine_matrix(coords)
        self._count=0
        
        if coord_type=='ras':
            self.coordinates['ras_coord'] = coords
            self.coordinates['voxel_coord'] = np.round(np.dot(self.ras2vox,affineM).T[:,:3])
        
        elif coord_type=='voxel':
            self.coordinates['voxel_coord'] = coords
            self.coordinates['ras_coord'] = np.dot(self.vox2ras,affineM).T[:,:3]
            
        else:
            raise ValueError('type should be either "ras" or "voxel"')
            
            
        ### to freesurfer coords
        
        rnum1 = np.random.randint(10**15,10**16) 
        rnum2 = np.random.randint(10**10,10**11)  
        rnum = '{rnum1}_{rnum2}'.format(rnum1=rnum1,rnum2=rnum2)       
        os.makedirs(os.path.join(os.path.abspath('.'),'temp_{rnum}'.format(rnum=rnum)))
        wf_dir = os.path.join(os.path.abspath('.'),'temp_{rnum}'.format(rnum=rnum))
        
        ## creating rawavg.mgz
        mc = MRIConvert()
        mc.inputs.in_file = img_file
        mc.inputs.out_file = os.path.join(wf_dir,'rawavg.mgz')
        mc.inputs.out_type = 'mgz'
        mc.run()

        ## creating orig.mgz

        mc = MRIConvert()
        mc.inputs.in_file = os.path.join(wf_dir,'rawavg.mgz')
        mc.inputs.out_file = os.path.join(wf_dir,'orig.mgz')
        mc.inputs.out_type = 'mgz'
        mc.inputs.conform = True
        mc.run()

    
        rawavg_file = os.path.join(wf_dir,'rawavg.mgz')            
        orig_file = os.path.join(wf_dir,'orig.mgz')
        
        
        
        ### loading 
        orig_img = nib.freesurfer.load(orig_file)
        self.ras2fsvox = orig_img.header.get_ras2vox()
        self.fsvox2ras_tkr = orig_img.header.get_vox2ras_tkr()
        self.ras2ras_tkr = np.dot(self.fsvox2ras_tkr,self.ras2fsvox)
            
        ras_affineM = self._to_affine_matrix(self.coordinates['ras_coord'])
        self.coordinates['ras_tkr_coord'] = np.dot(self.fsvox2ras_tkr,np.dot(self.ras2fsvox, ras_affineM)).T[:,:3]
        self.coordinates['fsvoxel_coord'] = np.round(np.dot(self.ras2fsvox,ras_affineM).T[:,:3])
        
        shutil.rmtree(wf_dir)  
          
        ##3 adding traits    
        self.traits_list = []
        for trait in traits:
            self.add_trait(trait,traits[trait])     
                  
    def _to_affine_matrix(self,coords):
        return np.hstack((coords, np.ones((coords.shape[0],1)))).T
        
    
    def add_trait(self,trait,value):
        
        """
        adds trait; 
        
        Parameters:
        ===========
        trait: trait name for instance "color" or "name" or "opacity" etc.
        value: value of the trait; the size of value should be equal to the number of npoints
        
        """
        
        if type(value)!=np.ndarray:
            try:
                value = np.atleast_1d(value)
            except:
                ValueError('{trait} should be numpy.ndarray!'.format(trait=trait))
        
        if value.shape[0]!=self.npoints:
            raise ValueError('{trait} shape should be equal to number of points'.format(trait=trait))    
            
        self.__setattr__(trait,value)
        if trait not in self.traits_list:
            self.traits_list.append(trait)
            
    def get_traits_dict(self):
        
        """
        return a dictionary of all the traits for the coordinates
        """
        traits={}
        for trait in self.traits_list:
            traits[trait] = self.__getattribute__(trait)
        return traits
        
        
    def get_coords_df(self, coord_types='all', subset_by=None, subset_vals=None):
        
        """
        returns a dataframe for the different coordinates of the points ( ras coords, voxel_coords,.. )
        """
        if coord_types=='all':
            coord_types = self.coordinates.keys()
    
        if subset_by: 
            targets = self.subset(subset_by,subset_vals)
        else:
            targets = self

        results = []    
        for coord_type in coord_types:
            coords = targets.coordinates[coord_type]
            if hasattr(self,'name'):
                coords_df = pd.DataFrame(coords,index=targets.name,columns=[coord_type+'_{axis}'.format(axis=a) for a in ['X','Y','Z']])
            else:
                coords_df = pd.DataFrame(coords,index=['coordinate_{i}'.format(i=i) for i in np.arange(self.npoints)], columns=[coord_type+'_{axis}'.format(axis=a) for a in ['X','Y','Z']])
            
            results.append(coords_df)
            
        return pd.concat(results,axis=1)  
        
        
    def map_to_surface(self, surface):
        
        """
        maps the points to a surface ( either pial or white) or an instance of Surf class. 
        """
        
        if not isinstance(surface,Surf):
            raise ("surface should be an instance of Surf")
        
        coords_ras_tkr = self.coordinates['ras_tkr_coord']
        indices = np.argmin(cdist(surface.vertices, coords_ras_tkr), axis=0)
        mapped_coords_ras_tkr = surface.vertices[indices,:]

        mapped_coords_ras_tkr_affineM = np.hstack((mapped_coords_ras_tkr,np.ones((mapped_coords_ras_tkr.shape[0],1))))
        mapped_coords_ras = np.dot(np.linalg.inv(self.ras2ras_tkr),mapped_coords_ras_tkr_affineM.T).T
        mapped_coords_ras = mapped_coords_ras[:,0:3]
        results = {'vertices': indices, 'ras_coord': mapped_coords_ras, 'ras_tkr_coord':mapped_coords_ras_tkr}
        
        return results      
    
    
    def img2imgcoord(self, ref_img, ref_name=None, method='linear', input_reorient2std=False, ref_reorient2std=False, wf_base_dir=None, wf_name='register',
                     linear_reg_file=None, warp_field_file = None, return_as_array=False):
        
        """
        registers the coordinates to another volume. 
        
        Parameters:
        ======================                         
        
        1. ref_img: destinateion volume
                     
        2. ref_name: name of the reference subject
        
        3. method: 'linear or nonlinear'
        
        4. input_reorient2std: reorient the input image to standard orientation using fslreorient2std

        5. ref_reorient2std: reorient the reference image to standard orientation using fslreorient2std
        
        6. wf_base_dir: the base directory where the results will be stored; if not specified, it will look for working_dir of Coords class; 
           if not found, it will use the currect directory as the base directory
        
        7. wf_name: name of the working process; the results will be saved under <wf_base_dir>/<wf_name>
        
        8. linear_reg_file:  registeration file if exists
        
        9. warp_field_file: nonlinear warp field if exits
                     
        10. return_as_array: if True, returns the new coordinates in a numpy ndarray; otherwise, it will return an instance of Coords using ref_img as the image volume
        
        """
                     
        ### wf_base_dir              
        if wf_base_dir is None and self.working_dir is not None:
            wf_base_dir = self.working_dir
            
        elif wf_base_dir is None and self.working_dir is None:
            print('Working dir has not been specified, results will be stored in:  ', os.path.abspath('.'))             
            wf_base_dir = os.path.abspath('.')
            
                
        img_file = self.img_file
        ras_coords = self.coordinates['ras_coord']
        new_coords = img2img_coord_register(ras_coords, img_file, ref_img, wf_base_dir, method=method, input_reorient2std=input_reorient2std, ref_reorient2std=ref_reorient2std,
        wf_name=wf_name,linear_reg_file=linear_reg_file, warp_field_file = warp_field_file)
        
        if return_as_array is False:
            traits={}
            for trait in self.traits_list:
                traits[trait] = self.__getattribute__(trait)
            
            new_coords = Coords(coords=new_coords, img_file=ref_img, subject=ref_name,**traits)

        return new_coords
            
   
    def subset(self,by,vals):
        
        """
        subsets the coordinates
        
        Parameters:
        =============
        
        1. by (string): which trait to use for subseting the coordinates; for instance you can use "name" to subset the coordinates
        
        2. vals (list): what values to use; for instance you can provide a list of the names
        
        """
        
        idx = []
        for val in vals:
            x = np.where(self.__getattribute__(by)==val)[0][0]
            idx.append(x)
            
        idx = np.array(idx)
                
        coords =  self.coordinates['ras_coord'][idx,:]
        
        traits={}
        for trait in self.traits_list:
            traits[trait] = self.__getattribute__(trait)[idx]
        
        return Coords(coords, self.img_file , coord_type='ras', working_dir=self.working_dir, **traits)
        
        
    def __iter__(self):
        return self
    
    def __next__(self):
        
        if self._count>=self.npoints:
            self._count = 0
            raise StopIteration
            
        else:
            kwargs = {}
            for trait in self.traits_list:
                if hasattr(self,trait):
                    kwargs[trait] = self.__getattribute__(trait)[self._count]
            ras_coord = self.coordinates['ras_coord'][self._count,:]
            voxel_coord = self.coordinates['voxel_coord'][self._count,:]
            c = _Coord(ras_coord,voxel_coord,**kwargs)
            self._count+=1
            return c
            
    def __getitem__(self,idx):
        
        kwargs = {}
        for trait in self.traits_list:
            if hasattr(self,trait):
                kwargs[trait] = self.__getattribute__(trait)[self._count]
                
        ras_coord = self.coordinates['ras_coord'][idx,:]
        voxel_coord = self.coordinates['voxel_coord'][idx,:]
        c = _Coord(ras_coord,voxel_coord,**kwargs)
        return c
                  
                  
class _Coord(object):
    def __init__(self, ras_coord, voxel_coord, **kwargs):
        self.ras_coord = ras_coord
        self.voxel_coord = voxel_coord
        self.traits_list = []
        for trait,value in kwargs.items():
            self.__setattr__(trait,value)
            self.traits_list.append(trait)
        
        
class MNICoords(Coords):
    
    def __init__(self,coords, mni_template='MNI152_T1_2mm.nii.gz',mni_directory=os.environ['FSLDIR']):
        
        mni_file = os.path.join(mni_directory,mni_tmeplate)
        Coords.__init__(self,coords,mni_file)
        
            
#############################################################################################################################
#                                                          FreesurferCoords class 
#############################################################################################################################           
class FreesurferCoords(Coords):

    def __init__(self, coords, subject, freesurfer_dir, guess_hemi=True, working_dir=None, coord_type='ras', **traits):

        """
        This class implements methods to transform between coordinates in the Freesurfer space.
        
        Parameters
        ==========
        1. coords: numpy array (n x 3). Coords are RAS coords defined in the native T1 space (rawavg). 
        
        2. subject: Freesurfer subject ID
        
        3. freesurfer_dir: Freesurfer freesurfer_dir
        
        4 .guess_hemi: uses Freesurfer processed volumes to guess which hemisphere each point belongs to. 
        
        5. working_dir : the directory where the results will be written to.
        
        6. **traits: dictionary containing other traits
        
        """
        self.freesurfer_dir = freesurfer_dir
        self.subject = subject
        self.working_dir = working_dir
        coords = np.atleast_2d(coords)
        self.coord_type = coord_type
        
        
        ## setting image file names
        rawavg_file = '{freesurfer_dir}/{subject}/mri/rawavg.mgz'.format(freesurfer_dir=freesurfer_dir,subject=subject)            
        orig_file = '{freesurfer_dir}/{subject}/mri/orig.mgz'.format(freesurfer_dir=freesurfer_dir,subject=subject)
        self.img_file = rawavg_file
        
        
        ### 
        self.img = nib.load(self.img_file)
        self.orig_img = nib.freesurfer.load(orig_file)
        
        ## transformations
        self.vox2ras = self.img.affine
        self.ras2vox = np.linalg.inv(self.vox2ras)
        self.ras2fsvox = self.orig_img.header.get_ras2vox()
        self.fsvox2ras_tkr = self.orig_img.header.get_vox2ras_tkr()
        self.ras2ras_tkr = np.dot(self.fsvox2ras_tkr,self.ras2fsvox)
    
    
        ## populating coordinates
        self.npoints = coords.shape[0]
        self.coordinates = {}
        self._count=0
        
        
        affineM = self._to_affine_matrix(coords)
        
        if coord_type=='ras':
            self.coordinates['ras_coord'] = coords
            self.coordinates['voxel_coord'] = np.round(np.dot(self.ras2vox,affineM).T[:,:3])
        
        elif coord_type=='voxel':
            self.coordinates['voxel_coord'] = coords
            self.coordinates['ras_coord'] = np.dot(self.vox2ras,affineM).T[:,:3]
        
        else:
            raise ValueError('type should be either "ras" or "voxel"')
            
       
        ras_affineM = self._to_affine_matrix(self.coordinates['ras_coord'])
        self.coordinates['ras_tkr_coord'] = np.dot(self.fsvox2ras_tkr,np.dot(self.ras2fsvox, ras_affineM)).T[:,:3]
        self.coordinates['fsvoxel_coord'] = np.round(np.dot(self.ras2fsvox,ras_affineM).T[:,:3])
        self.coordinates['talairach_coord'] = self._get_talairach_coords()
       
       
       ## guessing hemisphere
        if guess_hemi:
            self._guess_hemi() 
       
       ## adding traits
        self.traits_list = []
        for trait in traits:
            self.add_trait(trait,traits[trait])
            
            
   
    def _guess_hemi(self):
        """
        uses Freesurfer voxel coordinate to guess hemisphere.
        """
        
        self.hemi = []
        self.hemi_not_determined = []
        for s in np.arange(self.npoints):
            if self.coordinates['fsvoxel_coord'][s,0]> 128: 
                self.hemi.append('lh')   
            elif self.coordinates['fsvoxel_coord'][s,0] < 128:
                self.hemi.append('rh')   
            else:
                w = """Could not determine hemiphere for point {x},{y},{z}. Right hemiphere has been chosen arbitrarily for this point.
                 Manually set the hemiphere for this point by calling set_hemi_manually!"""
                w = w.format(x=self.coordinates['ras_coord'][s,0], y=self.coordinates['ras_coord'][s,1], z=self.coordinates['ras_coord'][s,2])
                warnings.warn(w)
                
                self.hemi_not_determined.append(s)
                self.hemi.append('rh')
                 
        self.hemi = np.array(self.hemi)
        
        
    def set_hemi_manually(self, n, hemi):
        """
        sets hemisphere manually for point n
        """
        self.hemi[n] = hemi
        if n in self.hemi_not_determined:
            self.hemi_not_determined.remove(n)
        
        
    def _read_talaraich_transformation(self):
        """ read talairach transformation from freesurfer talairach.xfm output"""
        
        fname = '{freesurfer_dir}/{subject}/mri/transforms/talairach.xfm'.format(freesurfer_dir=self.freesurfer_dir,
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
        return np.dot(talairach_tr,self._to_affine_matrix(self.coordinates['ras_coord'])).T[:,:3]
            
    def map_to_annot(self, annot, map_surface='white', inplace=True):
        """ map each point to specified annotation 
        
        Parameters
        -----------
        annot: string
            which annotation to use 
        
        map_surface: string
            the surface that points are projected into to get the vertices
            
        Returns
        ---------
        structures: numpy.ndarray
            a numpy array of structures ( annotations)
        
        color: numpy.ndarray
            a numpy array (npoints x 3) that specifies the color based on the annotation provided
            
        """
        
        if len(self.hemi_not_determined)>0:
            raise ValueError('Use set_hemi_manually to assign hemiphere to these points: %s'%(','.join(self.hemi_not_determined)))
        
     
        lh_annot = Annot('lh', annot, self.subject, self.freesurfer_dir)
        rh_annot = Annot('rh', annot, self.subject, self.freesurfer_dir)
        
        colors = np.zeros((self.npoints,3))
        structures = np.empty(self.npoints,dtype='object')
            
        mapped_vertices_indices = self.map_to_surface(surface=map_surface)['vertices']
        lh_mapped_vertices_indices = mapped_vertices_indices[self.hemi=='lh']
        rh_mapped_vertices_indices = mapped_vertices_indices[self.hemi=='rh']
        
        
        if np.sum(self.hemi=='lh')>0:
            lh_colors, lh_structures = lh_annot.get_vertices_colors(lh_mapped_vertices_indices),lh_annot.get_vertices_names(lh_mapped_vertices_indices)
            colors[self.hemi=='lh',:] = lh_colors
            structures[self.hemi=='lh'] = ['lh_' + x.decode('UTF-8') for x in lh_structures]

        if np.sum(self.hemi=='rh')>0:
            rh_colors, rh_structures = rh_annot.get_vertices_colors(rh_mapped_vertices_indices),rh_annot.get_vertices_names(rh_mapped_vertices_indices)
            colors[self.hemi=='rh',:] = rh_colors
            structures[self.hemi=='rh'] = ['rh_' + x.decode('UTF-8') for x in rh_structures]
        
        if inplace:
            self.add_trait('color',colors)
            self.add_trait('name', structures)
        
        return structures, colors
        
        
    def map_to_surface(self, surface='white'):
        
        """
        maps the points to a surface ( either pial or white) or an instance of Surf class. 
        """
        
        if len(self.hemi_not_determined)>0:
            raise ValueError('Use set_hemi_manually to assign hemiphere to these points: %s'%(','.join(self.hemi_not_determined)))
        
        lh_coords_ras_tkr = self.coordinates['ras_tkr_coord'][self.hemi=='lh',:]
        rh_coords_ras_tkr = self.coordinates['ras_tkr_coord'][self.hemi=='rh',:]
        
        if surface in ['white','pial']:
            lh_surf = FreesurferSurf('lh', surface,self.subject, self.freesurfer_dir)
            rh_surf = FreesurferSurf('rh', surface, self.subject, self.freesurfer_dir)
            lh_indices = np.argmin(cdist(lh_surf.vertices, lh_coords_ras_tkr), axis=0)
            rh_indices = np.argmin(cdist(rh_surf.vertices, rh_coords_ras_tkr), axis=0)
            lh_mapped_coords_ras_tkr= lh_surf.vertices[lh_indices,:]
            rh_mapped_coords_ras_tkr= rh_surf.vertices[rh_indices,:]
            
        elif isinstance(surface,Surf):
            lh_indices = np.argmin(cdist(surface.vertices, lh_coords_ras_tkr), axis=0)
            rh_indices = np.argmin(cdist(surface.vertices, rh_coords_ras_tkr), axis=0)
            lh_mapped_coords_ras_tkr = surface.vertices[lh_indices,:]
            rh_mapped_coords_ras_tkr = surface.vertices[rh_indices,:]
        
        
        mapped_vertices = np.empty(self.npoints, dtype='int')
        mapped_vertices[self.hemi=='lh']= lh_indices
        mapped_vertices[self.hemi=='rh'] = rh_indices
        
        mapped_coords_ras_tkr = np.zeros((self.npoints,3))
        mapped_coords_ras_tkr[self.hemi=='lh',:] = lh_mapped_coords_ras_tkr
        mapped_coords_ras_tkr[self.hemi=='rh',:] = rh_mapped_coords_ras_tkr
        

        mapped_coords_ras_tkr_affineM = np.hstack((mapped_coords_ras_tkr,np.ones((mapped_coords_ras_tkr.shape[0],1))))
        mapped_coords_ras = np.dot(np.linalg.inv(self.ras2ras_tkr),mapped_coords_ras_tkr_affineM.T).T
        mapped_coords_ras = mapped_coords_ras[:,0:3]
        results = {'vertices': mapped_vertices, 'ras_coord': mapped_coords_ras, 'ras_tkr_coord':mapped_coords_ras_tkr}
        
        return results

        
    def create_surf_roi(self, extents, surface='white', map_surface='white', annot=None, wf_base_dir=None,
      wf_name='surf_roi', add_vertex_to_name=True, label2vol=False,label2vol_tidy_up=False, **kwargs):
        
        """ creates surface ROIs for each stimulation target
        
        Parameters
        ----------
        extents: float or numpy.ndarray
            specifies the raidus of the growing ROI. Either one single number for all the points or a numpy array containing 
            radius for each point
            
        surface: string
            specifies which surface to use for growing the ROIs ( white or pial)
            
        annot: string
            specifies which annotation to use to assign name and color to each ROI
        
        out_dir: string
            output directory
        
            
        Returns
        -------
        
        roi: numpy.ndarray
            resulting ROIs
        """
        ## wf_base_dir
        if wf_base_dir is None and self.working_dir is not None:
            wf_base_dir = self.working_dir
        
        elif wf_base_dir is None and self.working_dir is None:
            print('Working dir has not been specified, results will be stored in:  ', os.path.abspath('.'))             
            wf_base_dir = os.path.abspath('.')
            
        if len(self.hemi_not_determined)>0:
            raise ValueError('Use set_hemi_manually to assign hemiphere to these points: %s'%(','.join(self.hemi_not_determined)))
            

        results = self.map_to_surface(map_surface)
        mapped_vertices, mapped_coords_ras_tkr, mapped_coords_ras = results['vertices'], results['ras_coord'], results['ras_tkr_coord']
             
        ### extents can be one number or an array, make it an array if it is a number
        if type(extents)==list or type(extents)==np.ndarray:
            assert(len(extents)==self.npoints,'extents can be either one number or a list where len(extents) is equal to number of points')
        else:
            extents = [extents]*self.npoints
        
        hemi = [0 if hemi=='lh' else 1 for hemi in self.hemi]
        rois = grow_labels(self.subject, mapped_vertices, extents, hemi, self.freesurfer_dir, surface=surface)
        
       
        ### get structures and color for labels according to annotation
        if annot:
            structures, colors = self.map_to_annot(annot, map_surface=map_surface)
            for i in range(self.npoints):
                rois[i].color = colors[i]
                vertex = mapped_vertices[i]
                if add_vertex_to_name:
                    rois[i].name = structures[i]+'_{r}mm_{surf}_{vertex}'.format(r=extents[i],surf=surface,vertex=vertex)
                else:
                    rois[i].name = structures[i]+'_{r}mm_{surf}'.format(r=extents[i],surf=surface)
                
        elif hasattr(self,'name') or hasattr(self,'color'):
            for i in range(self.npoints):
                if hasattr(self,'name'):
                    vertex = mapped_vertices[i]
                    if add_vertex_to_name:
                        rois[i].name = self.name[i]+'_{r}mm_{surf}_{vertex}'.format(r=extents[i],surf=surface,vertex=vertex)
                    else:
                        rois[i].name = self.name[i]+'_{r}mm_{surf}'.format(r=extents[i],surf=surface,vertex=vertex)
                if hasattr(self,'color'):
                    rois[i].color = self.color[i]

        else:
            for i in range(self.npoints):
                vertex = mapped_vertices[i]
                if add_vertex_to_name:
                    rois[i].name = 'coor_id_{i}_{r}mm_{surf}_{vertex}'.format(r=extents[i],surf=surface,vertex=vertex,i=i)
                else:
                    rois[i].name = 'coor_id_{i}_{r}mm_{surf}_{vertex}'.format(r=extents[i],surf=surface,vertex=vertex,i=i)
            

        #### saving ROI labels

        rois_path = []
        if not os.path.exists(os.path.join(wf_base_dir,wf_name)):
            os.makedirs(os.path.join(wf_base_dir,wf_name))
            
    
        for i,roi in enumerate(rois):
            os.environ['SUBJECTS_DIR'] = self.freesurfer_dir
            
            ### saving ROI label
            roi_path = '{wf_base_dir}/{wf_name}/{roi_name}-{hemi}.label'.format(wf_base_dir=wf_base_dir, wf_name=wf_name, roi_name=roi.name, hemi=roi.hemi)
            rois_path.append(roi_path)
            roi.save(roi_path)
            
            if label2vol:
                label2vol_wf_name = '{roi_name}-{hemi}'.format(roi_name=roi.name,hemi=roi.hemi)+'_vol'
                mri_label2vol(roi,subject=self.subject, freesurfer_dir=self.freesurfer_dir,
                wf_base_dir= os.path.join(wf_base_dir,wf_name), wf_name=label2vol_wf_name, tidy_up=label2vol_tidy_up,**kwargs)
                if label2vol_tidy_up:
                    shutil.move(os.path.join(wf_base_dir,wf_name,label2vol_wf_name,'mask_with_gm','{roi_name}-{hemi}.nii.gz'.format(roi_name=roi.name,hemi=roi.hemi)),
                os.path.join(wf_base_dir,wf_name))
                else:
                    shutil.move(os.path.join(wf_base_dir,wf_name,label2vol_wf_name,'label2vol','{roi_name}-{hemi}.nii.gz'.format(roi_name=roi.name,hemi=roi.hemi)),
                    os.path.join(wf_base_dir,wf_name))
                shutil.rmtree(os.path.join(wf_base_dir,wf_name,label2vol_wf_name))
        ### converting list to arrays
        self.add_trait('roi', np.array(rois))
        return self.roi,rois_path
        
    
    def img2imgcoord(self, ref_img, ref_name=None, method='linear', input_reorient2std=True, ref_reorient2std=False, wf_base_dir = None, wf_name='register', linear_reg_file=None, warp_field_file = None, return_as_array=False):
                     
        ## wf_base_dir
        if wf_base_dir is None and self.working_dir is not None:
            wf_base_dir = self.working_dir
        
        elif wf_base_dir is None and self.working_dir is None:
            print('Working dir has not been specified, results will be stored in:  ', os.path.abspath('.'))             
            wf_base_dir = os.path.abspath('.')
                 
        ## converting rawavg to nifti
        mc = MRIConvert()
        mc.inputs.in_file = self.img_file
        mc.inputs.out_file = 'rawavg.nii.gz'
        mc.inputs.out_type = 'niigz'
        mc_node = pe.Node(mc,name='rawavg_to_nifti')
        wf = pe.Workflow(name=wf_name,base_dir=wf_base_dir)
        wf.add_nodes([mc_node])
        wf.run()
        
        ras_coords = self.coordinates['ras_coord']
        
        new_coords = img2img_coord_register(ras_coords, os.path.join(wf_base_dir,wf_name,'rawavg_to_nifti','rawavg.nii.gz'), ref_img, wf_base_dir, method=method,
                                            input_reorient2std=input_reorient2std, ref_reorient2std=ref_reorient2std,
                                            wf_name=wf_name, linear_reg_file=linear_reg_file, warp_field_file = warp_field_file)
        
        if return_as_array is False:
            ## traits
            traits={}
            for trait in self.traits_list:
                traits[trait] = self.__getattribute__(trait)
        
            new_coords = Coords(coords=new_coords, img_file=ref_img, subject=ref_name,**traits)
        
        return new_coords
    
    def img2imgcoord_by_surf(self, target_subject, wf_base_dir=None, source_surface = 'pial', source_map_surface='pial', target_surface='pial'):
        
        if wf_base_dir is None and self.working_dir is not None:
            wf_base_dir = self.working_dir
            
        elif wf_base_dir is None and self.working_dir is None:
            print('Working dir has not been specified, results will be stored in:  ', os.path.abspath('.'))             
            wf_base_dir = os.path.abspath('.')
        
        rois,rois_paths = self.create_surf_roi(extents=2, wf_base_dir= wf_base_dir, wf_name='img2imgcoord_by_surf_roi', surface=source_surface, map_surface=source_map_surface, label2vol=False)
        
        wf = pe.Workflow(name='label2label',base_dir=wf_base_dir)
        for i in range(self.npoints):
            l2l = Label2Label()
            l2l.inputs.hemisphere = self.hemi[i]
            l2l.inputs.subject_id = target_subject
            l2l.inputs.sphere_reg = os.path.join(self.freesurfer_dir, target_subject, 'surf', self.hemi[i]+'.'+'sphere.reg')
            l2l.inputs.white = os.path.join(self.freesurfer_dir, target_subject, 'surf', self.hemi[i]+'.'+'white')
            
            l2l.inputs.source_subject = self.subject
            l2l.inputs.source_label = rois_paths[i]
            l2l.inputs.source_white = os.path.join(self.freesurfer_dir, self.subject, 'surf', self.hemi[i]+'.'+'white')
            l2l.inputs.source_sphere_reg = os.path.join(self.freesurfer_dir, self.subject, 'surf', self.hemi[i]+'.'+'sphere.reg')
            l2l.subjects_dir = self.freesurfer_dir
            l2l_node = pe.Node(l2l,'label2label_{i}'.format(i=i))
            wf.add_nodes([l2l_node])
        try:
            wf.run()
        except RuntimeError:
            pass
      
        for i in range(self.npoints):
            out_label_file = os.path.join(self.freesurfer_dir, target_subject, 'label', os.path.basename(rois_paths[i]).split('.label')[0]+'_converted'+'.label')
            shutil.move(out_label_file, os.path.join(wf_base_dir, 'label2label','label2label_{i}'.format(i=i)))
         
        new_coords = np.zeros((self.npoints,3))    
        for i in range(self.npoints):
            label_file = os.path.join(wf_base_dir, 'label2label','label2label_{i}'.format(i=i),os.path.basename(rois_paths[i]).split('.label')[0]+'_converted'+'.label')
            label_vertices = read_label(label_file)
            label_vertices.sort()
            label = Label(label_vertices,hemi=self.hemi[i],subject=target_subject)
            vertex = label.center_of_mass()
            targ_surf = FreesurferSurf(hemi=label.hemi, surf=target_surface, subject = target_subject, subjects_dir=self.freesurfer_dir)
            new_coords[i,:] = targ_surf.get_coords(vertex)
        return new_coords
            
            
        
    def __iter__(self):
        return self
        
    def __next__(self):
        
        if self._count>=self.npoints:
            self._count = 0
            raise StopIteration
        else:
            kwargs = {}
            for trait in self.traits_list:
                if hasattr(self,trait):
                    kwargs[trait] = self.__getattribute__(trait)[self._count]

            t = _FreesurferCoord(ras_coord = self.coordinates['ras_coord'][self._count,:],
                           voxel_coord = self.coordinates['voxel_coord'][self._count],
                           ras_tkr_coord = self.coordinates['ras_tkr_coord'][self._count,:],
                           fsvoxel_coord = self.coordinates['fsvoxel_coord'][self._count],
                           talairach_coord = self.coordinates['talairach_coord'][self._count],
                           hemi = self.hemi[self._count],
                           **kwargs)
                       
            self._count +=1
            return t
            
    def __getitem__(self,idx):
        
        kwargs = {}
        for trait in self.traits_list:
            if hasattr(self,trait):
                kwargs[trait] = self.__getattribute__(trait)[idx]
                
        t = _FreesurferCoord(ras_coord = self.coordinates['ras_coord'][idx,:],
                        voxel_coord = self.coordinates['voxel_coord'][idx],
                        ras_tkr_coord = self.coordinates['ras_tkr_coord'][idx,:],
                        fsvoxel_coord = self.coordinates['fsvoxel_coord'][idx],
                        talairach_coord = self.coordinates['talairach_coord'][idx],
                        hemi = self.hemi[idx],
                        **kwargs)
        return t
        
    def subset(self,by,vals):
        
        idx = []
        for val in vals:
            x = np.where(self.__getattribute__(by)==val)[0][0]
            idx.append(x)
            
        idx = np.array(idx)
                
        coords =  self.coordinates['ras_coord'][idx,:]
        
        traits={}
        for trait in self.traits_list:
            traits[trait] = self.__getattribute__(trait)[idx]
        
        return FreesurferCoords(coords, self.subject, self.freesurfer_dir, guess_hemi=True, working_dir=self.working_dir, **traits)

        



class FsaverageCoords(FreesurferCoords):
    
    def __init__(self, coords, subject='fsaverage', freesurfer_dir=os.environ['SUBJECTS_DIR'], guess_hemi=True, working_dir=None, **traits):

        """
        This class implements methods to transform between coordinates in the Freesurfer space.
        
        Parameters
        ==========
        coords: numpy array (n x 3). Coords are RAS coords defined in the native T1 space (rawavg). 
        subject: Freesurfer subject ID
        freesurfer_dir: Freesurfer freesurfer_dir
        guess_hemi: uses Freesurfer processed volumes to guess which hemisphere each point belongs to. 
        **traits: dictionary containing other traits
        
        """
        self.freesurfer_dir = freesurfer_dir
        self.subject = subject
        self.working_dir = working_dir
        coords = np.atleast_2d(coords)

        
        
        ## setting image file names
            
        orig_file = '{freesurfer_dir}/{subject}/mri/orig.mgz'.format(freesurfer_dir=freesurfer_dir,subject=subject)
        
        
        ## transformations
        self.orig_img = nib.freesurfer.load(orig_file)
        self.ras2fsvox = self.orig_img.header.get_ras2vox()
        self.fsvox2ras_tkr = self.orig_img.header.get_vox2ras_tkr()
        self.ras2ras_tkr = np.dot(self.fsvox2ras_tkr,self.ras2fsvox)
        
    
        ## populating coordinates
        self.npoints = coords.shape[0]
        self.coordinates = {}
        self.coordinates['ras_coord'] = coords
        self._count=0
        
        affineM = self._to_affine_matrix(coords)
        
        ### initiating Coords class.         
        self.coordinates['ras_tkr_coord'] = np.dot(self.fsvox2ras_tkr,np.dot(self.ras2fsvox, affineM)).T[:,:3]
        self.coordinates['fsvoxel_coord'] = np.round(np.dot(self.ras2fsvox,affineM).T[:,:3])
        self.coordinates['talairach_coord'] = self.coordinates['ras_coord']
        self.coordinates['voxel_coord'] = self.coordinates['fsvoxel_coord']
        
        ## guessing hemisphere
        if guess_hemi:
            self._guess_hemi()
            
        ## adding traits
        self.traits_list = []
        for trait in traits:
            self.add_trait(trait,traits[trait])

class _FreesurferCoord(object):
    def __init__(self,ras_coord, voxel_coord, ras_tkr_coord, fsvoxel_coord, talairach_coord, hemi, **kwargs):
        
        
        self.ras_coord = ras_coord
        self.voxel_coord = voxel_coord
        self.ras_tkr_coord = ras_tkr_coord
        self.fsvoxel_coord = fsvoxel_coord
        self.talairach_coord = talairach_coord
        self.hemi = hemi
        
        self.traits_list = []
        for trait,value in kwargs.items():
            self.__setattr__(trait, value)
            self.traits_list.append(trait)
            
    def add_trait(self, trait, value):
        
        self.__setattr__(trait,value)
        if trait not in self.traits_list:
            self.traits_list.append(trait)
            
            




        
