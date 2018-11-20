"""
classes to work with coordinates in native, freesurfer and MNI space
Author: Ehsan Tadayon, M.D. [sunny.tadayon@gmail.com / stadayon@bidcm.harvard.edu]
"""

import nibabel as nib
import numpy as np
import os
from nipype.interfaces.fsl import WarpPoints
import warnings
from mne.label import grow_labels
from nipype.interfaces.freesurfer import Label2Vol,Binarize,MRIsCalc
from nipype import Node, Workflow
from surface import Surf, FreesurferSurf
from annotation import Annot



class Coords(object):
    
    def __init__(self,coords, img_file, coord_type='ras', **traits):
        """ coordinate class"""
        
        self.img_file = img_file
        self.coord_type = coord_type
        self.img = nib.load(img_file)
        self.vox2ras = self.img.affine
        self.ras2vox = np.linalg.inv(self.vox2ras)
        self.npoints = coords.shape[0]
        self.coordinates = {}
        self._affineM = np.hstack((coords, np.ones((coords.shape[0],1)))).T
        self._count=0
        
        
        if coord_type=='ras':
            self.coordinates['ras_coord'] = coords
            self.coordinates['voxel_coord'] = np.round(np.dot(self.ras2vox,self._affineM).T[:,:3])
        elif coord_type=='voxel':
            self.coordinates['voxel_coord'] = coords
            self.coordinates['ras_coord'] = np.dot(self.vox2ras,self._affineM).T[:,:3]
            
        else:
            raise ValueError('type should be either "ras" or "voxel"')
            
        
        
        self.traits_list = []
        for trait in traits:
            self.add_trait(trait,traits[trait])     
            
    def add_trait(self,trait,value):
        
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
     
    
    def img2imgcoord(self, dest_img, reg_file, type):
        
       
        np.savetxt('./temp_coords.txt',self.coordinates['ras_coord'])
        warppoints = WarpPoints()
        warppoints.inputs.in_coords = './temp_coords.txt'
        warppoints.inputs.src_file = self.img_file
        warppoints.inputs.dest_file = dest_img
        
        if type=='xfm':
            warppoints.inputs.xfm_file = reg_file
        elif type=='warp':
            warppoints.inputs.warp_file = reg_file
        else:
            raise ValueError('type should be either xfm or warp')
            
        warppoints.inputs.coord_mm = True
        res = warppoints.run()
        res = np.loadtxt('./temp_coords_warped.txt')
        res = res[0:res.shape[0]-1,:]
        
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
    def __init__(self, ras_coord, voxel_coord,**kwargs):
        self.ras_coord = ras_coord
        self.voxel_coord = voxel_coord
        self.traits_list = []
        for trait,value in kwargs.items():
            self.__setattr__(trait,value)
            self.traits_list.append(trait)
        
    
    
class MNICoords(Coords):
    
    def __init__(self,coords,mni_template='MNI152_T1_2mm.nii.gz',mni_directory=os.environ['FSLDIR']):
        
        mni_file = os.path.join(mni_directory,mni_tmeplate)
        Coords.__init__(self,coords,mni_file)
        
            
           
class FreesurferCoords(Coords):

    def __init__(self, coords, subject, subjects_dir, guess_hemi=True, **traits):

        
        self.subjects_dir = subjects_dir
        self.subject = subject
        rawavg_file = '{subjects_dir}/{subject}/mri/rawavg.mgz'.format(subjects_dir=subjects_dir,subject=subject)
        orig_file = '{subjects_dir}/{subject}/mri/orig.mgz'.format(subjects_dir=subjects_dir,subject=subject)
        self.orig_img = nib.freesurfer.load(orig_file)
        self.ras2fsvox = self.orig_img.header.get_ras2vox()
        self.fsvox2ras_tkr = self.orig_img.header.get_vox2ras_tkr()
        
        self.ras2ras_tkr = np.dot(self.fsvox2ras_tkr,self.ras2fsvox)
        
        Coords.__init__(self, coords, rawavg_file, coord_type='ras', **traits)        
        self.coordinates['ras_tkr_coord'] = np.dot(self.fsvox2ras_tkr,np.dot(self.ras2fsvox, self._affineM)).T[:,:3]
        self.coordinates['fsvoxel_coord'] = np.round(np.dot(self.ras2fsvox,self._affineM).T[:,:3])
        self.coordinates['talairach_coord'] = self._get_talairach_coords()
        
        if guess_hemi:
            self._guess_hemi()
            
           
    def _guess_hemi(self):
        self.hemi = []
        self.hemi_not_determined = []
        for s in np.arange(self.npoints):
            if self.coordinates['fsvoxel_coord'][s,0]> 128: 
                self.hemi.append('lh')   
            elif self.coordinates['fsvoxel_coord'][s,0] < 128:
                self.hemi.append('rh')   
            else:
                warnings.warn('Could not determine hemiphere for point {x},{y},{z}. Right hemiphere has been chosen arbitrarily for this point. Manually set the hemiphere for this point by calling set_hemi_manually!'.format(x=self.coordinates['ras_coord'][s,0], y=self.coordinates['ras_coord'][s,1], z=self.coordinates['ras_coord'][s,2]))
                
                self._hemis_not_determined.append(s)
                self.hemi.append('rh')
                 
        self.hemi = np.array(self.hemi)
        
        
    def set_hemi_manually(self, n, hemi):
        
        self.hemi[n] = hemi
        if n in self._hemis_not_determined:
            self.hemis_not_determined.remove(n)
        
        
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
        
     
        lh_annot = Annot('lh', annot, self.subject, self.subjects_dir)
        rh_annot = Annot('rh', annot, self.subject, self.subjects_dir)
        
        colors = np.zeros((self.npoints,3))
        structures = np.empty(self.npoints,dtype='object')
        
        
        lh_coords = self.coordinates['ras_tkr_coord'][self.hemi=='lh',:]
        rh_coords = self.coordinates['ras_tkr_coord'][self.hemi=='rh',:]

        if np.sum(self.hemi=='lh')>0:
            colors[self.hemi=='lh',:] = lh_colors
            structures[self.hemi=='lh'] = ['lh_' + x.decode('UTF-8') for x in lh_structures]

        if np.sum(self.hemi=='rh')>0:
            colors[self.hemi=='rh',:] = rh_colors
            structures[self.hemi=='rh'] = ['rh_' + x.decode('UTF-8') for x in rh_structures]
        
    
        colors[self.hemi=='lh',:] = lh_colors
        colors[self.hemi=='rh',:] = rh_colors
    
        
        if inplace:
            self.add_trait('color',colors)
            self.add_trait('name', structures)
        
        return structures, colors
        
        
    def map_to_surface(self, surface='white'):
        
        if len(self.hemi_not_determined)>0:
            raise ValueError('Use set_hemi_manually to assign hemiphere to these points: %s'%(','.join(self.hemi_not_determined)))
        
        lh_coords = self.coordinates['ras_tkr_coord'][self.hemi=='lh',:]
        rh_coords = self.coordinates['ras_tkr_coord'][self.hemi=='rh',:]
        
        if surface in ['white','pial']:
            lh_surf = FreesurferSurf('lh', surface,self.subject, self.subjects_dir)
            rh_surf = FreesurferSurf('rh', surface, self.subject, self.subjects_dir)
            lh_mapped_vertices,lh_mapped_coords = lh_surf.project_coords(lh_coords)
            rh_mapped_vertices, rh_mapped_coords = rh_surf.project_coords(rh_coords)
            
        elif isinstance(surface,Surf):
            lh_mapped_vertices, lh_mapped_coords = surface.project_coords(lh_coords)
            rh_mapped_vertices, rh_mapped_coords = surface.project_coords(rh_coords)
        
        
        mapped_vertices = np.empty(self.npoints, dtype='int')
        mapped_vertices[self.hemi=='lh']= lh_mapped_vertices
        mapped_vertices[self.hemi=='rh'] = rh_mapped_vertices
        
        mapped_coords = np.zeros((self.npoints,3))
        mapped_coords[self.hemi=='lh',:] = lh_mapped_coords
        mapped_coords[self.hemi=='rh',:] = rh_mapped_coords
        
        return mapped_vertices, mapped_coords

        
    def create_surf_roi(self, extents, surface='white', annot=None, out_dir=None):
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
        
        if len(self.hemi_not_determined)>0:
            raise ValueError('Use set_hemi_manually to assign hemiphere to these points: %s'%(','.join(self.hemi_not_determined)))
        
        mapped_vertices, mapped_coords = self.map_to_surface(surface)
        
        hemi = [0 if hemi=='lh' else 1 for hemi in self.hemi]
        rois = grow_labels(self.subject, mapped_vertices, extents, hemi, self.subjects_dir)
        
        ### extents can be one number or an array, make it an array if it is a number
        try:
            len(extents)
        except:
            extents = [extents]*self.npoints
        
        ### get structures and color for labels according to annotation
        if annot:
            structures, colors = self.map_to_annot(annot, map_surface=surface)
            
            for i in range(self.npoints):
                rois[i].color = colors[i]
                vertex = mapped_vertices[i]
                rois[i].name = structures[i]+'_{r}mm_{surf}_{vertex}'.format(r=extents[i],surf=surface,vertex=vertex)
        
        #### saving ROI labels
        if out_dir:
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
        
            for i,roi in enumerate(rois):
                os.environ['SUBJECTS_DIR'] = self.subjects_dir
                
                ### saving ROI label
                roi_path = '{out_dir}/{roi_name}.label'.format(out_dir=out_dir,roi_name=roi.name)
                roi.save(roi_path)
                
                ### save volume
                wf = Workflow(name=roi.name+'_vol', base_dir=out_dir)
                
                label2vol = Node(Label2Vol(label_file=roi_path, template_file='{subjects_dir}/{subject}/mri/T1.mgz'.format(subjects_dir=self.subjects_dir, subject=self.subject),
                        hemi=roi.hemi, proj=(u'frac',0,1,0.01), identity=True, subject_id=self.subject), name='label2vol')
                        
                mask_dilate = Node(Binarize(dilate=1,erode=1,min=1),name='dilate_label_vol')
                mris_calc = Node(MRIsCalc(),name='mask_with_gm')
                mris_calc.inputs.in_file2='{subjects_dir}/{subject}/mri/{hemi}.ribbon.mgz'.format(subjects_dir=self.subjects_dir, subject=self.subject, hemi=roi.hemi)
                mris_calc.inputs.action='mul'
                mris_calc.inputs.out_file=roi.name+'.nii.gz'
                wf.connect([
                            (label2vol,mask_dilate,[("vol_label_file","in_file")]),
                            (mask_dilate,mris_calc,[('binary_file','in_file1')]),
                            ])
                    
                wf.run()
           
        ### converting list to arrays
        self.add_trait('roi', np.array(rois))
        return self.roi    
        
    def __iter__(self):
        return self
        
    def next(self):
        
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
                           hemi = self.hemi[self._count],
                           **kwargs)
                       
            self._count +=1
            return t
            
    def __getitem__(self,idx):
        
        kwargs = {}
        for trait in self.traits_list:
            if hasattr(self,trait):
                kwargs[trait] = self.__getattribute__(trait)[self._count]
                
        t = _FreesurferCoord(ras_coord = self.coordinates['ras_coord'][idx,:],
                        voxel_coord = self.coordinates['voxel_coord'][idx],
                        ras_tkr_coord = self.coordinates['ras_tkr_coord'][idx,:],
                        fsvoxel_coord = self.coordinates['fsvoxel_coord'][idx],
                        hemi = self.hemi[idx],
                        **kwargs)
        return t

        
class _FreesurferCoord(object):
    def __init__(self,ras_coord, voxel_coord, ras_tkr_coord, fsvoxel_coord, hemi, **kwargs):
        self.ras_coord = ras_coord
        self.voxel_coord = voxel_coord
        self.ras_tkr_coord = ras_tkr_coord
        self.fsvoxel_coord = fsvoxel_coord
        self.hemi = hemi
        
        self.traits_list = []
        for trait,value in kwargs.items():
            self.__setattr__(trait, value)
            self.traits_list.append(trait)

        
        
        
        
        

        