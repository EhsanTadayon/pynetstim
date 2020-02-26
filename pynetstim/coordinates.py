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
from nipype import Node, Workflow
from .freesurfer_files import Surf, FreesurferSurf, Annot
from nipype.interfaces.fsl import FLIRT,FNIRT
from .image_manipulation import img2img_register, mri_label2vol
from scipy.spatial.distance import cdist


class Coords(object):
    
    def __init__(self,coords, img_file, subject=None, coord_type='ras', working_dir=None, **traits):
        """ coordinate class"""
        
        self.img_file = img_file
        self.subject=subject
        self.coord_type = coord_type
        self.img = nib.load(img_file)
        self.working_dir = working_dir
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
            
    def get_traits_dict(self):
        
        traits={}
        for trait in self.traits_list:
            traits[trait] = self.__getattribute__(trait)
        return traits
    
    
    def img2imgcoord(self, dest_img, method='linear', wf_base_dir = os.path.abspath('.'), wf_name='register', reg_file=None):
        
        if method not in ['linear','nonlinear']:
            raise('method should be either linear or nonlinear')
            
        np.savetxt('./temp_coords.txt',self.coordinates['ras_coord'])
        warppoints = WarpPoints()
        warppoints.inputs.in_coords = './temp_coords.txt'
        warppoints.inputs.src_file = self.img_file
        warppoints.inputs.dest_file = dest_img
        
        if reg_file is None:
            
            img2img_register(img_file = self.img_file, ref_file = dest_img, wf_base_dir = wf_base_dir, wf_name=wf_name, method=method,
        flirt_out_reg_file = 'linear_reg.mat',flirt_out_file = 'img2img_linear.nii.gz',
                    fnirt_out_file = 'img2img_nonlinear.nii.gz')
                    
            if method=='linear':
                warppoints.inputs.xfm_file =  os.path.join(wf_base_dir,wf_name, method, 'linear_reg.mat')
            elif method=='nonlinear':
                warppoints.inputs.warp_file = os.path.join(wf_base_dir,wf_name, method, '')
                   
        else:
            if method=='linear':
                warppoints.inputs.xfm_file =  reg_file
            elif method=='nonlinear':
                warppoints.inputs.warp_file = reg_file        
            
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
        
        return Coords(coords, self.img_file, coord_type='ras', **traits)
        
    def get_coords_df(self, coord_types='all', by=None, vals=None, to_df=True):
    
        if coord_types=='all':
            coord_types = self.coordinates.keys()
    
        if by: 
            targets = self.subset(by,vals)
        else:
            targets = self

        results = []    
        for coord_type in coord_types:
            coords = targets.coordinates[coord_type]
            coords_df = pd.DataFrame(coords,index=targets.name,columns=[coord_type+'_{axis}'.format(axis=a) for a in ['X','Y','Z']])
            print(coords_df.head())
            results.append(coords_df)
        return pd.concat(results,axis=1)
                
        
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

    def __init__(self, coords, subject, freesurfer_dir, guess_hemi=True, working_dir=None, **traits):

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
        
        ## setting image file names
        rawavg_file = '{freesurfer_dir}/{subject}/mri/rawavg.mgz'.format(freesurfer_dir=freesurfer_dir,subject=subject)            
        orig_file = '{freesurfer_dir}/{subject}/mri/orig.mgz'.format(freesurfer_dir=freesurfer_dir,subject=subject)
        
        ### loading 
        self.orig_img = nib.freesurfer.load(orig_file)
        self.ras2fsvox = self.orig_img.header.get_ras2vox()
        self.fsvox2ras_tkr = self.orig_img.header.get_vox2ras_tkr()
        self.ras2ras_tkr = np.dot(self.fsvox2ras_tkr,self.ras2fsvox)
        
        
        ### initiating Coords class. 
        Coords.__init__(self, coords, rawavg_file, subject=self.subject, coord_type='ras', working_dir=working_dir, **traits)
                
        self.coordinates['ras_tkr_coord'] = np.dot(self.fsvox2ras_tkr,np.dot(self.ras2fsvox, self._affineM)).T[:,:3]
        self.coordinates['fsvoxel_coord'] = np.round(np.dot(self.ras2fsvox,self._affineM).T[:,:3])
        self.coordinates['talairach_coord'] = self._get_talairach_coords()
        
        ## guessing hemisphere
        if guess_hemi:
            self._guess_hemi()
            
           
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

        
    def create_surf_roi(self, extents, surface='white', map_surface='white', annot=None, label2vol=True, out_dir=None, tidy_up=True ):
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
            

        mapped_vertices, mapped_coords_ras_tkr, mapped_coords_ras = self.map_to_surface(map_surface)
             
        ### extents can be one number or an array, make it an array if it is a number
        try:
            len(extents)
        except:
            extents = [extents]*self.npoints
        
        hemi = [0 if hemi=='lh' else 1 for hemi in self.hemi]
        rois = grow_labels(self.subject, mapped_vertices, extents, hemi, self.freesurfer_dir, surface=surface)
        
       
        ### get structures and color for labels according to annotation
        if annot:
            structures, colors = self.map_to_annot(annot, map_surface=map_surface)
            for i in range(self.npoints):
                rois[i].color = colors[i]
                vertex = mapped_vertices[i]
                rois[i].name = structures[i]+'_{r}mm_{surf}_{vertex}'.format(r=extents[i],surf=surface,vertex=vertex)
                
        elif hasattr(self,'name') or hasattr(self,'color'):
            for i in range(self.npoints):
                if hasattr(self,'name'):
                    vertex = mapped_vertices[i]
                    rois[i].name = self.name[i]+'_{r}mm_{surf}_{vertex}'.format(r=extents[i],surf=surface,vertex=vertex)
                if hasattr(self,'color'):
                    rois[i].color = self.color[i]

        else:
            for i in range(self.npoints):
                vertex = mapped_vertices[i]
                rois[i].name = 'coor_id_{i}_{r}mm_{surf}_{vertex}'.format(r=extents[i],surf=surface,vertex=vertex,i=i)
            

        #### saving ROI labels

        if out_dir:
            rois_path = []
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
        
            for i,roi in enumerate(rois):
                os.environ['SUBJECTS_DIR'] = self.freesurfer_dir
                
                ### saving ROI label
                roi_path = '{out_dir}/{roi_name}-{hemi}.label'.format(out_dir=out_dir,roi_name=roi.name,hemi=roi.hemi)
                rois_path.append(roi_path)
                roi.save(roi_path)
                
                if label2vol:
                    wf_name = '{roi_name}-{hemi}'.format(roi_name=roi.name,hemi=roi.hemi)+'_vol'
                    mri_label2vol(roi,subject=self.subject, freesurfer_dir=self.freesurfer_dir,
                    wf_base_dir=out_dir, wf_name=wf_name, tidy_up=tidy_up)
            
           
        ### converting list to arrays
        self.add_trait('roi', np.array(rois))
        return self.roi,rois_path
        
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
        
        return FreesurferCoords(coords, self.subject, self.freesurfer_dir, guess_hemi=True, **traits)

        



class FSaverage(FreesurferCoords):
    
    def __init__(self, coords, subject, freesurfer_dir, guess_hemi=True, working_dir=None, **traits):

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
        
        ## setting image file names           
        orig_file = '{freesurfer_dir}/{subject}/mri/orig.mgz'.format(freesurfer_dir=freesurfer_dir,subject=self.subject)
        
        ### loading 
        self.orig_img = nib.freesurfer.load(orig_file)
        self.ras2fsvox = self.orig_img.header.get_ras2vox()
        self.fsvox2ras_tkr = self.orig_img.header.get_vox2ras_tkr()
        self.ras2ras_tkr = np.dot(self.fsvox2ras_tkr,self.ras2fsvox)
        
        
        ### initiating Coords class. 
        Coords.__init__(self, coords, orig_file, subject=self.subject, coord_type='ras', working_dir=working_dir, **traits)
                
        self.coordinates['ras_tkr_coord'] = np.dot(self.fsvox2ras_tkr,np.dot(self.ras2fsvox, self._affineM)).T[:,:3]
        self.coordinates['fsvoxel_coord'] = np.round(np.dot(self.ras2fsvox,self._affineM).T[:,:3])
        self.coordinates['talairach_coord'] = self.coordinates['ras_coord']
        
        ## guessing hemisphere
        if guess_hemi:
            self._guess_hemi()

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


        