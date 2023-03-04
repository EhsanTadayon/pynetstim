"""
classes to work with coordinates in native, freesurfer and MNI space
"""

import nibabel as nib
import numpy as np
import os
import nipype.pipeline.engine as pe
from nipype import Node, Workflow
from nipype.interfaces.fsl import WarpPoints,  Reorient2Std, FLIRT,FNIRT
from nipype.interfaces.fsl.maths import MathsCommand
from nipype.interfaces.freesurfer import MRIConvert,Label2Label
import warnings
from mne.label import grow_labels
from .freesurfer_files import Surf, FreesurferSurf, Annot
from .image_manipulation import  mri_label2vol, Registeration
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

    def __init__(self,coords, img_file, subject, working_dir, coord_type='ras', **traits):
        """ General class to work with coordinates

        Parameters
        -----------

        coords: numpy array
                x,y,z coordinates matrix (npoints x 3)

        img_file: str
                path to image file


        subject: str
                subject name

        coord_type: str, {'ras', 'voxel'}
                coordinate system

        working_dir:
                the path to working directory


        **traits: other traits of the coordinates;
                 the traits size should be equal to the number of npoints
                 for instance, one can add "name" or "color" as extra traits for each coordinate


        Returns
        --------

        An object of Coords class

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
        self._count=0  ## for __next__ to iterate


        ##adding traits
        self.traits_list = []
        self.traits_type = {}
        for trait in traits:
            self.add_trait(trait,traits[trait],trait_type='transferable')

        self.add_trait('index',np.arange(1,self.npoints+1),'nontransferable')


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
        os.makedirs(os.path.join(self.working_dir,'temp_{rnum}'.format(rnum=rnum)))
        wf_dir = os.path.join(self.working_dir,'temp_{rnum}'.format(rnum=rnum))

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


    def _to_affine_matrix(self,coords):

        """
        returns an affine matrix for the specified coordinates
        """
        return np.hstack((coords, np.ones((coords.shape[0],1)))).T


    def add_trait(self,trait,value,trait_type='transferable'):

        """
        Adds trait to the object

        Parameters:
        ------------

        trait: str
               name of trait such as "name" or "color" or "opacity" or "network" etc.

        value: list or array
               a list or array of values for the trait


        Examples:
        -----------
        coords = Coords(np.array([23,42,12],[12,45,12]),
                        img_file = anat_img)
        coords.add_trait('name',['ldlpfc','lipl'])
        coords.add_trait('color',[(0,1,0),(1,0,0)])

        print('coordinates names: ', coords.name)
        print('coordinates colors: ', coords.color)

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
            self.traits_type[trait] = trait_type


    def get_traits_dict(self,trait_type='all'):

        """
        return a dictionary of all the traits for the coordinates
        """

        traits={}
        for trait in self.traits_list:
            if trait_type=='all':
                traits[trait] = self.__getattribute__(trait)
            else:
                if self.traits_type[trait]==trait_type:
                    traits[trait] = self.__getattribute__(trait)
        return traits


    def get_coords_df(self, coord_types=['ras_coord','voxel_coord'], subset_by=None, subset_vals=None, adding_traits_to_columns=True):

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
                name = ['coordinate_{i}'.format(i=i) for i in np.arange(self.npoints)]
                coords_df = pd.DataFrame(coords,index=name, columns=[coord_type+'_{axis}'.format(axis=a) for a in ['X','Y','Z']])

            results.append(coords_df)


        df = pd.concat(results,axis=1)
        df.index_name = 'name'

        ## add traits
        if adding_traits_to_columns:

            traits = self.get_traits_dict()
            if 'name' in traits.keys():
                traits.pop('name')
            for trait in traits:
                df[trait] = self.__getattribute__(trait).tolist()

        return df


    def map_to_surface(self, surface):

        """
        maps the points to a surface

        Parameters:
        -----------

        surface: either {'pial','white'} or an instance of Surf

        Returns:
        ----------
        A dictionary with the following keys:

            vertices: the vertices numbers

            ras_coord: the mapped ras coordinates

            ras_tkr_coord: the mapped ras tkr coordinates

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


    def img2imgcoord(self, ref_img, ref_name, input_brain=None, ref_brain=None, method='linear', input_reorient2std=False, ref_reorient2std=False, wf_name='register',
                     linear_reg_file=None, warp_field_file = None, return_as_array=False):

        """
        registers the coordinates to another volume.

        Parameters:
        -----------

        ref_img:
               path to reference image

        ref_name: str
                reference subject name

        method: str, {'linear', 'nonlinear'}
                     FSL registration method (FLIRT for 'linear' or FNIRT for 'nonlinear')

        input_reorient2std: boolean
                     reorient moving volume (img_file) to standard space using fslreorient2std

        ref_reorient2std: boolean
                     reorient the reference image to standard orientation using fslreorient2std

        wf_name: str
                      name of the working process; the results will be saved under <wf_base_dir>/<wf_name>

        linear_reg_file: str
                      registeration file if exists

        warp_field_file: str
                      nonlinear warp field if exits

        return_as_array: boolean
                      if True, returns the new coordinates in a numpy ndarray;
                      otherwise, it will return an instance of Coords using ref_img as the image volume

        """




        input_file = self.img_file


        ras_coords = self.coordinates['ras_coord']

        reg = Registeration(self.img_file,ref_img,input_brain=input_brain,ref_brain=ref_brain)
        reg.change_config('input_reorient2std',input_reorient2std)
        reg.change_config('ref_reorient2std',ref_reorient2std)
        reg.run(method,wf_name,self.working_dir)
        new_coords = reg.warp_points(coords = ras_coords, method=method)


        if return_as_array is False:

            new_working_dir = os.path.join(self.working_dir,'coords_registered_to_{ref_name}'.format(ref_name=ref_name))
            new_coords = Coords(coords=new_coords, img_file=ref_img, subject=ref_name, working_dir=new_working_dir, **self.get_traits_dict(trait_type='transferable'))

        return new_coords


    def subset(self,by,vals):

        """
        subsets the coordinates

        Parameters:
        -----------

        by: str
             which trait to use for subseting the coordinates; for instance you can use "name" to subset the coordinates

        vals: list
             what values to use; for instance you can provide a list of the names

        Returns:
        ---------
        An object of class Coords with all the traits from the original coords instance


        Examples:
        ---------
        coords = Coords(np.array([1,2,3],[3,4,5],[6,7,8]),img_file='anat.nii.gz',name=['c1','c2','c3'],network=['network1','network2','network1'])
        coords_sub = coords.subset('network',['network1'])

        """

        idx = []
        for val in vals:
            x = np.where(self.__getattribute__(by)==val)[0].tolist()
            idx.extend(x)

        idx = np.array(idx)

        coords =  self.coordinates['ras_coord'][idx,:]

        traits={}
        for trait in self.traits_list:
            traits[trait] = self.__getattribute__(trait)[idx]

        return Coords(coords, self.img_file , coord_type='ras', working_dir=self.working_dir, **traits)


    def create_vol_roi(self, kernel_args = 'kernel sphere 5', wf_name='volume_roi'):

        if not os.path.exists(os.path.join(self.working_dir,'volume_roi')):
            os.makedirs(os.path.join(self.working_dir,'volume_roi'))

        ## roi names
        if hasattr(self,'name'):
            names = self.name
        else:
            names = ['coord_{i}'.format(i=i+1) for i in range(self.npoints)]


        rois_paths = []

        for i in range(self.npoints):

            vol = np.zeros(self.img.shape)
            indices = self.coordinates['voxel_coord'][i,:].astype(int)
            vol[indices[0],indices[1],indices[2]] = 1
            volImg = nib.Nifti1Image(vol,self.img.affine,self.img.header)
            output_name = '{name}-kernel_{kernel}-volumeroi.nii.gz'.format(name=names[i],kernel=kernel_args.replace(' ','_'))
            output = os.path.join(self.working_dir,wf_name,output_name)
            nib.save(volImg,output)
            
            #kernel
            kernel = MathsCommand()
            kernel.inputs.in_file = output
            kernel.inputs.out_file = output
            kernel.inputs.args = '-{kernel_args} -fmean'.format(kernel_args=kernel_args)
            print(kernel.cmdline)
            kernel.run()
            
            #binarize
            binarize = MathsCommand()
            binarize.inputs.in_file = output
            binarize.inputs.out_file = output
            binarize.inputs.args = '-bin'
            print(binarize.cmdline)
            binarize.run()
                        
            #append to the list
            rois_paths.append(output_name)

        return rois_paths

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





#############################################################################################################################
#                                                          MNICoords class
#############################################################################################################################

class MNICoords(Coords):

    def __init__(self,coords, mni_template='MNI152_T1_2mm.nii.gz',mni_directory=os.environ['FSLDIR']):

        mni_file = os.path.join(mni_directory,mni_template)
        Coords.__init__(self,coords,mni_file)


#############################################################################################################################
#                                                          FreesurferCoords class
#############################################################################################################################
class FreesurferCoords(Coords):

    fsaverage_files = ['fsaverage','fsaverage4','fsaverage6','fsaverage5']

    def __init__(self, coords, subject, freesurfer_dir, working_dir, guess_hemi=True, coord_type='ras', **traits):

        """
        Coords class when the freesurfer recon-all exists.

        Parameters:
        ----------

        coords: numpy array
                x,y,z coordinates (npoints x 3)

        subject: str
                subject name

        freesurfer_dir: str
                Freesurfer SUBJECTS_DIR

        guess_hemi: boolean
                for each coordinate, guesses hemisphere

        working_dir : str
                 the directory where the results will be written to.

        **traits:
                other traits for each coordinate such as "color" or "name"; look at the Coords class

        """

        if subject in FreesurferCoords.fsaverage_files:

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

            ## adding traits
            self.traits_list = []
            self.traits_type = {}
            for trait in traits:
                self.add_trait(trait,traits[trait],trait_type='transferable')

            self.add_trait('index',np.arange(1,self.npoints+1),'nontransferable')

            ### initiating Coords class.
            self.coordinates['ras_tkr_coord'] = np.dot(self.fsvox2ras_tkr,np.dot(self.ras2fsvox, affineM)).T[:,:3]
            self.coordinates['fsvoxel_coord'] = np.round(np.dot(self.ras2fsvox,affineM).T[:,:3])
            self.coordinates['talairach_coord'] = self.coordinates['ras_coord']
            self.coordinates['voxel_coord'] = self.coordinates['fsvoxel_coord']

            ## guessing hemisphere
            if guess_hemi:
                self._guess_hemi()

        else:
            self.freesurfer_dir = freesurfer_dir
            self.subject = subject
            self.working_dir = working_dir
            coords = np.atleast_2d(coords)
            self.coord_type = coord_type
            self.npoints = coords.shape[0]

            ## adding traits
            self.traits_list = []
            self.traits_type = {}
            for trait in traits:
                self.add_trait(trait,traits[trait],trait_type='transferable')
            self.add_trait('index',np.arange(1,self.npoints+1),'nontransferable')



            ## setting image file names
            rawavg_file = '{freesurfer_dir}/{subject}/mri/rawavg.mgz'.format(freesurfer_dir=freesurfer_dir,subject=subject)
            orig_file = '{freesurfer_dir}/{subject}/mri/orig.mgz'.format(freesurfer_dir=freesurfer_dir,subject=subject)
            self.img_file = rawavg_file
            self.brain_file  = os.path.join(freesurfer_dir,subject,'mri','brain.mgz')



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

        #self.hemi = np.array(self.hemi)
        self.add_trait('hemi',np.array(self.hemi))

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

        Parameters:
        -----------
        annot: str
            which annotation to use

        map_surface: str, {'pial','white'}
            the surface that points are projected into to get the vertices

        Returns:
        ---------
        structures: numpy array
            a numpy array of structures ( annotations)

        color: numpy array
            a numpy array (npoints x 3) that specifies the color based on the annotation provided

        """

        if len(self.hemi_not_determined)>0:
             raise ValueError('Use set_hemi_manually to assign hemiphere to these points: %s'%(','.join([str(i) for i in self.hemi_not_determined])))

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
            raise ValueError('Use set_hemi_manually to assign hemiphere to these points: %s'%(','.join([str(i) for i in self.hemi_not_determined])))

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


    def create_surf_roi(self, extents, surface='white', map_surface='white', map_to_annot=None,
      wf_name='surf_roi', add_vertex_to_name=True, colors = None, add_as_trait=True):

        """ creates surface ROIs for each coordinate point

        Parameters
        ----------
        extents: int or list or np.array
            specifies the raidus of the growing ROI. Either one single number for all the points or a list or numpy array containing
            radius for each point

        surface: str
            specifies which surface to use for growing the ROIs ( white or pial)

        map_to_annot: str
            specifies which annotation to use to assign name and color to use for ROI labeling

        wf_name: str
              workflow name , the results will be saved under <wf_base_dir>/<wf_name>

        add_vertex_to_name: boolean
              if True, adds the projected vertex number to the ROI name ( in case the ROIs might have similar naming, adding vertex number can help get rid of this issue)

        Returns
        -------

        rois: instances of class of Label

        rois_paths: path to the saved labels

        """

        if len(self.hemi_not_determined)>0:
             raise ValueError('Use set_hemi_manually to assign hemiphere to these points: %s'%(','.join([str(i) for i in self.hemi_not_determined])))


        results = self.map_to_surface(map_surface)
        mapped_vertices, mapped_coords_ras_tkr, mapped_coords_ras = results['vertices'], results['ras_coord'], results['ras_tkr_coord']

        ### extents can be one number or an array, make it an array if it is a number
        if type(extents)==list or type(extents)==np.ndarray:
            assert(len(extents)==self.npoints,'extents can be either one number or a list where len(extents) is equal to number of points')
            extents = map(int,extents)
        else:
            extents = [int(extents)]*self.npoints

        hemi = [0 if hemi=='lh' else 1 for hemi in self.hemi]
        rois = grow_labels(self.subject, mapped_vertices, extents, hemi, self.freesurfer_dir, surface=surface, colors= colors)


        ### get structures and color for labels according to annotation
        if map_to_annot:
            structures, colors = self.map_to_annot(map_to_annot, map_surface=map_surface)
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
                        rois[i].name = self.hemi[i]+ '.' + self.name[i]+'_{r}mm_{surf}_{vertex}'.format(r=extents[i],surf=surface,vertex=vertex)
                    else:
                        rois[i].name = self.hemi[i]+ '.'+ self.name[i]+'_{r}mm_{surf}'.format(r=extents[i],surf=surface,vertex=vertex)
                if hasattr(self,'color'):
                    rois[i].color = self.color[i]

        else:
            for i in range(self.npoints):
                vertex = mapped_vertices[i]
                if add_vertex_to_name:
                    rois[i].name = self.hemi[i]+ '.' + 'coor_id_{i}_{r}mm_{surf}_{vertex}'.format(r=extents[i],surf=surface,vertex=vertex,i=i)
                else:
                    rois[i].name = self.hemi[i]+ '.' + 'coor_id_{i}_{r}mm_{surf}_{vertex}'.format(r=extents[i],surf=surface,vertex=vertex,i=i)


        #### saving ROI labels

        rois_path = []
        if not os.path.exists(os.path.join(self.working_dir,wf_name)):
            os.makedirs(os.path.join(self.working_dir,wf_name))


        for i,roi in enumerate(rois):
            os.environ['SUBJECTS_DIR'] = self.freesurfer_dir

            ### saving ROI label
            roi_path = '{working_dir}/{wf_name}/{roi_name}.label'.format(working_dir=self.working_dir, wf_name=wf_name, roi_name=roi.name)
            rois_path.append(roi_path)
            roi.save(roi_path)

        ### converting list to arrays
        if add_as_trait:
            self.add_trait('roi', np.array(rois),'nontransferable')

        return np.array(rois),rois_path


    def img2imgcoord(self, ref_img, ref_name, ref_brain=None, method='linear', input_reorient2std=True, ref_reorient2std=False, wf_name='register', linear_reg_file=None, warp_field_file = None, return_as_array=False):

        ## converting rawavg to nifti
        mc = MRIConvert()
        mc.inputs.in_file = self.img_file
        mc.inputs.out_file = 'rawavg.nii.gz'
        mc.inputs.out_type = 'niigz'
        mc_node = pe.Node(mc,name='rawavg_to_nifti')
        wf = pe.Workflow(name=wf_name,base_dir=self.working_dir)
        wf.add_nodes([mc_node])
        wf.run()
        input_file = os.path.join(self.working_dir,wf_name,'rawavg_to_nifti','rawavg.nii.gz')

        ##converting brain.mgz to Nifti1Image
        mc = MRIConvert()
        mc.inputs.in_file = self.brain_file
        mc.inputs.out_file = 'brain.nii.gz'
        mc.inputs.out_type = 'niigz'
        mc_node = pe.Node(mc,name='brain_to_nifti')
        wf = pe.Workflow(name=wf_name,base_dir=self.working_dir)
        wf.add_nodes([mc_node])
        wf.run()
        input_brain = os.path.join(self.working_dir,wf_name,'brain_to_nifti','brain.nii.gz')

        ras_coords = self.coordinates['ras_coord']



        reg = Registeration(input_file, ref_img, input_brain,ref_brain)
        reg.change_config('input_reorient2std',input_reorient2std)
        reg.change_config('ref_reorient2std',ref_reorient2std)
        reg.run(method, wf_name, self.working_dir)


        new_coords = reg.warp_points(ras_coords, method=method)

        if return_as_array is False:
            new_working_dir = os.path.join(self.working_dir,'coords_registered_to_{ref_name}'.format(ref_name=ref_name))
            new_coords = Coords(coords=new_coords, img_file=ref_img, subject=ref_name, working_dir=new_working_dir, **self.get_traits_dict('transferable'))

        return new_coords

    def img2imgcoord_by_surf(self, target_subject, source_surface = 'pial', source_map_surface='pial', target_surface='pial'):


        rois,rois_paths = self.create_surf_roi(extents=2, wf_name='img2imgcoord_by_surf_roi', surface=source_surface, map_surface=source_map_surface,add_as_trait=False)

        wf = pe.Workflow(name='label2label',base_dir=os.path.join(self.working_dir,'img2imgcoord_by_surf_roi'))

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
            shutil.move(out_label_file, os.path.join(self.working_dir,'img2imgcoord_by_surf_roi', 'label2label','label2label_{i}'.format(i=i)))

        new_coords = np.zeros((self.npoints,3))
        for i in range(self.npoints):
            label_file = os.path.join(self.working_dir,'img2imgcoord_by_surf_roi', 'label2label','label2label_{i}'.format(i=i),os.path.basename(rois_paths[i]).split('.label')[0]+'_converted'+'.label')
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
                        **kwargs)
        return t

    def subset(self,by,vals):

        idx = []
        for val in vals:
            x = np.where(self.__getattribute__(by)==val)[0].tolist()
            idx.extend(x)

        idx = np.array(idx)

        coords =  self.coordinates['ras_coord'][idx,:]

        traits={}
        for trait in self.traits_list:
            traits[trait] = self.__getattribute__(trait)[idx]


        return FreesurferCoords(coords, self.subject, self.freesurfer_dir, guess_hemi=True, working_dir=self.working_dir, **traits)





class _FreesurferCoord(object):
    def __init__(self,ras_coord, voxel_coord, ras_tkr_coord, fsvoxel_coord, talairach_coord, **kwargs):


        self.ras_coord = ras_coord
        self.voxel_coord = voxel_coord
        self.ras_tkr_coord = ras_tkr_coord
        self.fsvoxel_coord = fsvoxel_coord
        self.talairach_coord = talairach_coord


        self.traits_list = []
        for trait,value in kwargs.items():
            self.__setattr__(trait, value)
            self.traits_list.append(trait)

    def add_trait(self, trait, value):

        self.__setattr__(trait,value)
        if trait not in self.traits_list:
            self.traits_list.append(trait)
