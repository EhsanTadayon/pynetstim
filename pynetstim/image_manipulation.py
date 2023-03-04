from nipype.interfaces import fsl
import nipype.pipeline.engine as pe
from nipype.interfaces.freesurfer import Label2Vol,Binarize,MRIsCalc,MRIConvert
import os
import numpy as np




class Registeration(object):

    def __init__(self, input_file, ref_file,
                input_brain=None, ref_brain=None,
                linear_reg_file=None, nonlinear_reg_file=None):

        self.input_file_orig = input_file
        self.ref_file_orig = ref_file
        self.input_brain_orig = input_brain
        self.ref_brain_orig = ref_brain

        self.input_file = input_file
        self.ref_file = ref_file
        self.input_brain = input_brain
        self.ref_brain = ref_brain
        self.linear_reg_file = nonlinear_reg_file
        self.nonlinear_reg_file = nonlinear_reg_file


        self.config = {'input_reorient2std':False,
                        'ref_reorient2std':False,
                        'flirt_out_reg_file':'linear_reg.mat',
                        'flirt_out_file':'img2img_linear.nii.gz',
                        'fnirt_out_file':'img2img_nonlinear.nii.gz',
                        'fnirt_fieldcoeff_file':'img2img_nonlinear_fieldcoeff.nii.gz'}

    def change_config(self,trait,value):
        assert trait in self.config ,'trait should be one of {traits}'.format(traits=','.join(list(self.config.keys)))
        self.config[trait] = value

    def _prepare(self,wf_base_dir):

        if self.config['input_reorient2std']:
            reorient = fsl.Reorient2Std()
            reorient.inputs.in_file = self.input_file
            reorient.inputs.out_file = 'input_img_reorient2std.nii.gz'
            reorient_Node = pe.Node(reorient,'input_file')
            wf = pe.Workflow(name='reorient2std',base_dir=wf_base_dir)
            wf.add_nodes([reorient_Node])
            wf.run()
            self.input_file = os.path.join(wf_base_dir,'reorient2std','input_file','input_img_reorient2std.nii.gz')

        if self.config['ref_reorient2std']:
            reorient = fsl.Reorient2Std()
            reorient.inputs.in_file = self.ref_file
            reorient.inputs.out_file = 'ref_img_reorient2std.nii.gz'
            reorient_Node = pe.Node(reorient,'ref_file')
            wf = pe.Workflow(name='reorient2std',base_dir=wf_base_dir)
            wf.add_nodes([reorient_Node])
            wf.run()
            self.ref_file = os.path.join(wf_base_dir,'reorient2std','ref_file','ref_img_reorient2std.nii.gz')

        if self.input_brain is None:
            bet = pe.Node(interface=fsl.BET(in_file=self.input_file,out_file='input_brain.nii.gz'),name='bet')
            wf = pe.Workflow(name='input_bet',base_dir=wf_base_dir)
            wf.add_nodes([bet])
            wf.run()
            self.input_brain = os.path.join(wf_base_dir,'input_bet','bet','input_brain.nii.gz')

        if self.ref_brain is None:
            bet = pe.Node(interface=fsl.BET(in_file=self.ref_file,out_file='ref_brain.nii.gz'),name='bet')
            wf = pe.Workflow(name='ref_bet',base_dir=wf_base_dir)
            wf.add_nodes([bet])
            wf.run()
            self.ref_brain = os.path.join(wf_base_dir,'ref_bet','bet','ref_brain.nii.gz')



    def _linear_reg(self,wf_name,wf_base_dir):

        flirt = pe.Node(interface=fsl.FLIRT(in_file = self.input_brain, reference=self.ref_brain,
                                        out_matrix_file=self.config['flirt_out_reg_file'],
                                        out_file = self.config['flirt_out_file']), name='linear')

        wf = pe.Workflow(name=wf_name, base_dir=wf_base_dir)

        wf.add_nodes([flirt])
        wf.run()

        self.linear_reg_file = os.path.join(wf_base_dir,wf_name, 'linear', self.config['flirt_out_reg_file'])

    def _nonlinear_reg(self,wf_name,wf_base_dir):

        flirt = pe.Node(interface=fsl.FLIRT(in_file = self.input_brain, reference=self.ref_brain,
                                        out_matrix_file=self.config['flirt_out_reg_file'],
                                        out_file = self.config['flirt_out_file']), name='linear')


        fnirt = pe.Node(interface=fsl.FNIRT(in_file=self.input_file, ref_file = self.ref_file,
                                       warped_file = self.config['fnirt_out_file'], fieldcoeff_file = self.config['fnirt_fieldcoeff_file']), name='nonlinear')


        wf = pe.Workflow(name=wf_name, base_dir=wf_base_dir)


        ## adding nodes to workflows
        wf.add_nodes([flirt,fnirt])
        wf.connect([
                    (flirt,fnirt,[('out_matrix_file','affine_file')]),
                    ])
        wf.run()

        self.linear_reg_file = os.path.join(wf_base_dir,wf_name, 'linear', self.config['flirt_out_reg_file'])
        self.nonlinear_reg_file = os.path.join(wf_base_dir,wf_name,'nonlinear',self.config['fnirt_fieldcoeff_file'])

    def run(self, method, wf_name,wf_base_dir):

        ## prepare
        self._prepare(wf_base_dir)

        if method=='linear':
            self._linear_reg(wf_name,wf_base_dir)
        elif method=='nonlinear':
            self._nonlinear_reg(wf_name,wf_base_dir)


    def warp_points(self, coords, method, coord_type='ras'):

        rnum = create_random_fname()
        fname = '{rnum}.txt'.format(rnum=rnum)
        np.savetxt(fname,coords)
        warppoints = fsl.WarpPoints()
        warppoints.inputs.in_coords = fname
        if coord_type=='ras':
            warppoints.inputs.coord_mm = True

        if method=='linear':
            assert self.linear_reg_file is not None, 'The linear registeration file does not exist. Use run() to do the registeration!'
            warppoints.inputs.xfm_file =  self.linear_reg_file

        elif method=='nonlinear':
            assert self.nonlinear_reg_file is not None, 'The nonlinear registeration file does not exist. Use run() to do the registeration!'
            warppoints.inputs.warp_file = self.nonlinear_reg_file

        warppoints.inputs.src_file = self.input_file
        warppoints.inputs.dest_file = self.ref_file

        res = warppoints.run()
        res = np.loadtxt('{rnum}_warped.txt'.format(rnum=rnum))

        ## removing the files
        os.remove(fname)
        os.remove('{rnum}_warped.txt'.format(rnum=rnum))

        return res



def mri_label2vol(label_file, hemi, subject, freesurfer_dir,  wf_name, wf_base_dir,
                    reg_file=None , reg_header=None , identity=None , surface='white',
                    proj_type='frac', start=0, end=1, delta=0.01,
                    fill_thresh = None, tidy_up=False , **kwargs):

    wf = pe.Workflow(name=wf_name, base_dir=wf_base_dir)

    ## args
    args = {}

    ## proj
    proj = (proj_type,start,end,delta)

    if reg_file is not None:
        args['reg_file'] = reg_file
    if reg_header is not None:
        args['reg_header'] = reg_header
    if identity is not None:
        args['identity'] = identity
    if fill_thresh is not None:
        args['fill_thresh'] = fill_thresh

    if ((reg_file is None) and (reg_header is None) and (identity is None)):
        raise('Should specify registeration method: reg_file, reg_header or identity')



    print('args: ',args)

    ## label2vol
    label2vol = pe.Node(Label2Vol(label_file=label_file,
     template_file=os.path.join(freesurfer_dir,subject,'mri/T1.mgz'),
     hemi=hemi, proj=proj, subject_id=subject,
     subjects_dir = freesurfer_dir,
     surface=surface,
    **args),
     name='label2vol',
)

    label_name = os.path.basename(label_file)
    label_name = os.path.splitext(label_name)[0]
    if tidy_up:

        mask_dilate = pe.Node(Binarize(dilate=1,erode=1,min=1),name='dilate_label_vol')
        mris_calc = pe.Node(MRIsCalc(),name='tidy_up')
        mris_calc.inputs.in_file2=os.path.join(freesurfer_dir,subject,'mri/{hemi}.ribbon.mgz'.format(hemi=hemi))
        mris_calc.inputs.action='mul'
        mris_calc.inputs.out_file=label_name+'.nii.gz'
        wf.add_nodes([label2vol, mask_dilate, mris_calc])
        wf.connect([
                    (label2vol,mask_dilate,[("vol_label_file","in_file")]),
                    (mask_dilate,mris_calc,[('binary_file','in_file1')]),
                    ])

        out_vol = '{wf_base_dir}/{wf_name}/tidy_up/{output}'.format(wf_base_dir=wf_base_dir,wf_name=wf_name,output=label_name+'.nii.gz')
    else:
        wf.add_nodes([label2vol])
        out_vol = '{wf_base_dir}/{wf_name}/label2vol/{output}'.format(wf_base_dir=wf_base_dir,wf_name=wf_name,output=label_name+'.nii.gz')

    wf.run()
    return out_vol


def make_head_model(anat_img, out_dir):

    """ create head models including skin and skull.

    Note that rawavg.mgz should be specified ( not orig.mgz) to create head models for the purpose of visualizaiton"""

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if not os.path.exists(os.path.join(out_dir,'outer_skin_surface')):

        cmd ='cd {out_dir}; mri_watershed -surf surf {anat_img} brain.mgz'.format(out_dir=out_dir, anat_img=anat_img)
        os.system(cmd)

        for f in ['lh.surf_brain_surface','lh.surf_inner_skull_surface','lh.surf_outer_skin_surface','lh.surf_outer_skull_surface']:
            cmd = 'mv {out_dir}/{f} {out_dir}/{f2}'.format(f=f,f2=f.split('lh.surf_')[1],out_dir=out_dir)
            os.system(cmd)
    else:
        print('head model exists!')


def create_random_fname():
    import numpy as np
    ## creating temporary file name
    rnum1 = np.random.randint(10**15,10**16)
    rnum2 = np.random.randint(10**10,10**11)
    rnum = '{rnum1}{rnum2}'.format(rnum1=rnum1,rnum2=rnum2)
    return rnum
