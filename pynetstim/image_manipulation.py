from nipype.interfaces import fsl
import nipype.pipeline.engine as pe
from nipype.interfaces.freesurfer import Label2Vol,Binarize,MRIsCalc
import os

def img2img_register(img_file, ref_file, wf_base_dir, wf_name, method='linear',
                    flirt_out_reg_file = 'linear_reg.mat', flirt_out_file = 'img2img_linear.nii.gz',
                    fnirt_out_file = 'img2img_nonlinear.nii.gz'):
    
    
    if method=='linear':
        
        flirt = pe.Node(interface=fsl.FLIRT(in_file=img_file, reference=ref_file,
                                        out_matrix_file=flirt_out_reg_file, 
                                        out_file = flirt_out_file), name='linear')
        
        wf = pe.Workflow(name=wf_name, base_dir=wf_base_dir)
     
        wf.add_nodes([flirt])
        wf.run()
        
    if method=='nonlinear':
        flirt = pe.Node(interface=fsl.FLIRT(in_file=img_file, reference=ref_file,
                                        out_matrix_file=flirt_out_reg_file, 
                                        out_file = flirt_out_file), name='linear')
        fnirt = pe.Node(interface=fsl.FNIRT(in_file=img_file, ref_file = ref_file, 
                                       warped_file = fnirt_out_file), name='nonlinear')

        wf = pe.Workflow(name=wf_name, base_dir=wf_base_dir)


        ## adding nodes to workflows
        wf.add_nodes([flirt,fnirt])
        wf.connect([(flirt,fnirt,[('out_matrix_file','affine_file')])])
        wf.run()


def mri_label2vol(label, subject, freesurfer_dir, wf_base_dir, wf_name, proj=(u'frac',0,1,0.01), identity=True, tidy_up=True):
    
    wf = pe.Workflow(name=wf_name, base_dir=wf_base_dir)
    
    # save label first
    label_file = '{wf_base_dir}/{label_name}-{hemi}.label'.format(wf_base_dir=wf_base_dir,
                                                                  label_name=label.name,
                                                                  hemi=label.hemi)
    label.save(label_file)
    
    label2vol = pe.Node(Label2Vol(label_file=label_file,
     template_file=os.path.join(freesurfer_dir,subject,'mri/T1.mgz'),
     hemi=label.hemi, proj=proj, identity=identity, subject_id=subject), name='label2vol')
    
    if tidy_up:    
        
        mask_dilate = pe.Node(Binarize(dilate=1,erode=1,min=1),name='dilate_label_vol')
        mris_calc = pe.Node(MRIsCalc(),name='mask_with_gm')
        mris_calc.inputs.in_file2=os.path.join(freesurfer_dir,subject,'mri/{hemi}.ribbon.mgz'.format(hemi=label.hemi))
        mris_calc.inputs.action='mul'
        mris_calc.inputs.out_file=label.name+'.nii.gz'
        wf.add_nodes([label2vol, mask_dilate, mris_calc])
        wf.connect([
                    (label2vol,mask_dilate,[("vol_label_file","in_file")]),
                    (mask_dilate,mris_calc,[('binary_file','in_file1')]),
                    ])
    else:
        wf.add_nodes([label2vol])
    
    wf.run()