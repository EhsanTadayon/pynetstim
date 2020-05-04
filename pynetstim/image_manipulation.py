from nipype.interfaces import fsl
import nipype.pipeline.engine as pe
from nipype.interfaces.freesurfer import Label2Vol,Binarize,MRIsCalc,MRIConvert
import os
import numpy as np

def img2img_register(img_file, ref_file, wf_base_dir, wf_name,
                    input_reorient2std=False, ref_reorient2std=False,
                    method='linear',
                    flirt_out_reg_file = 'linear_reg.mat', flirt_out_file = 'img2img_linear.nii.gz',
                    fnirt_out_file = 'img2img_nonlinear.nii.gz', fnirt_fieldcoeff_file = 'img2img_nonlinear_fieldcoeff.nii.gz'):
    
    img_file = os.path.abspath(img_file)
    ref_file = os.path.abspath(ref_file)
    wf_base_dir = os.path.abspath(wf_base_dir)
    
    if not img_file.endswith (('.nii','.nii.gz')):
        raise("input file should be in nifti format (.nii or .nii.gz)!")
    if not ref_file.endswith(('.nii','.nii.gz')):
        raise("destination file (ref_file) should be in nifti format (.nii or .nii.gz)!")
    
    
    ### reorient2std
    
    ### reorient img_file
    if input_reorient2std:
        reorient = fsl.Reorient2Std()
        reorient.inputs.in_file = img_file
        reorient.inputs.out_file = 'input_img_reorient2std.nii.gz'
        reorient_Node = pe.Node(reorient,'input_file')
        wf = pe.Workflow(name='reorient2std',base_dir=wf_base_dir)
        wf.add_nodes([reorient_Node])
        wf.run()
        img_file = os.path.join(wf_base_dir,'reorient2std','input_file','input_img_reorient2std.nii.gz')
    
    ### reorient ref_file if config is not MNI
    if ref_reorient2std and config!='MNI':
        reorient = fsl.Reorient2Std()
        reorient.inputs.in_file = ref_file
        reorient.inputs.out_file = 'ref_img_reorient2std.nii.gz'
        reorient_Node = pe.Node(reorient,'ref_file')
        wf = pe.Workflow(name='reorient2std',base_dir=wf_base_dir)
        wf.add_nodes([reorient_Node])
        wf.run()
        ref_file = os.path.join(wf_base_dir,'reorient2std','ref_file','ref_img_reorient2std.nii.gz')

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
                                       warped_file = fnirt_out_file, fieldcoeff_file = fnirt_fieldcoeff_file), name='nonlinear')


        
        
        wf = pe.Workflow(name=wf_name, base_dir=wf_base_dir)


        ## adding nodes to workflows
        wf.add_nodes([flirt,fnirt])
        wf.connect([(flirt,fnirt,[('out_matrix_file','affine_file')])])
        wf.run()
        
    linear_reg_file = os.path.join(wf_base_dir,wf_name, method, flirt_out_reg_file)
    nonlinear_warp_field_file = os.path.join(wf_base_dir,wf_name,method,fnirt_fieldcoeff_file)
        
    return {'img_file':img_file,
            'ref_file':ref_file,
            'linear_reg_file':linear_reg_file,
            'warp_field_file':nonlinear_warp_field_file}
            
            
            
            
def img2img_coord_register(ras_coords, img_file, dest_img, wf_base_dir, method='linear',
                     input_reorient2std=False, ref_reorient2std=False, wf_name='register',
                     linear_reg_file=None, warp_field_file = None):
        
        if method not in ['linear','nonlinear']:
            raise('method should be either linear or nonlinear')      
        
        rnum1 = np.random.randint(10**15,10**16) 
        rnum2 = np.random.randint(10**10,10**11)  
        rnum = '{rnum1}_{rnum2}'.format(rnum1=rnum1,rnum2=rnum2)       
        fname='./temp_{rnum}_coords.txt'.format(rnum=rnum)
            
        np.savetxt(fname,ras_coords)
        warppoints = fsl.WarpPoints()
        warppoints.inputs.in_coords = fname
        
        if linear_reg_file is None and warp_field_file is None:
            results = img2img_register(img_file = img_file, ref_file = dest_img, wf_base_dir = wf_base_dir, wf_name=wf_name, method=method,
                            input_reorient2std=input_reorient2std, ref_reorient2std=ref_reorient2std, 
                            flirt_out_reg_file = 'linear_reg.mat',flirt_out_file = 'img2img_linear.nii.gz', 
                            fnirt_out_file = 'img2img_nonlinear.nii.gz')
                    
            if method=='linear':
                warppoints.inputs.xfm_file =  results['linear_reg_file']
            elif method=='nonlinear':
                warppoints.inputs.warp_file = results['warp_field_file']
                
            warppoints.inputs.src_file = results['img_file']
            warppoints.inputs.dest_file = results['ref_file']
                  
        else:
            
            if method=='linear':
                warppoints.inputs.xfm_file =  linear_reg_file     
            elif method=='nonlinear':
                warppoints.inputs.warp_file = warp_field_file
            warppoints.inputs.src_file = img_file
            warppoints.inputs.dest_file = dest_img
            
        warppoints.inputs.coord_mm = True
        res = warppoints.run()
        res = np.loadtxt('./{fname}_warped.txt'.format(fname=fname.split('.txt')[0]))
        
        ## removing the files
        os.remove(fname)
        os.remove('./{fname}_warped.txt'.format(fname=fname.split('.txt')[0]))
        
        return res


def mri_label2vol(label, subject, freesurfer_dir, wf_base_dir, wf_name, proj=(u'frac',0,1,0.01), identity=True, tidy_up=True,**kwargs):
    
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
        mris_calc.inputs.out_file=label.name+'-'+label.hemi+'.nii.gz'
        wf.add_nodes([label2vol, mask_dilate, mris_calc])
        wf.connect([
                    (label2vol,mask_dilate,[("vol_label_file","in_file")]),
                    (mask_dilate,mris_calc,[('binary_file','in_file1')]),
                    ])
                    
        out_vol = '{wf_base_dir}/{wf_name}/mask_with_gm/{output}'.format(wf_base_dir=wf_base_dir,wf_name=wf_name,output=label.name+'-'+label.hemi+'.nii.gz')
    else:
        wf.add_nodes([label2vol])
        out_vol = '{wf_base_dir}/{wf_name}/label2vol/{output}'.format(wf_base_dir=wf_base_dir,wf_name=wf_name,output=label.name+'-'+label.hemi+'.nii.gz')
    
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
        