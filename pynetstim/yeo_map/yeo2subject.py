import numpy as np
from nipype. interfaces.freesurfer import SurfaceTransform
import os




def yeo2subject(subject, freesurfer_dir,nNet=7):

    os.environ['SUBJECTS_DIR'] = freesurfer_dir

    if not os.path.exists('{freesurfer_dir}/fsaverage'.format(freesurfer_dir=freesurfer_dir)):
        raise('fsaverage does not exist under freesurfer_dir!')


    for hemi in ['lh','rh']:
        sxfm = SurfaceTransform()
        annot = os.path.join(freesurfer_dir,'fsaverage','label',
                            '{hemi}.Yeo2011_{nNet}Networks_N1000.annot'.format(hemi=hemi,nNet=nNet))
        sxfm.inputs.source_annot_file = annot
        sxfm.inputs.source_subject = 'fsaverage'
        sxfm.inputs.target_subject =  subject
        sxfm.inputs.hemi = hemi
        out_file = os.path.join(freesurfer_dir,subject,'label','{hemi}.Yeo2011_{nNet}Networks_N1000.annot'.format(hemi=hemi,nNet=nNet))
        sxfm.inputs.out_file = out_file
        sxfm.inputs.subjects_dir = freesurfer_dir
        sxfm.run()


    cmd = "for hemi in lh rh;do mri_annotation2label --subject {subject} --hemi $hemi --annotation Yeo2011_{nNet}Networks_N1000 --outdir {freesurfer_dir}/{subject}/label; done"
    os.system(cmd.format(subject=subject, nNet=nNet, freesurfer_dir=freesurfer_dir))


    ## confidence
    for hemi in ['lh','rh']:
        cmd = 'mri_surf2surf --srcsubject fsaverage --sval {freesurfer_dir}/fsaverage/label/{hemi}.Yeo2011_{nNet}NetworksConfidence_N1000.mgz --hemi {hemi} --trgsubject {subject} --tval {freesurfer_dir}/{subject}/{hemi}.Yeo2011_{nNet}NetworksConfidence_N1000.mgh'
        os.system(cmd.format(freesurfer_dir=freesurfer_dir,hemi=hemi,nNet=nNet,subject=subject))
