### Utils
import os

def make_head_model(subject, freesurfer_dir):
    try:
        os.makedirs('{freesurfer_dir}/{subject}/bem'.format(freesurfer_dir=freesurfer_dir, subject=subject))
        cmd ='cd {freesurfer_dir}/{subject}/bem; mri_watershed -surf surf {freesurfer_dir}/{subject}/mri/rawavg.mgz brain.mgz'.format(freesurfer_dir=freesurfer_dir, subject = subject)
        os.system(cmd)

        for f in ['lh.surf_brain_surface','lh.surf_inner_skull_surface','lh.surf_outer_skin_surface','lh.surf_outer_skull_surface']:
            cmd = 'mv {freesurfer_dir}/{subject}/bem/{f} {freesurfer_dir}/{subject}/bem/{f2}'.format(f=f,f2=f.split('lh.surf_')[1],subject=subject, freesurfer_dir=freesurfer_dir)
            os.system(cmd)
    except OSError:
        print 'head model exists, if you want to recreate the head models, remove the directory:\n{d}'.format(d='{freesurfer_dir}/{subject}/bem'.format(freesurfer_dir=freesurfer_dir, subject=subject))
    
