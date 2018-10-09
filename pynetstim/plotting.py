"""
plotting and visualization functions and classes
Author: Ehsan Tadayon, M.D. [sunny.tadayon@gmail.com / stadayon@bidmc.harvard.edu]
"""

from surfer import Brain
from mayavi import mlab
from utils import Surf,FreesurferSurf
from nilearn.plotting import plot_anat
import matplotlib.pyplot as plt
import numpy as np
import os

def plot_targets(targets, hemi, surf, map_surface= None, annot=None, map_to_annot = None, show_rois=False, scale_factor=.8, show_names=False, opacity=1, background='black',
                img_basename=None, skin=True):
    
    # adding brain and annotations
    brain = Brain(targets.subject, surf=surf, hemi=hemi, subjects_dir=targets.subjects_dir, offset=False, background=background)
    
    ## adding skin
    if skin:
        skin_surf = Surf('{subjects_dir}/{subject}/bem/lh.watershed_outer_skin_surface'.format(subjects_dir=targets.subjects_dir, subject=targets.subject))
        mlab.triangular_mesh(skin_surf.vertices[:,0], skin_surf.vertices[:,1], skin_surf.vertices[:,2], skin_surf.faces, opacity=0.2, color=(1,1,0))

    ### annotation
    if annot:
        if hemi in ['lh','rh']:
            brain.add_annotation(annot, hemi=hemi, borders=False)
            
        elif hemi=='both':
            brain.add_annotation(annot, hemi='lh', borders=False)
            brain.add_annotation(annot, hemi='rh', borders=False, remove_existing=False)
            
    ### map to annot
    if map_to_annot:
        targets.names, targets.colors = targets.map_to_annot(map_to_annot) 

    ### show targets
    
    for target in targets:
        
        if target.hemi==hemi or hemi=='both':
            brain.add_foci(target.ras_tkr_coord, hemi = target.hemi, color = target.color, scale_factor = scale_factor, alpha = opacity, map_surface = map_surface)
            if show_rois:
                brain.add_label(target.roi, hemi=target.hemi, color=target.roi.color)
            if show_names:
                mlab.text3d(target.ras_tkr_coord[0],target.ras_tkr_coord[1],target.ras_tkr_coord[2],target.name,scale=4)


    if img_basename:
        brain.save_imageset(prefix='{figures_dir}/img_basename'.format(figures_dir=targets.figures_dir),views = ['lat','med','caud','dor'],filetype='png')
    
    mlab.show()
    
    
def plot_samples(samples, hemi='both', surf='pial', annot=None,annot_alpha=1, map_surface=None, show_average=False, background='white', out_dir=None, prefix=''):
    """
    plot individual pulses. 
    """
    
    # adding brain and annotations
    brain = Brain(samples.subject, surf='pial', hemi='both', subjects_dir=samples.subjects_dir, offset=False, background=background)
    if annot:
        if hemi=='both':
            brain.add_annotation(annot,hemi='lh',borders=False, alpha=annot_alpha)
            brain.add_annotation(annot,hemi='rh',borders=False,remove_existing=False, alpha=annot_alpha)
        if hemi in ['lh','rh']:
            brain.add_annotation(annot,hemi=hemi,borders=False, alpha=annot_alpha)
            
    skin_surf = Surf('{subjects_dir}/{subject}/bem/lh.watershed_outer_skin_surface'.format(subjects_dir=samples.subjects_dir, subject=samples.subject))
    mlab.triangular_mesh(skin_surf.vertices[:,0], skin_surf.vertices[:,1], skin_surf.vertices[:,2], skin_surf.faces, opacity=0.2, color=(1,1,0))
    
    
    
    mlab.points3d(samples.coordinates['ras_tkr_coords'][:,0],samples.coordinates['ras_tkr_coords'][:,1],samples.coordinates['ras_tkr_coords'][:,2],scale_factor=7,reset_zoom=False,
    color=(1,0,0),opacity=.2)
    
    
    #### show the average coordinate
    if show_average:
        
        avg_coord = np.mean(samples.coordinates['ras_tkr_coords'],axis=0)
        mlab.points3d(avg_coord[0],avg_coord[1],avg_coord[2],color=(0,0,1),scale_factor=10)
    
    
    ### project to surface
    if map_surface is not None:
        lh_map_surf = FreesurferSurf('lh', map_surface, samples.subject, samples.subjects_dir)
        rh_map_surf = FreesurferSurf('rh', map_surface, samples.subject, samples.subjects_dir)
        
        lh_samples = samples.coordinates['ras_tkr_coords'][samples.hemis=='lh']
        rh_samples = samples.coordinates['ras_tkr_coords'][samples.hemis=='rh']
        
        if lh_samples.shape[0]>0:
            lh_mapped_coords = lh_map_surf.project_coords(lh_samples)[1]
            mlab.points3d(lh_mapped_coords[:,0],lh_mapped_coords[:,1],lh_mapped_coords[:,2],scale_factor=7,color=(0,1,1),reset_zoom=False)
            
        if rh_samples.shape[0]>0:    
            rh_mapped_coords = rh_map_surf.project_coords(rh_samples)[1]
            mlab.points3d(rh_mapped_coords[:,0],rh_mapped_coords[:,1],rh_mapped_coords[:,2],scale_factor=7,color=(0,1,1), reset_zoom=False)
        
    ### save
    if out_dir is None:
        mlab.show()
    else:
        fname = os.path.join(out_dir,prefix)
        brain.save_imageset(fname,['med', 'lat', 'dor','caud','ros'], 'png')
        saved_images = [fname+'_'+view+'.png' for view in ['med','lat','dor','caud','ros']]
    
            
        fig,ax = plt.subplots(1,len(saved_images),figsize=(12,6))
        for i, img in enumerate(saved_images):
            img = plt.imread(img)
            ax[i].imshow(img)
            ax[i].set_axis_off()
        fig.savefig(fname+'.png',dpi=600,bbox_inches='tight')
        plt.close(fig)
        for img in saved_images:
            os.remove(img)   
        mlab.close()

        
        
        
       
    
    
    
    
    


        
             

        

    
        