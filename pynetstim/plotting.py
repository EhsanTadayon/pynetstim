#### plotting and visualization functions and classes
from surfer import Brain
from mayavi import mlab
from utils import Surf


def plot_points(points, hemi, surf, map_surface= None, annot=None, map_to_annot = None, show_rois=False, scale_factor=.8, show_names=False, opacity=1, background='black', img_basename=None, skin=True):
    
    # adding brain and annotations
    brain = Brain(points.subject, surf=surf, hemi=hemi, subjects_dir=points.subjects_dir, offset=False, background=background)
    
    ## adding skin
    if skin:
        skin_surf = Surf('{subjects_dir}/{subject}/bem/lh.watershed_outer_skin_surface'.format(subjects_dir=points.subjects_dir, subject=points.subject))
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
        points.names, points.colors = points.map_to_annot(map_to_annot) 

    ### show points
    
    for point in points:
        if point.hemi==hemi or hemi=='both':
            brain.add_foci(point.ras_tkr_coord, hemi = point.hemi, color = point.color, scale_factor = scale_factor, alpha = opacity, map_surface = map_surface)
            if show_rois:
                brain.add_label(point.roi, hemi=point.hemi, color=point.roi.color)
            if show_names:
                mlab.text3d(point.ras_tkr_coord[0],point.ras_tkr_coord[1],point.ras_tkr_coord[2],point.name,scale=4)


    if img_basename:
        brain.save_imageset(prefix='{figures_dir}/img_basename'.format(figures_dir=points.figures_dir),views = ['lat','med','caud','dor'],filetype='png')
    
    mlab.show()
             

        

    
        