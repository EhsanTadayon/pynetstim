"""
plotting and visualization functions and classes
Author: Ehsan Tadayon, M.D. [sunny.tadayon@gmail.com / stadayon@bidmc.harvard.edu]
"""

from surfer import Brain
from mayavi import mlab
from .freesurfer_files import Surf, FreesurferSurf, Annot
from nilearn.plotting import plot_anat
import matplotlib.pyplot as plt
import numpy as np
import os
import copy
from .image_manipulation import make_head_model


class plotting_points(object):
    
    def __init__(self, points, hemi='both', surf='pial', map_surface=None, annot=None, map_to_annot=None, show_skin=True, show_roi=False,
        show_name=False, name_scale=1, name_color=(0,0,0), show_directions=False, background='white',
        show_average=False, opacity=1, scale_factor=.5,color=np.array([(1,0,0)]),
        use_default=True):
     
        self.points = copy.deepcopy(points)
        self.subject = points.subject
        self.freesurfer_dir = points.freesurfer_dir
        self.hemi = hemi
        self.surf = surf
        self.background = background
        self.map_surface=map_surface
        self.annot = annot
        self.map_to_annot = map_to_annot
        self.show_skin = show_skin
        self.show_roi = show_roi
        self.scale_factor = scale_factor
        self.show_directions = show_directions
        self.show_name = show_name
        self.show_average = show_average
        self.opacity = opacity
        self.color=color
        self.scale_factor = scale_factor
        self.use_default = use_default
        self.name_scale = name_scale
        self.name_color = name_color

        
        self._set_config()
        self._plot()
        
        
    def _set_config(self):
    
        
        if not hasattr(self.points,'color'):
            if self.color.shape[0]==1:
                self.points.add_trait('color', np.repeat(self.color,self.points.npoints,axis=0))
                
            elif self.color.shape[0]==self.points.npoints:
                self.points.add_trait('color',self.color)
            else:
                raise 'Color size is incorrect'
            
        if not hasattr(self.points, 'scale_factor'):
            self.points.add_trait('scale_factor', np.ones(self.points.npoints)*self.scale_factor)
              
        if not hasattr(self.points,'opacity'):
            self.points.add_trait('opacity',np.ones(self.points.npoints)*self.opacity)
        
        
    def _plot(self):
        
        ## plot
        self.brain = Brain(hemi=self.hemi, surf=self.surf, subject_id=self.subject, subjects_dir=self.freesurfer_dir, background = self.background,offset=False)
         
        if self.show_skin:
            self._show_skin()
        
        if self.annot:
            self._show_annot()
  
        self._add_points()
        
        if self.show_average:
            self._show_average()
            
        
    def show(self):
        mlab.show()
    
         
         
    def _add_points(self):
                    
        if self.map_surface:
            temp = self.points.map_to_surface(self.map_surface)
            mapped_vertices, mapped_coords_ras_tkr = temp['vertices'], temp['ras_tkr_coord']
            mapped_coords = mapped_coords_ras_tkr
        
        ### FreesurferCoords.map_to_annot requires map_surface to be set ( set to white by default). Here, we use map_surface if it has been set, otherwise use 'white' by default. 
        
        if self.map_to_annot:            
            if self.map_surface:
                self.points.name, self.points.color = self.points.map_to_annot(self.map_to_annot, map_surface=self.map_surface)
            else:
                self.points.name, self.points.color = self.points.map_to_annot(self.map_to_annot, map_surface='white')
                
        
        for i in range(self.points.npoints):
            point = self.points[i]
        
            if point.hemi==self.hemi or self.hemi=='both':
               
                if self.map_surface:
                    self.brain.add_foci(mapped_coords[i,:], hemi=point.hemi, color=point.color, scale_factor=point.scale_factor, alpha=point.opacity)
                else:
                    self.brain.add_foci(point.ras_tkr_coord, hemi=point.hemi, color=point.color, scale_factor=point.scale_factor, alpha=point.opacity)
                
                if self.show_roi and hasattr(point,'roi'):
                    self.brain.add_label(point.roi, hemi=point.hemi, color=point.roi.color)
                    
                if self.show_name and hasattr(point,'name'):
                    mlab.text3d(point.ras_tkr_coord[0], point.ras_tkr_coord[1], point.ras_tkr_coord[2], point.name, scale=self.name_scale, color=self.name_color)
                    
                if self.show_directions and hasattr(point,'direction'):
                    origin = point.ras_tkr_coord.flatten().tolist()
                    X,Y,Z = zip(origin,origin,origin)
                    p0 = point.direction[0:3]
                    p1 = point.direction[3:6]
                    p2 = point.direction[6:9]
                    U,W,V = zip(p0,p1,p2)
                    plot_directions(X,Y,Z,U,W,V)
                            
                
    def _show_skin(self):
        
        ## create head model under freesurfer_dir/{subject}/bem (if none existing!)
        anat_img = os.path.join(self.points.freesurfer_dir,self.points.subject,'mri','rawavg.mgz')
        out_dir = os.path.join(self.points.freesurfer_dir,self.points.subject,'bem')
        make_head_model(anat_img,out_dir)
        
        skin_surf = Surf('{freesurfer_dir}/{subject}/bem/outer_skin_surface'.format(freesurfer_dir=self.freesurfer_dir, subject=self.subject))
        mlab.triangular_mesh(skin_surf.vertices[:,0], skin_surf.vertices[:,1], skin_surf.vertices[:,2], skin_surf.faces, opacity=0.2, color=(1,1,0))
        
    def _show_annot(self):
        
        if self.hemi in ['lh','rh']:
            self.brain.add_annotation(self.annot, hemi=self.hemi, borders=False)
            
        elif self.hemi=='both':
            self.brain.add_annotation(self.annot, hemi='lh', borders=False)
            self.brain.add_annotation(self.annot, hemi='rh', borders=False, remove_existing=False)
            
    def _show_average(self,scale_factor=6, color=(0,0,1)):
        
        avg_coord = np.mean(self.points.coordinates['ras_tkr_coord'],axis=0)
        mlab.points3d(avg_coord[0],avg_coord[1], avg_coord[2], color=color, scale_factor=scale_factor)
        
    def save_image(self, out_dir, views=['lat','med','dor','ros','caud'], prefix='',  filetype='jpg'):
        
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
            
        prefix = out_dir+'/'+ prefix
        
        self.brain.save_imageset(prefix=prefix,views=views,filetype=filetype)
        mlab.close()
        
        fig,ax = plt.subplots(1,len(views),figsize=(15,8))
        for i,view in enumerate(views):
            img_name = prefix+'_'+view+'.'+filetype
            img = plt.imread(img_name)
            ax[i].imshow(img)
            ax[i].set_axis_off()
            os.remove(img_name)
            
        fig.savefig(prefix+'.png',bbox_inches='tight',dpi=600)
        plt.close()
        
        ## return saved image path
        return prefix+'.png'



class plotting_points_fast(plotting_points):
    
    
    def _add_points(self):
        
        if self.map_to_annot:
            self.points.name, self.points.color = self.points.map_to_annot(self.map_to_annot) 
            
        if self.map_surface:
            mapped_vertices, mapped_coords_ras_tkr, mapped_coords_ras = self.points.map_to_surface(self.map_surface)
            mapped_coords = mapped_coords_ras_tkr
            mlab.points3d(mapped_coords[:,0],mapped_coords[:,1],mapped_coords[:,2],
                      scale_factor = self.scale_factor*8,reset_zoom=False, color=(1,0,0), opacity = self.opacity)
        else:
            mlab.points3d(self.points.coordinates['ras_tkr_coord'][:,0],self.points.coordinates['ras_tkr_coord'][:,1],self.points.coordinates['ras_tkr_coord'][:,2],
                      scale_factor = self.scale_factor*8,reset_zoom=False, color=(1,0,0), opacity = self.opacity)
            
        
    
        
def plot_directions(x,y,z,u,v,w):
    scales=[20,20,25]
    lines_width=[2,3,1]
    color=[(1,0,0),(0,1,0),(0,0,1)]
    
    obj = mlab.quiver3d(x, y, z, u, v, w,
            line_width=3, colormap='hsv', 
            scale_factor=0.8, mode='arrow')

    for i in range(len(x)):
        r,g,b = color[i]
        #print("R: {}, G: {}, B: {}".format(r,g,b))
        obj = mlab.quiver3d(x[i], y[i], z[i], u[i], v[i], w[i],
                line_width=lines_width[i], color=(r,g,b), colormap='hsv', 
                scale_factor=scales[i], mode='arrow')        
       
    
    
    
    
    


        
             

        

    
        