"""
Main stimulation project functions and classes
"""

from utils import Coords, FreesurferCoords, Surf, FreesurferSurf
from brainsight import BrainsightSessionFile, BrainsightSamples, BrainsightTargets
import os
from pymisc.htmlreport import HtmlDoc
from mayavi import mlab
from surfer import Brain



## TO DO:
#### define a Target class ( and brainsight targets will inherit from it) 
### define a Sample class ( and BrainsightSamples will inherit from it)
#### for plot_points, you can speicify map_surface to map your points to a surface and also, you can include ROI for each point 

class StimProject(object):
    
    """ main stimulation project class"""
    def __init__(self, subject, freesurfer_dir, project_dir, brainsight_file = None):
        
        self.subject = subject

        if not os.path.isabs(project_dir):
            raise ValueError('please provide absolute path for project_dir')
        if not os.path.isabs(freesurfer_dir):
            raise ValueError('please provide absolute path for freesurfer_dir')
    
        self.project_dir = project_dir
        self.subject_dir = os.path.join(self.project_dir,subject)
        self.freesurfer_dir = freesurfer_dir
        self.targets = []
        
        
        
        
        
        self._start_project()
        self.brainsight_file = brainsight_file
        self._read_brainsight_file()
        
    def _start_project(self):
        
        # TO DO: add more folders
        if not os.path.exists(self.subject_dir):
            os.makedirs('{subject_dir}/brainsight'.format(subject_dir=self.subject_dir))
            os.makedirs('{subject_dir}/figures'.format(subject_dir = self.subject_dir))
        
        self.brainsight_dir='{subject_dir}/brainsight'.format(subject_dir=self.subject_dir)
        self.figures_dir = '{subject_dir}/figures'.format(subject_dir=self.subject_dir)  
    
    def _read_brainsight_file(self):
        
        if self.brainsight_file is not None:
            bs = BrainsightSessionFile(self.brainsight_file,out_dir='{subject_dir}/brainsight'.format(subject_dir=self.subject_dir))
            self.brainsight_samples = BrainsightSamples('{subject_dir}/brainsight/samples.txt'.format(subject_dir=self.subject_dir))
            self.brainsight_targets = BrainsightTargets('{subject_dir}/brainsight/targets.txt'.format(subject_dir=self.subject_dir))
            
     
    def plot_points(self, points, surf, hemi, annot, scale_factor=6, opacity=.5, background='black', img_basename=None, skin=True):
        
        # adding brain and annotations
        brain = Brain(self.subject, surf=surf, hemi=hemi, subjects_dir=self.freesurfer_dir, offset=False, background=background)
        
        ## adding skin
        skin_surf = Surf('{freesurfer_dir}/{subject}/bem/lh.watershed_outer_skin_surface'.format(freesurfer_dir=self.freesurfer_dir, subject=self.subject))
        mlab.triangular_mesh(skin_surf.vertices[:,0], skin_surf.vertices[:,1], skin_surf.vertices[:,2],skin_surf.faces,opacity=0.2,color=(1,1,0))
        
        if hemi in ['lh','rh']:
            brain.add_annotation(annot, hemi=hemi, borders=False)
        elif hemi=='both':
            brain.add_annotation(annot, hemi='lh', borders=False)
            brain.add_annotation(annot, hemi='rh', borders=False, remove_existing=False)
        
        mlab.points3d(points.coordinates['ras_tkr_coords'][:,0], points.coordinates['ras_tkr_coords'][:,1],
                                                                 points.coordinates['ras_tkr_coords'][:,2], scale_factor=scale_factor, color=(1,0,0),reset_zoom=False, opacity=opacity)
        
        
        
        if img_basename:
            brain.save_imageset(prefix='{figures_dir}/img_basename'.format(figures_dir=self.figures_dir),views = ['lat','med','caud','dor'],filetype='png')
        
        mlab.show()
        
     
     
     
     
     
     
            
                  
    
    def summary(self):
       
       if self.brainsight_file:
           self._brainsight_summary()
   
   
   
   
   
    def _brainsight_summary(self):
   
       """ write summary of samples """
        

       html = HtmlDoc('samples summary') 
       html.add_image('../docs/stim_fig.png',200,800,'middle')
       html.add_header('h1','Brainsight sessions summary')
       html.add_header('h3','Number of session: %d'%len(self.brainsight_samples.sessions_names))

       for session in self.brainsight_samples.sessions_names:
   
           html.add_header('h2',session)
           targets = self.brainsight_samples.get_targets(session)
           html.add_paragraph('<b>Targets:</b> '+', '.join(targets)+'<br><br>')
   
           stims_sequences = self.brainsight_samples.get_stims_sequence(session)
           stims_sequences_sorted = sorted(stims_sequences)
        
           for seq in stims_sequences_sorted:
               t = stims_sequences[seq]
               html.add_paragraph('&nbsp&nbsp&nbsp -> Target: '+ t + ' (start:' + str(seq[0])+ ', end:' + str(seq[1]) + ', stimulations:' + str(seq[1]-seq[0]+1)+')<br><br>')
   
           html.write('{subject_dir}/brainsight/summary.html'.format(subject_dir=self.subject_dir))
        

        
            
            

    
    
        
    
        
        
        