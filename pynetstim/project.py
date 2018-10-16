"""
Main stimulation project functions and classes
Author: Ehsan Tadayon, M.D. [sunny.tadayon@gmail.com / stadayon@bidmc.harvard.edu]
"""

from coordinates import Coords, FreesurferCoords
from surface import Surf, FreesurferSurf
from brainsight import BrainsightSessionFile, BrainsightSamples, BrainsightTargets
import os
from pymisc.htmlreport import HtmlDoc
from nipype.interfaces.fsl import FLIRT,FNIRT
from plotting import plot_samples
import numpy as np
from scipy.spatial.distance import cdist

## TO DO:
#### define a Target class ( and brainsight targets will inherit from it) 
### define a Sample class ( and BrainsightSamples will inherit from it)
#### register the subject to MNI using flirt and fnirt. 


class StimProject(object):
    
    """ main stimulation project class"""
    def __init__(self, subject, project_dir, anatomical = None, freesurfer_dir = None, brainsight_file = None, to_mni = False,
     mni_template='MNI152_T1_2mm.nii.gz', mni_directory = os.path.join(os.environ['FSLDIR'],'data/standard')):
        
        self.subject = subject

        if not os.path.isabs(project_dir):
            raise ValueError('please provide absolute path for project_dir')
        if not os.path.isabs(freesurfer_dir):
            raise ValueError('please provide absolute path for freesurfer_dir')
            
        self.anatomical = anatomical
        self.project_dir = project_dir
        self.subject_dir = os.path.join(self.project_dir,subject)
        self.freesurfer_dir = freesurfer_dir
        self.brainsight_file = brainsight_file
        self.targets = []
    
        self._start_project()
        
        if to_mni:
            self._to_mni(mni_template)    
            
        self._read_brainsight_file()
        
    def _start_project(self):
        
        # TO DO: add more folders
        if not os.path.exists(self.subject_dir):
            
            
            os.makedirs('{subject_dir}/brainsight'.format(subject_dir=self.subject_dir))
            os.makedirs('{subject_dir}/figures'.format(subject_dir = self.subject_dir))
            os.makedirs('{subject_dir}/mni'.format(subject_dir = self.subject_dir))
            os.makedirs('{subject_dir}/anatomical'.format(subject_dir=self.subject_dir))
        
        self.brainsight_dir='{subject_dir}/brainsight'.format(subject_dir=self.subject_dir)
        self.figures_dir = '{subject_dir}/figures'.format(subject_dir=self.subject_dir)  
        self.mni_nonlinear = '{subject_dir}/mni'.format(subject_dir=self.subject_dir)


    def _to_mni(self,to_mni,mni_template):
        pass #### TO DO    
        
        
    def _read_brainsight_file(self):
        
        if self.brainsight_file is not None:
            bs = BrainsightSessionFile(self.brainsight_file,out_dir='{subject_dir}/brainsight'.format(subject_dir=self.subject_dir))
            self.brainsight_samples = BrainsightSamples('{subject_dir}/brainsight/samples.txt'.format(subject_dir=self.subject_dir))
            self.brainsight_targets = BrainsightTargets('{subject_dir}/brainsight/targets.txt'.format(subject_dir=self.subject_dir), self.subject, self.freesurfer_dir)
            

    def summary(self,plot_pulses=False, overwrite=False):
       
       if self.brainsight_file:
           if os.path.exists('{subject_dir}/brainsight/summary.html'.format(subject_dir=self.subject_dir))==False or overwrite==True:
               self._brainsight_summary(plot_pulses)
           else:
               print 'brainsight summary exists.If you wish to overwrite, change overwrite to True!'

    def _brainsight_summary(self,plot_pulses):
   
       """ write summary of samples """

       html = HtmlDoc('samples summary') 
       logo = os.path.abspath('../docs/logo.png')
       html.add_image(logo,200,800,'middle')
       html.add_header('h1','Brainsight sessions summary')
       html.add_header('h3','Number of session: %d'%len(self.brainsight_samples.sessions_names))

       for session in self.brainsight_samples.sessions_names:
   
           html.add_header('h2',session)
           targets = self.brainsight_samples.get_targets(session)
           html.add_paragraph('<b>Targets:</b> '+', '.join(targets)+'<br><br>')
   
           stims_sequences = self.brainsight_samples.get_stims_sequence(session)
           stims_sequences_sorted = sorted(stims_sequences)
        
           for j,seq in enumerate(stims_sequences_sorted):
               t = stims_sequences[seq]
               
               if plot_pulses:
                   pulses_= self.brainsight_samples.get_samples(session).iloc[seq[0]:seq[1]]
                   pulses_coords = pulses_[['loc_x','loc_y','loc_z']]
                   avg_coord = np.mean(pulses_coords,axis=0)
                   d = cdist(np.atleast_2d(avg_coord),pulses_coords)
                   idx = d<=np.mean(d)+3*np.std(d)
                   idx = idx.flatten()
                   pulses_coords = pulses_coords[idx]
                   pulses = FreesurferCoords(pulses_coords,self.subject,self.freesurfer_dir)
                   prefix=session+'_'+t+'_'+str(j)
                   plot_samples(pulses,out_dir=self.brainsight_dir, map_surface='pial', prefix=prefix,annot='aparc', annot_alpha=1,show_average=True)
                   html.add_paragraph('&nbsp&nbsp&nbsp -> Target: '+ t + ' (start:' + str(seq[0])+ ', end:' + str(seq[1]) + ', stimulations:' + str(seq[1]-seq[0]+1)+')<br><br>')
                   html.add_image(prefix+'.png',200,900,'middle')
                   html.add_paragraph('<br><br>')
               
               else:
                   html.add_paragraph('&nbsp&nbsp&nbsp -> Target: '+ t + ' (start:' + str(seq[0])+ ', end:' + str(seq[1]) + ', stimulations:' + str(seq[1]-seq[0]+1)+')<br><br>')
               
   
           html.write('{subject_dir}/brainsight/summary.html'.format(subject_dir=self.subject_dir))
        

        
            
            

    
    
        
    
        
        
