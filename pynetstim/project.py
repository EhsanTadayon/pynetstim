"""
Main stimulation project functions and classes
"""

from utils import Coords, FreesurferCoords, Surf, FreesurferSurf
from brainsight import BrainsightSessionFile, BrainsightSamples, BrainsightTargets
import os
from pymisc.htmlreport import HtmlDoc
from mayavi import mlab
from surfer import Brain
from nipype.interfaces.fsl import FLIRT,FNIRT


## TO DO:
#### define a Target class ( and brainsight targets will inherit from it) 
### define a Sample class ( and BrainsightSamples will inherit from it)
#### register the subject to MNI using flirt and fnirt. 


class StimProject(object):
    
    """ main stimulation project class"""
    def __init__(self, subject, freesurfer_dir, project_dir, brainsight_file = None, to_mni=True,
     mni_template='MNI152_T1_2mm.nii.gz', mni_directory = os.path.join(os.environ['FSLDIR'],'data/standard')):
        
        self.subject = subject

        if not os.path.isabs(project_dir):
            raise ValueError('please provide absolute path for project_dir')
        if not os.path.isabs(freesurfer_dir):
            raise ValueError('please provide absolute path for freesurfer_dir')
    
        self.project_dir = project_dir
        self.subject_dir = os.path.join(self.project_dir,subject)
        self.freesurfer_dir = freesurfer_dir
        self.brainsight_file = brainsight_file
        self.targets = []
    
        self._start_project()
        self._to_mni()
        self._read_brainsight_file()
        
    def _start_project(self):
        
        # TO DO: add more folders
        if not os.path.exists(self.subject_dir):
            
            os.makedirs('{subject_dir}/brainsight'.format(subject_dir=self.subject_dir))
            os.makedirs('{subject_dir}/figures'.format(subject_dir = self.subject_dir))
            os.makedirs('{subject_dir}/mni'.format(subject_dir = self.subject_dir))
        
        self.brainsight_dir='{subject_dir}/brainsight'.format(subject_dir=self.subject_dir)
        self.figures_dir = '{subject_dir}/figures'.format(subject_dir=self.subject_dir)  
        self.mni_nonlinear = '{subject_dir}/mni'.format(subject_dir=self.subject_dir)


    def _to_mni(self):
        pass ## TO DO
        
    def _read_brainsight_file(self):
        
        if self.brainsight_file is not None:
            bs = BrainsightSessionFile(self.brainsight_file,out_dir='{subject_dir}/brainsight'.format(subject_dir=self.subject_dir))
            self.brainsight_samples = BrainsightSamples('{subject_dir}/brainsight/samples.txt'.format(subject_dir=self.subject_dir))
            self.brainsight_targets = BrainsightTargets('{subject_dir}/brainsight/targets.txt'.format(subject_dir=self.subject_dir))
            

    def summary(self):
       
       if self.brainsight_file:
           self._brainsight_summary()

    def _brainsight_summary(self):
   
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
        
           for seq in stims_sequences_sorted:
               t = stims_sequences[seq]
               html.add_paragraph('&nbsp&nbsp&nbsp -> Target: '+ t + ' (start:' + str(seq[0])+ ', end:' + str(seq[1]) + ', stimulations:' + str(seq[1]-seq[0]+1)+')<br><br>')
   
           html.write('{subject_dir}/brainsight/summary.html'.format(subject_dir=self.subject_dir))
        

        
            
            

    
    
        
    
        
        
