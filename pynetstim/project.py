"""
Main stimulation project functions and classes
Author: Ehsan Tadayon, M.D. [sunny.tadayon@gmail.com / stadayon@bidmc.harvard.edu]
"""

from .coordinates import Coords, FreesurferCoords
from .surface import Surf, FreesurferSurf
from .brainsight import BrainsightSessionFile, BrainsightSamples, BrainsightTargets, BrainsightElectrodes
import os
from pymisc.htmlreport import HtmlDoc
from nipype.interfaces.fsl import FLIRT,FNIRT
from .plotting import plotting_points_fast
import numpy as np
from scipy.spatial.distance import cdist
from .utils import make_head_model
import shutil


## TO DO:


class BrainsightProject(object):
    
    """ main stimulation project class"""
    def __init__(self, subject, project_dir, brainsight_file, freesurfer_dir = None, simnibs_dir = None,to_mni = False,
     mni_template='MNI152_T1_2mm.nii.gz', mni_dir = os.path.join(os.environ['FSLDIR'],'data/standard'), overwrite=None):
        
        self.subject = subject

        self.simnibs_dir = simnibs_dir
        if self.simnibs_dir is not None:
            self.freesurfer_dir = self.simnibs_dir
        else:
            self.freesurfer_dir = freesurfer_dir
            
        self.project_dir = project_dir  
        self.brainsight_file = brainsight_file    
        self.subject_dir = os.path.join(self.project_dir,subject)
        self._to_mni = to_mni
        self._mni_template = mni_template
        self._mni_dir = mni_dir 
        self.overwrite = overwrite
        self.directories = {}
        
        
        self._start_project()
        
    def _start_project(self):
        
        self._start_log()
            
        ## read brainsight file
        self._read_brainsight_file()
        
        ## 
        if self._to_mni:
            self._to_mni(self.mni_template, self.mni_dir)
                
        ### create head models
        self._make_head_model()    
        
    def _start_log(self):
        
        try:
            previous_logs = open('{subject_dir}/logs/log.txt'.format(subject_dir = self.subject_dir),'r').read().split('\n')[0:-1]
            
            if self.overwrite is None:
                self.log = previous_logs
                
            else:
                if self.overwrite=='all':
                    self.log = []
                    f = open('{subject_dir}/logs/log.txt'.format(subject_dir = self.subject_dir),'w')
                    f.close()
                else:
                    for x in self.overwrite:
                        previous_logs.remove(x)
                    f = open('{subject_dir}/logs/log.txt'.format(subject_dir = self.subject_dir),'w')
                    for p in previous_logs:
                        f.write(p+'\n')
                    f.close()
                    self.log = previous_logs    
                    
        except FileNotFoundError:
            os.makedirs('{subject_dir}/logs'.format(subject_dir = self.subject_dir))
            f = open('{subject_dir}/logs/log.txt'.format(subject_dir = self.subject_dir),'w')
            f.close()
            f.close()
            self.log = []
            
                        
    def _add_to_log(self,name):
        
        f = open('{subject_dir}/logs/log.txt'.format(subject_dir = self.subject_dir),'a')
        f.write(name+'\n')
        f.close()
        self.log.append(name)
        
        
    def _add_dir(self,name):
        
        try:
           os.makedirs('{subject_dir}/{name}'.format(subject_dir=self.subject_dir,name=name))
           
        except:
           shutil.rmtree('{subject_dir}/{name}'.format(subject_dir=self.subject_dir,name=name),ignore_errors=True)
           os.makedirs('{subject_dir}/{name}'.format(subject_dir=self.subject_dir,name=name))
        
        self.directories[name] = '{subject_dir}/{name}'.format(subject_dir = self.subject_dir, name=name)
            
        
    def _make_head_model(self):
        
        if 'make_head_model' not in self.log:
            make_head_model(self.subject,self.freesurfer_dir)
            self._add_to_log('make_head_model')
        
    def _to_mni(self,mni_template):
        if 'to_mni' not in self.log:
            self._add_dir('mni')
            pass
            ### TO DO
            
            self._add_to_log('to_mni')
            
    def _read_brainsight_file(self):
        
        if 'read_brainsight_file' not in self.log:
            self._add_dir('brainsight')
            self._add_to_log('read_brainsight_file')
            
        bs = BrainsightSessionFile(self.brainsight_file,out_dir='{subject_dir}/brainsight'.format(subject_dir=self.subject_dir))
        self.brainsight_samples = BrainsightSamples('{subject_dir}/brainsight/Sample.txt'.format(subject_dir=self.subject_dir))
        self.brainsight_targets = BrainsightTargets('{subject_dir}/brainsight/Target.txt'.format(subject_dir=self.subject_dir), self.subject, self.freesurfer_dir)
        
        if 'Electrode' in bs.tables_names:
            self.brainsight_electrodes = BrainsightElectrodes('{subject_dir}/brainsight/Electrode.txt'.format(subject_dir=self.subject_dir))


    def summary(self,plot_pulses=False):
       
       if 'summary' not in self.log:
            self._brainsight_summary(plot_pulses) 
            self._add_to_log('summary')

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
                   plotting_points_fast(pulses, map_surface='pial', annot='aparc',show_average=True, out_dir = os.path.join(self.directories['brainsight'], 'images'), prefix=prefix)
                   html.add_paragraph('&nbsp&nbsp&nbsp -> Target: '+ t + ' (start:' + str(seq[0])+ ', end:' + str(seq[1]) + ', stimulations:' + str(seq[1]-seq[0]+1)+')<br><br>')
                   html.add_image('./images/'+prefix+'.png',200,900,'middle')
                   html.add_paragraph('<br><br>')
               
               else:
                   html.add_paragraph('&nbsp&nbsp&nbsp -> Target: '+ t + ' (start:' + str(seq[0])+ ', end:' + str(seq[1]) + ', stimulations:' + str(seq[1]-seq[0]+1)+')<br><br>')
               
   
           html.write('{subject_dir}/brainsight/summary.html'.format(subject_dir=self.subject_dir))
        

        
            

        
                    

    
    
        
    
        
        
