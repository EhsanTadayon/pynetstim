"""
Main stimulation project functions and classes
Author: Ehsan Tadayon, M.D. [sunny.tadayon@gmail.com / stadayon@bidmc.harvard.edu]
"""

from coordinates import Coords, FreesurferCoords
from surface import Surf, FreesurferSurf
from brainsight import BrainsightSessionFile, BrainsightSamples, BrainsightTargets, BrainsightElectrodes
import os
from pymisc.htmlreport import HtmlDoc
from nipype.interfaces.fsl import FLIRT,FNIRT
from plotting import plotting_points_fast
import numpy as np
from scipy.spatial.distance import cdist


## TO DO:



class StimProject(object):
    
    """ main stimulation project class"""
    def __init__(self, subject, project_dir, anatomical = None, freesurfer_dir = None, simnibs_dir=None, brainsight_file = None, to_mni = False,
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
        self.simnibs_dir = simnibs_dir
        
        if self.simnibs_dir is not None:
            self.freesurfer_dir = self.simnibs_dir
            
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
            os.makedirs('{subject_dir}/')
        
        self.brainsight_dir='{subject_dir}/brainsight'.format(subject_dir=self.subject_dir)
        self.figures_dir = '{subject_dir}/figures'.format(subject_dir=self.subject_dir)  
        self.mni_nonlinear = '{subject_dir}/mni'.format(subject_dir=self.subject_dir)
        
        ### create head models
        self._make_head_model()
        
        
    def _make_head_model(self):
        
        if not os.path.isfile('{freesurfer_dir}/{subject}/bem/outer_skin_surface'.format(freesurfer_dir=self.freesurfer_dir, subject=self.subject)):
            if not os.path.exists('{freesurfer_dir}/{subject}/bem'.format(freesurfer_dir=self.freesurfer_dir, subject=self.subject)):
                os.makedirs('{freesurfer_dir}/{subject}/bem'.format(freesurfer_dir=self.freesurfer_dir, subject=self.subject))
        
            cmd ='cd {freesurfer_dir}/{subject}/bem; mri_watershed -surf surf {freesurfer_dir}/{subject}/mri/rawavg.mgz brain.mgz'.format(freesurfer_dir=self.freesurfer_dir, subject = self.subject)
            os.system(cmd)
        
            for f in ['lh.surf_brain_surface','lh.surf_inner_skull_surface','lh.surf_outer_skin_surface','lh.surf_outer_skull_surface']:
                cmd = 'mv {freesurfer_dir}/{subject}/bem/{f} {freesurfer_dir}/{subject}/bem/{f2}'.format(f=f,f2=f.split('lh.surf_')[1],freesurfer_dir=self.freesurfer_dir,subject=self.subject)
                os.system(cmd)
        

    def _to_mni(self,to_mni,mni_template):
        pass #### TO DO    
        
        
    def _read_brainsight_file(self):
        
        if self.brainsight_file is not None:
            bs = BrainsightSessionFile(self.brainsight_file,out_dir='{subject_dir}/brainsight'.format(subject_dir=self.subject_dir))
            self.brainsight_samples = BrainsightSamples('{subject_dir}/brainsight/Sample.txt'.format(subject_dir=self.subject_dir))
            self.brainsight_targets = BrainsightTargets('{subject_dir}/brainsight/Target.txt'.format(subject_dir=self.subject_dir), self.subject, self.freesurfer_dir)
            if 'Electrode' in bs.tables_names:
                self.brainsight_electrodes = BrainsightElectrodes('{subject_dir}/brainsight/Electrode.txt'.format(subject_dir=self.subject_dir))
            

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
                   plot_points_fast(pulses, map_surface='pial', annot='aparc', annot_alpha=1,show_average=True, out_dir = os.path.join(self.brainsight_dir, 'images'), prefix=prefix)
                   html.add_paragraph('&nbsp&nbsp&nbsp -> Target: '+ t + ' (start:' + str(seq[0])+ ', end:' + str(seq[1]) + ', stimulations:' + str(seq[1]-seq[0]+1)+')<br><br>')
                   html.add_image('./images/'+prefix+'.png',200,900,'middle')
                   html.add_paragraph('<br><br>')
               
               else:
                   html.add_paragraph('&nbsp&nbsp&nbsp -> Target: '+ t + ' (start:' + str(seq[0])+ ', end:' + str(seq[1]) + ', stimulations:' + str(seq[1]-seq[0]+1)+')<br><br>')
               
   
           html.write('{subject_dir}/brainsight/summary.html'.format(subject_dir=self.subject_dir))
        

        
            
            

    
    
        
    
        
        
