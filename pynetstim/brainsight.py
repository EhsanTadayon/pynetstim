"""
miscellaneous functions and classes to work with Brainsight output
Author: Ehsan Tadayon, M.D. [sunny.tadayon@gmail.com / stadayon@bidmc.harvard.edu]
"""

import pandas as pd
import numpy as np
from collections import defaultdict
import os
from .coordinates import Coords, FreesurferCoords
from datetime import datetime
import matplotlib.pyplot as plt
from .freesurfer_files import Surf, FreesurferSurf, Annot
from nipype.interfaces.fsl import FLIRT,FNIRT
from .plotting import plotting_points_fast,plotting_points
from scipy.spatial.distance import cdist
from .utils import make_head_model, HtmlDoc, clean_plot
import shutil


##### Classes in this module: 
##### ---- 1. BrainsightSessionFile
##### ---- 2. BrainsightTargets
##### -----3. BrainsightElectrodes
##### -----4. BrainsightProject

##### functions to use with BrainsightSamples: 
##### ---- 1.chunk_samples
##### ---- 2.plot_chunks


###################################################################################################    
##                                      BrainsightSessionFile
###################################################################################################     
class BrainsightSessionFile(object):
    
    """ reading and parsing brainsight session output file, currently outputs: targets, samples, electrodes, planned landmarks, and session landmarks"""
    
    def __init__(self, bs_session_file, base_name='', out_dir='.'):
    
        self.bs_session_file = bs_session_file
        self.base_name = base_name
        self.out_dir = os.path.abspath(out_dir)
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
         
        
        
        self.f = open(self.bs_session_file,'r').read().split('\n')
        self.tables = [(i,header) for i,header in enumerate(self.f) if '# ' in header]  #### finding where sub-tables start
        self.tables = self.tables + [(len(self.f),'# End')]
        self.tables = self.tables[6:]
        self.tables_names = [x[1].split(' ')[1] for x in self.tables][0:-1]
        
        for name in self.tables_names:
            print(name)
            if name=='Target':
                self._parse_table(name,self.base_name+name+'.txt', exclusion=None)
            else:
                self._parse_table(name,self.base_name+name+'.txt')
        
    
    def _find_start_and_end(self,tables,table):
        
        for i,t in enumerate(tables):
            if '# {table}'.format(table=table) in t[1]:
                s = t[0]
                e = tables[i+1][0]
                return s,e
                
    def _parse_table(self,table,outname,exclusion=None):

        s,e = self._find_start_and_end(self.tables,table)
        header = self._clean_header(self.f[s])
        rows = [header] + self.f[s+1:e]
        if exclusion:
            rows = [x for x in rows if exclusion not in x]
        
        g = open('{out_dir}/{outname}'.format(out_dir=self.out_dir,outname=outname),'w')
        g.write('\n'.join(rows))
        g.close()
        
        
    def _clean_header(self,header):
        
        header = header.replace('.','')
        header = header.replace('# ','')
        header = header.replace(' ','_')
        header = header.lower()
        return header
                


###################################################################################################    
##                                      BrainsightTargets
###################################################################################################     

class BrainsightTargets(object):
    
    def __init__(self, subject, project_dir, anat_img=None, freesurfer_dir=None):
        
        self.project_dir = project_dir
        self.subject = subject
        self.targets_file = os.path.join(project_dir,subject,'brainsight','Target.txt')
        
        self.freesurfer_dir = freesurfer_dir
        self.anat_img = anat_img
        self._df = pd.read_table(self.targets_file)
        self.working_dir = os.path.join(project_dir,subject)
        
            
    def get_coord(self):
        return self._df[['loc_x','loc_y','loc_z']].values
        
    def get_direction(self):
        return self._df[['m0n0','m0n1','m0n2','m1n0','m1n1','m1n2','m2n0','m2n1','m2n2']].values    
        
    def get_name(self):
        return self._df['target_name'].values
        
    def get_table(self):
        return self._df.copy()
        
    def to_freesurfer_coords(self):
        
        if self.freesurfer_dir is None:
            raise('freesurfer_dir should be provided!')
            
        fscoords = FreesurferCoords(self.get_coord(), subject=self.subject, freesurfer_dir=self.freesurfer_dir,
                                                    name=self.get_name(), direction=self.get_direction(), working_dir=self.working_dir)
        return fscoords
        
    def to_coords(self):
        
        coords = Coords(self.get_coord(), img_file=self.anat_img, subject=self.subject, name=self.get_name(), direction=self.get_direction(),working_dir=self.working_dir)
        return coords
         
    
        

###################################################################################################    
##                                      BrainsightSamples
###################################################################################################        
class BrainsightSamples(object):
    
    def __init__(self, samples_file, rename_sessions=True):
        
        self.samples_file = samples_file
        self._dirname = os.path.dirname(self.samples_file)
        self._load_dataframe(rename_sessions)
        self.sessions_names = np.unique(self._df.session_name)
        self._sessions = defaultdict(dict)
        
        
        for session in self.sessions_names:
            session_samples, stims_sequence, targets_stims, session_targets = self._load_session(session)
            self._sessions[session]['samples'] = session_samples
            self._sessions[session]['stims_sequence'] = stims_sequence
            self._sessions[session]['targets_stims'] = targets_stims
            self._sessions[session]['targets'] = session_targets
        
            
    def _load_dataframe(self, rename_sessions):
        self._df = pd.read_table(self.samples_file)
        numeric_columns = [u'loc_x',
       u'loc_y', u'loc_z', u'm0n0', u'm0n1', u'm0n2', u'm1n0', u'm1n1',
       u'm1n2', u'm2n0', u'm2n1', u'm2n2', u'dist_to_target', u'target_error',
       u'angular_error', u'twist_error']
        self._df[numeric_columns] = self._df[numeric_columns].apply(pd.to_numeric,errors='coerce')
        
        ## renaming sessions names
        if rename_sessions:
            temp = self._df['session_name'].values.tolist()
            sidx = np.unique(self._df['session_name'],return_index=True)[1]
            s = [temp[idx] for idx in sorted(sidx)]
            ss = {x:'Session {i}'.format(i=i+1) for i,x in enumerate(s)} ### renaming the sessions to Session 1, Session 2, ... 
            def change_name(x):
                return ss[x]
            
            self._df['session_name'] = self._df['session_name'].apply(change_name)

            
    def _load_session(self,session):
        
        session_samples = self._df[(self._df.session_name==session)].copy()
        session_samples.index = np.arange(session_samples.shape[0])
        session_targets = np.unique(session_samples.assoc_target)
        session_targets = [target for target in session_targets if 'Sample' not in target]
        stims_sequence = {}
        targets_stims = {}
        
        for target in session_targets:
            s = np.diff(np.double(session_samples.assoc_target==target))
            starts = np.where(s==1)[0] + 1
            starts = starts.tolist()
            ends = np.where(s==-1)[0]
            ends = ends.tolist()
            if len(starts)>len(ends):
                ends = ends+[session_samples.shape[0]]
            if len(ends)>len(starts):
                starts = [0] + starts
                
            for start,end in zip(starts,ends):
                stims_sequence[(start,end)] = target
                
            targets_stims[target] = zip(starts,ends)
            
        return session_samples, stims_sequence, targets_stims, session_targets
        
    
    def get_samples(self,session=None):
        if session is None:
            return self._df.copy()
            
        else:
            return self._sessions[session]['samples']
        
    def get_stims_sequence(self,session):
        return self._sessions[session]['stims_sequence']
        
    def get_session_targets(self, session):
        return self._sessions[session]['targets']
                
    def get_target_stims(self,target,session=None, chunk=None):
        
        if not session: 
            return self._df[self._df.assoc_target==target]
        
        else:
            session_samples = self.get_samples(session)
            
            if chunk is None:
                return session_samples[session_samples.assoc_target==target]
                
            else:
                #try:
                start,end = self._sessions[session]['targets_stims'][target][chunk]
                return session_samples.iloc[start:end]
                #except:
                    #raise ValueError('chunk exceeds the range')
                    
                    
###################################################################################################    
##                                      Functions to use with BriansightSamples
###################################################################################################     
def chunk_samples(samples_df, thr=50):
    samples_df = samples_df.copy()
  
    ### reading time 
    fmt = '%H:%M:%S.%f'
    t = [datetime.strptime(x,fmt) for x in samples_df.time]
    d = map(lambda i: t[i] - t[i-1], np.arange(1,len(t)))
    d = [x.seconds for x in d] # to seconds
    d = [0] + d # to account for the first sample
    d = np.array(d)
    idx = np.where(d>thr)[0].tolist()

    idx = [0] + idx + [len(d)-1]
    print(idx)

    ## chunks
    chunks = {}
    for i in range(1,len(idx[1:])+1):
        chunks[i] = samples_df.iloc[idx[i-1]:idx[i]]
    return chunks
        
    
def plot_chunks(chunks,col='target_error',figsize=(8,6)):
    
    fig,ax = plt.subplots(figsize=figsize)
    for chunk in chunks:
        ax.plot(chunks[chunk].index, chunks[chunk][col].values,'*-')
    
    ymax = ax.get_ylim()[1]
    for chunk in chunks:
        ax.fill_between(chunks[chunk].index,0,ymax,alpha=.2)

    ax.set_ylabel(col)
    ax.set_xlabel('pulse number')
    
    #ylim
    ax.set_ylim(0,ymax)
        
    fig,ax = clean_plot(fig,ax)
    return fig,ax
    


###################################################################################################    
##                                      BrainsightElectrodes
###################################################################################################     
class BrainsightElectrodes(object):
    
    def __init__(self,subject, project_dir, anat_img=None, freesurfer_dir=None):
        
        self.subject = subject
        self.project_dir = project_dir
        self.electrodes_file = os.path.join(project_dir,subject,'brainsight','Electrode.txt')
        
        self.anat_img = anat_img
        self.freesurfer_dir = freesurfer_dir
        self.working_dir = os.path.join(project_dir,subject)
        
        
        self._df = pd.read_table(self.electrodes_file,na_values='(null)')
        cols = ['electrode_name', 'electrode_type', 'session_name', 'loc_x', 'loc_y',
         'loc_z', 'm0n0', 'm0n1', 'm0n2', 'm1n0', 'm1n1', 'm1n2', 'm2n0', 'm2n1', 'm2n2']
        self._df = self._df[cols]
      
        
    def _get_session_df(self, session, exclude_null=True):
        df2 = self._df.copy()
        df2 = df2[(df2.electrode_type=='EEG')&(df2.session_name==session)]
        if exclude_null:
            df2.dropna(axis=0,inplace=True)
        return df2
        
    def get_name(self,session=None,exclude_null=True):
        if session:
            df = self._get_session_df(session,exclude_null)
        else:
            df = self._df
        return df['electrode_name'].values
        
    def get_coord(self,session=None, exclude_null=True):
        if session: 
            df = self._get_session_df(session, exclude_null)
        else:
            df = self._df
        return df[['loc_x','loc_y','loc_z']].values
    
    def get_direction(self,session=None,exclude_null=True):
        if session:
            df = self._get_session_df(session,exclude_null)
        else:
            df = self._df
            
        return df[['m0n0','m0n1','m0n2','m1n0','m1n1','m1n2','m2n0','m2n1','m2n2']].values
        
    def to_freesurfer_coords(self,session=None):

        if self.freesurfer_dir is None:
            raise('both subject and freesurfer_dir should be provided!')
            
        fscoords = FreesurferCoords(self.get_coord(session), subject=self.subject, freesurfer_dir=self.freesurfer_dir,
                                                    name=self.get_name(session), direction=self.get_direction(session), working_dir=self.working_dir)
        return fscoords 
        
    def to_coords(self):
        coords = Coords(self.get_coord(), img_file=self.anat_img, subject=self.subject, name=self.get_name(), direction=self.get_direction(),working_dir=self.working_dir)
        return coords   
        



###################################################################################################    
##                                      BrainsightProject
###################################################################################################     

class BrainsightProject(object):
    
    """ main stimulation project class"""
    def __init__(self, subject, brainsight_file, project_dir, anat_img=None, freesurfer_dir = None):
        
        self.subject = subject
        self.project_dir = project_dir  
        self.brainsight_file = brainsight_file
        self.anat_img = anat_img
        self.subject_dir = os.path.join(self.project_dir,subject)
        self.freesurfer_dir = freesurfer_dir
        
        ## initializing
        self._parse_brainsight_file()
         
    def _parse_brainsight_file(self):
        
        if not os.path.exists(os.path.join(self.subject_dir,'brainsight')):
            os.makedirs(os.path.join(self.subject_dir,'brainsight'))
            
        bs = BrainsightSessionFile(self.brainsight_file,out_dir='{subject_dir}/brainsight'.format(subject_dir=self.subject_dir))
        self.brainsight_samples = BrainsightSamples('{subject_dir}/brainsight/Sample.txt'.format(subject_dir=self.subject_dir))
        self.brainsight_targets = BrainsightTargets(self.subject, self.project_dir, anat_img = self.anat_img, freesurfer_dir=self.freesurfer_dir)
        
        if 'Electrode' in bs.tables_names:
            self.brainsight_electrodes = BrainsightElectrodes(subject=self.subject, project_dir=self.project_dir, anat_img = self.anat_img, freesurfer_dir=self.freesurfer_dir)

    def summary(self, plot_pulses=False, remove_outliers=True, overwrite=False, heightpx=200, widthpx=800):
        
        """ write summary of samples """
        if os.path.exists(os.path.join(self.subject_dir,'brainsight','summary.html'))==False or overwrite==True:
            
            html = HtmlDoc('Brainsight sessions summary') 
            html.add_header('h1','Brainsight sessions summary')
            html.add_header('h3','Number of session: %d'%len(self.brainsight_samples.sessions_names))

            for session in self.brainsight_samples.sessions_names:

                html.add_header('h2',session)
                session_targets = self.brainsight_samples.get_session_targets(session)
                html.add_paragraph('<b>Targets:</b> '+', '.join(session_targets)+'<br><br>')

                stims_sequences = self.brainsight_samples.get_stims_sequence(session)
                stims_sequences_sorted = sorted(stims_sequences)

                for j,seq in enumerate(stims_sequences_sorted):
                   t = stims_sequences[seq]
                   html.add_paragraph('&nbsp&nbsp&nbsp -> Target: '+ t + ' (start:' + str(seq[0])+ ', end:' + str(seq[1]) + ', stimulations:' + str(seq[1]-seq[0]+1)+')<br><br>')

                   if plot_pulses:
               
                       if self.freesurfer_dir is None:
                           raise('freesurfer_dir should be provided to plot pulses!')
           
                       ### create head models
                       anat_img = os.path.join(self.freesurfer_dir,'mri','rawavg.mgz')
                       out_dir = os.path.join(self.freesurfer_dir,self.subject, 'bem')
                       make_head_model(anat_img,out_dir)
               
                       ## pulses
                       pulses_= self.brainsight_samples.get_samples(session).iloc[seq[0]:seq[1]]
                       pulses_coords = pulses_[['loc_x','loc_y','loc_z']].values  #pulses ras coords

                       ## remove outlier pulses 
                       if remove_outliers:
                           avg_coord = np.mean(pulses_coords,axis=0)
                           d = cdist(np.atleast_2d(avg_coord),pulses_coords)
                           idx = d<=np.mean(d)+2*np.std(d)
                           idx = idx.flatten()
                           pulses_coords = pulses_coords[idx]

                       ## plotting pulses   
                       pulses = FreesurferCoords(pulses_coords, self.subject, self.freesurfer_dir)
                       #
                       prefix=session+'_'+t+'_'+str(j)
                       img_out_dir = os.path.join(self.subject_dir,'brainsight', 'images')
                       p = plotting_points_fast(pulses, show_average=True)
                       img_path = p.save_image(out_dir = img_out_dir, prefix=prefix)
                       html.add_image(img_path,heightpx,widthpx,'middle')
                       html.add_paragraph('<br><br>')

                       html.write('{subject_dir}/brainsight/summary.html'.format(subject_dir=self.subject_dir))
        






    