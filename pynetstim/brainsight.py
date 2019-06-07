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
from pymisc.plotting import clean_plot
import matplotlib.pyplot as plt
from .surface import Surf, FreesurferSurf
from pymisc.htmlreport import HtmlDoc
from nipype.interfaces.fsl import FLIRT,FNIRT
from .plotting import plotting_points_fast
from scipy.spatial.distance import cdist
from .utils import make_head_model
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
    
    def __init__(self,bs_session_file,base_name='', out_dir='.'):
    
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

class BrainsightTargets(FreesurferCoords):
    
    def __init__(self, targets_file, subject, freesurfer_dir):
        
        self.targets_file = targets_file 
        self._df = pd.read_table(targets_file)
            
        coords = self._df[['loc_x','loc_y','loc_z']].values
        guess_hemi=True
        name = self._df['target_name'].values
        direction = self._df[['m0n0','m0n1','m0n2','m1n0','m1n1','m1n2','m2n0','m2n1','m2n2']].values
        FreesurferCoords.__init__(self,coords, subject, freesurfer_dir, name=name, direction=direction)
        
    def get_coords(self,targets=None):
        if targets is None:
            targets = self.name
            
        return self._df[self._df.target_name.apply(lambda x: x in targets)][['loc_x','loc_y','loc_z']].values
        
    def get_table(self):
        return self._df.copy()
    
    

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
        
    def get_targets(self, session):
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
    
    def __init__(self,electrodes_file):
        self.electrodes_file = electrodes_file
        self._df = pd.read_table(electrodes_file,na_values='(null)')
        
    def get_electrodes(self,session,exclude_null=True):
        df2 = self._df.copy()
        df2 = df2[(df2.electrode_type=='EEG')&(df2.session_name==session)]
        df2.dropna(axis=0,inplace=True)
        return df2




###################################################################################################    
##                                      BrainsightProject
###################################################################################################     

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
            
    
    
     