"""
miscellaneous functions and classes to work with Brainsight output
Author: Ehsan Tadayon, M.D. [sunny.tadayon@gmail.com / stadayon@bidmc.harvard.edu]
"""

import pandas as pd
import numpy as np
from collections import defaultdict
import os
from coordinates import Coords, FreesurferCoords
from datetime import datetime
from pymisc.plotting import clean_plot
import matplotlib.pyplot as plt


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
            print name
            if name=='Target':
                self._parse_table(name,self.base_name+name+'.txt', exclusion='Sample')
            else:
                self._parse_table(name,self.base_name+name+'.txt')
        
        # target table 
        #self._parse_table('Target',self.base_name+'targets.txt',exclusion='Sample')
    
        # Samples
        #self._parse_table('Sample',self.base_name+'samples.txt')
        
        # Electrodes
        #self._parse_table('Electrode',self.base_name+'electrodes.txt')
        
        # planned landmarks
        #self._parse_table('Planned Landmark',self.base_name+'planned_landmarks.txt')
        
        # session landmark
        #self._parse_table('Session Landmark',self.base_name+'session_landmarks.txt')
        
    
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
            s = np.unique(self._df['session_name'])
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



class BrainsightElectrodes(object):
    
    def __init__(self,electrodes_file):
        self.electrodes_file = electrodes_file
        self._df = pd.read_table(electrodes_file,na_values='(null)')
        
    def get_electrodes(self,session,exclude_null=True):
        df2 = self._df.copy()
        df2 = df2[(df2.electrode_type=='EEG')&(df2.session_name==session)]
        df2.dropna(axis=0,inplace=True)
        return df2



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
    print idx

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
    


            
        
    
    


        
        
        
    
    
     