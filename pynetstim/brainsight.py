"""
miscellaneous functions and classes for pynetstim projects
"""

import pandas as pd
import numpy as np
from collections import defaultdict


class BrainsightSessionFile(object):
    
    """ reading and parsing brainsight session output file, currently outputs: targets, samples, electrodes, planned landmarks, and session landmarks"""
    
    def __init__(self,bs_session_file, out_dir):
    
        self.bs_session_file = bs_session_file
        self.out_dir = out_dir
        
        
        self.f = open(self.bs_session_file,'r').read().split('\n')
        self.tables = [(i,header) for i,header in enumerate(self.f) if '# ' in header]  #### finding where sub-tables start
        self.tables = self.tables + [(len(self.f),'# End')]
        
        # target table 
        self._parse_table('Target','targets.txt',exclusion='Sample')
    
        # Samples
        self._parse_table('Sample','samples.txt')
        
        # Electrodes
        self._parse_table('Electrode','electrodes.txt')
        
        # planned landmarks
        self._parse_table('Planned Landmark','planned_landmarks.txt')
        
        # session landmark
        self._parse_table('Session Landmark','session_landmarks.txt')
        
    
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
                

class Targets(object):
    
    def __init__(self, targets_file):
        self.targets_file = targets_file 
        self.targets_df = pd.read_table(targets_file)
        self.targets_names = np.unique(self.targets_df.target_name)
        
    def get_targets_coords(self,targets):
      
        return self.targets_df[self.targets_df.target_name.apply(lambda x: x in targets)][['loc_x','loc_y','loc_z']].values
        
        
class Samples(object):
    
    def __init__(self,samples_file, print_summary=True):
        
        self.samples_file = samples_file
        self.samples_df = pd.read_table(samples_file)
        
        self.sessions_names = np.unique(self.samples_df.session_name)
        self.sessions = defaultdict(dict)
        
        for session in self.sessions_names:
            session_samples, stims_sequence, targets_stims = self._load_session(session)
            self.sessions[session]['samples'] = session_samples
            self.sessions[session]['stims_sequence'] = stims_sequence
            self.sessions[session]['targets_stims'] = targets_stims
        
        if print_summary:
            self.summary()
            
    def _load_session(self,session):
        
        session_samples = self.samples_df[(self.samples_df.session_name==session)]
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
            
        return session_samples, stims_sequence, targets_stims
        
        
            
    def get_target(self,target,session=None, chunk=None):
        if not session: 
            return self.samples[self.samples.assoc_target==target]
        
        else:
            session_samples = self.sessions[session]['samples']
            if not chunk:
                return session_samples[session_samples.assoc_target==target]
            else:
                try:
                    start,end = self.sessions[session]['targets_stims'][target][chunk-1]
                    return session_samples[session]['samples'].iloc[start:end]
                except:
                    raise ValueError('chunk exceeds the range')
                
    def summary(self):
        
        """ write summary of samples """
        
        print 'Number of sessions: ', len(self.sessions_names)
        print '======================'
        for session in self.sessions_names:
            message=session+'\n' + '------------\n'
            targets = np.unique(self.sessions[session]['targets_stims'].keys())
            message+= 'Session targets: ' + ','.join(targets)+'\n\n'
            sequences = self.sessions[session]['stims_sequence'].keys()
            sequences = sorted(sequences)
            for s in sequences:
                t = self.sessions[session]['stims_sequence'][s]
                message+= '\t-> Target: '+ t + ' (start:' + str(s[0])+ ', end:' + str(s[1]) + ', stimulations:' + str(s[1]-s[0]+1)+')\n\n'
            print message

        
            
        
    
    


        
        
        
    
    
     