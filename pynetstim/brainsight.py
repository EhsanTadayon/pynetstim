"""
miscellaneous functions and classes for pynetstim projects
"""

import pandas as pd



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
                
    def get_table(self,table):
        df = pd.read_table('{out_dir}/{table}.txt'.format(out_dir=self.out_dir, table = table))
        return df

        
        
        
    
    
     