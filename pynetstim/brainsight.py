"""
miscellaneous functions and classes for pynetstim projects
"""


class BrainsightSessionFile(object):
    
    """ reading and parsing brainsight session output file, currently outputs: targets, samples, electrodes, planned landmarks, and session landmarks"""
    
    def __init__(self,bs_session_file, out_dir):
    
        self.bs_session_file = bs_session_file
        self.out_dir = out_dir
        
        
        f = open(self.bs_session_file,'r').read().split('\n')
        tables = [(i,header) for i,header in enumerate(f) if '# ' in header]  #### finding where sub-tables start
        tables = tables + [(len(f),'# End')]
        
        # target table 
        s,e = self._find_start_and_end(tables,'Target')
        targets = f[s:e]
        targets = [x for x in targets if 'Sample' not in x]
        g = open('{out_dir}/Targets.txt'.format(out_dir=out_dir),'w')
        g.write('\n'.join(targets))
        g.close()
        
        # Samples
        s,e = self._find_start_and_end(tables,'Sample')
        samples = f[s:e]
        g = open('{out_dir}/Samples.txt'.format(out_dir=out_dir),'w')
        g.write('\n'.join(samples))
        g.close()
        
        # Electrodes
        s,e = self._find_start_and_end(tables,'Electrode')
        electrodes = f[s:e]
        g = open('{out_dir}/Electrodes.txt'.format(out_dir=out_dir),'w')
        g.write('\n'.join(electrodes))
        g.close()
        
        # planned landmarks
        s,e = self._find_start_and_end(tables,'Planned Landmark')
        planned_landmarks = f[s:e]
        g = open('{out_dir}/PlannedLandmarks.txt'.format(out_dir=out_dir),'w')
        g.write('\n'.join(planned_landmarks))
        g.close()
        
        # session landmark
        s,e = self._find_start_and_end(tables,'Session Landmark')
        session_landmarks = f[s:e]
        g = open('{out_dir}/SessionLandmarks.txt'.format(out_dir=out_dir),'w')
        g.write('\n'.join(session_landmarks))
        g.close()
        
    
    def _find_start_and_end(self,tables,table):
        
        for i,t in enumerate(tables):
            if '# {table}'.format(table=table) in t[1]:
                s = t[0]
                e = tables[i+1][0]
                return s,e
                
    def get_table(self,table):
        df = pd.read_table()

        
        
        
    
    
     