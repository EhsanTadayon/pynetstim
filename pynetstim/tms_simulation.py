#### registration module 
from nipype.interfaces.fsl import FLIRT
from nipype import Workflow, Node
import os
from targets import Targets
from coordinates import Coords
from simnibs import sim_struct, run_simulation
from surface import Surf
from utils import make_head_model

class Simnibs(object):
    
    SIMNIBSDIR='/Users/stadayon/simnibs_2.1.1'
    def __init__(self, subject, simnibs_dir, out_dir):
        self.subject = subject
        self.simnibs_dir = simnibs_dir
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        self.out_dir = out_dir
        make_head_model('fs_'+subject, self.simnibs_dir)
        
        
    def register_T1_to_simnibs(self):
        
        flirt = Node(FLIRT(),name='flirt')
        flirt.inputs.in_file = os.path.join(self.simnibs_dir,'m2m_'+self.subject,'T1fs.nii.gz')
        flirt.inputs.reference = os.path.join(self.simnibs_dir, 'm2m_'+self.subject,'T1fs_conform.nii.gz')
        flirt.inputs.out_file = 'T1_in_Simnibs.nii.gz'
        flirt.inputs.out_matrix_file = 'T12Simnibs.mat'
        flirt.inputs.searchr_x = [-180,180]
        flirt.inputs.searchr_y = [-180,180]
        flirt.inputs.searchr_z = [-180,180]
    
        wf = Workflow(name='T1_to_simnibs_registration',base_dir=self.out_dir)
        wf.add_nodes([flirt])
        wf.run()
        
    def add_targets(self, targets, csv_fname, t12simnibs=False, project2skin=True, distance=4):
        
        ras_coords = targets.coordinates['ras_coord']
        
        if t12simnibs:
            if not os.path.exists(os.path.join(self.out_dir,'T1_to_simnibs_registration')):
                self.register_T1_to_simnibs()
            
                
            coordsObj = Coords(targets.coordinates['ras_coord'],os.path.join(self.simnibs_dir,'m2m_'+self.subject,'T1fs.nii.gz'))
            t12simnibs_reg = os.path.join(self.out_dir,'T1_to_simnibs_registration','flirt','T12Simnibs.mat')
            ras_coords = coordsObj.img2imgcoord(os.path.join(self.simnibs_dir, 'm2m_'+self.subject,'T1fs_conform.nii.gz'),t12simnibs_reg,type='xfm')
            
        if project2skin is True:
            if not os.path.isfile('{subjects_dir}/{subject}/bem/outer_skin_surface'.format(subjects_dir=self.simnibs_dir,subject='fs_'+self.subject)):
                make_head_model('fs_'+self.subject,self.simnibs_dir)
                
            skin_surf = Surf('{simnibs_dir}/{subject}/bem/outer_skin_surface'.format(simnibs_dir=self.simnibs_dir, subject='fs_'+self.subject))
            mapped_vertices, ras_coords = skin_surf.project_coords(ras_coords)
            
        self.targets = Targets(ras_coords, 'fs_'+self.subject, self.simnibs_dir, direction = targets.direction)
        
        self.csv_fname = csv_fname   
        self._targets_to_csv(csv_fname)
            
            
    def _targets_to_csv(self, fname, dist=4):
        
        if not os.path.exists(os.path.join(self.out_dir,'simnibs_targets')):
            os.makedirs(os.path.join(self.out_dir,'simnibs_targets'))
            
        f = open(os.path.join(self.out_dir,'simnibs_targets',fname),'w')
            
        for target in self.targets:
            if not hasattr(target,'name'):
                target.name='tms_target'
                
            row = ['CoilPos']+ target.ras_coord.tolist() + target.direction.tolist()[6:] + target.direction.tolist()[3:6] + [dist] + [target.name]
            row = [str(x) for x in row]
            f.write(','.join(row)+'\n')
        
        
    
    def run_simnibs(self, sim_output_name,fnamecoil='MagVenture_MC_B70.ccd',didt=1e6, overwrite=False):
        
        s = sim_struct.SESSION()
        s.fnamehead = os.path.join(self.simnibs_dir,self.subject+'.msh')
        try:
            if overwrite:
               os.rmdir(os.path.join(self.out_dir,'simnibs_simulations',sim_output_name)) 
            os.makedirs(os.path.join(self.out_dir,'simnibs_simulations',sim_output_name))
        except:
            raise 
        s.pathfem = os.path.join(self.out_dir,'simnibs_simulations',sim_output_name)
        s.map_to_vol = True
        s.map_to_surf = True
        tms = s.add_tmslist()
        tms.fnamecoil = os.path.join(Simnibs.SIMNIBSDIR,'ccd-files',fnamecoil)
        tms.add_positions_from_csv(os.path.join(self.out_dir,'simnibs_targets',self.csv_fname))
        for p in tms.pos:
            p.didt = didt
        
        run_simulation(s)
    
            
    
        
    
            
        

    