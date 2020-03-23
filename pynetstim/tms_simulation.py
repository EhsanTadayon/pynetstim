#### registration module
from nipype.interfaces.fsl import FLIRT
from nipype import Workflow, Node
import os
from .coordinates import Coords,FreesurferCoords
from simnibs import sim_struct, run_simnibs
from .freesurfer_files import Surf
from .image_manipulations import img2img_coord_register, make_head_model

## finding simnibs_paths
paths = os.environ['PATH'].split(':')
simnibs_path = [x for x in paths if 'install_simnibs' in x][0]
simnibs_path = os.path.dirname(simnibs_path)

class Simnibs(object):

    SIMNIBSDIR= simnibs_path
    
    def __init__(self, subject, anat_img, mesh_dir, wf_base_dir):
        
        self.subject = subject
        self.anat_img = anat_img
        self.mesh_dir = mesh_dir
        self.wf_base_dir = wf_base_dir
                
        if not os.path.exists(wf_base_dir):
            os.makedirs(wf_base_dir)
            
            
    def add_targets(self, targets, csv_fname, t12simnibs=False, project2skin=True, distance=4):
        
        targets_ras = targets.coordinates['ras_coord']

        
        if t12simnibs:
            
                t12simnibs_reg = self.register_T1_to_simnibs()
                
                targets_ras = img2img_coord_register(targets_ras, img_file=self.anat_img, ref_img=os.path.join(self.mesh_dir, 'm2m_'+self.subject,'T1fs_conform.nii.gz'), method='linear',
                        linear_reg_file=t12simnibs_reg,return_as_array=True)
                

        if project2skin is True:
            
            input_img = os.path.join(self.mesh_dir, 'm2m_'+self.subject,'T1fs_conform.nii.gz')
            make_head_model(input_img, os.path.join(self.mesh_dir,'bem'))
            
            skin_surface = Surf(os.path.join(self.mesh_dir,'bem','outer_skin_surface'))
            temp = Coords(targets_ras,img_file=input_img)
            results = temp.map_to_surface(skin_surface)
            targets_ras = results['ras_coord']
            
        targets = Coords(targets_ras, img_file=self.anat_img, working_dir=targets.working_dir, **targets.get_traits_dict())

        self.targets = targets
        self.csv_fname = csv_fname
        self._targets_to_csv(csv_fname)
    
            
        
    def register_T1_to_simnibs(self):

        ### run flirt registartion if it has not been run before
        dest_img = os.path.join(self.mesh_dir, 'm2m_'+self.subject,'T1fs_conform.nii.gz')
        if not os.path.exists(os.path.join(self.wf_base_dir,'T1_to_simnibs_registration')):
            
            flirt = Node(FLIRT(),name='flirt')
            flirt.inputs.in_file = os.path.join(self.mesh_dir,'m2m_'+self.subject,'T1fs.nii.gz')
            flirt.inputs.reference = dest_img
            flirt.inputs.out_file = 'T1_in_Simnibs.nii.gz'
            flirt.inputs.out_matrix_file = 'T12Simnibs.mat'
            flirt.inputs.searchr_x = [-180,180]
            flirt.inputs.searchr_y = [-180,180]
            flirt.inputs.searchr_z = [-180,180]

            wf = Workflow(name='T1_to_simnibs_registration',base_dir=self.wf_base_dir)
            wf.add_nodes([flirt])
            wf.run()

        ## path to registration file
        t12simnibs_reg = os.path.join(self.wf_base_dir,'T1_to_simnibs_registration','flirt','T12Simnibs.mat')
        
        return t12simnibs_reg

   

    def _targets_to_csv(self, fname, dist=4):

        if not os.path.exists(os.path.join(self.wf_base_dir,'simnibs_targets')):
            os.makedirs(os.path.join(self.wf_base_dir,'simnibs_targets'))

        f = open(os.path.join(self.wf_base_dir,'simnibs_targets',fname),'w')

        for target in self.targets:
            if not hasattr(target,'name'):
                target.name='tms_target'

            row = ['CoilPos']+ target.ras_coord.tolist() + target.direction.tolist()[6:] + target.direction.tolist()[3:6] + [dist] + [target.name]
            row = [str(x) for x in row]
            f.write(','.join(row)+'\n')


    def run_simulation(self, sim_output_name,fnamecoil='MagVenture_MC_B70.ccd',didt=1e6, overwrite=False):

        s = sim_struct.SESSION()
        s.fnamehead = os.path.join(self.mesh_dir,self.subject+'.msh')
        print('output dir: ',self.wf_base_dir)
        
        try:
            if overwrite and os.path.exists(os.path.join(self.wf_base_dir),'simnibs_simulations',sim_output_name):
               os.rmdir(os.path.join(self.wf_base_dir,'simnibs_simulations',sim_output_name))
            os.makedirs(os.path.join(self.wf_base_dir,'simnibs_simulations',sim_output_name))
        except:
            raise OSError('simulation directory exists, overwrite=True to overwrite the results!')

        s.pathfem = os.path.join(self.wf_base_dir,'simnibs_simulations',sim_output_name)
        s.map_to_vol = True
        s.map_to_surf = True
        tms = s.add_tmslist()
        tms.fnamecoil = os.path.join(Simnibs.SIMNIBSDIR,'simnibs','ccd-files',fnamecoil)
        tms.add_positions_from_csv(os.path.join(self.wf_base_dir,'simnibs_targets',self.csv_fname))

        for p in tms.pos:
            p.didt = didt

        run_simnibs(s)
