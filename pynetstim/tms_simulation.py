#### registration module
from nipype.interfaces.fsl import FLIRT
from nipype import Workflow, Node
import os
from .coordinates import Coords,FreesurferCoords
from simnibs import sim_struct, run_simnibs
from .freesurfer_files import Surf
from .utils import make_head_model

class Simnibs(object):

    SIMNIBSDIR='/Users/davidemomi/Applications/SimNIBS/'
    def __init__(self, subject, freesurfer_dir, mri2mesh_dir, out_dir):
        self.subject = subject
        self.mri2mesh_dir = mri2mesh_dir
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        self.out_dir = out_dir
        self.freesurfer_dir = freesurfer_dir



    def register_T1_to_simnibs(self):

        ### run flirt registartion if it has not been run before
        if not os.path.exists(os.path.join(self.out_dir,'T1_to_simnibs_registration')):

            flirt = Node(FLIRT(),name='flirt')
            flirt.inputs.in_file = os.path.join(self.mri2mesh_dir,'m2m_'+self.subject,'T1fs.nii.gz')
            flirt.inputs.reference = os.path.join(self.mri2mesh_dir, 'm2m_'+self.subject,'T1fs_conform.nii.gz')
            flirt.inputs.out_file = 'T1_in_Simnibs.nii.gz'
            flirt.inputs.out_matrix_file = 'T12Simnibs.mat'
            flirt.inputs.searchr_x = [-180,180]
            flirt.inputs.searchr_y = [-180,180]
            flirt.inputs.searchr_z = [-180,180]

            wf = Workflow(name='T1_to_simnibs_registration',base_dir=self.out_dir)
            wf.add_nodes([flirt])
            wf.run()

        ## path to registration file
        t12simnibs_reg = os.path.join(self.out_dir,'T1_to_simnibs_registration','flirt','T12Simnibs.mat')

        return t12simnibs_reg

    def add_targets(self, targets, csv_fname, t12simnibs=False, project2skin=True, distance=4):

        ras_coords = targets.coordinates['ras_coord']

        if t12simnibs:
                t12simnibs_reg = self.register_T1_to_simnibs()
                targets = targets.img2imgcoord(os.path.join(self.mri2mesh_dir, 'm2m_'+self.subject,'T1fs_conform.nii.gz'),t12simnibs_reg,type='xfm')

        if project2skin is True:
            if os.path.exists(os.path.join(self.freesurfer_dir,self.subject,'bem','outer_skin_surface')):
                pass
            else:
                make_head_model(os.path.join(self.freesurfer_dir,self.subject,'mri','rawavg.mgz'), os.path.join(self.freesurfer_dir,self.subject,'bem'))
            skin_surface = Surf(os.path.join(self.freesurfer_dir,self.subject,'bem','outer_skin_surface'))
            results = targets.map_to_surface(skin_surface)

            targets = FreesurferCoords(results['ras_coord'], self.subject, self.freesurfer_dir, direction = targets.direction)

        self.targets = targets
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


    def run_simulation(self, sim_output_name,fnamecoil='MagVenture_MC_B70.ccd',didt=1e6, overwrite=False):

        s = sim_struct.SESSION()
        s.fnamehead = os.path.join(self.mri2mesh_dir,self.subject+'.msh')
        print('output dir: ',self.out_dir)
        try:
            if overwrite and os.path.exists(os.path.join(self.out_dir),'simnibs_simulations',sim_output_name):
               os.rmdir(os.path.join(self.out_dir,'simnibs_simulations',sim_output_name))
            os.makedirs(os.path.join(self.out_dir,'simnibs_simulations',sim_output_name))
        except:
            raise OSError('simulation directory exists, overwrite=True to overwrite the results!')

        s.pathfem = os.path.join(self.out_dir,'simnibs_simulations',sim_output_name)
        s.map_to_vol = True
        s.map_to_surf = True
        tms = s.add_tmslist()
        tms.fnamecoil = os.path.join(Simnibs.SIMNIBSDIR,'simnibs','ccd-files',fnamecoil)
        tms.add_positions_from_csv(os.path.join(self.out_dir,'simnibs_targets',self.csv_fname))

        for p in tms.pos:
            p.didt = didt

        run_simnibs(s)
