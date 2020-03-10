Node: linear (fsl)
==================


 Hierarchy : register.linear
 Exec ID : linear


Original Inputs
---------------


* angle_rep : <undefined>
* apply_isoxfm : <undefined>
* apply_xfm : <undefined>
* args : <undefined>
* bbrslope : <undefined>
* bbrtype : <undefined>
* bgvalue : <undefined>
* bins : <undefined>
* coarse_search : <undefined>
* cost : <undefined>
* cost_func : <undefined>
* datatype : <undefined>
* display_init : <undefined>
* dof : <undefined>
* echospacing : <undefined>
* environ : {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
* fieldmap : <undefined>
* fieldmapmask : <undefined>
* fine_search : <undefined>
* force_scaling : <undefined>
* in_file : /Users/ehsantadayon/Dropbox/Projects/Broad/projects/pynetstim/examples/working_with_brainsight/registration/reorient2std/input_file/input_img_reorient2std.nii.gz
* in_matrix_file : <undefined>
* in_weight : <undefined>
* interp : <undefined>
* min_sampling : <undefined>
* no_clamp : <undefined>
* no_resample : <undefined>
* no_resample_blur : <undefined>
* no_search : <undefined>
* out_file : img2img_linear.nii.gz
* out_log : <undefined>
* out_matrix_file : linear_reg.mat
* output_type : NIFTI_GZ
* padding_size : <undefined>
* pedir : <undefined>
* ref_weight : <undefined>
* reference : /usr/local/fsl/data/standard/MNI152_T1_2mm.nii.gz
* rigid2D : <undefined>
* save_log : <undefined>
* schedule : <undefined>
* searchr_x : <undefined>
* searchr_y : <undefined>
* searchr_z : <undefined>
* sinc_width : <undefined>
* sinc_window : <undefined>
* uses_qform : <undefined>
* verbose : <undefined>
* wm_seg : <undefined>
* wmcoords : <undefined>
* wmnorms : <undefined>


Execution Inputs
----------------


* angle_rep : <undefined>
* apply_isoxfm : <undefined>
* apply_xfm : <undefined>
* args : <undefined>
* bbrslope : <undefined>
* bbrtype : <undefined>
* bgvalue : <undefined>
* bins : <undefined>
* coarse_search : <undefined>
* cost : <undefined>
* cost_func : <undefined>
* datatype : <undefined>
* display_init : <undefined>
* dof : <undefined>
* echospacing : <undefined>
* environ : {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
* fieldmap : <undefined>
* fieldmapmask : <undefined>
* fine_search : <undefined>
* force_scaling : <undefined>
* in_file : /Users/ehsantadayon/Dropbox/Projects/Broad/projects/pynetstim/examples/working_with_brainsight/registration/reorient2std/input_file/input_img_reorient2std.nii.gz
* in_matrix_file : <undefined>
* in_weight : <undefined>
* interp : <undefined>
* min_sampling : <undefined>
* no_clamp : <undefined>
* no_resample : <undefined>
* no_resample_blur : <undefined>
* no_search : <undefined>
* out_file : img2img_linear.nii.gz
* out_log : <undefined>
* out_matrix_file : linear_reg.mat
* output_type : NIFTI_GZ
* padding_size : <undefined>
* pedir : <undefined>
* ref_weight : <undefined>
* reference : /usr/local/fsl/data/standard/MNI152_T1_2mm.nii.gz
* rigid2D : <undefined>
* save_log : <undefined>
* schedule : <undefined>
* searchr_x : <undefined>
* searchr_y : <undefined>
* searchr_z : <undefined>
* sinc_width : <undefined>
* sinc_window : <undefined>
* uses_qform : <undefined>
* verbose : <undefined>
* wm_seg : <undefined>
* wmcoords : <undefined>
* wmnorms : <undefined>


Execution Outputs
-----------------


* out_file : <undefined>
* out_log : <undefined>
* out_matrix_file : /Users/ehsantadayon/Dropbox/Projects/Broad/projects/pynetstim/examples/working_with_brainsight/registration/register/linear/linear_reg.mat


Runtime info
------------


* cmdline : flirt -in /Users/ehsantadayon/Dropbox/Projects/Broad/projects/pynetstim/examples/working_with_brainsight/registration/reorient2std/input_file/input_img_reorient2std.nii.gz -ref /usr/local/fsl/data/standard/MNI152_T1_2mm.nii.gz -out img2img_linear.nii.gz -omat linear_reg.mat
* duration : 185.220039
* hostname : ehsans-mbp.bidmc.harvard.edu
* prev_wd : /Users/ehsantadayon/Dropbox/Projects/Broad/projects/pynetstim/examples/working_with_brainsight
* working_dir : /Users/ehsantadayon/Dropbox/Projects/Broad/projects/pynetstim/examples/working_with_brainsight/registration/register/linear


Terminal output
~~~~~~~~~~~~~~~





Terminal - standard output
~~~~~~~~~~~~~~~~~~~~~~~~~~





Terminal - standard error
~~~~~~~~~~~~~~~~~~~~~~~~~





Environment
~~~~~~~~~~~


* Apple_PubSub_Socket_Render : /private/tmp/com.apple.launchd.NZEWn5cF1t/Render
* CLICOLOR : 1
* COLORFGBG : 7;0
* COLORTERM : truecolor
* CONDA_DEFAULT_ENV : py3.6
* CONDA_EXE : /Users/ehsantadayon/anaconda2/bin/conda
* CONDA_PREFIX : /Users/ehsantadayon/anaconda2/envs/py3.6
* CONDA_PROMPT_MODIFIER : (py3.6) 
* CONDA_PYTHON_EXE : /Users/ehsantadayon/anaconda2/bin/python
* CONDA_SHLVL : 1
* DISPLAY : /private/tmp/com.apple.launchd.cwtdRtY39p/org.macosforge.xquartz:0
* FIX_VERTEX_AREA : 
* FMRI_ANALYSIS_DIR : /Applications/freesurfer/fsfast
* FREESURFER_HOME : /Applications/freesurfer
* FSFAST_HOME : /Applications/freesurfer/fsfast
* FSF_OUTPUT_FORMAT : nii.gz
* FSLDIR : /usr/local/fsl
* FSLGECUDAQ : cuda.q
* FSLLOCKDIR : 
* FSLMACHINELIST : 
* FSLMULTIFILEQUIT : TRUE
* FSLOUTPUTTYPE : NIFTI_GZ
* FSLREMOTECALL : 
* FSLTCLSH : /usr/local/fsl/bin/fsltclsh
* FSLWISH : /usr/local/fsl/bin/fslwish
* FSL_BIN : /usr/local/fsl/bin
* FSL_DIR : /usr/local/fsl
* FS_OVERRIDE : 0
* FUNCTIONALS_DIR : /Applications/freesurfer/sessions
* GIT_PAGER : cat
* GROUP : staff
* HOME : /Users/ehsantadayon
* HOST : ehsans-mbp.bidmc.harvard.edu
* HOSTTYPE : unknown
* ITERM_PROFILE : Default
* ITERM_SESSION_ID : w0t0p4:6FA3CB23-02E6-42D5-8A32-BEB172BA8133
* JPY_PARENT_PID : 88302
* KERNEL_LAUNCH_TIMEOUT : 40
* KMP_DUPLICATE_LIB_OK : True
* KMP_INIT_AT_FORK : FALSE
* LANG : en_US.UTF-8
* LC_TERMINAL : iTerm2
* LC_TERMINAL_VERSION : 3.3.8
* LOCAL_DIR : /Applications/freesurfer/local
* LOGNAME : ehsantadayon
* MACHTYPE : x86_64
* MINC_BIN_DIR : /Applications/freesurfer/mni/bin
* MINC_LIB_DIR : /Applications/freesurfer/mni/lib
* MNI_DATAPATH : /Applications/freesurfer/mni/data
* MNI_DIR : /Applications/freesurfer/mni
* MNI_PERL5LIB : /Applications/freesurfer/mni/lib/../System/Library/Perl/5.8.6
* MPLBACKEND : module://ipykernel.pylab.backend_inline
* OS : Darwin
* OSTYPE : darwin
* PAGER : cat
* PATH : /Users/ehsantadayon/anaconda2/envs/py3.6/bin:/Users/ehsantadayon/anaconda2/bin:/Applications/freesurfer/bin:/Applications/freesurfer/fsfast/bin:/Applications/freesurfer/tktools:/usr/local/fsl/bin:/Applications/freesurfer/bin/freeview.app/Contents/MacOS/:/Applications/freesurfer/mni/bin:/usr/local/fsl/bin:/Library/Frameworks/Python.framework/Versions/2.7/bin:/Library/Frameworks/Python.framework/Versions/Current/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/Users/ehsantadayon/Library/Enthought/Canopy_64bit/User/bin:/Library/TeX/texbin:/opt/X11/bin:/usr/local/git/bin:/Applications/FSLeyes.app/Contents/MacOS
* PERL5LIB : /Applications/freesurfer/mni/lib/../System/Library/Perl/5.8.6
* PWD : /Users/ehsantadayon/Dropbox/Projects/Broad/projects/pynetstim/examples/working_with_brainsight
* SHELL : /bin/tcsh
* SHLVL : 2
* SSH_AUTH_SOCK : /private/tmp/com.apple.launchd.9LyJ8e0jgb/Listeners
* SUBJECTS_DIR : /Users/ehsantadayon/Dropbox/Projects/Broad/projects/pynetstim/examples/working_with_brainsight/data/example_data/freesurfer
* TERM : xterm-color
* TERM_PROGRAM : iTerm.app
* TERM_PROGRAM_VERSION : 3.3.8
* TERM_SESSION_ID : w0t0p4:6FA3CB23-02E6-42D5-8A32-BEB172BA8133
* TMPDIR : /var/folders/nw/z3t9lfvj6ls9mqzdrbk6b4l80000gn/T/
* USER : ehsantadayon
* VENDOR : apple
* XPC_FLAGS : 0x0
* XPC_SERVICE_NAME : 0
* _ : /Users/ehsantadayon/anaconda2/envs/py3.6/bin/jupyter
* __CF_USER_TEXT_ENCODING : 0x1F5:0x0:0x0

