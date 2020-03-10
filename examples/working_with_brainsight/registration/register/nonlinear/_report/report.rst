Node: nonlinear (fsl)
=====================


 Hierarchy : register.nonlinear
 Exec ID : nonlinear


Original Inputs
---------------


* affine_file : /Users/ehsantadayon/Dropbox/Projects/Broad/projects/pynetstim/examples/working_with_brainsight/registration/register/linear/linear_reg.mat
* apply_inmask : <undefined>
* apply_intensity_mapping : <undefined>
* apply_refmask : <undefined>
* args : <undefined>
* bias_regularization_lambda : <undefined>
* biasfield_resolution : <undefined>
* config_file : <undefined>
* derive_from_ref : <undefined>
* environ : {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
* field_file : <undefined>
* fieldcoeff_file : img2img_nonlinear_fieldcoeff.nii.gz
* hessian_precision : <undefined>
* in_file : /Users/ehsantadayon/Dropbox/Projects/Broad/projects/pynetstim/examples/working_with_brainsight/registration/reorient2std/input_file/input_img_reorient2std.nii.gz
* in_fwhm : <undefined>
* in_intensitymap_file : <undefined>
* inmask_file : <undefined>
* inmask_val : <undefined>
* intensity_mapping_model : <undefined>
* intensity_mapping_order : <undefined>
* inwarp_file : <undefined>
* jacobian_file : <undefined>
* jacobian_range : <undefined>
* log_file : <undefined>
* max_nonlin_iter : <undefined>
* modulatedref_file : <undefined>
* out_intensitymap_file : <undefined>
* output_type : NIFTI_GZ
* ref_file : /usr/local/fsl/data/standard/MNI152_T1_2mm.nii.gz
* ref_fwhm : <undefined>
* refmask_file : <undefined>
* refmask_val : <undefined>
* regularization_lambda : <undefined>
* regularization_model : <undefined>
* skip_implicit_in_masking : <undefined>
* skip_implicit_ref_masking : <undefined>
* skip_inmask : <undefined>
* skip_intensity_mapping : <undefined>
* skip_lambda_ssq : <undefined>
* skip_refmask : <undefined>
* spline_order : <undefined>
* subsampling_scheme : <undefined>
* warp_resolution : <undefined>
* warped_file : img2img_nonlinear.nii.gz


Execution Inputs
----------------


* affine_file : /Users/ehsantadayon/Dropbox/Projects/Broad/projects/pynetstim/examples/working_with_brainsight/registration/register/linear/linear_reg.mat
* apply_inmask : <undefined>
* apply_intensity_mapping : <undefined>
* apply_refmask : <undefined>
* args : <undefined>
* bias_regularization_lambda : <undefined>
* biasfield_resolution : <undefined>
* config_file : <undefined>
* derive_from_ref : <undefined>
* environ : {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
* field_file : <undefined>
* fieldcoeff_file : img2img_nonlinear_fieldcoeff.nii.gz
* hessian_precision : <undefined>
* in_file : /Users/ehsantadayon/Dropbox/Projects/Broad/projects/pynetstim/examples/working_with_brainsight/registration/reorient2std/input_file/input_img_reorient2std.nii.gz
* in_fwhm : <undefined>
* in_intensitymap_file : <undefined>
* inmask_file : <undefined>
* inmask_val : <undefined>
* intensity_mapping_model : <undefined>
* intensity_mapping_order : <undefined>
* inwarp_file : <undefined>
* jacobian_file : <undefined>
* jacobian_range : <undefined>
* log_file : <undefined>
* max_nonlin_iter : <undefined>
* modulatedref_file : <undefined>
* out_intensitymap_file : <undefined>
* output_type : NIFTI_GZ
* ref_file : /usr/local/fsl/data/standard/MNI152_T1_2mm.nii.gz
* ref_fwhm : <undefined>
* refmask_file : <undefined>
* refmask_val : <undefined>
* regularization_lambda : <undefined>
* regularization_model : <undefined>
* skip_implicit_in_masking : <undefined>
* skip_implicit_ref_masking : <undefined>
* skip_inmask : <undefined>
* skip_intensity_mapping : <undefined>
* skip_lambda_ssq : <undefined>
* skip_refmask : <undefined>
* spline_order : <undefined>
* subsampling_scheme : <undefined>
* warp_resolution : <undefined>
* warped_file : img2img_nonlinear.nii.gz


Execution Outputs
-----------------


* field_file : <undefined>
* fieldcoeff_file : /Users/ehsantadayon/Dropbox/Projects/Broad/projects/pynetstim/examples/working_with_brainsight/registration/register/nonlinear/img2img_nonlinear_fieldcoeff.nii.gz
* jacobian_file : <undefined>
* log_file : /Users/ehsantadayon/Dropbox/Projects/Broad/projects/pynetstim/examples/working_with_brainsight/registration/register/nonlinear/input_img_reorient2std_log.txt
* modulatedref_file : <undefined>
* out_intensitymap_file : <undefined>
* warped_file : /Users/ehsantadayon/Dropbox/Projects/Broad/projects/pynetstim/examples/working_with_brainsight/registration/register/nonlinear/img2img_nonlinear.nii.gz


Runtime info
------------


* cmdline : fnirt --aff=/Users/ehsantadayon/Dropbox/Projects/Broad/projects/pynetstim/examples/working_with_brainsight/registration/register/linear/linear_reg.mat --cout=/Users/ehsantadayon/Dropbox/Projects/Broad/projects/pynetstim/examples/working_with_brainsight/registration/register/nonlinear/img2img_nonlinear_fieldcoeff.nii.gz --in=/Users/ehsantadayon/Dropbox/Projects/Broad/projects/pynetstim/examples/working_with_brainsight/registration/reorient2std/input_file/input_img_reorient2std.nii.gz --logout=/Users/ehsantadayon/Dropbox/Projects/Broad/projects/pynetstim/examples/working_with_brainsight/registration/register/nonlinear/input_img_reorient2std_log.txt --ref=/usr/local/fsl/data/standard/MNI152_T1_2mm.nii.gz --iout=/Users/ehsantadayon/Dropbox/Projects/Broad/projects/pynetstim/examples/working_with_brainsight/registration/register/nonlinear/img2img_nonlinear.nii.gz
* duration : 1789.232479
* hostname : ehsans-mbp.bidmc.harvard.edu
* prev_wd : /Users/ehsantadayon/Dropbox/Projects/Broad/projects/pynetstim/examples/working_with_brainsight
* working_dir : /Users/ehsantadayon/Dropbox/Projects/Broad/projects/pynetstim/examples/working_with_brainsight/registration/register/nonlinear


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

