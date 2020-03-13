

#---------------------------------
# New invocation of recon-all Tue Feb  7 01:09:27 EST 2017 

 mri_convert /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_Raw/MNI152_T1_2mm.nii.gz /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_FS2/mri/orig/001.mgz 

#--------------------------------------------
#@# MotionCor Tue Feb  7 01:09:31 EST 2017

 cp /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_FS2/mri/orig/001.mgz /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_FS2/mri/rawavg.mgz 


 mri_convert /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_FS2/mri/rawavg.mgz /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_FS2/mri/orig.mgz --conform 


 mri_add_xform_to_header -c /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_FS2/mri/transforms/talairach.xfm /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_FS2/mri/orig.mgz /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_FS2/mri/orig.mgz 

#--------------------------------------------
#@# Talairach Tue Feb  7 01:09:41 EST 2017

 mri_nu_correct.mni --no-rescale --i orig.mgz --o orig_nu.mgz --n 1 --proto-iters 1000 --distance 50 


 talairach_avi --i orig_nu.mgz --xfm transforms/talairach.auto.xfm 

talairach_avi log file is transforms/talairach_avi.log...

 cp transforms/talairach.auto.xfm transforms/talairach.xfm 

#--------------------------------------------
#@# Talairach Failure Detection Tue Feb  7 01:13:24 EST 2017

 talairach_afd -T 0.005 -xfm transforms/talairach.xfm 


 awk -f /ncf/nrg/sw/apps/freesurfer/6_2015_04_21/bin/extract_talairach_avi_QA.awk /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_FS2/mri/transforms/talairach_avi.log 


 tal_QC_AZS /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_FS2/mri/transforms/talairach_avi.log 

#--------------------------------------------
#@# Nu Intensity Correction Tue Feb  7 01:13:24 EST 2017

 mri_nu_correct.mni --i orig.mgz --o nu.mgz --uchar transforms/talairach.xfm --n 2 


 mri_add_xform_to_header -c /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_FS2/mri/transforms/talairach.xfm nu.mgz nu.mgz 

#--------------------------------------------
#@# Intensity Normalization Tue Feb  7 01:17:00 EST 2017

 mri_normalize -g 1 -mprage nu.mgz T1.mgz 

#--------------------------------------------
#@# Skull Stripping Tue Feb  7 01:23:49 EST 2017

 mri_em_register -rusage /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_FS2/touch/rusage.mri_em_register.skull.dat -skull nu.mgz /ncf/nrg/sw/apps/freesurfer/6_2015_04_21/average/RB_all_withskull_2016-03-21.gca transforms/talairach_with_skull.lta 


 mri_watershed -rusage /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_FS2/touch/rusage.mri_watershed.dat -T1 -brain_atlas /ncf/nrg/sw/apps/freesurfer/6_2015_04_21/average/RB_all_withskull_2016-03-21.gca transforms/talairach_with_skull.lta T1.mgz brainmask.auto.mgz 


 cp brainmask.auto.mgz brainmask.mgz 

#-------------------------------------
#@# EM Registration Tue Feb  7 02:06:00 EST 2017

 mri_em_register -rusage /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_FS2/touch/rusage.mri_em_register.dat -uns 3 -mask brainmask.mgz nu.mgz /ncf/nrg/sw/apps/freesurfer/6_2015_04_21/average/RB_all_2016-03-21.gca transforms/talairach.lta 

#--------------------------------------
#@# CA Normalize Tue Feb  7 02:41:15 EST 2017

 mri_ca_normalize -c ctrl_pts.mgz -mask brainmask.mgz nu.mgz /ncf/nrg/sw/apps/freesurfer/6_2015_04_21/average/RB_all_2016-03-21.gca transforms/talairach.lta norm.mgz 

#--------------------------------------
#@# CA Reg Tue Feb  7 02:44:20 EST 2017

 mri_ca_register -rusage /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_FS2/touch/rusage.mri_ca_register.dat -nobigventricles -T transforms/talairach.lta -align-after -mask brainmask.mgz norm.mgz /ncf/nrg/sw/apps/freesurfer/6_2015_04_21/average/RB_all_2016-03-21.gca transforms/talairach.m3z 

#--------------------------------------
#@# Remove Neck Tue Feb  7 05:22:48 EST 2017

 mri_remove_neck -radius 25 nu.mgz transforms/talairach.m3z /ncf/nrg/sw/apps/freesurfer/6_2015_04_21/average/RB_all_2016-03-21.gca nu_noneck.mgz 

#--------------------------------------
#@# SubCort Seg Tue Feb  7 05:24:40 EST 2017

 mri_ca_label -rusage /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_FS2/touch/rusage.mri_ca_label.dat -relabel_unlikely 9 .3 -prior 0.5 -align norm.mgz transforms/talairach.m3z /ncf/nrg/sw/apps/freesurfer/6_2015_04_21/average/RB_all_2016-03-21.gca aseg.auto_noCCseg.mgz 


 mri_cc -aseg aseg.auto_noCCseg.mgz -o aseg.auto.mgz -lta /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_FS2/mri/transforms/cc_up.lta MNI_FS2 

#--------------------------------------
#@# Merge ASeg Tue Feb  7 07:06:47 EST 2017

 cp aseg.auto.mgz aseg.presurf.mgz 

#--------------------------------------------
#@# Intensity Normalization2 Tue Feb  7 07:06:47 EST 2017

 mri_normalize -mprage -aseg aseg.presurf.mgz -mask brainmask.mgz norm.mgz brain.mgz 

#--------------------------------------------
#@# Mask BFS Tue Feb  7 07:13:24 EST 2017

 mri_mask -T 5 brain.mgz brainmask.mgz brain.finalsurfs.mgz 

#--------------------------------------------
#@# WM Segmentation Tue Feb  7 07:13:28 EST 2017

 mri_segment -mprage brain.mgz wm.seg.mgz 


 mri_edit_wm_with_aseg -keep-in wm.seg.mgz brain.mgz aseg.presurf.mgz wm.asegedit.mgz 


 mri_pretess wm.asegedit.mgz wm norm.mgz wm.mgz 

#--------------------------------------------
#@# Fill Tue Feb  7 07:18:36 EST 2017

 mri_fill -a ../scripts/ponscc.cut.log -xform transforms/talairach.lta -segmentation aseg.auto_noCCseg.mgz wm.mgz filled.mgz 

#--------------------------------------------
#@# Tessellate lh Tue Feb  7 07:19:52 EST 2017

 mri_pretess ../mri/filled.mgz 255 ../mri/norm.mgz ../mri/filled-pretess255.mgz 


 mri_tessellate ../mri/filled-pretess255.mgz 255 ../surf/lh.orig.nofix 


 rm -f ../mri/filled-pretess255.mgz 


 mris_extract_main_component ../surf/lh.orig.nofix ../surf/lh.orig.nofix 

#--------------------------------------------
#@# Tessellate rh Tue Feb  7 07:20:06 EST 2017

 mri_pretess ../mri/filled.mgz 127 ../mri/norm.mgz ../mri/filled-pretess127.mgz 


 mri_tessellate ../mri/filled-pretess127.mgz 127 ../surf/rh.orig.nofix 


 rm -f ../mri/filled-pretess127.mgz 


 mris_extract_main_component ../surf/rh.orig.nofix ../surf/rh.orig.nofix 

#--------------------------------------------
#@# Smooth1 lh Tue Feb  7 07:20:19 EST 2017

 mris_smooth -nw -seed 1234 ../surf/lh.orig.nofix ../surf/lh.smoothwm.nofix 

#--------------------------------------------
#@# Smooth1 rh Tue Feb  7 07:20:30 EST 2017

 mris_smooth -nw -seed 1234 ../surf/rh.orig.nofix ../surf/rh.smoothwm.nofix 

#--------------------------------------------
#@# Inflation1 lh Tue Feb  7 07:20:41 EST 2017

 mris_inflate -no-save-sulc ../surf/lh.smoothwm.nofix ../surf/lh.inflated.nofix 

#--------------------------------------------
#@# Inflation1 rh Tue Feb  7 07:21:39 EST 2017

 mris_inflate -no-save-sulc ../surf/rh.smoothwm.nofix ../surf/rh.inflated.nofix 

#--------------------------------------------
#@# QSphere lh Tue Feb  7 07:22:37 EST 2017

 mris_sphere -q -seed 1234 ../surf/lh.inflated.nofix ../surf/lh.qsphere.nofix 

#--------------------------------------------
#@# QSphere rh Tue Feb  7 07:28:33 EST 2017

 mris_sphere -q -seed 1234 ../surf/rh.inflated.nofix ../surf/rh.qsphere.nofix 

#--------------------------------------------
#@# Fix Topology Copy lh Tue Feb  7 07:34:26 EST 2017

 cp ../surf/lh.orig.nofix ../surf/lh.orig 


 cp ../surf/lh.inflated.nofix ../surf/lh.inflated 

#--------------------------------------------
#@# Fix Topology Copy rh Tue Feb  7 07:34:26 EST 2017

 cp ../surf/rh.orig.nofix ../surf/rh.orig 


 cp ../surf/rh.inflated.nofix ../surf/rh.inflated 

#@# Fix Topology lh Tue Feb  7 07:34:27 EST 2017

 mris_fix_topology -rusage /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_FS2/touch/rusage.mris_fix_topology.lh.dat -mgz -sphere qsphere.nofix -ga -seed 1234 MNI_FS2 lh 

#@# Fix Topology rh Tue Feb  7 08:21:57 EST 2017

 mris_fix_topology -rusage /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_FS2/touch/rusage.mris_fix_topology.rh.dat -mgz -sphere qsphere.nofix -ga -seed 1234 MNI_FS2 rh 


 mris_euler_number ../surf/lh.orig 


 mris_euler_number ../surf/rh.orig 


 mris_remove_intersection ../surf/lh.orig ../surf/lh.orig 


 rm ../surf/lh.inflated 


 mris_remove_intersection ../surf/rh.orig ../surf/rh.orig 


 rm ../surf/rh.inflated 

#--------------------------------------------
#@# Make White Surf lh Tue Feb  7 09:12:26 EST 2017

 mris_make_surfaces -aseg ../mri/aseg.presurf -noaparc -whiteonly -mgz -T1 brain.finalsurfs MNI_FS2 lh 

#--------------------------------------------
#@# Make White Surf rh Tue Feb  7 09:20:21 EST 2017

 mris_make_surfaces -aseg ../mri/aseg.presurf -noaparc -whiteonly -mgz -T1 brain.finalsurfs MNI_FS2 rh 

#--------------------------------------------
#@# Smooth2 lh Tue Feb  7 09:28:19 EST 2017

 mris_smooth -n 3 -nw -seed 1234 ../surf/lh.white ../surf/lh.smoothwm 

#--------------------------------------------
#@# Smooth2 rh Tue Feb  7 09:28:30 EST 2017

 mris_smooth -n 3 -nw -seed 1234 ../surf/rh.white ../surf/rh.smoothwm 

#--------------------------------------------
#@# Inflation2 lh Tue Feb  7 09:28:41 EST 2017

 mris_inflate -rusage /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_FS2/touch/rusage.mris_inflate.lh.dat ../surf/lh.smoothwm ../surf/lh.inflated 

#--------------------------------------------
#@# Inflation2 rh Tue Feb  7 09:29:37 EST 2017

 mris_inflate -rusage /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_FS2/touch/rusage.mris_inflate.rh.dat ../surf/rh.smoothwm ../surf/rh.inflated 

#--------------------------------------------
#@# Curv .H and .K lh Tue Feb  7 09:30:32 EST 2017

 mris_curvature -w lh.white 


 mris_curvature -thresh .999 -n -a 5 -w -distances 10 10 lh.inflated 

#--------------------------------------------
#@# Curv .H and .K rh Tue Feb  7 09:33:05 EST 2017

 mris_curvature -w rh.white 


 mris_curvature -thresh .999 -n -a 5 -w -distances 10 10 rh.inflated 


#-----------------------------------------
#@# Curvature Stats lh Tue Feb  7 09:35:38 EST 2017

 mris_curvature_stats -m --writeCurvatureFiles -G -o ../stats/lh.curv.stats -F smoothwm MNI_FS2 lh curv sulc 


#-----------------------------------------
#@# Curvature Stats rh Tue Feb  7 09:35:45 EST 2017

 mris_curvature_stats -m --writeCurvatureFiles -G -o ../stats/rh.curv.stats -F smoothwm MNI_FS2 rh curv sulc 

#--------------------------------------------
#@# Sphere lh Tue Feb  7 09:35:53 EST 2017

 mris_sphere -rusage /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_FS2/touch/rusage.mris_sphere.lh.dat -seed 1234 ../surf/lh.inflated ../surf/lh.sphere 

#--------------------------------------------
#@# Sphere rh Tue Feb  7 10:50:27 EST 2017

 mris_sphere -rusage /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_FS2/touch/rusage.mris_sphere.rh.dat -seed 1234 ../surf/rh.inflated ../surf/rh.sphere 

#--------------------------------------------
#@# Surf Reg lh Tue Feb  7 12:19:16 EST 2017

 mris_register -curv -rusage /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_FS2/touch/rusage.mris_register.lh.dat ../surf/lh.sphere /ncf/nrg/sw/apps/freesurfer/6_2015_04_21/average/lh.curvature.buckner40.2016-03-20.tif ../surf/lh.sphere.reg 

#--------------------------------------------
#@# Surf Reg rh Tue Feb  7 13:37:24 EST 2017

 mris_register -curv -rusage /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/MNI_FS2/touch/rusage.mris_register.rh.dat ../surf/rh.sphere /ncf/nrg/sw/apps/freesurfer/6_2015_04_21/average/rh.curvature.buckner40.2016-03-20.tif ../surf/rh.sphere.reg 

#--------------------------------------------
#@# Jacobian white lh Tue Feb  7 15:07:28 EST 2017

 mris_jacobian ../surf/lh.white ../surf/lh.sphere.reg ../surf/lh.jacobian_white 

#--------------------------------------------
#@# Jacobian white rh Tue Feb  7 15:07:31 EST 2017

 mris_jacobian ../surf/rh.white ../surf/rh.sphere.reg ../surf/rh.jacobian_white 

#--------------------------------------------
#@# AvgCurv lh Tue Feb  7 15:07:34 EST 2017

 mrisp_paint -a 5 /ncf/nrg/sw/apps/freesurfer/6_2015_04_21/average/lh.curvature.buckner40.2016-03-20.tif#6 ../surf/lh.sphere.reg ../surf/lh.avg_curv 

#--------------------------------------------
#@# AvgCurv rh Tue Feb  7 15:07:37 EST 2017

 mrisp_paint -a 5 /ncf/nrg/sw/apps/freesurfer/6_2015_04_21/average/rh.curvature.buckner40.2016-03-20.tif#6 ../surf/rh.sphere.reg ../surf/rh.avg_curv 

#-----------------------------------------
#@# Cortical Parc lh Tue Feb  7 15:07:40 EST 2017

 mris_ca_label -l ../label/lh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 MNI_FS2 lh ../surf/lh.sphere.reg /ncf/nrg/sw/apps/freesurfer/6_2015_04_21/average/lh.DKatlas.2016-03-20.gcs ../label/lh.aparc.annot 

#-----------------------------------------
#@# Cortical Parc rh Tue Feb  7 15:08:04 EST 2017

 mris_ca_label -l ../label/rh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 MNI_FS2 rh ../surf/rh.sphere.reg /ncf/nrg/sw/apps/freesurfer/6_2015_04_21/average/rh.DKatlas.2016-03-20.gcs ../label/rh.aparc.annot 

#--------------------------------------------
#@# Make Pial Surf lh Tue Feb  7 15:08:29 EST 2017

 mris_make_surfaces -orig_white white -orig_pial white -aseg ../mri/aseg.presurf -nowhite -mgz -T1 brain.finalsurfs MNI_FS2 lh 

#--------------------------------------------
#@# Make Pial Surf rh Tue Feb  7 15:23:15 EST 2017

 mris_make_surfaces -orig_white white -orig_pial white -aseg ../mri/aseg.presurf -nowhite -mgz -T1 brain.finalsurfs MNI_FS2 rh 

#--------------------------------------------
#@# Surf Volume lh Tue Feb  7 15:38:02 EST 2017
#--------------------------------------------
#@# Surf Volume rh Tue Feb  7 15:38:07 EST 2017
#--------------------------------------------
#@# Cortical ribbon mask Tue Feb  7 15:38:13 EST 2017

 mris_volmask --aseg_name aseg.presurf --label_left_white 2 --label_left_ribbon 3 --label_right_white 41 --label_right_ribbon 42 --save_ribbon MNI_FS2 

#-----------------------------------------
#@# Parcellation Stats lh Tue Feb  7 15:47:43 EST 2017

 mris_anatomical_stats -th3 -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.stats -b -a ../label/lh.aparc.annot -c ../label/aparc.annot.ctab MNI_FS2 lh white 


 mris_anatomical_stats -th3 -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.pial.stats -b -a ../label/lh.aparc.annot -c ../label/aparc.annot.ctab MNI_FS2 lh pial 

#-----------------------------------------
#@# Parcellation Stats rh Tue Feb  7 15:50:01 EST 2017

 mris_anatomical_stats -th3 -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.stats -b -a ../label/rh.aparc.annot -c ../label/aparc.annot.ctab MNI_FS2 rh white 


 mris_anatomical_stats -th3 -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.pial.stats -b -a ../label/rh.aparc.annot -c ../label/aparc.annot.ctab MNI_FS2 rh pial 

#-----------------------------------------
#@# Cortical Parc 2 lh Tue Feb  7 15:52:19 EST 2017

 mris_ca_label -l ../label/lh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 MNI_FS2 lh ../surf/lh.sphere.reg /ncf/nrg/sw/apps/freesurfer/6_2015_04_21/average/lh.CDatlas.2016-03-20.gcs ../label/lh.aparc.a2009s.annot 

#-----------------------------------------
#@# Cortical Parc 2 rh Tue Feb  7 15:52:51 EST 2017

 mris_ca_label -l ../label/rh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 MNI_FS2 rh ../surf/rh.sphere.reg /ncf/nrg/sw/apps/freesurfer/6_2015_04_21/average/rh.CDatlas.2016-03-20.gcs ../label/rh.aparc.a2009s.annot 

#-----------------------------------------
#@# Parcellation Stats 2 lh Tue Feb  7 15:53:24 EST 2017

 mris_anatomical_stats -th3 -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.a2009s.stats -b -a ../label/lh.aparc.a2009s.annot -c ../label/aparc.annot.a2009s.ctab MNI_FS2 lh white 

#-----------------------------------------
#@# Parcellation Stats 2 rh Tue Feb  7 15:54:34 EST 2017

 mris_anatomical_stats -th3 -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.a2009s.stats -b -a ../label/rh.aparc.a2009s.annot -c ../label/aparc.annot.a2009s.ctab MNI_FS2 rh white 

#-----------------------------------------
#@# Cortical Parc 3 lh Tue Feb  7 15:55:44 EST 2017

 mris_ca_label -l ../label/lh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 MNI_FS2 lh ../surf/lh.sphere.reg /ncf/nrg/sw/apps/freesurfer/6_2015_04_21/average/lh.DKTatlas.2016-03-20.gcs ../label/lh.aparc.DKTatlas.annot 

#-----------------------------------------
#@# Cortical Parc 3 rh Tue Feb  7 15:56:12 EST 2017

 mris_ca_label -l ../label/rh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 MNI_FS2 rh ../surf/rh.sphere.reg /ncf/nrg/sw/apps/freesurfer/6_2015_04_21/average/rh.DKTatlas.2016-03-20.gcs ../label/rh.aparc.DKTatlas.annot 

#-----------------------------------------
#@# Parcellation Stats 3 lh Tue Feb  7 15:56:38 EST 2017

 mris_anatomical_stats -th3 -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.DKTatlas.stats -b -a ../label/lh.aparc.DKTatlas.annot -c ../label/aparc.annot.DKTatlas.ctab MNI_FS2 lh white 

#-----------------------------------------
#@# Parcellation Stats 3 rh Tue Feb  7 15:57:47 EST 2017

 mris_anatomical_stats -th3 -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.DKTatlas.stats -b -a ../label/rh.aparc.DKTatlas.annot -c ../label/aparc.annot.DKTatlas.ctab MNI_FS2 rh white 

#-----------------------------------------
#@# WM/GM Contrast lh Tue Feb  7 15:58:56 EST 2017

 pctsurfcon --s MNI_FS2 --lh-only 

#-----------------------------------------
#@# WM/GM Contrast rh Tue Feb  7 15:59:04 EST 2017

 pctsurfcon --s MNI_FS2 --rh-only 

#-----------------------------------------
#@# Relabel Hypointensities Tue Feb  7 15:59:13 EST 2017

 mri_relabel_hypointensities aseg.presurf.mgz ../surf aseg.presurf.hypos.mgz 

#-----------------------------------------
#@# AParc-to-ASeg aparc Tue Feb  7 15:59:59 EST 2017

 mri_aparc2aseg --s MNI_FS2 --volmask --aseg aseg.presurf.hypos 

#-----------------------------------------
#@# AParc-to-ASeg a2009s Tue Feb  7 16:02:03 EST 2017

 mri_aparc2aseg --s MNI_FS2 --volmask --aseg aseg.presurf.hypos --annot aparc.a2009s 

#-----------------------------------------
#@# AParc-to-ASeg DKTatlas Tue Feb  7 16:04:08 EST 2017

 mri_aparc2aseg --s MNI_FS2 --volmask --aseg aseg.presurf.hypos --annot aparc.DKTatlas 

#-----------------------------------------
#@# APas-to-ASeg Tue Feb  7 16:06:13 EST 2017

 apas2aseg --i aparc+aseg.mgz --o aseg.mgz 

#--------------------------------------------
#@# ASeg Stats Tue Feb  7 16:06:23 EST 2017

 mri_segstats --seg mri/aseg.mgz --sum stats/aseg.stats --pv mri/norm.mgz --empty --brainmask mri/brainmask.mgz --brain-vol-from-seg --excludeid 0 --excl-ctxgmwm --supratent --subcortgray --in mri/norm.mgz --in-intensity-name norm --in-intensity-units MR --etiv --surf-wm-vol --surf-ctx-vol --totalgray --euler --ctab /ncf/nrg/sw/apps/freesurfer/6_2015_04_21/ASegStatsLUT.txt --subject MNI_FS2 

#-----------------------------------------
#@# WMParc Tue Feb  7 16:11:14 EST 2017

 mri_aparc2aseg --s MNI_FS2 --labelwm --hypo-as-wm --rip-unknown --volmask --o mri/wmparc.mgz --ctxseg aparc+aseg.mgz 


 mri_segstats --seg mri/wmparc.mgz --sum stats/wmparc.stats --pv mri/norm.mgz --excludeid 0 --brainmask mri/brainmask.mgz --in mri/norm.mgz --in-intensity-name norm --in-intensity-units MR --subject MNI_FS2 --surf-wm-vol --ctab /ncf/nrg/sw/apps/freesurfer/6_2015_04_21/WMParcStatsLUT.txt --etiv 

#--------------------------------------------
#@# BA_exvivo Labels lh Tue Feb  7 16:26:42 EST 2017

 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.BA1_exvivo.label --trgsubject MNI_FS2 --trglabel ./lh.BA1_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.BA2_exvivo.label --trgsubject MNI_FS2 --trglabel ./lh.BA2_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.BA3a_exvivo.label --trgsubject MNI_FS2 --trglabel ./lh.BA3a_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.BA3b_exvivo.label --trgsubject MNI_FS2 --trglabel ./lh.BA3b_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.BA4a_exvivo.label --trgsubject MNI_FS2 --trglabel ./lh.BA4a_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.BA4p_exvivo.label --trgsubject MNI_FS2 --trglabel ./lh.BA4p_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.BA6_exvivo.label --trgsubject MNI_FS2 --trglabel ./lh.BA6_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.BA44_exvivo.label --trgsubject MNI_FS2 --trglabel ./lh.BA44_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.BA45_exvivo.label --trgsubject MNI_FS2 --trglabel ./lh.BA45_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.V1_exvivo.label --trgsubject MNI_FS2 --trglabel ./lh.V1_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.V2_exvivo.label --trgsubject MNI_FS2 --trglabel ./lh.V2_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.MT_exvivo.label --trgsubject MNI_FS2 --trglabel ./lh.MT_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.entorhinal_exvivo.label --trgsubject MNI_FS2 --trglabel ./lh.entorhinal_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.perirhinal_exvivo.label --trgsubject MNI_FS2 --trglabel ./lh.perirhinal_exvivo.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.BA1_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./lh.BA1_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.BA2_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./lh.BA2_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.BA3a_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./lh.BA3a_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.BA3b_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./lh.BA3b_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.BA4a_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./lh.BA4a_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.BA4p_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./lh.BA4p_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.BA6_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./lh.BA6_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.BA44_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./lh.BA44_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.BA45_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./lh.BA45_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.V1_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./lh.V1_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.V2_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./lh.V2_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.MT_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./lh.MT_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.entorhinal_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./lh.entorhinal_exvivo.thresh.label --hemi lh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/lh.perirhinal_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./lh.perirhinal_exvivo.thresh.label --hemi lh --regmethod surface 


 mris_label2annot --s MNI_FS2 --hemi lh --ctab /ncf/nrg/sw/apps/freesurfer/6_2015_04_21/average/colortable_BA.txt --l lh.BA1_exvivo.label --l lh.BA2_exvivo.label --l lh.BA3a_exvivo.label --l lh.BA3b_exvivo.label --l lh.BA4a_exvivo.label --l lh.BA4p_exvivo.label --l lh.BA6_exvivo.label --l lh.BA44_exvivo.label --l lh.BA45_exvivo.label --l lh.V1_exvivo.label --l lh.V2_exvivo.label --l lh.MT_exvivo.label --l lh.entorhinal_exvivo.label --l lh.perirhinal_exvivo.label --a BA_exvivo --maxstatwinner --noverbose 


 mris_label2annot --s MNI_FS2 --hemi lh --ctab /ncf/nrg/sw/apps/freesurfer/6_2015_04_21/average/colortable_BA.txt --l lh.BA1_exvivo.thresh.label --l lh.BA2_exvivo.thresh.label --l lh.BA3a_exvivo.thresh.label --l lh.BA3b_exvivo.thresh.label --l lh.BA4a_exvivo.thresh.label --l lh.BA4p_exvivo.thresh.label --l lh.BA6_exvivo.thresh.label --l lh.BA44_exvivo.thresh.label --l lh.BA45_exvivo.thresh.label --l lh.V1_exvivo.thresh.label --l lh.V2_exvivo.thresh.label --l lh.MT_exvivo.thresh.label --l lh.entorhinal_exvivo.thresh.label --l lh.perirhinal_exvivo.thresh.label --a BA_exvivo.thresh --maxstatwinner --noverbose 


 mris_anatomical_stats -th3 -mgz -f ../stats/lh.BA_exvivo.stats -b -a ./lh.BA_exvivo.annot -c ./BA_exvivo.ctab MNI_FS2 lh white 


 mris_anatomical_stats -th3 -mgz -f ../stats/lh.BA_exvivo.thresh.stats -b -a ./lh.BA_exvivo.thresh.annot -c ./BA_exvivo.thresh.ctab MNI_FS2 lh white 

#--------------------------------------------
#@# BA_exvivo Labels rh Tue Feb  7 16:34:13 EST 2017

 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.BA1_exvivo.label --trgsubject MNI_FS2 --trglabel ./rh.BA1_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.BA2_exvivo.label --trgsubject MNI_FS2 --trglabel ./rh.BA2_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.BA3a_exvivo.label --trgsubject MNI_FS2 --trglabel ./rh.BA3a_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.BA3b_exvivo.label --trgsubject MNI_FS2 --trglabel ./rh.BA3b_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.BA4a_exvivo.label --trgsubject MNI_FS2 --trglabel ./rh.BA4a_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.BA4p_exvivo.label --trgsubject MNI_FS2 --trglabel ./rh.BA4p_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.BA6_exvivo.label --trgsubject MNI_FS2 --trglabel ./rh.BA6_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.BA44_exvivo.label --trgsubject MNI_FS2 --trglabel ./rh.BA44_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.BA45_exvivo.label --trgsubject MNI_FS2 --trglabel ./rh.BA45_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.V1_exvivo.label --trgsubject MNI_FS2 --trglabel ./rh.V1_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.V2_exvivo.label --trgsubject MNI_FS2 --trglabel ./rh.V2_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.MT_exvivo.label --trgsubject MNI_FS2 --trglabel ./rh.MT_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.entorhinal_exvivo.label --trgsubject MNI_FS2 --trglabel ./rh.entorhinal_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.perirhinal_exvivo.label --trgsubject MNI_FS2 --trglabel ./rh.perirhinal_exvivo.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.BA1_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./rh.BA1_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.BA2_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./rh.BA2_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.BA3a_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./rh.BA3a_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.BA3b_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./rh.BA3b_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.BA4a_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./rh.BA4a_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.BA4p_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./rh.BA4p_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.BA6_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./rh.BA6_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.BA44_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./rh.BA44_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.BA45_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./rh.BA45_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.V1_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./rh.V1_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.V2_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./rh.V2_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.MT_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./rh.MT_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.entorhinal_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./rh.entorhinal_exvivo.thresh.label --hemi rh --regmethod surface 


 mri_label2label --srcsubject fsaverage --srclabel /users/stadayon/Projects/Allen_Brain_Atlas/Data/MNI_reconall/fsaverage/label/rh.perirhinal_exvivo.thresh.label --trgsubject MNI_FS2 --trglabel ./rh.perirhinal_exvivo.thresh.label --hemi rh --regmethod surface 


 mris_label2annot --s MNI_FS2 --hemi rh --ctab /ncf/nrg/sw/apps/freesurfer/6_2015_04_21/average/colortable_BA.txt --l rh.BA1_exvivo.label --l rh.BA2_exvivo.label --l rh.BA3a_exvivo.label --l rh.BA3b_exvivo.label --l rh.BA4a_exvivo.label --l rh.BA4p_exvivo.label --l rh.BA6_exvivo.label --l rh.BA44_exvivo.label --l rh.BA45_exvivo.label --l rh.V1_exvivo.label --l rh.V2_exvivo.label --l rh.MT_exvivo.label --l rh.entorhinal_exvivo.label --l rh.perirhinal_exvivo.label --a BA_exvivo --maxstatwinner --noverbose 


 mris_label2annot --s MNI_FS2 --hemi rh --ctab /ncf/nrg/sw/apps/freesurfer/6_2015_04_21/average/colortable_BA.txt --l rh.BA1_exvivo.thresh.label --l rh.BA2_exvivo.thresh.label --l rh.BA3a_exvivo.thresh.label --l rh.BA3b_exvivo.thresh.label --l rh.BA4a_exvivo.thresh.label --l rh.BA4p_exvivo.thresh.label --l rh.BA6_exvivo.thresh.label --l rh.BA44_exvivo.thresh.label --l rh.BA45_exvivo.thresh.label --l rh.V1_exvivo.thresh.label --l rh.V2_exvivo.thresh.label --l rh.MT_exvivo.thresh.label --l rh.entorhinal_exvivo.thresh.label --l rh.perirhinal_exvivo.thresh.label --a BA_exvivo.thresh --maxstatwinner --noverbose 


 mris_anatomical_stats -th3 -mgz -f ../stats/rh.BA_exvivo.stats -b -a ./rh.BA_exvivo.annot -c ./BA_exvivo.ctab MNI_FS2 rh white 


 mris_anatomical_stats -th3 -mgz -f ../stats/rh.BA_exvivo.thresh.stats -b -a ./rh.BA_exvivo.thresh.annot -c ./BA_exvivo.thresh.ctab MNI_FS2 rh white 

