

#---------------------------------
# New invocation of recon-all Wed Aug  1 23:25:02 EDT 2018 

 mri_convert /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/RAW/Broad_82_YZ/3/DICOM/1.2.840.113619.6.408.331269118540224283087655827658524882091-3-100-11iemj5.dcm /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects/Broad_82_YZ/mri/orig/001.mgz 

#--------------------------------------------
#@# MotionCor Wed Aug  1 23:25:16 EDT 2018

 cp /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects/Broad_82_YZ/mri/orig/001.mgz /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects/Broad_82_YZ/mri/rawavg.mgz 


 mri_convert /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects/Broad_82_YZ/mri/rawavg.mgz /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects/Broad_82_YZ/mri/orig.mgz --conform 


 mri_add_xform_to_header -c /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects/Broad_82_YZ/mri/transforms/talairach.xfm /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects/Broad_82_YZ/mri/orig.mgz /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects/Broad_82_YZ/mri/orig.mgz 

#--------------------------------------------
#@# Talairach Wed Aug  1 23:25:31 EDT 2018

 mri_nu_correct.mni --no-rescale --i orig.mgz --o orig_nu.mgz --n 1 --proto-iters 1000 --distance 50 


 talairach_avi --i orig_nu.mgz --xfm transforms/talairach.auto.xfm 

talairach_avi log file is transforms/talairach_avi.log...

 cp transforms/talairach.auto.xfm transforms/talairach.xfm 

#--------------------------------------------
#@# Talairach Failure Detection Wed Aug  1 23:28:49 EDT 2018

 talairach_afd -T 0.005 -xfm transforms/talairach.xfm 


 awk -f /ncf/nrg/sw/apps/freesurfer/6.0.0/bin/extract_talairach_avi_QA.awk /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects/Broad_82_YZ/mri/transforms/talairach_avi.log 


 tal_QC_AZS /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects/Broad_82_YZ/mri/transforms/talairach_avi.log 

#--------------------------------------------
#@# Nu Intensity Correction Wed Aug  1 23:28:50 EDT 2018

 mri_nu_correct.mni --i orig.mgz --o nu.mgz --uchar transforms/talairach.xfm --n 2 


 mri_add_xform_to_header -c /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects/Broad_82_YZ/mri/transforms/talairach.xfm nu.mgz nu.mgz 

#--------------------------------------------
#@# Intensity Normalization Wed Aug  1 23:32:02 EDT 2018

 mri_normalize -g 1 -mprage nu.mgz T1.mgz 

#--------------------------------------------
#@# Skull Stripping Wed Aug  1 23:35:51 EDT 2018

 mri_em_register -rusage /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects/Broad_82_YZ/touch/rusage.mri_em_register.skull.dat -skull nu.mgz /ncf/nrg/sw/apps/freesurfer/6.0.0/average/RB_all_withskull_2016-05-10.vc700.gca transforms/talairach_with_skull.lta 


 mri_watershed -rusage /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects/Broad_82_YZ/touch/rusage.mri_watershed.dat -T1 -brain_atlas /ncf/nrg/sw/apps/freesurfer/6.0.0/average/RB_all_withskull_2016-05-10.vc700.gca transforms/talairach_with_skull.lta T1.mgz brainmask.auto.mgz 


 cp brainmask.auto.mgz brainmask.mgz 

#-------------------------------------
#@# EM Registration Thu Aug  2 00:03:54 EDT 2018

 mri_em_register -rusage /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects/Broad_82_YZ/touch/rusage.mri_em_register.dat -uns 3 -mask brainmask.mgz nu.mgz /ncf/nrg/sw/apps/freesurfer/6.0.0/average/RB_all_2016-05-10.vc700.gca transforms/talairach.lta 

#--------------------------------------
#@# CA Normalize Thu Aug  2 00:30:45 EDT 2018

 mri_ca_normalize -c ctrl_pts.mgz -mask brainmask.mgz nu.mgz /ncf/nrg/sw/apps/freesurfer/6.0.0/average/RB_all_2016-05-10.vc700.gca transforms/talairach.lta norm.mgz 

#--------------------------------------
#@# CA Reg Thu Aug  2 00:33:37 EDT 2018

 mri_ca_register -rusage /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects/Broad_82_YZ/touch/rusage.mri_ca_register.dat -nobigventricles -T transforms/talairach.lta -align-after -mask brainmask.mgz norm.mgz /ncf/nrg/sw/apps/freesurfer/6.0.0/average/RB_all_2016-05-10.vc700.gca transforms/talairach.m3z 

#--------------------------------------
#@# SubCort Seg Thu Aug  2 04:51:30 EDT 2018

 mri_ca_label -relabel_unlikely 9 .3 -prior 0.5 -align norm.mgz transforms/talairach.m3z /ncf/nrg/sw/apps/freesurfer/6.0.0/average/RB_all_2016-05-10.vc700.gca aseg.auto_noCCseg.mgz 


 mri_cc -aseg aseg.auto_noCCseg.mgz -o aseg.auto.mgz -lta /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects/Broad_82_YZ/mri/transforms/cc_up.lta Broad_82_YZ 

#--------------------------------------
#@# Merge ASeg Thu Aug  2 06:27:46 EDT 2018

 cp aseg.auto.mgz aseg.presurf.mgz 

#--------------------------------------------
#@# Intensity Normalization2 Thu Aug  2 06:27:46 EDT 2018

 mri_normalize -mprage -aseg aseg.presurf.mgz -mask brainmask.mgz norm.mgz brain.mgz 

#--------------------------------------------
#@# Mask BFS Thu Aug  2 06:33:45 EDT 2018

 mri_mask -T 5 brain.mgz brainmask.mgz brain.finalsurfs.mgz 

#--------------------------------------------
#@# WM Segmentation Thu Aug  2 06:33:48 EDT 2018

 mri_segment -mprage brain.mgz wm.seg.mgz 


 mri_edit_wm_with_aseg -keep-in wm.seg.mgz brain.mgz aseg.presurf.mgz wm.asegedit.mgz 


 mri_pretess wm.asegedit.mgz wm norm.mgz wm.mgz 

#--------------------------------------------
#@# Fill Thu Aug  2 06:38:02 EDT 2018

 mri_fill -a ../scripts/ponscc.cut.log -xform transforms/talairach.lta -segmentation aseg.auto_noCCseg.mgz wm.mgz filled.mgz 

#--------------------------------------------
#@# Tessellate lh Thu Aug  2 06:39:14 EDT 2018

 mri_pretess ../mri/filled.mgz 255 ../mri/norm.mgz ../mri/filled-pretess255.mgz 


 mri_tessellate ../mri/filled-pretess255.mgz 255 ../surf/lh.orig.nofix 


 rm -f ../mri/filled-pretess255.mgz 


 mris_extract_main_component ../surf/lh.orig.nofix ../surf/lh.orig.nofix 

#--------------------------------------------
#@# Tessellate rh Thu Aug  2 06:39:26 EDT 2018

 mri_pretess ../mri/filled.mgz 127 ../mri/norm.mgz ../mri/filled-pretess127.mgz 


 mri_tessellate ../mri/filled-pretess127.mgz 127 ../surf/rh.orig.nofix 


 rm -f ../mri/filled-pretess127.mgz 


 mris_extract_main_component ../surf/rh.orig.nofix ../surf/rh.orig.nofix 

#--------------------------------------------
#@# Smooth1 lh Thu Aug  2 06:39:39 EDT 2018

 mris_smooth -nw -seed 1234 ../surf/lh.orig.nofix ../surf/lh.smoothwm.nofix 

#--------------------------------------------
#@# Smooth1 rh Thu Aug  2 06:39:49 EDT 2018

 mris_smooth -nw -seed 1234 ../surf/rh.orig.nofix ../surf/rh.smoothwm.nofix 

#--------------------------------------------
#@# Inflation1 lh Thu Aug  2 06:39:59 EDT 2018

 mris_inflate -no-save-sulc ../surf/lh.smoothwm.nofix ../surf/lh.inflated.nofix 

#--------------------------------------------
#@# Inflation1 rh Thu Aug  2 06:41:03 EDT 2018

 mris_inflate -no-save-sulc ../surf/rh.smoothwm.nofix ../surf/rh.inflated.nofix 

#--------------------------------------------
#@# QSphere lh Thu Aug  2 06:42:04 EDT 2018

 mris_sphere -q -seed 1234 ../surf/lh.inflated.nofix ../surf/lh.qsphere.nofix 

#--------------------------------------------
#@# QSphere rh Thu Aug  2 06:48:23 EDT 2018

 mris_sphere -q -seed 1234 ../surf/rh.inflated.nofix ../surf/rh.qsphere.nofix 

#--------------------------------------------
#@# Fix Topology Copy lh Thu Aug  2 06:54:31 EDT 2018

 cp ../surf/lh.orig.nofix ../surf/lh.orig 


 cp ../surf/lh.inflated.nofix ../surf/lh.inflated 

#--------------------------------------------
#@# Fix Topology Copy rh Thu Aug  2 06:54:31 EDT 2018

 cp ../surf/rh.orig.nofix ../surf/rh.orig 


 cp ../surf/rh.inflated.nofix ../surf/rh.inflated 

#@# Fix Topology lh Thu Aug  2 06:54:31 EDT 2018

 mris_fix_topology -rusage /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects/Broad_82_YZ/touch/rusage.mris_fix_topology.lh.dat -mgz -sphere qsphere.nofix -ga -seed 1234 Broad_82_YZ lh 

#@# Fix Topology rh Thu Aug  2 07:50:57 EDT 2018

 mris_fix_topology -rusage /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects/Broad_82_YZ/touch/rusage.mris_fix_topology.rh.dat -mgz -sphere qsphere.nofix -ga -seed 1234 Broad_82_YZ rh 


 mris_euler_number ../surf/lh.orig 


 mris_euler_number ../surf/rh.orig 


 mris_remove_intersection ../surf/lh.orig ../surf/lh.orig 


 rm ../surf/lh.inflated 


 mris_remove_intersection ../surf/rh.orig ../surf/rh.orig 


 rm ../surf/rh.inflated 

#--------------------------------------------
#@# Make White Surf lh Thu Aug  2 08:46:32 EDT 2018

 mris_make_surfaces -aseg ../mri/aseg.presurf -white white.preaparc -noaparc -whiteonly -mgz -T1 brain.finalsurfs Broad_82_YZ lh 

#--------------------------------------------
#@# Make White Surf rh Thu Aug  2 08:54:25 EDT 2018

 mris_make_surfaces -aseg ../mri/aseg.presurf -white white.preaparc -noaparc -whiteonly -mgz -T1 brain.finalsurfs Broad_82_YZ rh 

#--------------------------------------------
#@# Smooth2 lh Thu Aug  2 09:01:55 EDT 2018

 mris_smooth -n 3 -nw -seed 1234 ../surf/lh.white.preaparc ../surf/lh.smoothwm 

#--------------------------------------------
#@# Smooth2 rh Thu Aug  2 09:02:04 EDT 2018

 mris_smooth -n 3 -nw -seed 1234 ../surf/rh.white.preaparc ../surf/rh.smoothwm 

#--------------------------------------------
#@# Inflation2 lh Thu Aug  2 09:02:13 EDT 2018

 mris_inflate -rusage /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects/Broad_82_YZ/touch/rusage.mris_inflate.lh.dat ../surf/lh.smoothwm ../surf/lh.inflated 

#--------------------------------------------
#@# Inflation2 rh Thu Aug  2 09:03:04 EDT 2018

 mris_inflate -rusage /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects/Broad_82_YZ/touch/rusage.mris_inflate.rh.dat ../surf/rh.smoothwm ../surf/rh.inflated 

#--------------------------------------------
#@# Curv .H and .K lh Thu Aug  2 09:03:54 EDT 2018

 mris_curvature -w lh.white.preaparc 


 mris_curvature -thresh .999 -n -a 5 -w -distances 10 10 lh.inflated 

#--------------------------------------------
#@# Curv .H and .K rh Thu Aug  2 09:06:08 EDT 2018

 mris_curvature -w rh.white.preaparc 


 mris_curvature -thresh .999 -n -a 5 -w -distances 10 10 rh.inflated 


#-----------------------------------------
#@# Curvature Stats lh Thu Aug  2 09:08:18 EDT 2018

 mris_curvature_stats -m --writeCurvatureFiles -G -o ../stats/lh.curv.stats -F smoothwm Broad_82_YZ lh curv sulc 


#-----------------------------------------
#@# Curvature Stats rh Thu Aug  2 09:08:24 EDT 2018

 mris_curvature_stats -m --writeCurvatureFiles -G -o ../stats/rh.curv.stats -F smoothwm Broad_82_YZ rh curv sulc 

#--------------------------------------------
#@# Sphere lh Thu Aug  2 09:08:30 EDT 2018

 mris_sphere -rusage /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects/Broad_82_YZ/touch/rusage.mris_sphere.lh.dat -seed 1234 ../surf/lh.inflated ../surf/lh.sphere 

#--------------------------------------------
#@# Sphere rh Thu Aug  2 10:07:29 EDT 2018

 mris_sphere -rusage /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects/Broad_82_YZ/touch/rusage.mris_sphere.rh.dat -seed 1234 ../surf/rh.inflated ../surf/rh.sphere 

#--------------------------------------------
#@# Surf Reg lh Thu Aug  2 11:14:08 EDT 2018

 mris_register -curv -rusage /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects/Broad_82_YZ/touch/rusage.mris_register.lh.dat ../surf/lh.sphere /ncf/nrg/sw/apps/freesurfer/6.0.0/average/lh.folding.atlas.acfb40.noaparc.i12.2016-08-02.tif ../surf/lh.sphere.reg 

#--------------------------------------------
#@# Surf Reg rh Thu Aug  2 12:33:08 EDT 2018

 mris_register -curv -rusage /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects/Broad_82_YZ/touch/rusage.mris_register.rh.dat ../surf/rh.sphere /ncf/nrg/sw/apps/freesurfer/6.0.0/average/rh.folding.atlas.acfb40.noaparc.i12.2016-08-02.tif ../surf/rh.sphere.reg 

#--------------------------------------------
#@# Jacobian white lh Thu Aug  2 13:55:02 EDT 2018

 mris_jacobian ../surf/lh.white.preaparc ../surf/lh.sphere.reg ../surf/lh.jacobian_white 

#--------------------------------------------
#@# Jacobian white rh Thu Aug  2 13:55:05 EDT 2018

 mris_jacobian ../surf/rh.white.preaparc ../surf/rh.sphere.reg ../surf/rh.jacobian_white 

#--------------------------------------------
#@# AvgCurv lh Thu Aug  2 13:55:08 EDT 2018

 mrisp_paint -a 5 /ncf/nrg/sw/apps/freesurfer/6.0.0/average/lh.folding.atlas.acfb40.noaparc.i12.2016-08-02.tif#6 ../surf/lh.sphere.reg ../surf/lh.avg_curv 

#--------------------------------------------
#@# AvgCurv rh Thu Aug  2 13:55:10 EDT 2018

 mrisp_paint -a 5 /ncf/nrg/sw/apps/freesurfer/6.0.0/average/rh.folding.atlas.acfb40.noaparc.i12.2016-08-02.tif#6 ../surf/rh.sphere.reg ../surf/rh.avg_curv 

#-----------------------------------------
#@# Cortical Parc lh Thu Aug  2 13:55:13 EDT 2018

 mris_ca_label -l ../label/lh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 Broad_82_YZ lh ../surf/lh.sphere.reg /ncf/nrg/sw/apps/freesurfer/6.0.0/average/lh.DKaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs ../label/lh.aparc.annot 

#-----------------------------------------
#@# Cortical Parc rh Thu Aug  2 13:55:36 EDT 2018

 mris_ca_label -l ../label/rh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 Broad_82_YZ rh ../surf/rh.sphere.reg /ncf/nrg/sw/apps/freesurfer/6.0.0/average/rh.DKaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs ../label/rh.aparc.annot 

#--------------------------------------------
#@# Make Pial Surf lh Thu Aug  2 13:55:58 EDT 2018

 mris_make_surfaces -orig_white white.preaparc -orig_pial white.preaparc -aseg ../mri/aseg.presurf -mgz -T1 brain.finalsurfs Broad_82_YZ lh 

#--------------------------------------------
#@# Make Pial Surf rh Thu Aug  2 14:19:49 EDT 2018

 mris_make_surfaces -orig_white white.preaparc -orig_pial white.preaparc -aseg ../mri/aseg.presurf -mgz -T1 brain.finalsurfs Broad_82_YZ rh 

#--------------------------------------------
#@# Surf Volume lh Thu Aug  2 14:42:53 EDT 2018
#--------------------------------------------
#@# Surf Volume rh Thu Aug  2 14:42:58 EDT 2018
#--------------------------------------------
#@# Cortical ribbon mask Thu Aug  2 14:43:03 EDT 2018

 mris_volmask --aseg_name aseg.presurf --label_left_white 2 --label_left_ribbon 3 --label_right_white 41 --label_right_ribbon 42 --save_ribbon Broad_82_YZ 

#-----------------------------------------
#@# Parcellation Stats lh Thu Aug  2 14:57:33 EDT 2018

 mris_anatomical_stats -th3 -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.stats -b -a ../label/lh.aparc.annot -c ../label/aparc.annot.ctab Broad_82_YZ lh white 


 mris_anatomical_stats -th3 -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.pial.stats -b -a ../label/lh.aparc.annot -c ../label/aparc.annot.ctab Broad_82_YZ lh pial 

#-----------------------------------------
#@# Parcellation Stats rh Thu Aug  2 14:59:50 EDT 2018

 mris_anatomical_stats -th3 -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.stats -b -a ../label/rh.aparc.annot -c ../label/aparc.annot.ctab Broad_82_YZ rh white 


 mris_anatomical_stats -th3 -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.pial.stats -b -a ../label/rh.aparc.annot -c ../label/aparc.annot.ctab Broad_82_YZ rh pial 

#-----------------------------------------
#@# Cortical Parc 2 lh Thu Aug  2 15:02:01 EDT 2018

 mris_ca_label -l ../label/lh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 Broad_82_YZ lh ../surf/lh.sphere.reg /ncf/nrg/sw/apps/freesurfer/6.0.0/average/lh.CDaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs ../label/lh.aparc.a2009s.annot 

#-----------------------------------------
#@# Cortical Parc 2 rh Thu Aug  2 15:02:31 EDT 2018

 mris_ca_label -l ../label/rh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 Broad_82_YZ rh ../surf/rh.sphere.reg /ncf/nrg/sw/apps/freesurfer/6.0.0/average/rh.CDaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs ../label/rh.aparc.a2009s.annot 

#-----------------------------------------
#@# Parcellation Stats 2 lh Thu Aug  2 15:03:01 EDT 2018

 mris_anatomical_stats -th3 -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.a2009s.stats -b -a ../label/lh.aparc.a2009s.annot -c ../label/aparc.annot.a2009s.ctab Broad_82_YZ lh white 

#-----------------------------------------
#@# Parcellation Stats 2 rh Thu Aug  2 15:04:09 EDT 2018

 mris_anatomical_stats -th3 -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.a2009s.stats -b -a ../label/rh.aparc.a2009s.annot -c ../label/aparc.annot.a2009s.ctab Broad_82_YZ rh white 

#-----------------------------------------
#@# Cortical Parc 3 lh Thu Aug  2 15:05:16 EDT 2018

 mris_ca_label -l ../label/lh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 Broad_82_YZ lh ../surf/lh.sphere.reg /ncf/nrg/sw/apps/freesurfer/6.0.0/average/lh.DKTaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs ../label/lh.aparc.DKTatlas.annot 

#-----------------------------------------
#@# Cortical Parc 3 rh Thu Aug  2 15:05:39 EDT 2018

 mris_ca_label -l ../label/rh.cortex.label -aseg ../mri/aseg.presurf.mgz -seed 1234 Broad_82_YZ rh ../surf/rh.sphere.reg /ncf/nrg/sw/apps/freesurfer/6.0.0/average/rh.DKTaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs ../label/rh.aparc.DKTatlas.annot 

#-----------------------------------------
#@# Parcellation Stats 3 lh Thu Aug  2 15:06:03 EDT 2018

 mris_anatomical_stats -th3 -mgz -cortex ../label/lh.cortex.label -f ../stats/lh.aparc.DKTatlas.stats -b -a ../label/lh.aparc.DKTatlas.annot -c ../label/aparc.annot.DKTatlas.ctab Broad_82_YZ lh white 

#-----------------------------------------
#@# Parcellation Stats 3 rh Thu Aug  2 15:07:11 EDT 2018

 mris_anatomical_stats -th3 -mgz -cortex ../label/rh.cortex.label -f ../stats/rh.aparc.DKTatlas.stats -b -a ../label/rh.aparc.DKTatlas.annot -c ../label/aparc.annot.DKTatlas.ctab Broad_82_YZ rh white 

#-----------------------------------------
#@# WM/GM Contrast lh Thu Aug  2 15:08:15 EDT 2018

 pctsurfcon --s Broad_82_YZ --lh-only 

#-----------------------------------------
#@# WM/GM Contrast rh Thu Aug  2 15:08:26 EDT 2018

 pctsurfcon --s Broad_82_YZ --rh-only 

#-----------------------------------------
#@# Relabel Hypointensities Thu Aug  2 15:08:35 EDT 2018

 mri_relabel_hypointensities aseg.presurf.mgz ../surf aseg.presurf.hypos.mgz 

#-----------------------------------------
#@# AParc-to-ASeg aparc Thu Aug  2 15:09:12 EDT 2018

 mri_aparc2aseg --s Broad_82_YZ --volmask --aseg aseg.presurf.hypos --relabel mri/norm.mgz mri/transforms/talairach.m3z /ncf/nrg/sw/apps/freesurfer/6.0.0/average/RB_all_2016-05-10.vc700.gca mri/aseg.auto_noCCseg.label_intensities.txt 

#-----------------------------------------
#@# AParc-to-ASeg a2009s Thu Aug  2 15:16:32 EDT 2018

 mri_aparc2aseg --s Broad_82_YZ --volmask --aseg aseg.presurf.hypos --relabel mri/norm.mgz mri/transforms/talairach.m3z /ncf/nrg/sw/apps/freesurfer/6.0.0/average/RB_all_2016-05-10.vc700.gca mri/aseg.auto_noCCseg.label_intensities.txt --a2009s 

#-----------------------------------------
#@# AParc-to-ASeg DKTatlas Thu Aug  2 15:23:52 EDT 2018

 mri_aparc2aseg --s Broad_82_YZ --volmask --aseg aseg.presurf.hypos --relabel mri/norm.mgz mri/transforms/talairach.m3z /ncf/nrg/sw/apps/freesurfer/6.0.0/average/RB_all_2016-05-10.vc700.gca mri/aseg.auto_noCCseg.label_intensities.txt --annot aparc.DKTatlas --o mri/aparc.DKTatlas+aseg.mgz 

#-----------------------------------------
#@# APas-to-ASeg Thu Aug  2 15:31:13 EDT 2018

 apas2aseg --i aparc+aseg.mgz --o aseg.mgz 

#--------------------------------------------
#@# ASeg Stats Thu Aug  2 15:31:22 EDT 2018

 mri_segstats --seg mri/aseg.mgz --sum stats/aseg.stats --pv mri/norm.mgz --empty --brainmask mri/brainmask.mgz --brain-vol-from-seg --excludeid 0 --excl-ctxgmwm --supratent --subcortgray --in mri/norm.mgz --in-intensity-name norm --in-intensity-units MR --etiv --surf-wm-vol --surf-ctx-vol --totalgray --euler --ctab /ncf/nrg/sw/apps/freesurfer/6.0.0/ASegStatsLUT.txt --subject Broad_82_YZ 

#-----------------------------------------
#@# WMParc Thu Aug  2 15:35:32 EDT 2018

 mri_aparc2aseg --s Broad_82_YZ --labelwm --hypo-as-wm --rip-unknown --volmask --o mri/wmparc.mgz --ctxseg aparc+aseg.mgz 


 mri_segstats --seg mri/wmparc.mgz --sum stats/wmparc.stats --pv mri/norm.mgz --excludeid 0 --brainmask mri/brainmask.mgz --in mri/norm.mgz --in-intensity-name norm --in-intensity-units MR --subject Broad_82_YZ --surf-wm-vol --ctab /ncf/nrg/sw/apps/freesurfer/6.0.0/WMParcStatsLUT.txt --etiv 

INFO: fsaverage subject does not have uptodate *rhinal labels!

 cd /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects; rm -Rf fsaverage; cd - 

INFO: Creating symlink to fsaverage subject...

 cd /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects; ln -s /ncf/nrg/sw/apps/freesurfer/6.0.0/subjects/fsaverage; cd - 

#--------------------------------------------
#@# BA_exvivo Labels lh Thu Aug  2 15:49:37 EDT 2018

 mri_label2label --srcsubject fsaverage --srclabel /net/rcss9/srv/export/ncf_pascual-leone_cnbs/share_root/lab/broad/baseline/subjects/fsaverage/label/lh.BA1_exvivo.label --trgsubject Broad_82_YZ --trglabel ./lh.BA1_exvivo.label --hemi lh --regmethod surface 

