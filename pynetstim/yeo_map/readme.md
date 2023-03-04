### Morphing Yeo2011 networks into subject

Yeo2011 maps are 7 or 17 networks in the fsaverage common template. One way to use these maps is to morph them back into subject space. We did this for a project that started in 2017. We particularly were interested in stimulating networks using TMS while recording the propagation of the perturbation by EEG. To do that, we morphed the Yeo2011 maps into subject space and used the created masks for source analysis of EEG data. 

Here is how you can use this code to morph Yeo2011 maps into individual space: 

1. Run Freesurfer recon-all for your subject
2. run this command: " bash yeo2subject.sh <subject> 7 ". 7 refers to 7 network maps. Note that SUBJECTS_DIR should be defined prior to running this code and the subject's recon-all results should be there. 
	
