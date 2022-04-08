# Paper Grid simulation code. 
Work in Progress

Scripts that belong to the paper


## MATLAB preparation
- SPM12 is required for some of the script, make sure SPM12 is added as a path and working
- Download the code from this github page


## Dataset preparation
The data should be organized according to the BIDS structure.
The processing scripts expect the BIDS-structure and the shared data from the article on OSF is structured accordingly.
To illustrate based on one participants from the open-source data (sub-09), we require:
/BIDS_GridSim/sub-09/func/sub-09_task-HandGesture_bold.nii
/BIDS_GridSim/sub-09/func/sub-09_task-HandGesture_bold.json
/BIDS_GridSim/sub-09/func/sub-09_task-HandGesture_events.tsv
/BIDS_GridSim/derivatives/freesurfer/sub-09/surf/lh.pial
/BIDS_GridSim/derivatives/freesurfer/sub-09/surf/rh.pial
/BIDS_GridSim/derivatives/freesurfer/sub-09/mri/orig.mgz

Some required files were derived from the freesurfer output and need to be organized accordingly:
/BIDS_GridSim/derivatives/freesurfer/sub-09/surf/lh.pial.gii
/BIDS_GridSim/derivatives/freesurfer/sub-09/surf/rh.pial.gii
/BIDS_GridSim/derivatives/freesurfer/sub-09/mri/orig.mat
/BIDS_GridSim/derivatives/freesurfer/sub-09/label/lh.aparc.Ext.annot
/BIDS_GridSim/derivatives/freesurfer/sub-09/label/rh.aparc.Ext.annot

Other required files need to be generated, such as the grey matter volume mask (resliced to functional, native space):
/BIDS_GridSim/derivatives/gm_masks/sub-09/func_fs_mask.nii

To reproduce the results in the article, we have provided the output of the first two steps, yielding the following files:
/BIDS_GridSim/derivatives/lh_simulations/sub-09/sampleSet.mat
/BIDS_GridSim/derivatives/lh_simulations/sub-09/lh_ext_hull.gii
/BIDS_GridSim/derivatives/lh_simulations/sub-09/lh_sensorymotor_hull.gii
