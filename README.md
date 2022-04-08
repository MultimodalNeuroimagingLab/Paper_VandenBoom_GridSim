# Functional MRI based simulations of ECoG grid configurations for optimal measurement of spatially distributed hand-gesture information
Work in Progress

Source-code that belong to article: Van Den Boom, M. A., Miller, K. J., Ramsey, N. F., & Hermes, D. (2021). Functional MRI based simulations of ECoG grid configurations for optimal measurement of spatially distributed hand-gesture information. Journal of neural engineering, 18(2), 026013.
Shared datasets can be found on OSF: https://osf.io/z6j3x/.

## Software requirements
- Matlab v2019 or higher
- SPM12 (make sure as path to Matlab and working)
- The code from this github page

## Data input
The data is expected to be organized according to the BIDS structure (https://bids-specification.readthedocs.io/).
The shared data from the article on OSF is structured accordingly.
To illustrate based on one participants from the open-source data (sub-09), we require:
```
/BIDS_GridSim/sub-09/func/sub-09_task-HandGesture_bold.nii
/BIDS_GridSim/sub-09/func/sub-09_task-HandGesture_bold.json
/BIDS_GridSim/sub-09/func/sub-09_task-HandGesture_events.tsv
/BIDS_GridSim/derivatives/freesurfer/sub-09/surf/lh.pial
/BIDS_GridSim/derivatives/freesurfer/sub-09/surf/rh.pial
/BIDS_GridSim/derivatives/freesurfer/sub-09/mri/orig.mgz
```

Some required files were derived from the freesurfer output and need to be organized accordingly:
```
/BIDS_GridSim/derivatives/freesurfer/sub-09/surf/lh.pial.gii
/BIDS_GridSim/derivatives/freesurfer/sub-09/surf/rh.pial.gii
/BIDS_GridSim/derivatives/freesurfer/sub-09/mri/orig.mat
/BIDS_GridSim/derivatives/freesurfer/sub-09/label/lh.aparc.Ext.annot
/BIDS_GridSim/derivatives/freesurfer/sub-09/label/rh.aparc.Ext.annot
```

Other required files need to be generated, such as the grey matter volume mask (resliced to functional, native space):
```
/BIDS_GridSim/derivatives/gm_masks/sub-09/func_fs_mask.nii
```

To reproduce the results in the article, we have provided the output of the first two steps, yielding the following files:
```
/BIDS_GridSim/derivatives/lh_simulations/sub-09/sampleSet.mat
/BIDS_GridSim/derivatives/lh_simulations/sub-09/lh_ext_hull.gii
/BIDS_GridSim/derivatives/lh_simulations/sub-09/lh_sensorymotor_hull.gii
```