%
%  Run all steps of the grid simulation on a single participant
%

% limit the number of threads to process data (less threads also means less memory)
numThreads = 1;

% set the path to the BIDS data directory
bids_rootPath   = 'Y:\OpenData\BIDS_GridSim';

% set the subject, task and hemisphere to simulate on
bids_sub        = '09';
bids_task       = 'HandGesture';
hemi            = 'lh';

% In order to reproduce the results in the article we start at step 2, and
% set the samples filename to the previously generated sample-set
hullFilename = 'lh_ext_hull.gii';               % <--- delete (and perform step 0 below) if working on own dataset
samplesFilename = 'sampleSet.mat';              % <--- delete (and perform step 1 below) if working on own dataset

% set the classification configuration
classConfig = [];
classConfig.isolateTrials = [];
classConfig.silent = 'yes';
classConfig.hrfMethod = 'new';
classConfig.earlyFeatureSelectionMethod = 'tmaphighest';
classConfig.skipStartVolumes = 5;
classConfig.minFeatures = 2;                                            % the minimum number of features
classConfig.minVoxelMeanForFeature = 15000;                             % voxels with a signal over time below this minimum will be removed (raw input values, pre-normalization and detrending) - the value under which functional voxel are considered to be outside of the brain
classConfig.timeLockActive = 'search_highestbold_hrfpeaks';             % search_highestbold_hrfpeaks, static
classConfig.timeLockActiveSearchScansBeforePeak = 1;
classConfig.timeLockActiveSearchScansAfterPeak = 1;
classConfig.timeLockActiveVolumeOffsets = [0];
classConfig.normAndDetrend = {'NormVoxelsPercSignal', 'SPMDetrend'};
classConfig.crossValidation = 'LeaveOneOut';
classConfig.lateFeatureSelectionMethod = 'none';
classConfig.featureTimeMethod = 'average';                              % how features over time are used, average or serialized
classConfig.classificationMethod = 'svm';




%%
% retrieve the root path and make sure dependencies can be found
%
gridSimRoot = gridSimRootPath;
addpath([gridSimRoot, filesep, 'scripts']);



%%
%  Execute steps
%

% step 0 - Create a hull based on freesurfer parcellation areas
% hullFilename = s0_createHull(bids_rootPath, bids_sub, hemi);
%
% - skip because each hull was created with different parameters for
%   each participant and manually corrected. The subject's hull that
%   was used is included with the data on OSF


% step 1 - Generate random 3D sample-points on a given hull
% samplesFilename = s1_generateSamples(bids_rootPath, bids_sub, hemi);
%
% - skip because the sample-points are generated at random, and
%   the previously generated sample-sets on which the rest of the analysis
%   build are unique and included in the data on OSF


% step 2 - Peform searchlight classification on each of the 3D sample-points with a specific radius
s2_searchlight(bids_rootPath, bids_sub, bids_task, hemi, samplesFilename, 7, classConfig, numThreads);

%... more

