%
%  Run all steps of the grid simulation on a single participant
%

% set the path to the BIDS data directory
bids_rootPath   = 'F:\HMI\OpenData\BIDS_GridSim';

% set the subject, task and hemisphere to simulate on
bids_sub        = '09';
bids_task       = 'HandGesture';
hemi            = 'lh';

% In order to reproduce the results in the article we start at step 2, and
% set the samples filename to the previously generated sample-set
samplesFilename = 'sampleSet.mat';



%%
% retrieve the root path and make sure dependencies can be found
%
gridSimRoot = gridSimRootPath;
addpath([gridSimRoot, filesep, 'scripts']);



%%
%  Execute steps
%

% step 0 - Create a hull based on freesurfer parcellation areas
% s0_createHull(bids_rootPath, bids_sub, hemi);
%
% - skip because each hull was created with different parameters for
%   each participant and manually corrected. The subject's hull that
%   was used is included with the data on OSF


% step 1 - Generate random 3D sample-points on a given hull
% [samplesFilename] = s1_generateSamples(bids_rootPath, bids_sub, hemi);
%
% - skip because the sample-points are generated at random, and
%   the previously generated sample-sets on which the rest of the analysis
%   build are unique and included in the data on OSF


% step 2 - Classify the 3D sample-points using a searchlight with a specific radius
% to be continued...



