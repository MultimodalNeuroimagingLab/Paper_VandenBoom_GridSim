%
%  Run all steps of the grid simulation on a single participant
%

% set the path to the BIDS data directory
bids_rootPath   = 'F:\HMI\OpenData\BIDS_GridSim';

% set the subject, task and hemisphere to simulate on
bids_sub        = '09';
bids_task       = 'HandGesture';
hemi            = 'lh';




%%
% retrieve the root path and make sure dependencies can be found
%
gridSimRoot = gridSimRootPath;
addpath([gridSimRoot, filesep, 'scripts']);



%%
%  Execute steps
%

% step 1 - Create a hull based on freesurfer parcellation areas
% s1_createHull(bids_rootPath, bids_sub, hemi)
%
% - skip because each hull was created with different parameters for
%   each participant and manually corrected. The subject's hull that
%   was used is included with the data on OSF

% step 2 - Generate

% to be continued...


