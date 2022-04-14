%
%  Run all steps of the grid simulation on a single participant
%

% limit the number of threads to process data (less threads also means less memory)
numThreads = 1;

% set the path to the BIDS data directory
bids_rootPath   = 'Y:\WorkData\BIDS_GridSim';

% set the subject, task and hemisphere to simulate on
bids_sub        = '09';
bids_task       = 'HandGesture';
hemi            = 'lh';

% output path for the simulations
bids_simPath = fullfile(bids_rootPath, 'derivatives', [hemi, '_simulations'], ['sub-' bids_sub]);
if ~exist(bids_simPath, 'dir'),     mkdir(bids_simPath);    end

% In order to reproduce the results in the article we start at step 2, and
% set the samples filename to the previously generated sample-set
hullFile = 'lh_ext_hull.gii';                                   % <--- delete (and perform step 0 below) if working on own dataset
samplesFile_s1_out = 'sampleSet';                               % <--- delete (and perform step 1 below) if working on own dataset
load(fullfile(bids_simPath, [samplesFile_s1_out, '.mat']));     % <--- delete (and perform step 1 below) if working on own dataset

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
addpath([gridSimRoot, filesep, 'functions']);



%%
% Step 0 - Create a hull based on freesurfer parcellation areas
%
%{

%
% - skip because each hull was created with different parameters for
%   each participant and manually corrected. The subject's hull that
%   was used is included with the data on OSF

%gHull = s0_createHull(  bids_rootPath, bids_sub, hemi, [hemi, '.aparc.annot'], { 'precentral', 'postcentral' });
gHull = s0_createHull(  bids_rootPath, bids_sub, hemi, [hemi, '.aparc.Ext.annot'], ...
                        { 'superiorfrontal', 'caudalmiddlefrontal', 'parsopercularis', 'precentral', ...
                          'postcentral', 'superiorparietal', 'supramarginal', 'inferiorparietal'  });

% save
hullFile = [hemi, '_ext_hull__.gii'];
save(gHull, fullfile(bids_simPath, hullFile))
                          
%}



%% 
%  Step 1 - Generate random 3D sample-points on a given hull
%
%{

% - skip because the sample-points are generated at random, and to reproduce the
%   output for this subject we need to use exactly the same input data, which
%   are the previously generated sample-sets, included in the data on OSF

hullPath = fullfile(bids_simPath, [hemi, '_ext_hull.gii']);
[SS, suffix] = s1_generateSamples(hullPath);

samplesFile_s1_out = ['sampleSet-', suffix];
save(fullfile(bids_simTaskPath, [samplesFile_s1_out, '.mat']), 'SS');

%}



%%
%  Step 2 - Peform searchlight classification on each of the 3D sample-points with a specific radius
%
[SS, suffix]                = s2_searchlight(   SS, bids_rootPath, bids_sub, bids_task, hemi, ...
                                                7, ...                                                      % <-- searchlight radius
                                                classConfig, numThreads);

% save results
bids_simTaskPath            = fullfile(bids_simPath, bids_task);
if ~exist(bids_simTaskPath, 'dir'),     mkdir(bids_simTaskPath);    end
samplesFile_s2_out          = [samplesFile_s1_out, '_', bids_task, '_', suffix];
save(fullfile(bids_simTaskPath, [samplesFile_s2_out, '.mat']), 'SS');
%samplesFile_s2_out         = 'sampleSet_HandGesture_search-rad7';                                          % <-- for debugging or to pick up after this step (make sure to load the file first)



%%
%  Step 3 - Select the 3D sample-points with a searchlight classification score above a given threshold
%
[SS, suffix]                = s3_sampleSelection(SS);

% save results
bids_simTaskPath            = fullfile(bids_simPath, bids_task);
samplesFile_s3_out          = fullfile(bids_simTaskPath, [samplesFile_s2_out, '_', suffix]);
if ~exist(samplesFile_s3_out, 'dir'),     mkdir(samplesFile_s3_out);    end
save(fullfile(samplesFile_s3_out, [samplesFile_s2_out, '_', suffix, '.mat']), 'SS');
%samplesFile_s3_out         = fullfile(bids_simTaskPath, 'sampleSet_HandGesture_search-rad7_P95nofoc');      % <-- for debugging or to pick up after this step (make sure to load the file first)



%%
%  Step 4 - Project all virtual grid-configurations on a hull using the 3D sample-points as the center, with their random rotations
%
samplesFiles_s4_out         = s4_projectGrids(  SS, ...
                                                fullfile(bids_simPath, [hemi, '_ext_hull.gii']), ...
                                                fullfile(samplesFile_s3_out, [samplesFile_s2_out, '_', suffix, '_proj']));



%... more

