%   
%   Step 2 - Perform searchlight classification on each of the 3D sample-points with a specific radius
%            Sample-points outside of the FOV of the task are removed beforehand
%
%   [SS, filenameSuffix] = s2_searchlight(SS, bids_rootPath, bids_sub, bids_task, hemi, searchLightRadius, classConfig, numThreads)
% 
%       SS                  = Structure that holds the sample-points to apply the searchlight classifications on
%       bids_rootPath       = path to the BIDS root-directory where the data is located
%       bids_sub            = the subject to perform the searchlight step on
%       bids_task           = the task to use for the searchlight classification
%       hemi                = the hemisphere that this step is applied on
%       searchLightRadius   = radius of the searchlight in mm
%       classConfig         = The classification configuration struct
%       numThreads          = The number of threads used to classifys (0 = set to #cores)
%
%
%   Returns: 
%       SS                  = The input structure with the searchlight classifications results added
%       filenameSuffix      = Depending on the searchlight radius, a suggestion for a filename suffix
%
%
%   Copyright 2020, Max van den Boom (Multimodal Neuroimaging Lab, Mayo Clinic, Rochester MN)

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [SS, filenameSuffix] = s2_searchlight(SS, bids_rootPath, bids_sub, bids_task, hemi, searchLightRadius, classConfig, numThreads)
    if exist('numThreads', 'var') == 0,  numThreads = 0;     end
    
    searchlightLimit = 200000;              % limit on the number of searchlights

    
    %%
    % retrieve the root path and make sure dependencies can be found
    %
    
    addpath('../');
    gridSimRoot = gridSimRootPath;
    addpath([gridSimRoot, filesep, 'functions']);
    addpath(genpath([gridSimRoot, filesep, 'external']));



    %%
    %  Build paths and read the files
    %

    % build paths to the data
    bids_simPath = fullfile(bids_rootPath, 'derivatives', [hemi, '_simulations'], ['sub-' bids_sub]);
    bids_fsPath = fullfile(bids_rootPath, 'derivatives', 'freesurfer', ['sub-' bids_sub]);
    bids_gmFilepath = fullfile(bids_rootPath, 'derivatives', 'gm_masks', ['sub-' bids_sub], 'func_fs_mask.nii');

    bids_func_eventFilepath = fullfile(bids_rootPath, ['sub-' bids_sub], 'func', ['sub-' bids_sub '_task-' bids_task '_events.tsv']);
    bids_func_dataFilepath = fullfile(bids_rootPath, ['sub-' bids_sub], 'func', ['sub-' bids_sub '_task-' bids_task '_bold.nii']);
    bids_func_jsonFilepath = fullfile(bids_rootPath, ['sub-' bids_sub], 'func', ['sub-' bids_sub '_task-' bids_task '_bold.json']);

    % load the grey matter volume (resliced to functional, native space)
    gmData = mx.nifti.readVol(bids_gmFilepath);

    % load the task volumes, event and JSON file
    [funcVolData, funcVol]      = mx.nifti.readVol(bids_func_dataFilepath);  funcVol = funcVol(1);
    events                      = readtable(bids_func_eventFilepath, 'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', {'N/A', 'n/a'}, 'ReadVariableNames', true);
    json                        = jsondecode(fileread(bids_func_jsonFilepath));

    % generate the real-word coordinates for the voxels in the functional images
    [R, C, P] = ndgrid(1:funcVol.dim(1), 1:funcVol.dim(2), 1:funcVol.dim(3));
    RCP       = [R(:)'; C(:)'; P(:)'];
    clear R C P
    RCP(4, :) = 1;
    funcXyzData = funcVol.mat(1:3, :) * RCP;
    funcXyzData = funcXyzData';
    clear RCP;

    % check whether the grey matter volume has the same dimensions as the functional data
    if  size(funcVolData, 1) ~= size(gmData, 1) || size(funcVolData, 2) ~= size(gmData, 2) || size(funcVolData, 3) ~= size(gmData, 3)
        error('Error: the grey matter volume and functional volume data dimensions mismatch');
    end

    % NaN out the non-grey-matter values in the functional data
    aa = reshape(funcVolData, [], size(funcVolData, 4));
    aa(gmData == 0, :) = NaN;
    funcVolData = reshape(aa, size(funcVolData));
    clear aa;



    %%
    %  Tranform everything related in freesurfer space (samplepoints) to match the "MRI/NIFTI" space
    %  Note: the freesurfer pial (surface) and the mgz/nii have a slight spatial difference, in the 
    %        header of each mgz there is information to correct for this
    %

    % load the orig.mgz header information and calculate the transformation matrix to bring it back to native space
    load(fullfile(bids_fsPath, 'mri', filesep, 'orig.mat'));
    transMatFsToNative      = orig.vox2ras * inv(orig.tkrvox2ras);

    % check if the sampleset positions have not yet been corrected (transformed)
    if isfield(SS, 'transformedToNiftiSpace') == 0 || isempty(SS.transformedToNiftiSpace)
        % not corrected

        % transform the positions
        SS.preTrans_samplePositions     = SS.samplePositions;
        SS.samplePositions              = mx.three_dimensional.transVerticesByAffineMat(SS.samplePositions, transMatFsToNative);
        if isfield(SS, 'roiBox')
            SS.preTrans_roiBox          = SS.roiBox;
            SS.roiBox.vertices          = mx.three_dimensional.transVerticesByAffineMat(SS.roiBox.vertices, transMatFsToNative);
        end

        % flag as corrected by storing the transformation
        SS.transformedToNiftiSpace      = transMatFsToNative;

    end



    %%
    %  Check the functional field of view and remove the samples which are outside of the fMRI FOV
    %

    % determine the corners of the fov
    numVoxelsInSlice = size(funcVolData, 1) * size(funcVolData, 2);
    fovCorners = [  funcXyzData(1, :); ...
                    funcXyzData(size(funcVolData, 1), :); ...
                    funcXyzData(numVoxelsInSlice, :); ...
                    funcXyzData(numVoxelsInSlice - size(funcVolData, 1) + 1, :); ...
                    funcXyzData(end, :); ...                            % end is same as numVoxelsInSlice * size(funcVolData, 3)
                    funcXyzData(end - size(funcVolData, 1) + 1, :); ...
                    funcXyzData(end - numVoxelsInSlice + 1, :); ...
                    funcXyzData(end - numVoxelsInSlice + size(funcVolData, 1), :)];

    % create a box based on the FOV (FOV might not be completely rectangular)
    [~, ~, fovTris] = mx.three_dimensional.rectPointsTo3D(fovCorners);
    fovBox = [];
    fovBox.vert = fovCorners;
    fovBox.tri = fovTris;
    fovBox = gifti(fovBox);

    % determine which sample grid points (electrodes) are inside of the fov-box
    fovRaw = [];
    fovRaw.vertices = double(fovBox.vertices);
    fovRaw.faces = fovTris;
    insideFOV = mx.three_dimensional.ext.inpolyhedron(fovRaw, SS.samplePositions);
    outsideFovIndices = find(insideFOV == 0);
    clear insideFOV;

    % check if there are samples outside of the fMRI FOV
    if ~isempty(outsideFovIndices)
        % one or more samples are out of FOV

        % remove the samples
        SS = removeFromSampleSetByIndices(SS, outsideFovIndices);

        % message remaining
        disp (['Samples remaining (within fMRI FOV): ', num2str(size(SS.samplePositions, 1))]);

    else
        % no samples outside FOV

        % message
        disp ('All samples are within the fMRI FOV');

    end
    clear outsideFovIndices;



    %%
    %  Limit the number of grids. The classification of each grid takes a lot of time
    %  so apply an upper cap on the number of grids/classifications
    %
    if searchlightLimit ~= 0 && size(SS.samplePositions, 1) > searchlightLimit

        % remove the samples
        SS = removeFromSampleSetByIndices(SS, searchlightLimit + 1:size(SS.samplePositions, 1));

    end



    %%
    %  Prepare variables for classification
    %

    % determine which voxels (indices) are grey matter
    gmIndices = find(gmData == 1);
    clear gmData;

    % select only the grey matter voxels to determine their distances to the samples
    %funcGmXyz = funcXyzData(:, :);
    funcGmXyz = funcXyzData(gmIndices, :);
    clear funcXyzData;

    % create variables to store the classification results in
    SS.sampleSearchFocalResult = nan(size(SS.samplePositions, 1), 1);
    SS.sampleSearchNonFocalResult = nan(size(SS.samplePositions, 1), 1);
    SS.sampleSearchNonFocalResultTotal = zeros(size(SS.samplePositions, 1), 1);
    SS.sampleSearchNonFocalResultCounter = zeros(size(SS.samplePositions, 1), 1);



    %%
    % Searchlight on the samples
    %

    % determine which voxels are within each sample point's searchlight distance
    splitConfig = [];
    splitConfig.numPoints = 2000;
    sampleVoxels = mx.three_dimensional.retrieveRadialProximatePoints_split(funcGmXyz, ...
                                                                            SS.samplePositions, ...
                                                                            searchLightRadius, ...
                                                                            splitConfig);
    clear funcGmXyz;

    % lookup the voxel indices in the original (now they are grey matter voxel indexed)
    % and clear out the feature voxels if less than the minimum features (don't
    % delete from the array because we want to keep the numbering (index)
    % intact)
    for iSample = 1:length(sampleVoxels)
        sampleVoxels{iSample} = gmIndices(sampleVoxels{iSample});

        if length(sampleVoxels{iSample}) < classConfig.minFeatures
            sampleVoxels{iSample} = [];
        end
    end
    clear gmIndices;

    % transfer the information
    sets = {};
    sets(1).volumes = funcVolData;
    sets(1).conditions = events.conditions;
    if any(sets(1).conditions == 0)
       error('Make sure none of the conditions in the events file are 0, conditions should be 1-based'); 
    end
    sets(1).units = 'seconds';
    sets(1).onsets = events.onset;
    sets(1).durations = events.duration;
    sets(1).tr = json.RepetitionTime;

    % classify
    splitConfig = [];
    splitConfig.numPoints = 2000;
    if searchlightLimit ~= 0
        if numThreads == 0
            splitConfig.numPoints = ceil(searchlightLimit / feature('numcores'));
        else
            splitConfig.numPoints = ceil(searchlightLimit / numThreads);
        end
    end
    splitConfig.threads = numThreads;
    [clsResults, success] = mx.class.classifyFmriDataSets_split(classConfig, ...
                                                                sampleVoxels, ...
                                                                sets, ...
                                                                'within', ...
                                                                splitConfig);

    % loop through the samples
    numSamples = size(SS.samplePositions, 1);
    for iSample = 1:numSamples

        % check if classification failed
        if ~success(iSample)

            % continue to the next
            continue;

        end

        % store the searchlight focal result
        SS.sampleSearchFocalResult(iSample) = clsResults(iSample).accTotal;                                    

        % calculate the distance from the sample to each other sample and
        % determine which samples are within the sample searchlight radius
        saDist = (SS.samplePositions(iSample, 1)' - SS.samplePositions(:, 1)) .^ 2 + ...
                 (SS.samplePositions(iSample, 2)' - SS.samplePositions(:, 2)) .^ 2 + ...
                 (SS.samplePositions(iSample, 3)' - SS.samplePositions(:, 3)) .^ 2;
        saDist = find(saDist < (searchLightRadius  .^ 2));

        % store the searchlight non-focal results
        for iNonFocalSample = 1:length(saDist)

            SS.sampleSearchNonFocalResultTotal(saDist(iNonFocalSample)) = SS.sampleSearchNonFocalResultTotal(saDist(iNonFocalSample)) + clsResults(iSample).accTotal;
            SS.sampleSearchNonFocalResultCounter(saDist(iNonFocalSample)) = SS.sampleSearchNonFocalResultCounter(saDist(iNonFocalSample)) + 1;

        end

    end     % end sample point loop

    % average the non-focal sample results
    for iSample = 1:length(SS.sampleSearchNonFocalResult)

        % calculate the mean based on the total and the counter
        SS.sampleSearchNonFocalResult(iSample) = SS.sampleSearchNonFocalResultTotal(iSample) / SS.sampleSearchNonFocalResultCounter(iSample);

    end

    % remove the total field from the structure (leave the counter, it might be
    % handy to see from how many values the sample average is made up)
    SS = rmfield(SS, 'sampleSearchNonFocalResultTotal');



    %%
    %  Undo the transformation of the samplepoint from before (to "MRI/NIFTI" space from freesurfer space)
    %

    % check if the sampleset positions have been corrected (translated)
    if isfield(SS, 'transformedToNiftiSpace') == 1 && ~isempty(SS.transformedToNiftiSpace)
        % corrected

        % transform the positions back
        %SS.samplePositions      = SS.samplePositions - SS.transformedToNiftiSpace;
        SS.samplePositions      = SS.preTrans_samplePositions;
        SS = rmfield(SS, 'preTrans_samplePositions');
        if isfield(SS, 'roiBox')
            %SS.roiBox.vertices  = SS.roiBox.vertices - SS.transformedToNiftiSpace;
            SS.roiBox.vertices  = SS.preTrans_roiBox.vertices;
            SS = rmfield(SS, 'preTrans_roiBox');
        end

        % flag as uncorrected by removing the field
        SS = rmfield(SS, 'transformedToNiftiSpace');

    end

    % return a file suffix
    filenameSuffix = ['search-rad', num2str(searchLightRadius)];
    
end