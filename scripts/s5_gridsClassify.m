%   
%   Step 5 - Classify each projected grid (3D sample point)
%   outputFiles = ...
% 
%       gridSampleSets              = A cell array with where each cell represents the filepath to a grid-configuration output file
%       ...
%       hullPath                    = path to the hull gifti file
%       classConfig                 = the classification configuration struct
%       numCylThreads               = the number of threads used to find the voxels in electrode cylinders (0 = set to #cores)
%       numClassThreads             = the number of threads used to classify (0 = set to #cores)
%
%
%   Returns: 
%       outputFiles                 = A cell array containing the paths to all grid-configurations output files (with the classifications)
%
%
%   Copyright 2020, Max van den Boom (Multimodal Neuroimaging Lab, Mayo Clinic, Rochester MN)

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function outputFiles = s5_gridsClassify(gridSampleSets, bids_rootPath, bids_sub, bids_task, hemi, hullPath, classConfig, numCylThreads, numClassThreads)


    gridLimit = 2000;                       % limit on the number of classifications
    allowGridPartialOutOfFOV = 0;           % flag whether a grid is allowed to partially stick outside the FOV (0 = not allowed, 1 = allowed)
    allowGridPartialOutOfROI = 1;           % flag whether a grid is allowed to partially stick outside the ROI (0 = not allowed, 1 = allowed)
    useVoxelsOutsideROI = 1;                % flag whether voxels outside of the ROI can be used (0 = voxels outside the ROI are excluded, 1 = voxels outside the ROI are included)


    % message
    disp(['- Hemisphere: ', hemi]);
    if ~isfield(classConfig, 'isolateTrials') || ~strcmpi(classConfig.isolateTrials, 'condition')
        disp('- No isolation of conditions');
    else
        disp('- Isolation of conditions');
    end
    disp(['- Minimum voxel mean for feature: ', num2str(classConfig.minVoxelMeanForFeature)]);
    if allowGridPartialOutOfFOV == 0
        disp('- Grids with electrodes outside of the FOV will be removed');
    else
        disp('- Grids with electrodes outside of the FOV will be allowed');
    end
    if allowGridPartialOutOfROI == 0
        disp('- Grids with electrodes outside of the ROI-box will be removed');
    else
        disp('- Grids with electrodes outside of the ROI-box will be allowed');
    end
    if useVoxelsOutsideROI == 0
        disp('- Only voxels inside of the ROI-box will be used for classifcation');
    else
        disp('- Voxels outside of the ROI-box will also be used for classifcation');
    end



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

    %%
    %  Load the roi hull, grey matter data and functional data
    %

    % load the grey matter volume (resliced to functional, native space)
    gmData = mx.nifti.readVol(bids_gmFilepath);

    % read the roi hull
    gHull = gifti(hullPath);

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
    %   Translate everything which is in pial space (pial, hull) to match the "MRI/NIFTI" space
    %   the freesurfer pial (surface) and the mgz/nii have a slight spatial difference.
    %   In the header of each mgz there is information to correct for this
    %

    % load the orig.mgz header information and calculate the transformation matrix to bring it back to native space
    load(fullfile(bids_fsPath, 'mri', filesep, 'orig.mat'));
    transMatFsToNative      = orig.vox2ras * inv(orig.tkrvox2ras);

    %gSurfPial.vertices      = mx.three_dimensional.transVerticesByAffineMat(gSurfPial.vertices, transMatFsToNative);
    gHull.vertices          = mx.three_dimensional.transVerticesByAffineMat(gHull.vertices, transMatFsToNative);
    
    
    %%
    %  Classify each grid configuration on the sample-set points
    %

    outputFiles = {};

    % loop through the grid sample-sets
    for iGridSampleSet = 1:length(gridSampleSets)

        % create a clear workspace, while keeping essential input variables
        clearvars -except iGridSampleSet gridSampleSets ...
            bids_rootPath bids_sub bids_task hemi classConfig samplesFile niiGmFilepath ...
            gHull transMatFsToNative events json ...
            funcXyzData funcVol funcVolData gmData ...
            gridLimit numCylThreads numClassThreads ...
            allowGridPartialOutOfFOV allowGridPartialOutOfROI useVoxelsOutsideROI outputFiles;

        % message
        disp(['file: ', gridSampleSets{iGridSampleSet}]);


        %%
        %  Load the samples
        %

        % load the samples
        load(gridSampleSets{iGridSampleSet});

        % transfer the classification settings in the config
        SS.gridConfig.allowGridPartialOutOfFOV = allowGridPartialOutOfFOV;    % flag whether a grid is allowed to partially stick outside the FOV (0 = not allowed, 1 = allowed)
        SS.gridConfig.allowGridPartialOutOfROI = allowGridPartialOutOfROI;    % flag whether a grid is allowed to partially stick outside the ROI (0 = not allowed, 1 = allowed)
        SS.gridConfig.useVoxelsOutsideROI = useVoxelsOutsideROI;              % flag whether voxels outside of the ROI can be used (0 = voxels outside the ROI are excluded, 1 = voxels outside the ROI are included)



        %%
        %  Tranform everything related in freesurfer space (samplepoints) to match the "MRI/NIFTI" space
        %  Note: the freesurfer pial (surface) and the mgz/nii have a slight spatial difference, in the 
        %        header of each mgz there is information to correct for this
        %

        % check if the sampleset positions have not yet been corrected (translated)
        if isfield(SS, 'transformedToNiftiSpace') == 0 || isempty(SS.transformedToNiftiSpace) || length(SS.transformedToNiftiSpace) ~= 3
            % not corrected

            % transform the positions
            SS.preTrans_samplePositions     = SS.samplePositions;
            SS.samplePositions              = mx.three_dimensional.transVerticesByAffineMat(SS.samplePositions, transMatFsToNative);
            SS.preTrans_sampleGrids         = SS.sampleGrids;
            SS.sampleGrids                  = mx.three_dimensional.transVerticesByAffineMat(SS.sampleGrids, transMatFsToNative);
            if isfield(SS, 'roiBox')
                SS.preTrans_roiBox          = SS.roiBox;
                SS.roiBox.vertices          = mx.three_dimensional.transVerticesByAffineMat(SS.roiBox.vertices, transMatFsToNative);
            end

            % flag as corrected by storing the transformation
            SS.transformedToNiftiSpace      = transMatFsToNative;

        end


        %%
        %  Check the functional field of view and remove the grids where one or more electrodes are outside of the fMRI FOV
        %
        if ~isfield(SS.gridConfig, 'allowGridPartialOutOfFOV') || SS.gridConfig.allowGridPartialOutOfFOV == 0

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
            [~, fovLineSegMatrix, fovTris] = mx.three_dimensional.rectPointsTo3D(fovCorners);
            fovBox = [];
            fovBox.vert = fovCorners;
            fovBox.tri = fovTris;
            fovBox = gifti(fovBox);

            % determine which sample grid points (electrodes) are inside of the fov-box
            fovRaw = [];
            fovRaw.vertices = double(fovBox.vertices);
            fovRaw.faces = fovTris;
            insideFOV = mx.three_dimensional.ext.inpolyhedron(fovRaw, SS.sampleGrids);
            outsideFovIndices = find(insideFOV == 0);

            % check if there are electrodes outside of the fMRI FOV
            if ~isempty(outsideFovIndices)
                % one or more electrodes are out of FOV

                % determine the start electrode of the grids where (at least one) electrode
                % is outside of the FOV, the samples which are to be removed
                remSamples = floor((outsideFovIndices - 1) / SS.gridConfig.elecTotal) + 1;
                remSamples = unique(remSamples);

                % message
                disp (['Sample grids to be removed (one or more electrodes outside of functional MRI FOV): ', num2str(length(remSamples))]);

                % remove the samples (and their grids)
                SS = removeFromSampleSetByIndices(SS, remSamples);
                clear remSampleGrids;

                % message remaining
                disp (['Sample grids remaining: ', num2str(size(SS.samplePositions, 1))]);

            else
                % no electrodes outside FOV

                % message
                disp ('All sample grids are within the fMRI FOV');

            end

        end


        %%
        %  Check the ROI and remove the grids where one or more electrodes are outside of the ROI box
        %
        if ~isfield(SS.gridConfig, 'allowGridPartialOutOfROI') || SS.gridConfig.allowGridPartialOutOfROI == 0

            % check if a ROI-box is available
            if ~isfield(SS, 'roiBox')

                % message and return
                error('Error: grids should be restricted to the ROI-box, however the samplefile does not specify a ROI-box');

            end

            % determine which sample grid points (electrodes) are inside of the ROI-box
            roiRaw = [];
            roiRaw.vertices = double(SS.roiBox.vertices);
            roiRaw.faces = SS.roiBox.faces;
            insideROI = mx.three_dimensional.ext.inpolyhedron(roiRaw, SS.sampleGrids);
            outsideROIIndices = find(insideROI == 0);

            % check if there are electrodes outside of the ROI-box
            if ~isempty(outsideROIIndices)
                % one or more electrodes are out of ROI-box

                % determine the start electrode of the grids where (at least one) electrode
                % is outside of the ROI-box, the samples which are to be removed
                remSamples = floor((outsideROIIndices - 1) / SS.gridConfig.elecTotal) + 1;
                remSamples = unique(remSamples);

                % message
                disp (['Sample grids to be removed (one or more electrodes outside of ROI-box): ', num2str(length(remSamples))]);

                % remove the samples (and their grids)
                SS = removeFromSampleSetByIndices(SS, remSamples);
                clear remSampleGrids;

                % message remaining
                disp (['Sample grids remaining: ', num2str(size(SS.samplePositions, 1))]);

            else
                % no electrodes outside ROI-box

                % message
                disp ('All sample grids are within the ROI-box');

            end

        end

    
        %%
        %   Limit the number of grids. The classification of each grid takes a lot of time
        %   so apply an upper cap on the number of grids/classifications
        %
        if gridLimit ~= 0 && size(SS.samplePositions, 1) > gridLimit

            % remove the samples
            SS = removeFromSampleSetByIndices(SS, gridLimit + 1:size(SS.samplePositions, 1));

        end
    

        %%
        %  Determine which voxels are available for the "electrode"-cylinders to use
        %

        % check if classification should only be based on the voxels inside of the ROI-box
        if isfield(SS.gridConfig, 'useVoxelsOutsideROI') && SS.gridConfig.useVoxelsOutsideROI == 0
            % voxels in ROI-box and grey voxels

            % check if a ROI-box is available
            if ~isfield(SS, 'roiBox')

                % message and return
                error('Error: voxels should be restricted to the ROI-box , however the samplefile does not specify a ROI-box');

            end

            % determine which sample-positions are in the ROI-Box
            roiRaw = [];
            roiRaw.vertices = double(SS.roiBox.vertices);
            roiRaw.faces = SS.roiBox.faces;
            roiVoxels = mx.three_dimensional.ext.inpolyhedron(roiRaw, funcXyzData);

            % set the grey matter voxels as available voxels
            availableVoxelIndices = find(gmData(:) == 1 & roiVoxels == 1);
            %availableVoxelIndices = 1:numel(roiVoxels == 1);            % debug, use all voxels, not just grey matter


        else
            % grey voxels only

            % set the grey matter voxels as available voxels
            availableVoxelIndices = find(gmData == 1);
            %availableVoxelIndices = 1:numel(gmData);            % debug, use all voxels, not just grey matter

        end

        % select only the available (grey matter/ROI-box) voxels to determine their distances to the electrodes
        funcAvailableXyz = funcXyzData(availableVoxelIndices, :);
    
        
        %%
        %  Determine the voxels for each "electrode"-cylinder

        % get the first edges of each triangle as unit vectors and calculate the normals (as unit vectors)
        e0 = double(gHull.vertices(gHull.faces(:, 2), :) - gHull.vertices(gHull.faces(:, 1), :));
        e1 = double(gHull.vertices(gHull.faces(:, 3), :) - gHull.vertices(gHull.faces(:, 1), :));
        e0 = e0 ./ vecnorm(e0')';
        normals = cross(e0, e1, 2);
        normals = normals ./ vecnorm(normals')';

        % convert the points to cylinder axis
        cylAxis = [SS.sampleGrids, SS.sampleGrids + normals(SS.sampleGridsTriangles, :) * 5];
        cylAxisOut = [SS.sampleGrids, SS.sampleGrids - normals(SS.sampleGridsTriangles, :) * 5];

        % determine which voxels are within which electrode radius distance
        % (both for inward cylinders and outward ones)
        splitConfig = [];
        splitConfig.numPoints = 2000;
        if gridLimit ~= 0
            if numCylThreads == 0
                splitConfig.numPoints = ceil(gridLimit / feature('numcores'));
            else
                splitConfig.numPoints = ceil(gridLimit / numCylThreads);
            end
        end
        splitConfig.threads = numCylThreads;
        [elecVoxels, distOnAxis, distFromAxis] = ...
            mx.three_dimensional.retrieveCylindricalProximatePoints_split(  funcAvailableXyz, ...
                                                                            cylAxis, ...
                                                                            SS.gridConfig.elecSizeRadius, ...
                                                                            splitConfig, ...
                                                                            1);
        [elecVoxelsCylOut, distOnAxisCylOut, distFromAxisCylOut] = ...
            mx.three_dimensional.retrieveCylindricalProximatePoints_split(  funcAvailableXyz, ...
                                                                            cylAxisOut, ...
                                                                            SS.gridConfig.elecSizeRadius, ...
                                                                            splitConfig, ...
                                                                            1);
    
        % lookup which voxels they were in the original (now they are "available" voxel indexed)
        n = size(SS.sampleGrids, 1);
        for i = 1:n
            elecVoxels{i} = availableVoxelIndices(elecVoxels{i});
            elecVoxelsCylOut{i} = availableVoxelIndices(elecVoxelsCylOut{i});
            elecVoxels{i} = unique([elecVoxels{i}; elecVoxelsCylOut{i}]);
        end
    
    
        %%
        %  Classify each grid
        %

        % a variable to the grids, their electrodes and the voxel per electrode
        gridElectrodeVoxels = cell(1, size(SS.samplePositions, 1));

        % loop through the grids
        for iGrid = 1:size(SS.samplePositions, 1)

            % retrieve the electrode voxels
            elecStart = ((iGrid - 1) * SS.gridConfig.elecTotal) + 1;
            voxelsPerElectrode = elecVoxels(elecStart:elecStart + SS.gridConfig.elecTotal - 1);

            gridElectrodeVoxels{iGrid} = voxelsPerElectrode;
        end       

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
        if gridLimit ~= 0
            if numClassThreads == 0
                splitConfig.numPoints = ceil(gridLimit / feature('numcores'));
            else
                splitConfig.numPoints = ceil(gridLimit / numClassThreads);
            end
        end
        splitConfig.threads = numClassThreads;
        [clsResults, success] = mx.class.classifyFmriDataSets_split(classConfig, ...
                                                                    gridElectrodeVoxels, ...
                                                                    sets, ...
                                                                    'within', ...
                                                                    splitConfig     );

        % make sure subjResults does not exists, this will cause it to be a struct array, rather than something else
        clear subjResults;

        % build the maskCode
        maskCode = 'Elec3D';
        if isfield(SS.gridConfig, 'useVoxelsOutsideROI') && SS.gridConfig.useVoxelsOutsideROI == 0
            maskCode = [maskCode, '_roiVoxelsOnly'];
        end

        % loop through the grids
        counter = 1;
        for iGrid = 1:size(SS.samplePositions, 1)

            % check if classification failed
            if ~success(iGrid)

                % continue to the next
                continue;

            end


            % store the result
            objResult = mx.class.createResult(  bids_sub, ...                       % subject
                                                bids_task, ...                      % run/task
                                                0, ...                              % numVoxels
                                                maskCode, ...                       % maskCode
                                                classConfig, ...                    % config
                                                clsResults(iGrid), ...              % results
                                                sets(1).conditions, ...             % conditions
                                                sets(1).onsets, ...                 % onsets (in seconds)
                                                sets(1).tr, ...                     % tr
                                                [], ...                             % tmaps
                                                hemi, ...
                                                iGrid);
            objResult.elecTotal         = SS.gridConfig.elecTotal;
            objResult.elecNumHorz       = SS.gridConfig.elecNumHorz;
            objResult.elecNumVert       = SS.gridConfig.elecNumVert;
            objResult.spacing           = SS.gridConfig.elecHorzSpacing;
            objResult.radius            = SS.gridConfig.elecSizeRadius;

            % add to the results
            subjResults(counter) = objResult;
            counter = counter + 1;

        end 

        if ~exist('subjResults', 'var')
            subjResults = [];
        end
    
        % message
        disp(['Number of grids classified: ', num2str(length(subjResults))]);

    

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

            %SS.sampleGrids      = SS.sampleGrids - SS.transformedToNiftiSpace;
            SS.sampleGrids      = SS.preTrans_sampleGrids;
            SS = rmfield(SS, 'preTrans_sampleGrids');

            if isfield(SS, 'roiBox')
                %SS.roiBox.vertices  = SS.roiBox.vertices - SS.transformedToNiftiSpace;
                SS.roiBox.vertices  = SS.preTrans_roiBox.vertices;
                SS = rmfield(SS, 'preTrans_roiBox');
            end

            % flag as uncorrected by removing the field
            SS = rmfield(SS, 'transformedToNiftiSpace');

        end

    
        %%
        % save the result
        %
        
        newFile = strrep(gridSampleSets{iGridSampleSet}, '.mat', '_Results.mat');
        save(newFile, 'SS', 'subjResults');

        % store the path to the grid-configuration sampleset with classifications
        outputFiles{end + 1} = newFile;
        
    end
    
end