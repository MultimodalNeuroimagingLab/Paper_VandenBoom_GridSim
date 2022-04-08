%
%   Perform classification based on one or more fMRI datasets, either per set (within), between sets (between) or concatenated sets (combine),
%   splitting up the input for speed optimization (multi-threading) or memory maximization
%
%   Usage: [results, success, signalSets] = classifyFmriDataSets(config, featureVoxels, sets, setMethod)
%
%   config                  = 
%   featureVoxels           = Indicates which voxels to use as features from the datasets, this 
%                             parameter can be used in several ways:
%                               [indices]    = one dimensional array that holds the linear indices of
%                                              the voxels in the input datasets to use as features
%                               [x-y-z]      = three-dimensional array (binary mask) to indicate which 
%                                              voxels to use as features for classification. The mask should have
%                                              the exact same dimension sizes as the first three dimensions
%                                              of the input datasets
%                               [cell array] = array which specifies that multiple classifications should
%                                              be performed using different feature-sets. Each cell
%                                              will result in classification using a specific set of 
%                                              voxels as features. The voxels to use for each cell (=featureset=classification)
%                                              are defined as the options above ([], [indices] or [x-y-z] mask)
%                               [cell array of cell arrays] = array which specifies that multiple classifications 
%                                              should be performed using different feature-sets (first array of cells).
%                                              Each cell in the first array of cells holds another cell array. Each cell
%                                              in this second cell array represents a feature in the featureset which 
%                                              consists of one or more voxels. The voxels to select can be defined
%                                              as above ([], [indices] or [x-y-z] mask) and will be averaged to
%                                              serve as one feature
%   sets                    = struct array of fmri datasets
%   sets(#).volumes         = The volume input data to pick the features from and classify on. The data 
%                             should hold 3D volume data in the first 3 dimensions, and time in the 4th
%   sets(#).conditions      = The experimental conditions per trial
%   sets(#).units           = The units the onsets and durations are in (can either be 'seconds' or 'scans')
%   sets(#).onsets          = A vector with the trial onsets
%   sets(#).durations       = A vector with the trial durations
%   sets(#).tr              = The repetition time used for the dataset
%   setMethod               = The dataset classification method, either:
%                               'within': classify per dataset, cross-validation
%                                         will be performed within each set
%                               'between': train on one dataset, classify on the
%                                          other. Exactly two input sets are expected
%                               'combine': concatinate all input sets and perform
%                                          cross-validation over the joint set
% 
%   splitConfig             = the configuration on splitting the input
%   splitConfig.numSets     = (optional) split input into x sets of 3D-points [0 = no splitting into sets]
%   splitConfig.numPoints   = (optional) splits the input in sets of x 3D-points [0 = no splitting into sets]
%   splitConfig.threads     = (optional) number of threads to run the sets on (requires the
%                             set to be split using either numSets or numPoints)
% 
% 
% 
%   Returns: 
%       results             = struct array with classification result(s). The first dimension (rows)
%                             represents the different feature-sets. The second dimension will represent
%                             the input datasets, the size of the second dimension depends on the classification
%                             method: 'within' will give one result per dataset (for each featureset); 'combined' and
%                             'between' will give one result per 2 or more datasets (for each featureset)
%                             Note: empty results structs can be returned on failure, use the success return
%                             variable to check
%       success             = vector which - for each selection of features (if featureVoxels is not cellarray
%                             then 1 flag; if featureVoxels is a cellarray then a flag for each cell) - indicates
%                             success (1) or failure (0) given that featureset
%       signalSets          = struct-array with the signal-sets and meta-info (onsets, durations, 
%                             conditions) that were used for classification. The first dimension (rows)
%                             represents the different feature-sets, the second dimension (columns) holds
%                             the input datasets.
%
%
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [results, success, signalSets] = classifyFmriDataSets_split(config, featureVoxels, sets, setMethod, splitConfig)
    warning('backtrace', 'off');
    
    % check if one featureset, then just consider it as one feature-set
    if ~iscell(featureVoxels)
        featureVoxels = {featureVoxels};
    end
    
    % determine the number of feature-sets
    numFeaturesSets = length(featureVoxels);
    
    % initialize the success values to zero
    success = zeros(1, numFeaturesSets);
    
    % loop through the sets
    for iSet = 1:length(sets)
        
        % create an empty struct
        % (and with that the struct array)
        signalSetsPrep(iSet) = mx.class.fmriSignalSet();
        
        % display
        if strcmpi(config.silent, 'yes') == 0,  disp(['set ', num2str(iSet)]); end
        
        % store the input volumes size
        signalSetsPrep(iSet).inputVolumesDims = size(sets(iSet).volumes);
        
        % save the input units, onsets, durations, conditions and tr. 
        % These initial input (units) are copied and adjusted further
        % down and will be used when the timelock method is set 
        % to 'highest bold' to generate a hrf.
        signalSetsPrep(iSet).inputUnits = sets(iSet).units;
        signalSetsPrep(iSet).inputOnsets = sets(iSet).onsets;
        signalSetsPrep(iSet).inputDurations = sets(iSet).durations;
        signalSetsPrep(iSet).inputConditions = sets(iSet).conditions;
        signalSetsPrep(iSet).inputTr = sets(iSet).tr;
        
        % create copies of the input onsets, durations and conditions which can be adjusted
        % (these are also used when the timelock method is set to 'highest bold' to generate a hrf)
        signalSetsPrep(iSet).units = signalSetsPrep(iSet).inputUnits;
        signalSetsPrep(iSet).onsets = signalSetsPrep(iSet).inputOnsets;
        signalSetsPrep(iSet).durations = signalSetsPrep(iSet).inputDurations;
        signalSetsPrep(iSet).conditions = signalSetsPrep(iSet).inputConditions;
        
        % make sure a version of the onsets is available in scans. This is used for
        % several checks and to epoch the data.
        if strcmpi(sets(iSet).units, 'seconds') == 1
            signalSetsPrep(iSet).onsetsScans = round(sets(iSet).onsets / sets(iSet).tr);
            if signalSetsPrep(iSet).onsetsScans(1) == 0
               signalSetsPrep(iSet).onsetsScans(1) = 1; 
            end
        elseif strcmpi(sets(iSet).units, 'scans') == 1
            signalSetsPrep(iSet).onsetsScans = sets(iSet).onsets;
        else
            fprintf(2, 'Error: unknown unit type\n');
            success(:) = 0;
            results = [];
            signalSets = []; 
            return;
        end
        
        % check if there are scan onsets equal or higher than the number of volumes
        excessOnsets = nnz(signalSetsPrep(iSet).onsetsScans >= signalSetsPrep(iSet).inputVolumesDims(4));
        if excessOnsets > 0
            warning(['Warning: One or more onsets are above the number of volumes. This can occur when scanner acquisition stopped early (giving a smaller amount of volumes to build a HRF on), while the task was still running (trials/conditions were still being added). Removing the last ', num2str(excessOnsets) ,' onsets']);
            signalSetsPrep(iSet).onsets = signalSetsPrep(iSet).onsets(1:end - excessOnsets);
            signalSetsPrep(iSet).onsetsScans = signalSetsPrep(iSet).onsetsScans(1:end - excessOnsets);
        end
        
        % check the number of conditions
        diffConditions = length(signalSetsPrep(iSet).conditions) - length(signalSetsPrep(iSet).onsetsScans);
        if diffConditions > 0
            warning(['Warning: more conditions than onsets, removing the last ', num2str(diffConditions) ,' conditions']);
            signalSetsPrep(iSet).conditions = signalSetsPrep(iSet).conditions(1:end - diffConditions);
        elseif diffConditions < 0
            fprintf(2, 'Error: less conditions than onsets\n');
            success(:) = 0;
            results = [];
            signalSets = [];
            return;
        end

        % check the number of durations
        diffDurations = length(signalSetsPrep(iSet).durations) - length(signalSetsPrep(iSet).onsetsScans);
        if diffDurations > 0
            warning(['Warning: more durations than onsets, removing the last ', num2str(diffDurations) ,' durations']);
            signalSetsPrep(iSet).durations = signalSetsPrep(iSet).durations(1:end - diffDurations);
        elseif diffDurations < 0
            fprintf(2, 'Error: less durations than onsets\n');
            success(:) = 0;
            results = [];
            signalSets = []; 
            return;
        end

    end


    %%%
    %% prepare and check input data(sets) and featuresets
    %%%
    
    % for each dataset, create version of input data where the first dimension represents
    % space (indexed linearly) and the second dimension time (this will serve as the base to pick the signal from)
	data = [];
    numVoxels = 0;
    for iSet = 1:length(sets)
        if length(size(sets(iSet).volumes)) == 4
            numVoxels = prod(size(sets(iSet).volumes(:, :, :, 1)));
            data(:, :, iSet) = reshape(sets(iSet).volumes, prod(size(sets(iSet).volumes(:, :, :, 1))), size(sets(iSet).volumes, 4));
        else
            %TODO: implement to accept more volume data format (voxels x time for example)
            fprintf(2, 'Error: unknown formatting of the volume data (expecting 4D)\n');
            success(:) = 0;
            results = [];
            signalSets = []; 
            return;
        end
    end
    
    % check the featureVoxels parameter(s) and
    % sure each featureSet (set of voxels) is based on linear voxel indices
    for iFeat = 1:numFeaturesSets
        
        % 
        if isempty(featureVoxels{iFeat})
            % empty
            % no voxels is also an option and will result in an empty struct later
            
        elseif ~iscell(featureVoxels{iFeat}) && isvector(featureVoxels{iFeat})
            % array of voxel indices
            
            % already in the format we want, leave as it is
            
        elseif ~iscell(featureVoxels{iFeat}) && length(size(featureVoxels{iFeat})) == 3
            % 3D-mask
            
            % retrieve the indices
            featureVoxels{iFeat} = find(featureVoxels{iFeat}(:));
            
            % TODO: allow mask where different voxel values indicate
            % different features, should be easy if converted to cell array
            % for now, just allow binary (1 and 0)
            
        elseif iscell(featureVoxels{iFeat})
            % cell array (where the voxels from each cell will be averaged to a feature)
            
            % loop throug the each cell (feature) and ensure voxel indices
            for iFeat2 = length(featureVoxels{iFeat}):-1:1
                if isempty(featureVoxels{iFeat}{iFeat2})
                    % empty, 
                    
                    % no voxels for the feature, remove the feature from
                    featureVoxels{iFeat}(iFeat2) = [];
                    
                elseif ~iscell(featureVoxels{iFeat}{iFeat2}) && isvector(featureVoxels{iFeat}{iFeat2})
                    % array of indices
                    
                    % already in the format we want, leave as it is

                elseif ~iscell(featureVoxels{iFeat}{iFeat2}) && length(size(featureVoxels{iFeat}{iFeat2})) == 3
                    % 3D-mask

                    % retrieve the indices
                    featureVoxels{iFeat}{iFeat2} = find(featureVoxels{iFeat}{iFeat2}(:));
                    
                else
                    fprintf(2, ['Error: unknown format to define the voxels to use as features in featureVoxels parameter, featureset (first array cell) ', num2str(iFeat), ' and feature (second array cell) ', num2str(iFeat2), ')\n']);
                    success(:) = 0;
                    results = [];
                    signalSets = []; 
                    return;
                end
                
            end
            
            % TODO: check featureVoxels input is valid instead of failing
            %       half-way through processing. Like if voxel indices are
            %       within array range

        else
            fprintf(2, ['Error: unknown format to define the voxels to use as features in featureVoxels parameter, featureset ', num2str(iFeat), '\n']);
            success(:) = 0;
            results = [];
            signalSets = []; 
            return;
        end
        
        % TODO: check featureVoxels input is valid instead of failing
        %       half-way through processing. Like if voxel indices are
        %       within array range
        
    end
    
    % duplicate the signal-sets for each featureset
    % (to store the signal in and to return later)
    if numFeaturesSets > 1
        signalSetsPrep = repmat(signalSetsPrep, numFeaturesSets, 1);
    end
    
    
    
    %%%
    %% split
    %%%
    [retVal, splitConfig, numSets, ~, ranges] = mx.misc.prepareSplitProcessing(numFeaturesSets, splitConfig);
    if retVal == 0
        success(:) = 0;
        results = [];
        signalSets = []; 
        return;
    end
    

    
    
    %%%
    %% split the datasets into signalSets and process
    %%%
    
    % create the variable to store the cumulated results
    setResults = cell(numSets, 1);
    setSuccess = cell(numSets, 1);
    if nargout > 2
        setSignalSets = cell(numSets, 1);
    end
    
    % set a starttime
    splitTic = tic;
    
    % check if multiple threads are required
    if splitConfig.threads > 1
        % multiple threads
        
        % loop through the sets
        numArguments = nargout;
        parfor i = 1:numSets
            subSignalSets = [];
            
            % message 
            if numSets > 1
                disp(['Started working on set ', num2str(i)]);
            end
            
            % process a subset of featuresets
            if numArguments > 2
                [subResults, subSuccess, subSignalSets] = classifyFmriDataSets_sub( config, ...
                                                                                    featureVoxels(ranges(i, 1):ranges(i, 2)), ...
                                                                                    signalSetsPrep(ranges(i, 1):ranges(i, 2), :), ...
                                                                                    setMethod, ...
                                                                                    data);
            else
                [subResults, subSuccess] = classifyFmriDataSets_sub(config, ...
                                                                    featureVoxels(ranges(i, 1):ranges(i, 2)), ...
                                                                    signalSetsPrep(ranges(i, 1):ranges(i, 2), :), ...
                                                                    setMethod, ...
                                                                    data);
            end
            
            % transfer to subset array
            setResults{i} = subResults;
            setSuccess{i} = subSuccess;
            if numArguments > 2
                setSignalSets{i} = subSignalSets;
            end
            
        end
        clear subResults subSuccess subSignalSets;
        
    else
        % no threading
        
        % loop through the sets
        for i = 1:numSets
            
            % message 
            if numSets > 1
                disp(['Started working on set ', num2str(i)]);
            end
            
            % process a subset of featuresets
            if nargout > 2
                [subResults, subSuccess, subSignalSets] = classifyFmriDataSets_sub( config, ...
                                                                                    featureVoxels(ranges(i, 1):ranges(i, 2)), ...
                                                                                    signalSetsPrep(ranges(i, 1):ranges(i, 2), :), ...
                                                                                    setMethod, ...
                                                                                    data);
            else
                [subResults, subSuccess] = classifyFmriDataSets_sub(config, ...
                                                                    featureVoxels(ranges(i, 1):ranges(i, 2)), ...
                                                                    signalSetsPrep(ranges(i, 1):ranges(i, 2), :), ...
                                                                    setMethod, ...
                                                                    data);
            end
            
            % transfer to subset array
            setResults{i} = subResults;
            setSuccess{i} = subSuccess;
            if nargout > 2
                setSignalSets{i} = subSignalSets;
            end
            clear subResults subSuccess subSignalSets;
            
        end
        
    end
    
    % message
    if strcmpi(splitConfig.silent, 'yes') == 0
        disp(['Processed ', num2str(numFeaturesSets), ' featuresets in ', num2str(toc(splitTic)), ' seconds']);
    end
    
    % loop through the sets
    for i = 1:numSets
        
        % only add a block if there is at least one success (because the first assignment will determine the
        % structure for all to use, and if these are all empty (even in a struct array) then a next block (with at
        % least one result and a different struct array) cannot be assigned
        
        % so check if there is a success in the block (the output result and signalset have a valid structure)
        if nnz(setSuccess{i}) > 0
           
            % set the datapoints in the return variable
            results(ranges(i, 1):ranges(i, 2), :) = setResults{i};
            success(ranges(i, 1):ranges(i, 2)) = setSuccess{i};
            if nargout > 2
                signalSets(ranges(i, 1):ranges(i, 2), :) = setSignalSets{i};
            end
            
        end
        
    end
    
    % check if there are empty result structs missing at the end (this happens when the last
    % subsets fail, so only the flags in the success array will get set, but no result 
    % will be added to the result array), and there is at least one result (required because 
    % else there are no fields defined for the struct and we need those to extend the array)
    if exist('results', 'var') && length(results) < length(success)
        
        % fill out the result- (and signal-)sets with empty structs so they end up have the 
        % same length as the successflag array (matlab does not allow for just extending struct-arrays, 
        % so we copy the first to the length of the successflag array + 1 and remove the last one)
        results(length(success) + 1) = results(1);
        results(length(success) + 1) = [];
        if nargout > 2
            signalSets(length(success) + 1) = signalSets(1);
            signalSets(length(success) + 1) = [];
        end
        
    end
    
    
    % check if there are no results or signalsets
    % and make sure it return something
    if ~exist('results', 'var')
        results = [];
    end
    if nargout > 2 && ~exist('signalSets', 'var')
        signalSets = [];
    end
    
end

% process a subset of featuresets
function [results, success, signalSets] = classifyFmriDataSets_sub(config, featureVoxels, signalSetsPrep, setMethod, data)

    % variables to remove voxels if their value over time is below a given mean
    useMinVoxelMean = isfield(config, 'minVoxelMeanForFeature');
    if useMinVoxelMean
        minVoxelMean = str2num(string(config.minVoxelMeanForFeature));
        if isempty(minVoxelMean)
            minVoxelMean = nan;
            useMinVoxelMean = 0;
        end
    end
    
    % variables to remove features if their value over time is below a given mean
    useMinFeatureMean = isfield(config, 'minFeatureMean');
    if useMinFeatureMean
        minFeatureMean = str2num(string(config.minFeatureMean));
        if isempty(minFeatureMean)
            minFeatureMean = nan;
            useMinFeatureMean = 0;
        end
    end
    
    % determine the number of feature-sets
    numFeaturesSets = length(featureVoxels);
    
    % initialize the success values to zero
    success = zeros(1, numFeaturesSets);

    % loop through the feature-sets
    for iFeat = 1:numFeaturesSets
        
        % store the number of features that are used featureset (before any removing)
        featureSetNumFeatures = length(featureVoxels{iFeat});
        
        % loop through the sets
        for iSet = 1:size(signalSetsPrep, 2)
                
            % check if the featureset is defined as a cell
            if iscell(featureVoxels{iFeat})
                % cell, every cell is a feature (average over the voxels defined in the cell)
                
                % allocate signaldata matrix
                signalSetsPrep(iFeat, iSet).signalData = nan(length(featureVoxels{iFeat}), size(data, 2));
                
                % loop through the features in the featureset
                for iFeat2 = 1:length(featureVoxels{iFeat})
                    
                    % retrieve the data for this feature (by voxel indices)
                    featureData = data(featureVoxels{iFeat}{iFeat2}, :, iSet);
                    
                    % check whether a minimum voxel mean applies
                    if (useMinVoxelMean)
                        
                        % determine if the mean signal of a voxel is below a given threshold or nan
                        emptyVoxels = nanmean(featureData, 2);
                        emptyVoxels = emptyVoxels < minVoxelMean;

                        % remove the empty voxels
                        featureData(emptyVoxels, :) = [];
                        
                    end
                    
                    % remove the rows (voxels) from the feature that are all 0 or nan
                    featureData = featureData(any(featureData, 2), :);
                    
                    % average the voxels to one feature
                    signalSetsPrep(iFeat, iSet).signalData(iFeat2, :) = mean(featureData, 1);
                    
                end
                

            else
                % not cell, every voxel is a seperate feature
                
                % extract the signal using indices of the voxels from the volumes
                signalSetsPrep(iFeat, iSet).signalData = data(featureVoxels{iFeat}, :, iSet);

                % check whether a minimum voxel mean applies
                if (useMinVoxelMean)

                    % determine if the mean signal of a voxel is below a given threshold or nan
                    emptyVoxels = nanmean(signalSetsPrep(iFeat, iSet).signalData, 2);
                    emptyVoxels = emptyVoxels < minVoxelMean;

                    % remove the empty voxels
                    signalSetsPrep(iFeat, iSet).signalData(emptyVoxels, :) = [];

                end


            end
            
            % remove the rows (features/voxels) from the signal that are all 0
            signalSetsPrep(iFeat, iSet).signalData = signalSetsPrep(iFeat, iSet).signalData(any(signalSetsPrep(iFeat, iSet).signalData, 2), :);
                
            % check whether a minimum feature mean applies
            % this applies to the features (which could be an average of voxels or could have a 
            % different minimum) instead of to voxels like above
            if (useMinFeatureMean)

                % determine if the mean signal of a feature is below a given threshold or nan
                emptyVoxels = nanmean(signalSetsPrep(iFeat, iSet).signalData, 2);
                emptyVoxels = emptyVoxels < minFeatureMean;

                % remove the empty voxels
                signalSetsPrep(iFeat, iSet).signalData(emptyVoxels, :) = [];

            end
            
            % check if the signal is completely empty, or whether the signal is less than the minimum number of features
            if isempty(signalSetsPrep(iFeat, iSet).signalData) || size(signalSetsPrep(iFeat, iSet).signalData, 1) < config.minFeatures

                % message
                if strcmpi(config.silent, 'yes') == 0

                    if isempty(signalSetsPrep(iFeat, iSet).signalData)
                        warning(['Warning: less features than the minimum of ', num2str(config.minFeatures)]);
                    else
                        warning('Warning: no signal');
                    end

                end

                % flag as failure and continue
                success(iFeat) = 0;
                continue;
                
            end

            % make the 0 values to nan
            %setsOut(iSet).signalData(setsOut(iSet).signalData == 0) = nan;
            
        end
        
        
        % classify the signal-sets using the feature set
        [resultsOut, successOut, signalSetsOut] = mx.class.classifyFmriSignalSets(config, signalSetsPrep(iFeat, :), setMethod);
        if successOut == 0
            % on failure
            
            % flag feature set as failure and continue
            success(iFeat) = 0;
            continue;
            
        else
            
            % save the results-sets to the return variable
            results(iFeat, :) = resultsOut;
            
            % check if the signals need to be returned
            % (more than 2 return arguments)
            if nargout > 2
                % more than 2 return arguments (returns signalSets)
                
                % save the updated signal-sets to the return variable
                signalSets(iFeat, :) = signalSetsOut;
                
            end
            
            % clear the input signals from the signalsets to save memory %while processing the rest
            for iSet = 1:size(signalSetsPrep, 2)
                signalSetsPrep(iFeat, iSet).signalData = [];
            end
            
            % store the number of input features in the result
            for iResultSet = 1:size(results, 2)
                results(iFeat, iResultSet).numFeatinput = featureSetNumFeatures;
            end
            
            % set featureset as success
            success(iFeat) = 1;
            
        end
        clear resultsOut successOut signalSetsOut;
        
    end     % end feat loop
    
    % check if there are empty result structs missing at the end (this happens when the last
    % featuresets in the set fail, so only the flag in the success array will get set, but no result 
    % will be added to the result array), and there is at least one result (required because 
    % else there are no fields defined for the struct and we need those to extend the array)
    if exist('results', 'var') && length(results) < length(success)
        
        % fill out the result- (and signal-)sets with empty structs so they end up have the 
        % same length as the successflag array (matlab does not allow for just extending struct-arrays, 
        % so we copy the first to the length of the successflag array + 1 and remove the last one)
        results(length(success) + 1, :) = results(1, :);
        results(length(success) + 1, :) = [];
        if nargout > 2
            signalSets(length(success) + 1, :) = signalSets(1, :);
            signalSets(length(success) + 1, :) = [];
        end
        
    end
    
    % check if there are no results or signalsets
    % and make sure it return something
    if ~exist('results', 'var')
        results = [];
    end
    if nargout > 2 && ~exist('signalSets', 'var')
        signalSets = [];
    end
    
end
