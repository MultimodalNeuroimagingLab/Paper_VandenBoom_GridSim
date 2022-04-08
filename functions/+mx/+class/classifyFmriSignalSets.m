%
%   Perform classification based on one or more fMRI signal-sets
%
%   Usage: [results, success, signalSets] = classifyFmriSignalSets(config, signalSets, setMethod)
%
%   config                      = 
%   signalSets                  = struct-array of fMRI signal-sets
%   signalSets(#).signalData    = The signal data to classify on. The first dimension (rows) represent the
%                                 features/voxels and the second dimension (columns) represent the volumes/timepoints
%   signalSets(#).conditions    = The experimental conditions per trial
%   signalSets(#).units         = The units the onsets and durations are in (can either be 'seconds' or 'scans')
%   signalSets(#).onsets        = A vector with the trial onsets (in the specified 'units')
%   signalSets(#).durations     = A vector with the trial durations (in the specified 'units')
%   signalSets(#).onsetsScans   = A vector with the trial onsets in scans
%   setMethod                   = The dataset classification method, either:
%                                   'within': classify per dataset, cross-validation
%                                             will be performed within each set
%                                   'between': train on one dataset, classify on the
%                                              other. Exactly two input sets are expected
%                                   'combine': concatinate all input sets and perform
%                                               cross-validation over the joint set
% 
%   Returns: 
%       results             = struct array with classification result(s). The number of results depends on
%                             the classification method: 'within' will give one result per input dataset;
%                             'combined' and 'between' will give one result per 2 or more datasets
%       success             = flag whether the classification was succesfull (1) or errors have occured (0)
%       signalSets          = struct-array with the input signal-sets, with the adjusted meta-info (onsets, durations,
%                             conditions) and epoched data that were used for classification added.
% 
%
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [results, success, signalSets] = classifyFmriSignalSets(config, signalSets, setMethod)
    results = [];
    
    % flag as classification as success until contrary
    success = 1;
    
    % loop through the sets
    for iSet = 1:length(signalSets)
        
        % make sure the onsets, durations and conditions are in row for
        if ~isrow(signalSets(iSet).conditions),     signalSets(iSet).conditions = signalSets(iSet).conditions';     end
        if ~isrow(signalSets(iSet).onsets),         signalSets(iSet).onsets = signalSets(iSet).onsets';             end
        if ~isrow(signalSets(iSet).onsetsScans),    signalSets(iSet).onsetsScans = signalSets(iSet).onsetsScans';   end
        if ~isrow(signalSets(iSet).durations),      signalSets(iSet).durations = signalSets(iSet).durations';       end

        % store the number of features at the start of the classification of the set
        signalSets(iSet).numFeatStart = size(signalSets(iSet).signalData, 1);
        
        % normalize / detrend the signal (make sure 0s and nans are out of the data)
        signalSets(iSet).signalData = mx.class.cls_normalizeAndDetrendVolumeData(config, signalSets(iSet).signalData, signalSets(iSet).inputTr);
        
        % again remove the rows (voxels) from the signal that are all 0 or nan
        signalSets(iSet).signalData = signalSets(iSet).signalData(any(signalSets(iSet).signalData, 2), :);
        signalSets(iSet).signalData = signalSets(iSet).signalData(any(~isnan(signalSets(iSet).signalData), 2),:);
        
        % store the number of features after detrending, after the removal of 0's and Nans
        signalSets(iSet).numFeatNo0OrNan = size(signalSets(iSet).signalData, 1);
        
        % check if the signal is completely empty, or whether the signal is less than the minimum number of features
        if isempty(signalSets(iSet).signalData) || size(signalSets(iSet).signalData, 1) < config.minFeatures
            
            % message
            if strcmpi(config.silent, 'yes') == 0
                
                if isempty(signalSets(iSet).signalData)
                    warning(['Warning: less features than the minimum of ', num2str(config.minFeatures)]);
                else
                    warning('Warning: no signal');
                end
                
            end
            
            % flag as fail
            success = 0;
            
            return;
            
        end
        
        % message trials
        if strcmpi(config.silent, 'yes') == 0,  disp(['Trials: ', num2str(length(signalSets(iSet).onsetsScans))]);  end
        
        % check if timelock method requires HRF and peak selection from these HRFs
        % also determine if there are trials that need to be skipped at the
        % start (due to config.skipStartVolumes)
        if strcmpi(config.timeLockActive, 'search_highestbold_hrfpeaks') == 1

            % generate the regressors. One where assuming a transient response (a duration of 1) and one assuming a sustained response (with the actual durations)
            if strcmpi(signalSets(iSet).units, 'seconds') == 1
                if isfield(config, 'hrfMethod') && strcmpi(config.hrfMethod, 'new')
                    signalSets(iSet).reg_1scan = mx.misc.genHRFRegressor(   signalSets(iSet).inputTr, ...
                                                                            signalSets(iSet).units, ...
                                                                            signalSets(iSet).onsets, ...
                                                                            repmat(signalSets(iSet).inputTr, 1, length(signalSets(iSet).onsets)), ...
                                                                            signalSets(iSet).inputVolumesDims(4)  );    
                else
                    signalSets(iSet).reg_1scan = mx.misc.genHRFRegressor_old(   signalSets(iSet).inputTr, ...
                                                                                signalSets(iSet).units, ...
                                                                                signalSets(iSet).onsets, ...
                                                                                repmat(signalSets(iSet).inputTr, 1, length(signalSets(iSet).onsets)), ...
                                                                                signalSets(iSet).inputVolumesDims(4)  );    
                end
            else
                if isfield(config, 'hrfMethod') && strcmpi(config.hrfMethod, 'new')
                    signalSets(iSet).reg_1scan = mx.misc.genHRFRegressor(   signalSets(iSet).inputTr, ...
                                                                            signalSets(iSet).units, ...
                                                                            signalSets(iSet).onsets, ...
                                                                            ones(1, length(signalSets(iSet).onsets)), ...
                                                                            signalSets(iSet).inputVolumesDims(4)  );
                else
                    signalSets(iSet).reg_1scan = mx.misc.genHRFRegressor_old(   signalSets(iSet).inputTr, ...
                                                                                signalSets(iSet).units, ...
                                                                                signalSets(iSet).onsets, ...
                                                                                ones(1, length(signalSets(iSet).onsets)), ...
                                                                                signalSets(iSet).inputVolumesDims(4)  );
                end
            end
            if isfield(config, 'hrfMethod') && strcmpi(config.hrfMethod, 'new')
                signalSets(iSet).reg_dur = mx.misc.genHRFRegressor( signalSets(iSet).inputTr, ...
                                                                    signalSets(iSet).units, ...
                                                                    signalSets(iSet).onsets, ...
                                                                    signalSets(iSet).durations, ...
                                                                    signalSets(iSet).inputVolumesDims(4)  );
            else
                signalSets(iSet).reg_dur = mx.misc.genHRFRegressor_old( signalSets(iSet).inputTr, ...
                                                                        signalSets(iSet).units, ...
                                                                        signalSets(iSet).onsets, ...
                                                                        signalSets(iSet).durations, ...
                                                                        signalSets(iSet).inputVolumesDims(4)  );
            end            
            
            % normale hrf between 0 and 1
            %sets(iSet).reg_1scan = (sets(iSet).reg_1scan * (1 / max(sets(iSet).reg_1scan)));
            %sets(iSet).reg_dur = (sets(iSet).reg_dur * (1 / max(sets(iSet).reg_dur)));

            % retrieve the peaks of each trial from the HRF (assuming a transient response = a duration of 1 scan)
            [~, signalSets(iSet).transientPeakIndices, ...
                signalSets(iSet).transientSkippedStartTrials] = mx.class.cls_getFeatureScansFromHrf(config, ...
                                                                                                    signalSets(iSet).conditions, ...
                                                                                                    signalSets(iSet).reg_1scan);

            % retrieve the peaks of each trial from the HRF (assuming a sustained response = based on the true durations of the trials)
            [~, signalSets(iSet).sustainedPeakIndices, ...
                signalSets(iSet).sustainedSkippedStartTrials] = mx.class.cls_getFeatureScansFromHrf(config, ...
                                                                                                    signalSets(iSet).conditions, ...
                                                                                                    signalSets(iSet).reg_dur);

            % remove trials at the start if needed (due to config.skipStartVolumes and cls_getFeatureScansFromHrf)
            signalSets(iSet).skippedStartTrials = max(signalSets(iSet).transientSkippedStartTrials, signalSets(iSet).sustainedSkippedStartTrials);
            if signalSets(iSet).skippedStartTrials > 0
                signalSets(iSet).onsets               = signalSets(iSet).onsets(signalSets(iSet).skippedStartTrials + 1:end);
                signalSets(iSet).onsetsScans          = signalSets(iSet).onsetsScans(signalSets(iSet).skippedStartTrials + 1:end);
                signalSets(iSet).conditions           = signalSets(iSet).conditions(signalSets(iSet).skippedStartTrials + 1:end);
                signalSets(iSet).durations            = signalSets(iSet).durations(signalSets(iSet).skippedStartTrials + 1:end);
            end

            % remove trials at the end if needed (can be because a peak and valley set was missing at the end of the hrf)
            signalSets(iSet).skippedEndTrials = max(length(signalSets(iSet).onsetsScans) - length(signalSets(iSet).transientPeakIndices), length(signalSets(iSet).onsetsScans) - length(signalSets(iSet).sustainedPeakIndices));
            if signalSets(iSet).skippedEndTrials > 0
                signalSets(iSet).onsets               = signalSets(iSet).onsets(1:end - signalSets(iSet).skippedEndTrials);
                signalSets(iSet).onsetsScans          = signalSets(iSet).onsetsScans(1:end - signalSets(iSet).skippedEndTrials);
                signalSets(iSet).conditions           = signalSets(iSet).conditions(1:end - signalSets(iSet).skippedEndTrials);
                signalSets(iSet).durations            = signalSets(iSet).durations(1:end - signalSets(iSet).skippedEndTrials);
                signalSets(iSet).transientPeakIndices = signalSets(iSet).transientPeakIndices(1:length(signalSets(iSet).onsetsScans));
                signalSets(iSet).sustainedPeakIndices = signalSets(iSet).sustainedPeakIndices(1:length(signalSets(iSet).onsetsScans));
            end

        elseif strcmpi(config.timeLockActive, 'static') == 1


            signalSets(iSet).skippedStartTrials = 0;
            % todo: config.skipStartVolumes -> skippedStartTrials

            % don't forget to adjust
            %sets(iSet).onsets
            %sets(iSet).onsetsScans
            %sets(iSet).conditions
            %sets(iSet).durations

        end

        % epoch the dataset
        signalSets(iSet).epochedSignalData = mx.class.cls_epochData(signalSets(iSet).signalData, signalSets(iSet).onsetsScans);

        % check if there are trials which are outliers in terms of number
        % of scans. These could be incomplete trials at the end of a run
        nansPerTrial = sum(isnan(squeeze(signalSets(iSet).epochedSignalData(:,1,:))), 2);
        meanNans = mean(nansPerTrial);
        outlierThreshold = meanNans * 3;
        outlierTrials = find(nansPerTrial > outlierThreshold);
        if length(outlierTrials) > 0
            
            % message
            warning(['Warning: found ', num2str(length(outlierTrials)), ' trials with an outlying length (', num2str(outlierTrials), '), removing those trials from the set']);
            
            % remove those trials
            signalSets(iSet).epochedSignalData(outlierTrials, :, :) = [];
            signalSets(iSet).onsets(outlierTrials) = [];
            signalSets(iSet).onsetsScans(outlierTrials) = [];
            signalSets(iSet).conditions(outlierTrials) = [];
            signalSets(iSet).durations(outlierTrials) = [];
            if strcmpi(config.timeLockActive, 'search_highestbold_hrfpeaks') == 1
                signalSets(iSet).transientPeakIndices(outlierTrials) = [];
                signalSets(iSet).sustainedPeakIndices(outlierTrials) = [];
            end
            
        end
        

        %{
        
        % TODO: remove trials with signal outliers
        % should be added but not now, as it will probably not give much difference and the article needs to be done
        
        % determine the mean and std over all the average signal (over features)
        trialOutlierMean = mean(mean(signalData, 1));
        trialOutlierStd = std(mean(signalData, 1));

        % remove trials that have signal outliers
        outlierTrials = [];
        for iTrial = 1:size(epochedSignalData, 1)
            trialSignal = nanmean(squeeze(epochedSignalData(iTrial,:,:)), 1);
            if any(find(trialSignal > trialOutlierMean + trialOutlierStd * 3)) | any(find(trialSignal < trialOutlierMean - trialOutlierStd * 3))
                % the mean (over features) of one or more timepoints in the
                % trial is outlying, remove the trial

                outlierTrials = [outlierTrials, iTrial];

            end
        end
        if length(outlierTrials) > 0

            % message
            warning(['Warning: found ', num2str(length(outlierTrials)), ' trials with an outlying amplitudes (', num2str(outlierTrials), '), removing those trials from the set']);

            epochedSignalData(outlierTrials, :, :) = [];
            epochedRegTransient(outlierTrials, :, :) = [];
            epochedRegSustained(outlierTrials, :, :) = [];
            onsetsScans(outlierTrials) = [];
            durations(outlierTrials) = [];

        end
        %}

        % check if certain trials need to be isolated (only ones included)
        if strcmpi(config.isolateTrials, 'duration') == 1
            
            % message
            if strcmpi(config.silent, 'yes') == 0,  disp(['Isolating by duration: ', num2str(config.isolateTrialsValue)]); end
            
            % round the durations and create an array of unique durations
            roundedDurations = round(signalSets(iSet).durations);
            
            % determine the trials that do not fit the isolation condition
            %excludedTrials = find(roundedDurations ~= config.isolateTrialsValue);
            excludedTrials = ~ismember(roundedDurations, config.isolateTrialsValue);
            excludedTrials = find(excludedTrials == 1);
            clear roundedDurations;
            
            % remove those trials
            signalSets(iSet).epochedSignalData(excludedTrials, :, :) = [];
            signalSets(iSet).onsets(excludedTrials) = [];
            signalSets(iSet).onsetsScans(excludedTrials) = [];
            signalSets(iSet).conditions(excludedTrials) = [];
            signalSets(iSet).durations(excludedTrials) = [];
            if strcmpi(config.timeLockActive, 'search_highestbold_hrfpeaks') == 1
                signalSets(iSet).transientPeakIndices(excludedTrials) = [];
                signalSets(iSet).sustainedPeakIndices(excludedTrials) = [];
            end
            
        elseif strcmpi(config.isolateTrials, 'condition') == 1

            % message
            if strcmpi(config.silent, 'yes') == 0,  disp(['Isolating by condition: ', num2str(config.isolateTrialsValue)]); end

            % retrieve the conditions from the set (these are updated)
            conditions = signalSets(iSet).conditions;
            
            % determine the trials that do not fit the isolation condition
            excludedTrials = ~ismember(conditions, config.isolateTrialsValue);
            excludedTrials = find(excludedTrials == 1);
            clear conditions;
            
            % remove those trials
            signalSets(iSet).epochedSignalData(excludedTrials, :, :) = [];
            signalSets(iSet).onsets(excludedTrials) = [];
            signalSets(iSet).onsetsScans(excludedTrials) = [];
            signalSets(iSet).conditions(excludedTrials) = [];
            signalSets(iSet).durations(excludedTrials) = [];
            if strcmpi(config.timeLockActive, 'search_highestbold_hrfpeaks') == 1
                signalSets(iSet).transientPeakIndices(excludedTrials) = [];
                signalSets(iSet).sustainedPeakIndices(excludedTrials) = [];
            end
            
        end
        
        % message trials
        if strcmpi(config.silent, 'yes') == 0,  disp(['Trials to classify on: ', num2str(length(signalSets(iSet).onsetsScans))]);    end
        
        % check if timelock method requires a search for the highest bold
        if strcmpi(config.timeLockActive, 'search_highestbold_hrfpeaks') == 1

            % convert the (transient and sustained) peak indices to be relative to
            % the task onset (this the epoch formatting)
            signalSets(iSet).epochedTransientPeakIndices = signalSets(iSet).transientPeakIndices - signalSets(iSet).onsetsScans + 1;
            signalSets(iSet).epochedSustainedPeakIndices = signalSets(iSet).sustainedPeakIndices - signalSets(iSet).onsetsScans + 1;


            % update the HRF regressor peaks to peaks based on the average signal (perform this both for the linear and epoched signals)
            [signalSets(iSet).linearScanIndices, ...
             signalSets(iSet).linearSignalMean] = mx.class.cls_getMeanSignalPeaks_Linear(   config, ...
                                                                                            signalSets(iSet).signalData, ...
                                                                                            signalSets(iSet).transientPeakIndices, ...
                                                                                            signalSets(iSet).sustainedPeakIndices);
            [signalSets(iSet).epochedScanIndices, ...
             signalSets(iSet).epochedSignalMean] = mx.class.cls_getMeanSignalPeaks_Epoch(   config, ...
                                                                                            signalSets(iSet).epochedSignalData, ...
                                                                                            signalSets(iSet).epochedTransientPeakIndices, ...
                                                                                            signalSets(iSet).epochedSustainedPeakIndices);

            % display signal + hrf regressor with peaks
            %mx.class.cls_displayFeatureSelection_Linear(sets(iSet).linearSignalMean, sets(iSet).linearScanIndices, sets(iSet).reg_1scan);
            %sets(iSet).epochedHRF = squeeze(cls_epochData(sets(iSet).reg_1scan, sets(iSet).onsetsScans));
            %mx.class.cls_displayFeatureSelection_Epoch(sets(iSet).epochedSignalMean, sets(iSet).epochedScanIndices, sets(iSet).epochedHRF);


        elseif strcmpi(config.timeLockActive, 'static') == 1

            signalSets(iSet).linearSignalMean = nanmean(signalSets(iSet).signalData);
            signalSets(iSet).linearScanIndices = signalSets(iSet).onsetsScans + config.timeLockActiveStaticVolume - 1;

            signalSets(iSet).epochedSignalMean = squeeze(nanmean(signalSets(iSet).epochedSignalData, 2));
            signalSets(iSet).epochedScanIndices = repmat(config.timeLockActiveStaticVolume, 1, length(signalSets(iSet).onsetsScans));

            %mx.class.cls_displayFeatureSelection_Linear(sets(iSet).linearSignalMean, sets(iSet).linearScanIndices, []);
            %mx.class.cls_displayFeatureSelection_Epoch(sets(iSet).epochedSignalMean, sets(iSet).epochedScanIndices, []);

        end
        

    end
    
    if strcmpi(setMethod, 'within') == 1
        
        % loop through the given sets and classify them seperately
        for iSet = 1:length(signalSets)
            
            % see if there are trials to classiffy on, if not continue
            if size(signalSets(iSet).epochedSignalData, 1) == 0
                
                % message
                warning(['No trials to classify on, skipping set ', num2str(iSet)]);
                
                % skip
                continue;
                
            end
            
            % recode conditions (to ensure they are 0 based)
            recoded = unique(signalSets(iSet).conditions);
            recoded = sort(recoded);
            reConditions = nan(1, length(signalSets(iSet).conditions));
            for r = 1:length(recoded)
                reConditions(signalSets(iSet).conditions == recoded(r)) = r - 1;
            end
            
            [result] = mx.class.classifyFmriSignalSet(  config, 1, ...
                                                        signalSets(iSet).epochedSignalData, signalSets(iSet).epochedScanIndices, ...
                                                        reConditions, signalSets(iSet).durations);
            %[result] = mx.class.classifyFmriSignalSet(   config, 0, ...
                %                                         sets(iSet).signalData, sets(iSet).linearScanIndices, ...
                %                                         reConditions, sets(iSet).durations);
            result.linearScanIndices = signalSets(iSet).linearScanIndices;             % store the linear scan indices in the results (despite using epoched for classification)
            
            % add to the struct-array of results
            if (iSet == 1)
                results = result;           % first assignment to the variable determines the struct fields
            else
                results(iSet) = result;
            end
            
            % store the number of features at the start of the classification 
            % of the set and the number of features after detrending, after the removal of 0's and Nans
            results(iSet).numFeatStart = signalSets(iSet).numFeatStart;
            results(iSet).numFeatNo0OrNan = signalSets(iSet).numFeatNo0OrNan;
            
        end
        
    elseif strcmpi(setMethod, 'combine') == 1
        
        %
        % combine the sets
        %
        
        % determine the dimensions
        maxTrialLength = 0;
        totalTrials = 0;
        maxNumFeatures = 0;
        for iSet = 1:length(signalSets)
            if size(signalSets(iSet).epochedSignalData, 3) > maxTrialLength
               maxTrialLength =  size(signalSets(iSet).epochedSignalData, 3);
            end
            if size(signalSets(iSet).epochedSignalData, 2) > maxNumFeatures
               maxNumFeatures =  size(signalSets(iSet).epochedSignalData, 2);
            end
            totalTrials = totalTrials + size(signalSets(iSet).epochedSignalData, 1);
        end
        
        % create the combined matrices
        epochedSignalData = nan(totalTrials, maxNumFeatures, maxTrialLength);
        epochedScanIndices = [];
        conditions = [];
        durations = [];

        % TODO: make sure the number of features over sets stays the same
        
        
        % fill and combine the matrices
        trialOffset = 1;
        for iSet = 1:length(signalSets)
            %signalData
            %linearScanIndices
            %indicesConditions
            
            % fill the signal data matrix
            epochedSignalData(trialOffset:trialOffset + size(signalSets(iSet).epochedSignalData, 1) - 1, 1:size(signalSets(iSet).epochedSignalData, 2), 1:size(signalSets(iSet).epochedSignalData, 3)) = signalSets(iSet).epochedSignalData;
            
            % raise the trial offset for the signal data matrix
            trialOffset = trialOffset + size(signalSets(iSet).epochedSignalData, 1);

            % other matrices simple append
            epochedScanIndices = [epochedScanIndices signalSets(iSet).epochedScanIndices];
            conditions = [conditions signalSets(iSet).conditions];
            durations = [durations signalSets(iSet).durations];
            
        end
        
        % recode conditions (to ensure they are 0 based)
        recoded = unique(conditions);
        recoded = sort(recoded);
        reConditions = nan(1, length(conditions));
        for r = 1:length(recoded)
            reConditions(conditions == recoded(r)) = r - 1;
        end 
        
        %result = mx.class.classifyFmriSignalSet(config, signalData, linearScanIndices, conditions, durations);
        result = mx.class.classifyFmriSignalSet(config, 1, epochedSignalData, epochedScanIndices, reConditions, durations);
        %result = mx.class.classifyFmriSignalSet(config, 0, epochedSignalData, epochedScanIndices, reConditions, durations);

        results = result;       % first assignment to the variable  determines the struct fields
            
        % store the number of features at the start of the classification 
        % of the set and the number of features after detrending, after the removal of 0's and Nans
        results.numFeatStart = signalSets(1).numFeatStart;
        results.numFeatNo0OrNan = signalSets(1).numFeatNo0OrNan;
        

    elseif strcmpi(setMethod, 'between') == 1
        
        %result = mx.class.classifyFmriSignalSetBetween(config, signalData, linearScanIndices, conditions, durations);
        result = mx.class.classifyFmriSignalSetBetween_Epoch(config, signalSets);
        
        results = result;       % first assignment to the variable determines the struct fields
        
        % store the number of features at the start of the classification 
        % of the set and the number of features after detrending, after the removal of 0's and Nans
        results.numFeatStart = signalSets(1).numFeatStart;
        results.numFeatNo0OrNan = signalSets(1).numFeatNo0OrNan;
        
    else
        fprintf(2, 'Error: classification method, use either ''within'', ''between'' or ''combine''\n');
        success = 0;
        return;
        
    end
    
end