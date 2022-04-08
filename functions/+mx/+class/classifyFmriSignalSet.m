%
%   Classify an epoched or linear signal-set of trials using leave-one-out or split validation
%   classifyFmriSignalSet(config, signal, scanIndices, conditions, durations) classifies an entire dataset
%
%   epoched                 = whether the incoming signal and scanIndices are epoched [1] or linear [0]
%   signal                  = If the epoched argument is true [1] then this parameter
%                             should be a matrix with the epoched signal data for
%                             classification (the first dimension should be the trials, the
%                             second dimension should be the features and the third dimension
%                             the scans/timepoints within each trial).
%                             If the epoched argument is false [0] then this parameter
%                             should be a matrix with the linear signal data (the first dimension
%                             should be the features and the second dimension the scans/timepoints)
%                        
%   scanIndices             = vector with the indices of the scan/timepoints to be used for classification
%                             (if epoched, then the index is the offset per trial; if linear then the offset
%                              is for the entire series of scans/timepoints)
%   conditions              = vector with the condition for the trials
%   durations               = vector with the durations for the trials
%
%   config                  = the configuration for the classification
%   config.crossValidation  = the method ['LeaveOneOut' or 'Split' (not implemented yet)]
%
%   Returns:
%       result              = the result of the classification as a signalset result struct
%       
%
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
% 
function [result] = classifyFmriSignalSet(config, epoched, signal, scanIndices, conditions, durations)
    if epoched ~= 1,    epoched = 0;    end

    % create a signal-set result struct to return
    result = mx.class.fmriSignalSetResult();
    
    % set epoched or linear
    if epoched == 1        
        result.signalFormatting = 'epoched';
        result.epochedScanIndices = scanIndices;
    else
        result.signalFormatting = 'linear';
        result.linearScanIndices = scanIndices;
    end
    result.conditions = conditions;
    result.durations = durations;
    
    % create a percentage correct counter for every condition
    correctPerCondition = zeros(max(conditions) + 1, max(conditions) + 1);
    
    % round the durations and create an array of unique durations
    roundedDurations = round(durations);
    uniqueDurations = unique(roundedDurations);
    uniqueDurations = sort(uniqueDurations, 'ascend');
    
    % create an array to store all the duration results (duration, #correct, #incorrect, boldwave)
    resDurations{length(uniqueDurations), 5} = [];
    
    % create an array to hold the indices of the unique values
    roundedDurationIndices = nan(1, length(roundedDurations));
    
    % loop through the unique durations
    for i = 1:length(uniqueDurations)
        
        % store the unique duration
        resDurations{i, 1} = uniqueDurations(i);
        resDurations{i, 2} = 0;
        resDurations{i, 3} = 0;
        resDurations{i, 4} = 0;
        resDurations{i, 5} = {};
        
        % store the indices of the duration in the unique duration array
        roundedDurationIndices(roundedDurations == uniqueDurations(i)) = i;
        
    end
    
    % create the feature matrix
    if epoched == 1  
        featureData = mx.class.cls_extractFeatureMatrix_Epoch(config, signal, scanIndices);
    else
        featureData = mx.class.cls_extractFeatureMatrix_Linear(config, signal, scanIndices);
    end
    
    % set the number of features used in the result
    result.numFeatUsed = size(featureData, 1);
    
    % select the cross validation method
    if strcmpi(config.crossValidation, 'LeaveOneOut') == 1
        
        % loop through the trials
        result.numCorrect = 0;
        result.numIncorrect = 0;
        result.arrCorrect = nan(1, size(scanIndices, 2));
        for i = 1:size(scanIndices, 2)
            
            % copy the trials
            % and remove the test trial from the matrix
            trainData = featureData;
            trainConditions = conditions;
            trainData(:,i) = [];
            trainConditions(:,i) = [];

            % retrieve the test trial
            testData = featureData(:,i);
            testCondition = conditions(i);


            % classify the single trial
            [predictedLabel, ~, ~, ~] = mx.class.cls_classifyTrialData(config, trainData, trainConditions, testData, testCondition);

            % display
            %disp(['trial ', num2str(i), '  -  value: ', num2str(testTrialCondition), '  -  classified as: ', num2str(predictedLabel)]);

            % check if correct
            if (testCondition == predictedLabel)
                
                % flag this trial as correct
                result.arrCorrect(i) = 1;

                % count correct for the run
                result.numCorrect = result.numCorrect + 1;

                % count correct for the duration
                resDurations{roundedDurationIndices(i), 2} = resDurations{roundedDurationIndices(i), 2} + 1;

            else

                % flag this trial as incorrect
                result.arrCorrect(i) = 0;
                
                % count incorrect for the run
                result.numIncorrect = result.numIncorrect + 1;

                % count incorrect for the duration
                resDurations{roundedDurationIndices(i), 3} = resDurations{roundedDurationIndices(i), 3} + 1;

            end

            % count prediction per condition
            correctPerCondition(testCondition + 1, predictedLabel + 1) = correctPerCondition(testCondition + 1, predictedLabel + 1) + 1;

        end

    elseif strcmpi(config.crossValidation, 'Split') == 1
        
        % todo: implement split cross validation

    else

        % message
        fprintf(2, ['Unknown cross validation method: ', config.crossValidation, '\n']);
        return;

    end
    
    % store the #correct and #incorrect in the results (per condition)
    result.correctPerCondition = correctPerCondition;
    
    % calculate and store the percentage correct
    result.accTotal = result.numCorrect / (result.numCorrect + result.numIncorrect) * 100;
    
    % calculate the percentage correct per condition (row)
    for i=1:size(correctPerCondition, 1)
        for j=1:size(correctPerCondition, 2)
            correctPerCondition(i, j) = correctPerCondition(i, j) / nnz(conditions == (i - 1)) * 100;
        end
    end
    
    % calculate the percentage correct per duration (row)
    for i=1:size(resDurations, 1)
        resDurations{i, 4} = resDurations{i, 2} / (resDurations{i, 2} + resDurations{i, 3}) * 100;
    end
    
    % store the percentage correct per condition
    result.accPerCondition = correctPerCondition;
        
    % store the duration result matrix
    result.durationResults = resDurations;
    
end
