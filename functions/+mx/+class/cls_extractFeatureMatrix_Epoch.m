%
%   Extract the feature matrix from epoched data according to the given configuration.
%   featureData = cls_extractFeatureMatrix_Epoch(config, epochedSignal, epochedScanIndices)
%
%   config 				= the configuration structure with the fields:
%       					timeLockActiveVolumeOffsets =
%       					featureTimeMethod           =
%	epochedSignal		=
%	epochedScanIndices  = 
%
%   Returns:
%       featureData 	= the feature data as a matrix where the first
%                     	  dimension (rows) represents the features and
%                         the second dimension (columns) the trials
%
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function featureData = cls_extractFeatureMatrix_Epoch(config, epochedSignal, epochedScanIndices)


    % todo: if needed create baseline feature matrix


    % create a active feature matrix where the rows represent the features and columns the trials
    featureData = [];
    for i = 1:length(epochedScanIndices)
        
        % get the trial
        trial = squeeze(epochedSignal(i,:,:));
        
        % check if there is just a single feature
        if size(epochedSignal, 2) == 1
            
            % transpose so the columns are time and rows are features
            trial = trial';
            
        end
        
        % remove the timepoints which are only nans
        trial = trial(:, any(~isnan(trial), 1));
        
        % determine the scans to use as active
        trialScans = config.timeLockActiveVolumeOffsets + epochedScanIndices(i);
        if nnz(trialScans <= 0) > 0
            warning('Warning: the current config.timeLockActiveVolumeOffsets results in scans before the trial, dismissing those scans');
            trialScans(trialScans <= 0) = [];
        end
        if nnz(trialScans > size(trial, 2))
            warning('Warning: the current config.timeLockActiveVolumeOffsets results in scans after the trial, dismissing those scans');
            trialScans(trialScans > size(trial, 2)) = [];
        end
        if length(trialScans) == 0
            fprintf(2, 'Error: No active scans for the trial\n');
            return; 
        end
        
        % retrieve the active trials
        trialData = trial(:, trialScans);
        
        % summarize the timepoint to one feature vector if needed
        if size(trialData, 2) > 1
            
            if strcmpi(config.featureTimeMethod, 'average')
                
                % average the timepoints
                trialData = mean(trialData, 2);
                
            else
                
                % make all values sequential
                trialData = trialData(:);
                
            end
        end
        
        % add to the feature matrix
        featureData = [featureData, trialData];
        
    end
    
end
