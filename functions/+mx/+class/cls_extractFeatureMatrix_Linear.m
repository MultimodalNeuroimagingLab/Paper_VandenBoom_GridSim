%
%   Extract the feature matrix from linear data  according to the given configuration.
%   featureData = cls_extractFeatureMatrix_Linear(config, signal, scanIndices)
%
%   config 			= the configuration structure with the fields:
%       				timeLockActiveVolumeOffsets =
%        				featureTimeMethod           =
%   signal			= 
%	scanIndices		= 
% 
%
%   Returns:
%       featureData = the feature data as a matrix where the first
%                     dimension (rows) represents the features and
%                     the second dimension (columns) the trials
%
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function featureData = cls_extractFeatureMatrix_Linear(config, signal, scanIndices)


    % todo: if needed create baseline feature matrix


    % create a active feature matrix where the rows represent the features and columns the trials
    featureData = [];
    
    % loop through the trials
    for i = 1:length(scanIndices)
        
        % determine the scans to use as active
        trialScans = config.timeLockActiveVolumeOffsets + scanIndices(i);
        if nnz(trialScans <= 0) > 0
            warning('Warning: the current config.timeLockActiveVolumeOffsets results in scans below 1, dismissing those scans');
            trialScans(trialScans <= 0) = [];
        end
        if nnz(trialScans > size(signal, 2))
            warning('Warning: the current config.timeLockActiveVolumeOffsets results in scans larger than the number of scans/timepoints, dismissing those scans');
            trialScans(trialScans > size(signal, 2)) = [];
        end
        if length(trialScans) == 0
            fprintf(2, 'Error: No active scans for the trial\n');
            return; 
        end
        
        % retrieve the active trials
        trialData = signal(:, trialScans);
        
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
