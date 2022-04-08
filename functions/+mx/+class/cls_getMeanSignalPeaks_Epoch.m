%
%   Search for the peak of the average signal for each trial (epoched signal input)
%   [signalPeakIndices, signalMean] = cls_getMeanSignalPeaks_Epoch(config, epochedSignal, epochedTransientPeakIndices, epochedSustainedPeakIndices)
%
%   	config 					= the configuration structure which holds the timeLockActiveSearchScansBeforePeak
%								  and timeLockActiveSearchScansAfterPeak fields
%   	epochedSignal 			= 
% 		transientPeakIndices	= 
%		sustainedPeakIndices	=
%
%
%   Returns:
%		signalPeakIndices 	= The indices of the average signal peaks within a certain range per trial. The range
%   						  is determined by the peak of a theoretical transient response and the peak
%   						  of a theoretical sustained response. In addition, the configuration
%   						  might allow for more scans before the transient peak and/or more scans
%   						  after the sustained peak
%		signalMean			= 
%
%
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [signalPeakIndices, signalMean] = cls_getMeanSignalPeaks_Epoch(config, epochedSignal, epochedTransientPeakIndices, epochedSustainedPeakIndices)
    
    %
    % transient-sustained peaks
    %
    
    % check if there are the same number
    if length(epochedTransientPeakIndices) ~= length(epochedSustainedPeakIndices)
        fprintf(2, 'Error: number of transient and sustained peaks should be the same\n');
        return;
    end
    
    % calculate the difference in scans between the sustained and transient
    diffIndices = epochedSustainedPeakIndices - epochedTransientPeakIndices;
    
    % check if there is no sustained peak before a transient peak
    if nnz(diffIndices < 0) ~= 0
        warning('Error: one or more trials have a sustained peak index before the transient peak index, switching those trials');
        
        % switch trials indices
        temp = epochedSustainedPeakIndices;
        epochedSustainedPeakIndices(diffIndices < 0) = epochedTransientPeakIndices(diffIndices < 0);
        epochedTransientPeakIndices(diffIndices < 0) = temp(diffIndices < 0);
        
    end
    
    % calculate the signal mean
    signalMean = squeeze(nanmean(epochedSignal, 2));
    
    % loop through the trials
    signalPeakIndices = [];
    for i = 1:size(signalMean, 1)
        if i > length(epochedTransientPeakIndices) || i > length(epochedTransientPeakIndices)
           continue; 
        end
        
        % extract the trial
        trial = signalMean(i, :);
        trial(isnan(trial)) = [];
        
        % determine the begin and endscan
        startScan = max((epochedTransientPeakIndices(i) - config.timeLockActiveSearchScansBeforePeak), 1);
        endScan = min(epochedSustainedPeakIndices(i) + config.timeLockActiveSearchScansAfterPeak, length(trial));
        
        % search the highest bold value from the transient peak minus x scans till the sustained response plus x scans
        [highestValue, highestValueIndex] = max(trial(startScan:endScan));      
        highestValueIndex = startScan + highestValueIndex - 1;

        % add to signal peak index and condition to the list 
        signalPeakIndices(size(signalPeakIndices, 2) + 1) = highestValueIndex;

    end
    clear highestValue highestValueIndex i;

end
