%
%   Search for the peak of the average signal for each trial (linear signal input)
%   [signalPeakIndices, signalMean] = cls_getMeanSignalPeaks_Linear(config, signal, transientPeakIndices, sustainedPeakIndices)
%
%   	config 					= the configuration structure which holds the timeLockActiveSearchScansBeforePeak
%   						      and timeLockActiveSearchScansAfterPeak fields
%   	signal 					= data matrix to search for the peak signal. The first dimension (rows) represents the 
%            					  features/voxels and the second dimension (columns) represent the volumes/timepoints
% 		transientPeakIndices	=
%		sustainedPeakIndices	=
%
%
%   Returns:
%		signalPeakIndices 	= The indices of the average signal peaks within a certain range per
%							  trial. The range is determined by the peak of a theoretical transient
%							  response and the peak of a theoretical sustained response. In addition, the
%							  configuration might allow for more scans before the transient peak and/or more
%   						  scans after the sustained peak
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
function [signalPeakIndices, signalMean] = cls_getMeanSignalPeaks_Linear(config, signal, transientPeakIndices, sustainedPeakIndices)
    
    %
    % transient-sustained peaks
    %
    
    % check if there are the same number
    if length(transientPeakIndices) ~= length(sustainedPeakIndices)
        fprintf(2, 'Error: number of transient and sustained peaks should be the same\n');
        return;
    end
    
    % calculate the difference in scans between the sustained and transient
    diffIndices = sustainedPeakIndices - transientPeakIndices;
    
    % check if there is no sustained peak before a transient peak
    if nnz(diffIndices < 0) ~= 0
        fprintf(2, 'Error: one or more trials have a sustained peak index before the transient peak index, this should not be possible\n');
        return;
    end
    
    % calculate the signal mean
    signalMean = mean(signal);

    signalPeakIndices = [];
    for i = 1:length(transientPeakIndices)
        scanIndex = transientPeakIndices(i);

        % search the highest bold value from the transient peak minus x scans till the sustained response plus x scans
        highestValue = nan;
        highestValueIndex = nan;
        for j = scanIndex - config.timeLockActiveSearchScansBeforePeak:scanIndex + diffIndices(i) + config.timeLockActiveSearchScansAfterPeak
           if isnan(highestValue) || (j < length(signalMean) && signalMean(j) > highestValue)
               highestValue = signalMean(j);
               highestValueIndex = j;
           end
        end

        % add to signal peak index and condition to the list 
        signalPeakIndices(size(signalPeakIndices, 2) + 1) = highestValueIndex;

    end
    clear highestValue highestValueIndex skiploop scanIndex i j;



end