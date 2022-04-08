%
%   Retrieve the peaks and valley indices from a hrf vector
%	[valleyIndices, peakIndices, skippedStartTrials] = cls_getFeatureScansFromHrf(config, conditions, hrf)
%   
%       config			 = the classification configuration struct
%       conditions       = the stimulus conditions
%       hrf           	 = the HRF vector to extract the peaks and valleys from
%
%
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [valleyIndices, peakIndices, skippedStartTrials] = cls_getFeatureScansFromHrf(config, conditions, hrf)
    
    % find the locations of peaks in the HRF model
    [~, peakIndices] = findpeaks(hrf);
    
    % check if there at least one peak
    if length(peakIndices) == 0
        fprintf(2, 'Error: No peaks were found, check hrf\n');
        return;
    end
    
    % find the locations of valleys in the HRF model
    [~, valleyIndices] = findpeaks(max(hrf) - hrf);
    
    % remove any valley that comes after the last peak
    valleyIndices(valleyIndices > peakIndices(end)) = [];
    
    % find the lowest point before the first peak
    % and add the lowest point before the first peak as the first valley
    lowestPointIndex = -1;
    lowestPointValue = 0;
    for i = peakIndices(1):-1:1
        if lowestPointIndex == -1 || hrf(i) < lowestPointValue
            lowestPointIndex = i;
            lowestPointValue = hrf(i);
        end
    end
    valleyIndices = horzcat(lowestPointIndex, valleyIndices);
    clear lowestPointIndex lowestPointValue;

    % check if the number of valleys and peaks match
    if length(peakIndices) ~= length(valleyIndices)
        disp(['Number of peaks: ', num2str(length(peakIndices))]);
        disp(['Number of valleys: ', num2str(length(valleyIndices))]);
        fprintf(2, 'Error: Number of peaks and valleys in hrf do not match\n');
        return;
    end
    
    % check if there are more peaks than conditions
    if length(peakIndices) > length(conditions)
        
        % error
        disp(['Number of peaks: ', num2str(length(peakIndices))]);
        disp(['Number of conditions: ', num2str(length(conditions))]);
        fprintf(2, 'Error: More peaks in hrf than there are number of conditions\n');
        return;
        
    end
    
    % check if there are more conditions than hrf peaks
    if length(peakIndices) < length(conditions)
        
        % warning
        disp(['Number of peaks: ', num2str(length(peakIndices))]);
        disp(['Number of conditions: ', num2str(length(conditions))]);
        warning('Warning: More conditions than peaks in hrf. This can occur when scanner acquisition stopped early (giving a smaller amount of volumes to build a HRF on), while the task was still running (trials/conditions were still being added)');
        
    end
    
    % if peak or valleys are inside of the skipped volumes then discard the peak and condition
    skippedStartTrials = nnz(peakIndices <= config.skipStartVolumes | valleyIndices <= config.skipStartVolumes);
    for i = 1:skippedStartTrials
       peakIndices = peakIndices(2:end);
       valleyIndices = valleyIndices(2:end);
    end

end

