%
%   Epoch a dataset based on volume onsets
%   epochedData = cls_epochData(data, onsets) 
% 
%   	data 			= data matrix to epoch. The rows represent the voxels and the
%  		        		  columns represent the volumes/timepoints
%   	onsets 			= vector holding the onsets (volume numbers, 1-based)
%
%
%   Returns:
%		epochedData		= The trials from a dataset based on the given volume onsets.
%						  The output matrix is in the format <trials> x <voxels> x <volumes>
%
%
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [epochedData] = cls_epochData(data, onsets)

    % remove onsets which are higher than the number of volumes in the data
    onsets(onsets > size(data, 2)) = [];
    
    % ensure the onsets are rows
    if ~isrow(onsets), onsets = onsets'; end
    
    % find the longest interval between the onsets
    if length(onsets) == 1
        longestInterval = size(data, 2) - onsets;
    else
        longestInterval = max([onsets(2:end) - onsets(1:end-1), size(data, 2) - onsets(end) + 1]);
    end
    if longestInterval == 0
        fprintf(2, 'Error: could not find an interval between the trials with more than 0 scans\n');
        return;
    end
    
    % create the output matrix
    epochedData = nan(length(onsets), size(data, 1), longestInterval);
    
    % loop over the trials
    for i = 1:length(onsets)
        if i == length(onsets)
            len = size(data, 2) - onsets(end) + 1;
        else
            len = onsets(i + 1) - onsets(i);
        end
        epochedData(i, :, 1:len) = data(:, onsets(i):onsets(i) + len - 1);
    end
    
end







