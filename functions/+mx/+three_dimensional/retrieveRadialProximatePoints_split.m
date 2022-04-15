%
%   Determine which points in 3D space are close to another set of points in the same space given a specific radius,
%   splitting up the input for speed optimization (multi-threading) or memory maximization
%
%   [proximatePoints] = retrieveRadialProximatePoints(retrievalPoints, searchPoints, radius, splitConfig)
%
%   	retrievalPoints       = the points in 3D space (n-by-3 matrix) which are being searched in (to be retrieved)
%   	searchPoints          = the points in 3D space (n-by-3 matrix), where each point is a source in space from which to search
%   	radius                = the radius within each search point where a retrieval point will be considered close    
% 
%   	splitConfig           = the configuration on splitting the input
%   	splitConfig.numSets   = (optional) split input into x sets of 3D-points [0 = no splitting into sets]
%   	splitConfig.numPoints = (optional) splits the input in sets of x 3D-points [0 = no splitting into sets]
%   	splitConfig.threads   = (optional) number of threads to run the sets on (requires the
%                               set to be split using either numSets or numPoints)
%
%   Returns: 
%       proximatePoints   	  = cell array. Each cell represents an (input) search point, the values in the cell are the indices
%                               of the (input) retrieval points which are within the search point's radius
%   
%
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [proximatePoints] = retrieveRadialProximatePoints_split(retrievalPoints, searchPoints, radius, splitConfig)

    % determine the number of input points
    totalPoints = size(searchPoints, 1);
        
    % initialize the cell array
    proximatePoints = cell(totalPoints, 1);
    
    % return if there are no points to process
    if totalPoints == 0, return; end
    
    %%
    % split
    %
    [retVal, splitConfig, numSets, ~, ranges] = mx.misc.prepareSplitProcessing(totalPoints, splitConfig);
    if retVal == 0
        return;
    end
    
    
    %%
    % function
    %
        
    % create the variable to store the cumulated results
    setPoints = cell(numSets, 1);
    
    % set a starttime
    splitTic = tic;

    % check if multiple threads are required
    if splitConfig.threads > 1
        % multiple threads

        % loop through the sets
        parfor i = 1:numSets

            % process a subset of points
            [subPoints] = mx.three_dimensional.retrieveRadialProximatePoints(retrievalPoints, searchPoints(ranges(i, 1):ranges(i, 2), :), radius);
            
            % transfer to set
            setPoints{i} = subPoints;    

        end
        clear subPoints;
        
    else
        % no threading
        
        % loop through the sets
        for i = 1:numSets

            % process a subset of points
            [subPoints] = mx.three_dimensional.retrieveRadialProximatePoints(retrievalPoints, searchPoints(ranges(i, 1):ranges(i, 2), :), radius);
            
            % transfer to set
            setPoints{i} = subPoints;
            clear subPoints;
            
        end
        
    end
    
    % message
    if strcmpi(splitConfig.silent, 'yes') == 0
        disp(['Processed ', num2str(totalPoints), ' points in ', num2str(toc(splitTic)), ' seconds']);
    end
    
    % loop through the sets and set the datapoints in the return variable
    for i = 1:numSets
        proximatePoints(ranges(i, 1):ranges(i, 2)) = setPoints{i};
    end
    
end
