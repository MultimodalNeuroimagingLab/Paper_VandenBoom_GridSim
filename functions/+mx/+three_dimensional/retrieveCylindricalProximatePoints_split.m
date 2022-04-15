%
%   Determine which points in 3D space are close to a set of cylinders in the same space given a specific radius,
%   splitting up the input for speed optimization (multi-threading) or memory maximization
%
%   [proximatePoints, distOnAxis, distFromAxis] = retrieveRadialProximatePoints(retrievalPoints, searchPoints, radius, splitConfig, [usePreExclusion])
%
%   	retrievalPoints 	  = the points in 3D space (n-by-3 matrix) which are being searched in (to be retrieved)
%   	searchCylinders 	  = the cylinders in 3D space (n-by-6 matrix), where each row defines a cylinder within which to search
%   	radius          	  = the radius of the cylinders from which a retrieval point will be considered close
%   	usePreExclusion 	  = [optional, 0 or 1; default = 1) perform pre-exclusion of retrieval-points for each cylinder based on
%       	              	    a search radius. This is beneficial when with a larger number of retrieval-points/search-cylinders, 
%           	          		and/or when limited memory is available; for smaller numbers of retrieval-points and search-cylinders, set
%               	      		this parameter to 0 
% 
%   	splitConfig           = the configuration on splitting the input
%   	splitConfig.numSets   = (optional) split input into x sets of 3D-points [0 = no splitting into sets]
%   	splitConfig.numPoints = (optional) splits the input in sets of x 3D-points [0 = no splitting into sets]
%   	splitConfig.threads   = (optional) number of threads to run the sets on (requires the set to be split using either numSets or numPoints)
%
%   
%   Returns: 
%       proximatePoints       = cell array. Each cell represents an (input) cylinder, the values in the cell are the indices
%                               of the (input) retrieval points which are within the search cylinder's radius
%       distOnAxis            = cell array. Each cell represents an (input) cylinder, the values in the cell are the distances
%                               of the retrieval points on the cylinder axis. So the distance to the start cap point (first three
%                               columns in the searchCylinder parameters) if the point would be projected by the shortest distance onto the axis
%       distFromAxis          = cell array. Each cell represents an (input) cylinder, the values in the cell are the distances
%                               of the retrieval points to the cylinder axis. So the shortest distance of the point to the axis
%
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [proximatePoints, distOnAxis, distFromAxis] = retrieveCylindricalProximatePoints_split(retrievalPoints, searchCylinders, radius, splitConfig, usePreExclusion)
    
% set defaults for optional parameters
    if ~exist('usePreExclusion',  'var') || isempty(usePreExclusion)
        usePreExclusion = 1;
    end
    
    % determine the number of input cylinders
    totalCylinder = size(searchCylinders, 1);
        
    % initialize the cell arrays
    proximatePoints = cell(totalCylinder, 1);
    distOnAxis = cell(totalCylinder, 1);
    distFromAxis = cell(totalCylinder, 1);
    
    % return if there are no cylinders to process
    if totalCylinder == 0, return; end
    
    %%
    % split
    %
    [retVal, splitConfig, numSets, ~, ranges] = mx.misc.prepareSplitProcessing(totalCylinder, splitConfig);
    if retVal == 0
        return;
    end
    
    
    %%
    % function
    %
        
    % create the variable to store the cumulated results
    setCylinders = cell(numSets, 1);
    setDistOnAxis = cell(numSets, 1);
    setDistFromAxis = cell(numSets, 1);
    
    % set a starttime
    splitTic = tic;

    % check if multiple threads are required
    if splitConfig.threads > 1
        % multiple threads

        % loop through the sets
        parfor i = 1:numSets

            % process a subset of cylinders
            [subCylinders, subDistOnAxis, subDistFromAxis] = mx.three_dimensional.retrieveCylindricalProximatePoints(retrievalPoints, searchCylinders(ranges(i, 1):ranges(i, 2), :), radius, usePreExclusion);
            
            % transfer to set
            setCylinders{i} = subCylinders;    
            setDistOnAxis{i} = subDistOnAxis;
            setDistFromAxis{i} = subDistFromAxis;
            
        end
        clear subPoints subDistOnAxis subDistFromAxis;
        
    else
        % no threading
        
        % loop through the sets
        for i = 1:numSets

            % process a subset of cylinders
            [subCylinders, subDistOnAxis, subDistFromAxis] = mx.three_dimensional.retrieveCylindricalProximatePoints(retrievalPoints, searchCylinders(ranges(i, 1):ranges(i, 2), :), radius, usePreExclusion);
            
            % transfer to set
            setCylinders{i} = subCylinders;
            setDistOnAxis{i} = subDistOnAxis;
            setDistFromAxis{i} = subDistFromAxis;
            clear subPoints subDistOnAxis subDistFromAxis;
            
        end
        
    end
    
    % message
    if strcmpi(splitConfig.silent, 'yes') == 0
        disp(['Processed ', num2str(totalCylinder), ' cylinders in ', num2str(toc(splitTic)), ' seconds']);
    end
    
    % loop through the sets and set the datapoints in the return variable
    for i = 1:numSets
        proximatePoints(ranges(i, 1):ranges(i, 2)) = setCylinders{i};
        distOnAxis(ranges(i, 1):ranges(i, 2)) = setDistOnAxis{i};
        distFromAxis(ranges(i, 1):ranges(i, 2)) = setDistFromAxis{i};        
    end
    
end
