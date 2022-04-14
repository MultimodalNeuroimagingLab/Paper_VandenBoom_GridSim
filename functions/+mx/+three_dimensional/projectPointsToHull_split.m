%
%   Project points onto a 3D object, splitting up the input for speed
%   optimization (multi-threading) or memory maximization
%
%   [points, pointsProjected, pointsTriangles] = projectPointsToHull_split(hull, points, normdist, intersval, splitConfig)
%
%   	hull                  = the hull onto which the vertices are projected
%   	points                = the 3D points to project onto the hull (n-by-3 matrix)
%   	normdist              = the radius at which, from the electrode, vertices are
%       	                    included to generate an average normal to determine the
%           	                intersection
%   	intersval             = the search distance for intersection
%
%   	splitConfig           = the configuration on splitting the input
%   	splitConfig.numSets   = (optional) split input into x sets of 3D-points [0 = no splitting into sets]
%   	splitConfig.numPoints = (optional) splits the input in sets of x 3D-points [0 = no splitting into sets]
%   	splitConfig.threads   = (optional) number of threads to run the sets on (requires the
%       	                    set to be split using either numSets or numPoints)
% 
%
%   Returns: 
%       points           	= the new 3D points projected onto the hull (n-by-3 matrix).
%                             Unprojectable points are left in their original input position
%       pointsProjected   	= a flag for each point whether it was projected
%                             succesfully (1 if successfull, 0 if not)
%       pointsTriangles   	= the index of the hull triangle onto which the point is projected
%                             Unprojectable points will have NaN
%
%
%   Copyright 2019, Max van den Boom
%

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [points, pointsProjected, pointsTriangles] = projectPointsToHull_split(hull, points, normdist, intersval, splitConfig)
    
    % determine the number of input points
    totalPoints = size(points, 1);

    % preallocate the output arrays
    pointsProjected = false(totalPoints, 1);
    pointsTriangles = NaN(totalPoints, 1);

    % return if there are no points to process
    if totalPoints == 0, return; end
    
    %%%
    %% split
    %%%
    [retVal, splitConfig, numSets, ~, ranges] = mx.misc.prepareSplitProcessing(totalPoints, splitConfig);
    if retVal == 0
        return;
    end
    
    
    %%%
    %% function
    %%%
    
    
    % create the variable to store the cumulated results
    setPoints = cell(numSets, 1);
    setPointsProjected = cell(numSets, 1);
    setPointsTriangles = cell(numSets, 1);

    % set a starttime
    splitTic = tic;

    % check if multiple threads are required
    if splitConfig.threads > 1
        % multiple threads

        % loop through the sets
        parfor i = 1:numSets

            % process a subset of points
            [subPoints, subPointsProjected, subPointsTriangles] = ...
                mx.three_dimensional.projectPointsToHull(   hull, ...
                                                            points(ranges(i, 1):ranges(i, 2), :), ...
                                                            normdist, ...
                                                            intersval);

            % transfer to set
            setPoints{i} = subPoints;
            setPointsProjected{i} = subPointsProjected;
            setPointsTriangles{i} = subPointsTriangles;
            
        end
        clear subPoints subPointsProjected subPointsTriangles;
        
    else
        % no threading
        
        % loop through the sets
        for i = 1:numSets

            % process a subset of points
            [subPoints, subPointsProjected, subPointsTriangles] = ...
                mx.three_dimensional.projectPointsToHull(   hull, ...
                                                            points(ranges(i, 1):ranges(i, 2), :), ...
                                                            normdist, ...
                                                            intersval);

            % transfer to set
            setPoints{i} = subPoints;
            setPointsProjected{i} = subPointsProjected;
            setPointsTriangles{i} = subPointsTriangles;
            clear subPoints subPointsProjected subPointsTriangles;
            
        end
        
    end
    
    % message
    disp(['Processed ', num2str(totalPoints), ' points in ', num2str(toc(splitTic)), ' seconds']);
    
    % loop through the sets and set outputs in the return variables
    for i = 1:numSets
        points(ranges(i, 1):ranges(i, 2), :) = setPoints{i};
        pointsProjected(ranges(i, 1):ranges(i, 2)) = setPointsProjected{i};
        pointsTriangles(ranges(i, 1):ranges(i, 2)) = setPointsTriangles{i};
    end
    
end
