% 	
%   Determine which points in 3D space are close to a set of cylinders in the same space given a specific radius
%   [proximatePoints, distOnAxis, distFromAxis] = retrieveRadialProximatePoints(retrievalPoints, searchPoints, radius, [usePreExclusion])
%
%   	retrievalPoints   = the points in 3D space (n-by-3 matrix) which are being searched in (to be retrieved)
%   	searchCylinders   = the cylinders in 3D space (n-by-6 matrix), where each row defines a cylinder within which to search
%   	radius            = the radius of the cylinders from which a retrieval point will be considered close
%   	usePreExclusion   = [optional, 0 or 1; default = 1) perform pre-exclusion of retrieval-points for each cylinder based on
%                    		a search radius. This is beneficial when with a larger number of retrieval-points/search-cylinders, 
%                     		and/or when limited memory is available; for smaller numbers of retrieval-points and search-cylinders, set
%                     		this parameter to 0
% 
%
%   Returns: 
%       proximatePoints   = cell array. Each cell represents an (input)
%                           cylinder, the values in the cell are the indices
%                           of the (input) retrieval points which are within
%                           the search cylinder's radius
%       distOnAxis        = cell array. Each cell represents an (input)
%                           cylinder, the values in the cell are the distances
%                           of the retrieval points on the cylinder axis. So
%                           the distance to the start cap point (first three
%                           columns in the searchCylinder parameters) if the point
%                           would be projected by the shortest distance onto the axis
%       distFromAxis      = cell array. Each cell represents an (input)
%                           cylinder, the values in the cell are the distances
%                           of the retrieval points to the cylinder axis. So
%                           the shortest distance of the point to the axis
%
%
%   Copyright, 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [proximatePoints, distOnAxis, distFromAxis] = retrieveCylindricalProximatePoints(retrievalPoints, searchCylinders, radius, usePreExclusion)
    
    % set defaults for optional parameters
    if ~exist('usePreExclusion',  'var') || isempty(usePreExclusion)
        usePreExclusion = 1;
    end
    
    % allocate the return cell array
    n = size(searchCylinders, 1);
    proximatePoints = cell(n, 1);
    distOnAxis = cell(n, 1);
    distFromAxis = cell(n, 1);
    
    % square the radius
    sqRadius = radius * radius; 
    
    % calculate the x, y, and z difference between the start and end of each cylinder (the cylinder axis line)
    cylAxis = searchCylinders(:, 4:6) - searchCylinders(:, 1:3);
    
    % calculate each cylinder axis line's sum of squares 
    ssqAxis = sum(cylAxis' .^ 2);
    
    % calculate the squared distance from the start of each cylinder to each retrieval point
    dist = (searchCylinders(:, 1)' - retrievalPoints(:, 1)) .^ 2 + ...
           (searchCylinders(:, 2)' - retrievalPoints(:, 2)) .^ 2 + ...
           (searchCylinders(:, 3)' - retrievalPoints(:, 3)) .^ 2;
    
    % check if we will perform pre-exclusion
    if usePreExclusion
        % using pre-exclusion and looping through the cylinders

        % determine which points are within the reach of the cylinder (assuming a sphere)
        close = dist < (ssqAxis + sqRadius);

        % loop through each cylinder
        for iCyl = 1:n

            % retrieve the indices relevant for the current cylinder
            % (points within the cylinder's search radius)
            indices = find(close(:, iCyl));

            % calculate the dot product of the difference between the start of the
            % current cylinder and the difference between the start of the current
            % cylinder and each point
            dot = sum((retrievalPoints(indices, :) - searchCylinders(iCyl, 1:3)) .* cylAxis(iCyl, :), 2);
            
            % calculate the distance of each point to current cylinder axis
            dist2 = dist(indices, iCyl) - (dot .^ 2) ./ ssqAxis(iCyl);
            
            % determine which points are within the current cylinder's caps and within the radius
            within = dot >= 0 & dot <= ssqAxis(iCyl) & dist2 < sqRadius;
            
            % store the indices of the points within the cylinder and their metrics
            proximatePoints{iCyl} = indices(within);
            distOnAxis{iCyl} = dot(within) / sqrt(ssqAxis(iCyl));
            distFromAxis{iCyl} = sqrt(dist2(within));

        end
        
    else
        % not using pre-exclusion, process all data in one go
        
        % calculate dot product of the difference between the start of each
        % cylinder and the difference between the start of each cylinder and each point
        %{
        % commented out method takes more memory
        pntDistX = retrievalPoints(:, 1) - searchCylinders(:, 1)';
        pntDistY = retrievalPoints(:, 2) - searchCylinders(:, 2)';
        pntDistZ = retrievalPoints(:, 3) - searchCylinders(:, 3)';
        dot = pntDistX .* cylAxis(:, 1)' + pntDistY .* cylAxis(:, 2)' + pntDistZ .* cylAxis(:, 3)';
        %}
        dot = (retrievalPoints(:, 1) - searchCylinders(:, 1)') .* cylAxis(:, 1)' + ...
              (retrievalPoints(:, 2) - searchCylinders(:, 2)') .* cylAxis(:, 2)' + ...
              (retrievalPoints(:, 3) - searchCylinders(:, 3)') .* cylAxis(:, 3)';
        
        % calculate the distance of each point to each cylinder's axis
        dist = dist - (dot .^ 2) ./ ssqAxis;
        
        % determine which points are within each's cylinder caps and within the radius
        within = dot >= 0 & dot <= ssqAxis & dist < sqRadius;
  
        % loop through each cylinder
        for iCyl = 1:n
            
            % store the points within and their metrics
            proximatePoints{iCyl} = find(within(:, iCyl));
            distOnAxis{iCyl} = dot(proximatePoints{iCyl}, iCyl) / sqrt(ssqAxis(iCyl));
            distFromAxis{iCyl} = sqrt(dist(proximatePoints{iCyl}, iCyl));
            
        end
        
    end
       
end
