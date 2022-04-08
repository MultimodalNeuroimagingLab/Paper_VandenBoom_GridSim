%
%   Retrieve metrics from a 3D object
%   [metrics] = get3DMetrics(input, ...)
%   [metrics] = get3DMetrics(vertices, faces, ...)
%
%   input             = the input 3D to extract the metrics from
%   ...               = [optional] only calculate specific metric(s):
%                       'edges'
%                       'averageEdge'
%                       'minEdge'
%                       'maxEdge'
%                       'areas'
%                       'averageArea'
%                       'minArea'
%                       'maxArea'
%                       'totalArea'
%                       
%
%
%   Returns: 
%       metrics.edgeLengths         = the edge lengths
%       metrics.averageEdgeLength   = the average edge length
%       metrics.minEdgeLength       = the smallest edge length
%       metrics.maxEdgeLength       = the largest edge length
%       metrics.areas               = the areas
%       metrics.averageArea         = the average area
%       metrics.minArea             = the smallest area
%       metrics.maxArea             = the largest area
%       metrics.totalArea           = the total area
%
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [metrics] = get3DMetrics(varargin)

    % 
    specifyMetrics = 0;         % whether metric arguments are specified
    mEdgeLengths = 0;           % flag to return the edge lengths
    mAverageEdgeLength = 0;     % flag to return the average edge length
    mMinEdgeLength = 0;         % flag to return the minimum edge length
    mMaxEdgeLength = 0;         % flag to return the maximum edge length
    mAreas = 0;                 % flag to return the areas
    mAverageArea = 0;           % flag to return the average area
    mMinArea = 0;               % flag to return the minimum area
    mMaxArea = 0;               % flag to return the maximum area
    mTotalArea = 0;             % flot to return the total area
    
    % 
    vertices = [];
    faces = [];
    
    % loop through the input arguments
    for i = 1:nargin
        if isstruct(varargin{i}) || isobject(varargin{i})
            % if struct or object
            
            if isa(varargin{i}, 'gifti')
                % gifti object
                
                vertices = varargin{i}.vertices;
                faces = varargin{i}.faces;
                
            elseif isstruct(varargin{i})
                % struct
                
                if isfield(varargin{i}, 'vertices')
                    vertices = varargin{i}.vertices;
                elseif isfield(varargin{i}, 'vert')
                    vertices = varargin{i}.vert;
                end
                
                if isfield(varargin{i}, 'faces')
                    faces = varargin{i}.faces;
                elseif isfield(varargin{i}, 'face')
                    faces = varargin{i}.face;
                elseif isfield(varargin{i}, 'tri')
                    faces = varargin{i}.tri;
                end
                
            end
            
        elseif size(varargin{i}, 1) == 3 || size(varargin{i}, 2) == 3
            % if vector/matrix where the size of one dimenion is 3
            
            if i == 1
                % first argument should be vertices
                
                vertices = varargin{i};
                
            elseif i == 2
                % second argument should be faces
                
                faces = varargin{i};
                
            end

        elseif isa(varargin{i}, 'char')
            % check if the argument is a string (char array)
            
            % edge metrics
            if strcmpi(varargin{i}, 'edges') || strcmpi(varargin{i}, 'edgeLengths')
                specifyMetrics = 1;    
                mEdgeLengths = 1;
            end
            if strcmpi(varargin{i}, 'averageEdge') || strcmpi(varargin{i}, 'averageEdgeLength')
                specifyMetrics = 1;    
                mAverageEdgeLength = 1;
            end
            if strcmpi(varargin{i}, 'minEdge') || strcmpi(varargin{i}, 'minEdgeLength')
                specifyMetrics = 1;    
                mMinEdgeLength = 1;
            end
            if strcmpi(varargin{i}, 'maxEdge') || strcmpi(varargin{i}, 'maxEdgeLength')
                specifyMetrics = 1;    
                mMaxEdgeLength = 1;
            end
            
            % area metrics
            if strcmpi(varargin{i}, 'areas')
                specifyMetrics = 1;    
                mAreas = 1;
            end
            if strcmpi(varargin{i}, 'averageArea')
                specifyMetrics = 1;
                mAverageArea = 1;
            end
            if strcmpi(varargin{i}, 'minArea')
                specifyMetrics = 1;
                mMinArea = 1;
            end
            if strcmpi(varargin{i}, 'maxArea')
                specifyMetrics = 1;
                mMaxArea = 1;
            end
            if strcmpi(varargin{i}, 'totalArea')
                specifyMetrics = 1;
                mTotalArea = 1;
            end
            
        end
        
    end
    
    % check if there are vertices and faces
    if isempty(vertices) || isempty(faces)
        fprintf(2, 'Error: vertices and faces found. Make sure either an 3D object, or vertices and faces are given as first arguments\n');
        return;
    end
    
    % convert the vertices to double for increased precision
    vertices = double(vertices);
    
    % create a empty struct
    metrics = struct;

    % check if edge metrics need to be calculated
    if  specifyMetrics == 0 || (mAverageEdgeLength || mMinEdgeLength || mMaxEdgeLength || mEdgeLengths)

        % list the vertices of the unique edges
        uniqueEdges = [ faces(:, 1), faces(:, 2); ...
                        faces(:, 1), faces(:, 3); ...
                        faces(:, 2), faces(:, 3)];
        uniqueEdges = unique(sort(uniqueEdges, 2), 'rows');
        edgeDiff = vertices(uniqueEdges(:, 1), :) - vertices(uniqueEdges(:, 2), :);
        edgeLengths = sqrt(sum(edgeDiff .^ 2, 2));

        % edge
        if specifyMetrics == 0 || mEdgeLengths
            metrics.edges = edgeLengths;
        end
        if specifyMetrics == 0 || mAverageEdgeLength
            metrics.averageEdgeLength = mean(edgeLengths);
        end
        if specifyMetrics == 0 || mMinEdgeLength
            metrics.minEdgeLength = min(edgeLengths);
        end
        if specifyMetrics == 0 || mMaxEdgeLength
            metrics.maxEdgeLength = max(edgeLengths);
        end

    end

    % check if area metrics need to be calculated
    if  specifyMetrics == 0 || (mAreas || mAverageArea || mMinArea || mMaxArea || mTotalArea)

        % calculate the areas
        % (this is different than the edges calculation, because here
        %  we need all edges, not just the unique ones)
        v1 = vertices(faces(:, 1), :);
        v2 = vertices(faces(:, 2), :);
        v3 = vertices(faces(:, 3), :);
        diff12 = v1 - v2;
        diff13 = v1 - v3;
        %diff23 = v2 - v3;
        areas = cross(diff12, diff13, 2);
        areas = 1/2 * sqrt(sum(areas .^ 2, 2));

        % area
        if specifyMetrics == 0 || mAreas
            metrics.areas = areas;
        end
        if specifyMetrics == 0 || mAverageArea
            metrics.averageArea = mean(areas);
        end
        if specifyMetrics == 0 || mMinArea
            metrics.minArea = min(areas);
        end
        if specifyMetrics == 0 || mMaxArea
            metrics.maxArea = max(areas);
        end
        if specifyMetrics == 0 || mTotalArea
            metrics.totalArea = sum(areas);
        end

    end
    
end
