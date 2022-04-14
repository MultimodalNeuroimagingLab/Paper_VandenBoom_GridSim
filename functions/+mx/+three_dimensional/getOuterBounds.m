%
%   Retrieve outer bounds from a 3D object
%   [metrics] = getOuterBounds(input, ...)
%   [metrics] = getOuterBounds(vertices, faces, ...)
%
%   	input         = the input 3D to extract the outer bounds from
%
%
%   Returns: 
%       edges         = vector with the edges that define the outer bounds
%       vertexIndices = vector with the indices of the vertices that define
%                       the outer bounds
%
%
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [edges, vertexIndices] = getOuterBounds(input)

    % 
    faces = [];
    
    % loop through the input arguments
    if isstruct(input) || isobject(input)
        % if struct or object

        if isa(input, 'gifti')
            % gifti object

            faces = input.faces;

        elseif isstruct(input)
            % struct
            
            if isfield(input, 'faces')
                faces = input.faces;
            elseif isfield(input, 'face')
                faces = input.face;
            elseif isfield(input, 'tri')
                faces = input.tri;
            end

        end
    end
    
    % check if there are vertices and faces
    if isempty(faces)
        fprintf(2, 'Error: no faces found. Make sure the input is an 3D object\n');
        return;
    end
    
    % TODO: below can probably go faster, replace loop
    
    % determine the vertices on the outside of the patch
    % these are the vertices which have no shared edges
    allEdges = [faces(:, 1), faces(:, 2); ...
                faces(:, 1), faces(:, 3); ...
                faces(:, 2), faces(:, 3)];
    allEdges = unique(sort(allEdges, 2), 'rows');
    edges = [];    
    for iEdge = 1:length(allEdges)
        numTrWithEdge = nnz(sum(ismember(faces, allEdges(iEdge, :)), 2) == 2);
        if numTrWithEdge < 2
            edges(end + 1, :) = allEdges(iEdge,:);
        end
    end
    vertexIndices = unique(edges(:));
    
end
