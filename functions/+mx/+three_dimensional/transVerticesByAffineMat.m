%   
%   Transform a matrix of N-x-3 vertices by a 4x4 affine matrix
%   transformedVertices = transVerticesByAffineMat(vertices, affineMat)
%
%       vertices            = the vertex (1-x-3) or matrix of vertices (N-x-3) to transform
%       affineMat           = the affine matrix (4-x-4) to transform the vertices with
%
%
%   Returns: 
%       The transformated vertices in a N-x-3 matrix
%
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function transformedVertices = transVerticesByAffineMat(vertices, affineMat)
    transformedVertices = double(([vertices ones(size(vertices, 1), 1)])');
    transformedVertices = affineMat * transformedVertices;
    transformedVertices(4,:) = [];
    transformedVertices = transformedVertices';
end