%
%   Generate random sample points on a given 3D object
%   [samplePoints, sampleTriangles] = genRandomSamplePoints3DObject(obj3D, numSamples)
%
%       obj3D           = the 3D object onto which the vertices are projected
%       numSamples      = the number of 3D sample pont to be generated
%
%
%   Returns: 
%       samplePoints    = The random 3D sample points
%       sampleTriangles = The triangles to which each 3D sample point belongs
%
%
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [samplePoints, sampleTriangles] = genRandomSamplePoints3DObject(obj3D, numSamples)

    %
    % Prepare and check the input
    %
    
    % input variables 
    obj3DVertices = [];
    obj3DFaces = [];
  
    % check the type of input 3d object
    if isstruct(obj3D) || isobject(obj3D)
        % if struct or object

        if isa(obj3D, 'gifti')
            % gifti object

            obj3DVertices = obj3D.vertices;
            obj3DFaces = obj3D.faces;

        elseif isstruct(obj3D)
            % struct

            if isfield(obj3D, 'vertices')
                obj3DVertices = obj3D.vertices;
            elseif isfield(obj3D, 'vert')
                obj3DVertices = obj3D.vert;
            end

            if isfield(obj3D, 'faces')
                obj3DFaces = obj3D.faces;
            elseif isfield(obj3D, 'face')
                obj3DFaces = obj3D.face;
            elseif isfield(obj3D, 'tri')
                obj3DFaces = obj3D.tri;
            end

        end
        
    end
    
    % check if there are vertices and faces
    if isempty(obj3DVertices) || isempty(obj3DFaces)
        fprintf(2, 'Error: no vertices and faces found. Make sure a 3D object is given as the first argument\n');
        return;
    end
    
    % convert the vertices to double for increased precision
    obj3DVertices = double(obj3DVertices);
    
    
    
    %
    % generate a sample-list with randomly selected triangles (weighted by the triangle area)
    %

    % retrieve the metrics for this 3d object (includes the areas of each triangle)
    metrics = mx.three_dimensional.get3DMetrics(obj3DVertices, obj3DFaces, 'areas');
    areas = metrics.areas;

    % generate a list with the cumulative areas
    cumAreas = cumsum(areas);

    % randomize the matlab number generator
    rng('shuffle')

    % create a list of (numsamples) random numbers between 0 and the cumulative total
    %
    % note 1: be aware of the cumulative floating point arithmetic error, where
    % the cumulative 'total' will not be the same as the total of the areas
    % because of floating point inprecision. To account for this, we use the
    % cumulative total
    %
    % note 2: be aware that rand generates an open interval and that the randomization
    % therefore will never include 0 nor the total area in it's randomization.
    % However, given the precision of rand, the bias on not-selection of the 
    % first or last triangle is neglectable)
    sampleTriangles = rand(1, numSamples) * cumAreas(end);
	
    % determine which triangle based on the random (area) value
    inTriangle = (sampleTriangles >= [0; cumAreas(1:end - 1)] & sampleTriangles < cumAreas);
    [sampleTriangles, ~] = find(inTriangle == 1);
    clear inTriangle;

    % check if indeed an equal amount of random triangles are found as the
    % number of examples required (inprecisions might incidentally cause
    % strange effects)
    if numSamples ~= length(sampleTriangles)
        fprintf(2, 'Warning: random drawn triangles not equal to number of expected samples, debug\n');
    end


    %
    % generate a sample-list with random barycentric coordinates
    %

    % generate a list of random barycentric coordinates
    randCoords = rand(numSamples, 2);
    randCoords(sum(randCoords, 2) > 1, :) = 1 - randCoords(sum(randCoords, 2) > 1, :);
    
    
    %
    % prepare barycentric conversion
    %

    % retrieve the first vertex of each face (we will use this as origin for the triangle)
    p0 = obj3DVertices(obj3DFaces(:, 1), :);

    % shift the origin to the first vertex of every triangle
    % and calculate the (x, y, z) difference between the origin and the hit point
    E1 = obj3DVertices(obj3DFaces(:, 2), :) - p0;
    E2 = obj3DVertices(obj3DFaces(:, 3), :) - p0;

    
    %
    % calculate for each sample the point in 3D space (based on the sample's triangle and barycentric coordinate)
    %

    samplePoints = p0(sampleTriangles, :) + E1(sampleTriangles, :) .* randCoords(:, 1) + E2(sampleTriangles, :) .* randCoords(:, 2);

end
