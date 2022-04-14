%
%   Project points onto a 3D object
%   [points, pointsProjected, pointsTriangles] = projectPointsToHull(hull, points, normdist, intersval)
%
%   	hull        	= the hull onto which the vertices are projected
%   	points      	= the 3D points to project onto the hull (n-by-3 matrix)
%   	normdist    	= the radius at which, from the electrode, vertices are	included to generate an average
%						  normal to determine the intersection
%   	intersval   	= the search distance for intersection
%
%
%   Returns: 
%       points            = the new 3D points projected onto the hull (n-by-3 matrix).
%                           Unprojectable points are left in their original input position
%       pointsProjected   = a flag for each point whether it was projected
%                           succesfully (1 if successfull, 0 if not)
%       pointsTriangles   = the index of the hull triangle onto which the point is projected
%                           Unprojectable points will have NaN
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
function [points, pointsProjected, pointsTriangles] = projectPointsToHull(hull, points, normdist, intersval)

    %
    % Prepare and check the input
    %
    
    % input variables 
    hullVertices = [];
    hullFaces = [];
    
    % check the type of input 3d object
    if isstruct(hull) || isobject(hull)
        % if struct or object

        if isa(hull, 'gifti')
            % gifti object

            hullVertices = hull.vertices;
            hullFaces = hull.faces;

        elseif isstruct(hull)
            % struct

            if isfield(hull, 'vertices')
                hullVertices = hull.vertices;
            elseif isfield(hull, 'vert')
                hullVertices = hull.vert;
            end

            if isfield(hull, 'faces')
                hullFaces = hull.faces;
            elseif isfield(hull, 'face')
                hullFaces = hull.face;
            elseif isfield(hull, 'tri')
                hullFaces = hull.tri;
            end

        end
        
    end
    
    % check if there are vertices and faces
    if isempty(hullVertices) || isempty(hullFaces)
        fprintf(2, 'Error: vertices and faces found. Make sure either an 3D object, or vertices and faces are given as first arguments\n');
        return;
    end
    
    % preallocate the output arrays
    pointsProjected = false(size(points, 1), 1);
    pointsTriangles = NaN(size(points, 1), 1);
    
    % convert the vertices to double for increased precision
    hullVertices = double(hullVertices);
    
    % calculate the normals for each triangle
    tr = triangulation(double(hullFaces), hullVertices);
    triangleNormals = faceNormal(tr);
    [~, MSGID] = lastwarn();
    warning('off', MSGID);
    vertexNormals = vertexNormal(tr);
    [~, MSGID] = lastwarn();
    warning('off', MSGID);
    
    %compute the distances of all vertices to all electrodes
    tic
    veDist = (points(:, 1)' - hullVertices(:, 1)) .^ 2 + ...
             (points(:, 2)' - hullVertices(:, 2)) .^ 2 + ...
             (points(:, 3)' - hullVertices(:, 3)) .^ 2;
    %toc

    % sort the distances per electrode
    [veDistSorted, veDistSortedIndices] = sort(veDist);

    % check if there are electrodes with no vertices to calculate the normal on
    tic
    veSelectU = veDistSorted < normdist^2;
    [~, veSelectU] = min(veSelectU);
    veSelectU = veSelectU - 1;
    noNormElec = find(veSelectU == 0);
    if ~isempty(noNormElec)
        disp(['Average normal(s) could not be calculated for electrode(s): ', num2str(noNormElec), '. No vertices within the range of ', num2str(normdist), ' from these electrodes.']);
        return;
    end
    %toc

    % compute the normal based on the vertices within a fixed normal distance
    tic
    normals2av = nan(size(veDist, 2), 3);
    veIdx = veDist < normdist^2;                % faster than having find in the loop
    for i = 1:size(veDist, 2)
        normals2av(i, :) = mean(vertexNormals(veIdx(:, i), :));
    end
    %toc

    % normal vector computation (computation of line points p1 and p2):
    points_ends = points + normals2av;

    % debug, display some of the normal vectors
    %viewgii(gRoiHull, [p1(1:27, :), p2(1:27, :)], 'wireframe');


    %% 

    % use a fixed radius to include vertices and faces with which normals can intersect
    intersdist2 = intersval^2; %fixed distance


    tic

    % determine for each electrode which vertices lie within the range for
    % intersection (more correctly, the intersection with these vertices their
    % triangles)
    veSelectU = veDistSorted < intersdist2;
    [~, veSelectU] = min(veSelectU);    % finds the 
    veSelectU = veSelectU - 1;

    % check if there are electrodes which have no vertices within range for intersection
    noInterElec = find(veSelectU == 0);
    if ~isempty(noInterElec)
        disp(['Could not find vertices for electrode(s): ', num2str(noInterElec), ' based on the range of ', num2str(intersval), '.']);
        return;
    end
    %toc

    tic
    % loop through the electrodes
    for i = 1:size(veDist, 2)

        % retrieve the vertices sorted by distance (closest first)
        sortVert = veDistSortedIndices(1:veSelectU(i), i);

        % method 1, unsorted (about 10s)
        %[withinFaces] = ismember(tri, sortVert);
        %withinFaces = any(withinFaces, 2);
        %withinFaceIndices = find(withinFaces);
        

        % method 2, sorted using ismember and sorts (about 60s, this method is faster when there are more faces)
        [~, withinFaces] = ismember(hullFaces, sortVert);
        [a, b] = sort(withinFaces, 'asc');
        a2 = a(:);
        b2 = b(:);
        b2 = b2(a2~=0);
        a2 = a2(a2~=0);
        [~, d] = sort(a2, 'asc');
        sorttri = unique(b2(d), 'stable');


        %{
        % method 3, sorted using loop (about 61s)
        withinFaces = tri; %tri will be altered
        sorttri=[];
        Lcv = length(sortVert);
        for k = 1 : Lcv,
            cv = sortVert(k);
            [rows, columns] = find(withinFaces == cv);
            sorttri = [sorttri; rows];
            withinFaces(rows, :) = 0; %assign zero so that they're not added again
        end
        %}

        Lct = length(sorttri);
        inside = 0;
        % loop through the ordered triangles
        for k = 1 : Lct
            tr = hullFaces(sorttri(k), :);
            t = [hullVertices(tr(1), :); hullVertices(tr(2), :); hullVertices(tr(3), :)]';
            triangleNormal = triangleNormals(sorttri(k), :);

            %
            %  Find the intersection of the plane and the line.
            %
            [ ival, pint ] = mx.three_dimensional.ext.plane_normal_line_exp_int_3d ( t(:, 1)', triangleNormal, points(i, :), points_ends(i, :) );

            if ( ival == 1 )
                % the line and plane intersect at a single point

                %
                %  Now, check that all three triangles made by two vertices and
                %  the intersection point have the same "clock sense" as the
                %  triangle's normal vector.
                %
                v1 = t(:, 2)  - t(:, 1);
                v2 = pint - t(:, 1)';

                normal2(1) = v1(2) * v2(3) - v1(3) * v2(2);
                normal2(2) = v1(3) * v2(1) - v1(1) * v2(3);
                normal2(3) = v1(1) * v2(2) - v1(2) * v2(1);

                if ( triangleNormal * normal2' < 0.0 )
                    continue;
                end

                v1 = t(:, 3)  - t(:, 2);
                v2 = pint - t(:, 2)';

                normal2(1) = v1(2) * v2(3) - v1(3) * v2(2);
                normal2(2) = v1(3) * v2(1) - v1(1) * v2(3);
                normal2(3) = v1(1) * v2(2) - v1(2) * v2(1);

                if ( triangleNormal * normal2' < 0.0 )
                    continue;
                end

                v1 = t(:, 1)  - t(:, 3);
                v2 = pint - t(:, 3)';

                normal2(1) = v1(2) * v2(3) - v1(3) * v2(2);
                normal2(2) = v1(3) * v2(1) - v1(1) * v2(3);
                normal2(3) = v1(1) * v2(2) - v1(2) * v2(1);

                if ( triangleNormal * normal2' < 0.0 )
                    continue;
                end

                % intersection was found
                inside = 1;
                break;

            end

        end     % end triangle loop

        % check if inside a triangle
        if inside
            
            % set the intersecting position
            points(i, :) = pint(:)';
            
            % flag as projected
            pointsProjected(i) = 1;
            
            % set the hu;; triangle onto which the point was projected
            pointsTriangles(i) = sorttri(k);
            
        end

    end
    %toc    
    
end
