%   
%   Interprets the points of a rectangle (defined by 8 corner points) to retrieve line segments and triangles
%
%   [lineSegIndices, lineSegCoords] = rectPointsTo3D(rectPoints)
%   
%   rectPoints           = (8-by-3 matrix) the eight corner points that define the rectangle in space
%   
%   Returns: 
%       lineSegIndices   = 12-by-2-matrix, the 12 lines that make up
%                          the rectangle by referring to the indices of the
%                          input (corner) points
%       lineSegCoords    = 12-by-6-matrix, the 12 lines that make up
%                          the rectangle by line coordinates, the first 3 colums(second
%                          dimension in the matrix) represent the x, y and z coordinates
%                          of the start of the line, whereas the last three columns 
%                          represent the x, y and z coordinates of the end of the line
%       triangleIndices  = 12-by-3-matrix, the 12 triangles that make up
%                          the rectangle by referring to the indices of the
%                          input (corner) points
% 
%   Example:
%   	rectPoints = [ [0, 0, 0]; [100, 0, 0]; [100, 200, 50]; [100, 200, 0]; ...
%                      [0, 200, 0]; [0, 0, 50]; [100, 0, 50]; [0, 200, 50]  ];
%       rectPoints = ((rectPoints * rotx(30)) * roty(40)) * rotz(50);
%       [lineSegIndices, lineSegCoords] = rectToLineSegments(fovCorners);
% 
%   Copyright (C) 2019 Max van den Boom
%
function [lineSegIndices, lineSegCoords, triangleIndices] = rectPointsTo3D(rectPoints)
    
    %{
    % debug, show initial input
    figure;
    scatter3(rectPoints(:, 1), rectPoints(:, 2), rectPoints(:, 3));
    figure;
    scatter3(trPts(:, 1), trPts(:, 2), trPts(:, 3));
    %}
    
    % translate input points to an origin of 0
    trPts = rectPoints - rectPoints(1, :);

    % find the points that are adjacent to the point at 0
    % 
    % all diagonals in 3d should be further away than adjacent, so we can
    % assume the two points which are the closest are adjacent
    distPts =   trPts(2:end, 1) .^ 2 + ...
                trPts(2:end, 2) .^ 2 + ...
                trPts(2:end, 3) .^ 2;
    [~, closestPts] = sort(distPts, 'asc');
    closestPts = closestPts(1:2) + 1;
    
    % TODO, if not exactly rectangular, then try to get as close as possible
    
    % generate the rotation matrix (based on the adjacent points) to remove the angulation
    e0 = double(trPts(closestPts(1), :));
    e1 = double(trPts(closestPts(2), :));
    e0 = e0 / vecnorm(e0);
    e1 = e1 / vecnorm(e1);
    normal = cross(e0, e1, 2);
    trRotMat = [cross(e0, normal, 2); e0; normal]';
    
    %{
    % debug, show the normals and adjacents
    figure;
    scatter3(trPts(:, 1), trPts(:, 2), trPts(:, 3));
    hold on;
    %daspect([1 1 1])
    scatter3(trPts(closestPts, 1), trPts(closestPts, 2), trPts(closestPts, 3), 'r');
    scatter3(trPts(1, 1), trPts(1, 2), trPts(1, 3), 'g');
    scatter3(e0(1) * 20, e0(2) * 20, e0(3) * 20, 'k');
    scatter3(e1(1) * 100, e1(2) * 100, e1(3) * 100, 'k');
    scatter3(normal(1) * 170, normal(2) * 170, normal(3) * 170, 'k');
    hold off;
    %}
    
    % rotate the rectangle to remove the angulation
    trPts = trPts * trRotMat;

    %{
    % debug, display rectangle
    figure;
    scatter3(trPts(:, 1), trPts(:, 2), trPts(:, 3));
    %}
    
    % determine the spatial middle of the rectangle
    rectMin = min(trPts);
    rectMax = max(trPts);
    middle = rectMin + abs(rectMin - rectMax) / 2;
    
    %{
    % debug, display rectangle with middle
    trPts = [trPts; middle];
    figure;
    scatter3(trPts(:, 1), trPts(:, 2), trPts(:, 3));
    %}
    
    % determine the corner coordinates
    bigger = trPts > middle;
    c0 = find( bigger(:, 1) &  bigger(:, 2) &  bigger(:, 3), 1);   % +x, +y, +z
    c1 = find(~bigger(:, 1) &  bigger(:, 2) &  bigger(:, 3), 1);   % -x; +y, +z
    c2 = find(~bigger(:, 1) & ~bigger(:, 2) &  bigger(:, 3), 1);   % -x; -y, +z
    c3 = find( bigger(:, 1) & ~bigger(:, 2) &  bigger(:, 3), 1);   % +x; -y, +z
    
    c4 = find( bigger(:, 1) &  bigger(:, 2) & ~bigger(:, 3), 1);   % +x, +y, -z
    c5 = find(~bigger(:, 1) &  bigger(:, 2) & ~bigger(:, 3), 1);   % -x; +y, -z
    c6 = find(~bigger(:, 1) & ~bigger(:, 2) & ~bigger(:, 3), 1);   % -x; -y, -z
    c7 = find( bigger(:, 1) & ~bigger(:, 2) & ~bigger(:, 3), 1);   % +x; -y, -z
    
    
    % define the line segments by indices
    lineSegIndices = [  [c0 c1]; [c1 c2]; [c2 c3]; [c3 c0]; ...
                        [c4 c5]; [c5 c6]; [c6 c7]; [c7 c4]; ...
                        [c0 c4]; [c1 c5]; [c2 c6]; [c3 c7]; ...
                     ];

    % define the triangles by indices
    triangleIndices = [  [c0 c1 c2];    [c0 c2 c3]; ...
                        [c1 c5 c6];    [c6 c2 c1]; ...
                        [c6 c5 c4];    [c7 c6 c4]; ...
                        [c0 c3 c4];    [c7 c4 c3]; ...
                        [c4 c1 c0];    [c1 c4 c5]; ...
                        [c6 c3 c2];    [c3 c6 c7]; ...
                     ];
                 
    % define the line segments by coordinates
    lineSegCoords = [   rectPoints(c0, :), rectPoints(c1, :); ...
                        rectPoints(c1, :), rectPoints(c2, :); ...
                        rectPoints(c2, :), rectPoints(c3, :); ...
                        rectPoints(c3, :), rectPoints(c0, :); ...
                        rectPoints(c4, :), rectPoints(c5, :); ...
                        rectPoints(c5, :), rectPoints(c6, :); ...
                        rectPoints(c6, :), rectPoints(c7, :); ...
                        rectPoints(c7, :), rectPoints(c4, :); ...
                        rectPoints(c0, :), rectPoints(c4, :); ...
                        rectPoints(c1, :), rectPoints(c5, :); ...
                        rectPoints(c2, :), rectPoints(c6, :); ...
                        rectPoints(c3, :), rectPoints(c7, :); ...
                    ];
            
    %{
    % debug, lines
    figure;
    scatter3(rectPoints(:, 1), rectPoints(:, 2), rectPoints(:, 3));
    hold on;
    for iLine = 1:size(lineSegCoords, 1)
        plot3(  [lineSegCoords(iLine, 1), lineSegCoords(iLine, 4)], ...
                [lineSegCoords(iLine, 2), lineSegCoords(iLine, 5)], ...
                [lineSegCoords(iLine, 3), lineSegCoords(iLine, 6)], ...
                'Color', [1 0 0], ...
                'LineWidth', 1);
    end
    
    hold off;
    %}
    
end
