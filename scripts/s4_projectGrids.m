%   
%   Step 4 - Project virtual grids on a hull using the 3D sample-points as the center, with their random rotations and the given
%            grid configurations. Output files will be written in the provided path.
%
%   outputFiles = s4_projectGrids(SS, hullPath, projectionsOutputDir)
% 
%       SS                          = Structure that holds the sample-points to project the different grid-configurations on
%       hullPath                    = path to the hull gifti file, this is where the grids will be projected on
%       projectionsOutputDir        = path a folder to output the files to
%
%
%   Returns: 
%       outputFiles                 = A cell array containing the paths to all grid-configurations output files
%
%
%   Copyright 2020, Max van den Boom (Multimodal Neuroimaging Lab, Mayo Clinic, Rochester MN)

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function outputFiles = s4_projectGrids(SS, hullPath, projectionsOutputDir)


    % define the grid configuration to loop through
    %   <numHorz>, <numVert>, <radius>, <spacing>
    %gridConfigurations = [3, 3, 3, 6];  
    gridConfigurations = [  2, 2, 1, 4;     2, 2, 2, 4;     2, 2, 3, 4; ...
                            2, 2, 1, 6;     2, 2, 2, 6;     2, 2, 3, 6; ...
                            2, 2, 1, 8;     2, 2, 2, 8;     2, 2, 3, 8; ...
                            2, 2, 1, 10;    2, 2, 2, 10;    2, 2, 3, 10; ...

                            3, 3, 1, 4;     3, 3, 2, 4;     3, 3, 3, 4; ...
                            3, 3, 1, 6;     3, 3, 2, 6;     3, 3, 3, 6; ...
                            3, 3, 1, 8;     3, 3, 2, 8;     3, 3, 3, 8; ...
                            3, 3, 1, 10;    3, 3, 2, 10;    3, 3, 3, 10; ...

                            4, 4, 1, 4;     4, 4, 2, 4;     4, 4, 3, 4; ...
                            4, 4, 1, 6;     4, 4, 2, 6;     4, 4, 3, 6; ...
                            4, 4, 1, 8;     4, 4, 2, 8;     4, 4, 3, 8; ...
                            4, 4, 1, 10;    4, 4, 2, 10;    4, 4, 3, 10; ...

                            5, 5, 1, 4;     5, 5, 2, 4;     5, 5, 3, 4; ...
                            5, 5, 1, 6;     5, 5, 2, 6;     5, 5, 3, 6; ...
                            5, 5, 1, 8;     5, 5, 2, 8;     5, 5, 3, 8; ...
                            5, 5, 1, 10;    5, 5, 2, 10;    5, 5, 3, 10; ...

                            6, 6, 1, 4;     6, 6, 2, 4;     6, 6, 3, 4; ...
                            6, 6, 1, 6;     6, 6, 2, 6;     6, 6, 3, 6; ...
                            6, 6, 1, 8;     6, 6, 2, 8;     6, 6, 3, 8; ...
                            6, 6, 1, 10;    6, 6, 2, 10;    6, 6, 3, 10; ...

                            7, 7, 1, 4;     7, 7, 2, 4;     7, 7, 3, 4; ...
                            7, 7, 1, 6;     7, 7, 2, 6;     7, 7, 3, 6; ...
                            7, 7, 1, 8;     7, 7, 2, 8;     7, 7, 3, 8; ...
                            7, 7, 1, 10;    7, 7, 2, 10;    7, 7, 3, 10; ...

                            8, 8, 1, 4;     8, 8, 2, 4;     8, 8, 3, 4; ...
                            8, 8, 1, 6;     8, 8, 2, 6;     8, 8, 3, 6; ...
                            8, 8, 1, 8;     8, 8, 2, 8;     8, 8, 3, 8; ...
                            8, 8, 1, 10;    8, 8, 2, 10;    8, 8, 3, 10 ];

                    
    %%
    %  Retrieve the root path and make sure dependencies can be found
    %
    addpath('../');
    gridSimRoot = gridSimRootPath;
    addpath([gridSimRoot, filesep, 'functions']);
    addpath(genpath([gridSimRoot, filesep, 'external']));

    

    %%
    %  Load the hull to project on
    %
	gHull = gifti(hullPath);
    

    %%%
    %%  Warn if the samplepoints in the set have been transformed to match the "MRI/NIFTI" space
    %   it is possible that the pial and hull are not, this will severly mess up the projection
    %%%

    % check if the sampleset positions have been corrected (transformed)
    if isfield(SS, 'transformedToNiftiSpace') == 1 && ~isempty(SS.transformedToNiftiSpace)
        % corrected

        % message
        fprintf(2, 'Warning!!: the samples in the set have been transformed to match the "MRI/NIFTI" space, this step/function\nmight expect the samples to be in the FS pial space\n');

    end


    % create an copy of the original input sample-set
    origSS = SS;

    
    %%
    % project the grid for each grid configuration on the sample-set points
    %

    outputFiles = {};
    
    % loop through the different grid configurations
    for iGridConf = 1:size(gridConfigurations, 1)

        % create a clear workspace, while keeping essential input variables
        clearvars -except iGridConf gridConfigurations origSS bids_simTaskSelectionDir gHull projectionsOutputDir outputFiles;

        % reset the sample-set to the original input
        SS = origSS;

        % retrieve the grid configuration
        numHorz = gridConfigurations(iGridConf, 1);
        numVert = gridConfigurations(iGridConf, 2);
        radius = gridConfigurations(iGridConf, 3);
        spacing = gridConfigurations(iGridConf, 4);

        % message
        disp([  'numHorz: ', num2str(numHorz), ' - ' ...
                'numVert: ', num2str(numVert), ' - ' ...
                'radius: ', num2str(radius), ' - ' ...
                'spacing: ', num2str(spacing)]);

        %%%
        %% define the grid configuration in the samples file and calculate some variables
        %%%
        clear gridConfig;
        SS.gridConfig.elecNumHorz = numHorz;            % number of horizontal electrodes
        SS.gridConfig.elecNumVert = numVert;            % number of vertical electrodes
        SS.gridConfig.elecHorzSpacing = spacing;        % spacing in mm between electrode centers
        SS.gridConfig.elecVertSpacing = spacing;        % spacing in mm between electrode centers
        SS.gridConfig.elecSizeRadius = radius;          % radius of the electrods in mm
        SS.gridConfig.elecTotal = SS.gridConfig.elecNumHorz * SS.gridConfig.elecNumVert;

        SS.gridConfig.sampleHullEdgeExclusionMargin = 0.5; % the safety margin (to correct for brain curvature) by which a sample is removed when to close to the edge of the hull
                                                        % lower is safer, but could later turn out to be computationally more expensive for at electrode projection
                                                        % (1 = when the shortest distance between a samples center and it's border is higher
                                                        % than the distance of the center to the closest hull's edge; 0.5 then half of the
                                                        % shortest distance between a sample center and it's border is higher than ...)

        % calculate the horizontal and vertical size of the grid
        gridHSizeWithoutElec = (SS.gridConfig.elecNumHorz - 1) * SS.gridConfig.elecHorzSpacing;
        gridVSizeWithoutElec = (SS.gridConfig.elecNumVert - 1) * SS.gridConfig.elecVertSpacing;
        gridHSizeWithElec = gridHSizeWithoutElec + (SS.gridConfig.elecSizeRadius * 2);
        gridVSizeWithElec = gridVSizeWithoutElec + (SS.gridConfig.elecSizeRadius * 2);



        %%%
        %% remove points that are roughly too close to the border of the hull (before projection)
        %
        %  grid-samples where the shortest distance between the center and it's borders is 
        %  larger than the sample's center to any of the RIO edges can be removed beforehand
        %  as they will be sure to extend over the borders of the hull. We pre-remove the
        %  samples here because projecting computationally demanding, so it is
        %  better to remove the samples here
        %
        %  note: because the distance to the border does not take into account the
        %        curveture of the brain, we use a margin. Which means only samples
        %        half the distance (between their center and the closest border)
        %        will be removed
        %%%


        % message
        disp ('----');
        disp ('Identifying samples too close to the edge...');

        % determine the exclusion distance as the smallest distance between the center of the grid and it's border
        if SS.gridConfig.elecNumHorz * SS.gridConfig.elecHorzSpacing < SS.gridConfig.elecNumVert * SS.gridConfig.elecVertSpacing
            exclusionDist = gridHSizeWithElec / 2;
        else
            exclusionDist = gridVSizeWithElec / 2;
        end

        % shorten the exclusion distance by the margin factor (to correct for the brain curveture)
        exclusionDist = exclusionDist * SS.gridConfig.sampleHullEdgeExclusionMargin;

        % determine the outer bounds 
        [outerEdges, outerVertexIndices] = mx.three_dimensional.getOuterBounds(gHull);
        outerVertices = gHull.vertices(outerVertexIndices, :);

        % calculate the differences and sum of the squared x, y and z between each edge's end and start-point
        edgesDiff = gHull.vertices(outerEdges(:, 2), :) - gHull.vertices(outerEdges(:, 1), :);
        edgesSSQ = sum(edgesDiff .^ 2, 2);

        % variable to store the sample in that need to be removed (because to close)
        remSamples = [];

        % reset timer
        tic

        % loop through the samples
        % TODO, might be faster without loop
        for iSample = 1:size(SS.samplePositions, 1)

            % calculate the differences and sum of the squared x, y and z between the sample and the start point of each edge
            d12 = SS.samplePositions(iSample, :) - gHull.vertices(outerEdges(:, 1), :);
            S1 = sum(edgesDiff .* d12, 2);

            % calculate the t and cap it
            t = S1 ./ edgesSSQ;
            t(t < 0) = 0;   t(t > 1) = 1;

            % calculate the distance between the sample and the closest point on each edge
            dist = sqrt(sum((edgesDiff .* t - d12) .^ 2, 2));

            % determine the index of the edge which is the closest
            [dist, iEdge] = min(dist);

            % determine if the sample point is too close to the closest edge
            if (dist < exclusionDist)

                % add to the list of samples to remove
                remSamples(end + 1) = iSample;

            end

        end

        % message
        disp (['Sample edge removal: ', num2str(toc), ' ms']);
        disp (['Samples to be removed (too close to hull edge): ', num2str(length(remSamples))]);

        % remove the samples
        SS = removeFromSampleSetByIndices(SS, remSamples);

        % message remaining
        disp (['Samples remaining: ', num2str(size(SS.samplePositions, 1))]);

        

        %%%
        %% project the grids on the brain
        %%%

        % set starttime and display number of samples
        projGridTic = tic;
        disp(['samples to be projected: ', num2str(size(SS.samplePositions, 1))]);

        %
        % generate the grids and apply the random rotation
        %

        % define the grid as electrode positions in (local) 2d space
        [X, Y] = meshgrid(  0:SS.gridConfig.elecHorzSpacing:(SS.gridConfig.elecNumHorz - 1) * SS.gridConfig.elecHorzSpacing, ...
                            0:SS.gridConfig.elecVertSpacing:(SS.gridConfig.elecNumVert - 1) * SS.gridConfig.elecVertSpacing);
        X = X - gridHSizeWithoutElec / 2;
        Y = Y - gridVSizeWithoutElec / 2;
        sampleGridLocal = [X(:), Y(:)];

        % reset timer
        tic

        % create an (concatenated) array of sample grids based on the local grid
        SS.sampleGrids = repmat(sampleGridLocal, size(SS.samplePositions, 1), 1);

        % generate the rotation cos and sin for each sample
        rCos = cosd(SS.sampleGridRotations);
        rSin = sind(SS.sampleGridRotations);

        % (replicate to match the number of electrodes in the concatenated array of sample grids)
        rCos = repmat(rCos', SS.gridConfig.elecTotal, 1);
        rCos = rCos(:);
        rSin = repmat(rSin', SS.gridConfig.elecTotal, 1);
        rSin = rSin(:);

        % rotate each grid in the concatenated array of sample grids and add a third coordinate
        SS.sampleGrids = [ rCos .* SS.sampleGrids(:, 1) + rSin .* SS.sampleGrids(:, 2), ...
                          -rSin .* SS.sampleGrids(:, 1) + rCos .* SS.sampleGrids(:, 2), ...
                           zeros(size(SS.sampleGrids, 1), 1)];

        % message elapsed and reset timer
        disp (['Local rotation: ', num2str(toc), ' ms']);
        tic


        %
        % determine the rotation matrices for each triangle
        %

        % get the first edges of each triangle as unit vectors and calculate the normals (as unit vectors)
        e0 = double(gHull.vertices(gHull.faces(:, 2), :) - gHull.vertices(gHull.faces(:, 1), :));
        e1 = double(gHull.vertices(gHull.faces(:, 3), :) - gHull.vertices(gHull.faces(:, 1), :));
        e0 = e0 ./ vecnorm(e0')';
        normals = cross(e0, e1, 2);
        normals = normals ./ vecnorm(normals')';

        % the rotation matrix will be the triangle's coordinate system, using the
        % the normal, first edge and their orthogonal as z, y and x axis
        %trRotMat = [cross(e0, normals, 2), ...     % x, the orthagonal of the first edge and the normal (since both orthogonal, will result in unit vectors)
        %            e0, ...                        % y, the first edge
        %            normals];                      % z, the normals
        trRotMat = cat(3, cross(e0, normals, 2), e0, normals);
        trRotMat = permute(trRotMat, [1 3 2]);


        %
        % rotations of the grids to their triangles coordinate systems
        % TODO: perhaps could be faster, but couldn't get element-wise multiplication
        %       of matrix x matrix, with matrix being unequal sizes to work as quick
        %
        for i = 1:size(trRotMat, 1)

            % determine which grid point (indices) are affected by the current triangle rotation matrix
            iRotPoints = find(SS.sampleTriangles == i);
            iRotPoints = ((iRotPoints - 1) * SS.gridConfig.elecTotal) + 1;
            [~, gridI] = meshgrid(1:length(iRotPoints), 0:SS.gridConfig.elecTotal - 1);
            iRotPoints = gridI + iRotPoints';
            iRotPoints = iRotPoints(:);

            % rotate the grid points by the corrent triangle rotation matrix
            SS.sampleGrids(iRotPoints, :) = SS.sampleGrids(iRotPoints, :) * squeeze(trRotMat(i, :, :));

        end

        % message elapsed and reset timer
        disp (['Rotations to triangle: ', num2str(toc), ' ms']);
        tic


        %
        % translation of the grids to the sample positions
        %

        % calculate the x, y and z offsets of each sample
        xRep = repmat(SS.samplePositions(:,1)', SS.gridConfig.elecTotal, 1);
        yRep = repmat(SS.samplePositions(:,2)', SS.gridConfig.elecTotal, 1);
        zRep = repmat(SS.samplePositions(:,3)', SS.gridConfig.elecTotal, 1);
        trans = double([xRep(:), yRep(:), zRep(:)]);

        % transpose each grid in the concatenated array of sample grids
        SS.sampleGrids = SS.sampleGrids + trans;

        % message elapsed
        disp (['Transposing: ', num2str(toc), ' ms']);


        
        %
        % project the grid points (electrodes) onto the hull
        %

        %% 
        normdist = 25;
        intersval = 20;

        % Project all grid points (electrodes) onto the hull
        splitConfig = [];
        splitConfig.numPoints = 7000;
        [SS.sampleGrids, SS.sampleGridsProjected, SS.sampleGridsTriangles] = mx.three_dimensional.projectPointsToHull_split(gHull, SS.sampleGrids, normdist, intersval, splitConfig);

        % check whether some electrodes of the grids could not be projected
        if nnz(~SS.sampleGridsProjected) > 0

            % message
            fprintf(2, 'One or more grids electrodes could not be projected\n');

            % TODO, try again with larger distances?

        end

        % check if there are electrodes which could not be projected
        if nnz(~SS.sampleGridsProjected) > 0

            % message
            disp ('----');
            disp ('Identifying electrodes outside of hull...');

            % determine which grids were projected successfully (the grids of which all electrodes were projected correctly)
            gridsProjected = reshape(SS.sampleGridsProjected, SS.gridConfig.elecTotal, [])';
            gridsProjected = all(gridsProjected, 2);

            % determine the indices of the grids (samples) which were not projected successfully
            remSamples = find(gridsProjected == 0);

            % message
            disp (['Sample grids to be removed (one or more electrodes outside of hull): ', num2str(length(remSamples))]);

            % remove the samples (and their grids)
            SS = removeFromSampleSetByIndices(SS, remSamples);

            % message remaining
            disp (['Sample grids remaining: ', num2str(size(SS.samplePositions, 1))]);

        end

        % display elapsed time for projection
        disp(['Time elapsed on projection of samples: ', num2str(toc(projGridTic))]);
        
        % save the grid-configuration sampleset
        saveFolder =    fullfile(projectionsOutputDir, ['n', num2str(SS.gridConfig.elecTotal)]);
        if ~exist(saveFolder, 'dir'),   mkdir(saveFolder);  end
        outFilename = ['projections_n', num2str(SS.gridConfig.elecTotal), '_sp', num2str(SS.gridConfig.elecHorzSpacing), '_ra', num2str(SS.gridConfig.elecSizeRadius), '.mat'];
        save(fullfile(saveFolder, outFilename), 'SS');
        
        % store the path to the grid-configuration sampleset
        outputFiles{end + 1} = fullfile(saveFolder, outFilename);
        
    end     % end grid configuration loop

end
