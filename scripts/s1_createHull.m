%   
%   Step 1 - create a hull based on freesurfer parcellation areas
%   s1_createHull(bids_rootPath, bids_sub, hemi)
% 
%      bids_rootPath    = path to the BIDS root-directory where the data is located
%      bids_sub         = the subject to create a hull for
%      hemi             = which hemisphere to simulate on
%
%   Copyright (C) 2019 Max van den Boom  (Multimodal Neuroimaging Lab, Mayo Clinic, Rochester MN)

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function s1_createHull(bids_rootPath, bids_sub, hemi)

    % specify the freesurfer parcallation areas to base the hull on
    outputPrefix    = 'ext';
    fsROIAreas 	    = { 'superiorfrontal', 'caudalmiddlefrontal', 'parsopercularis', ...
                        'precentral', 'postcentral', 'superiorparietal', ...
                        'supramarginal', 'inferiorparietal'  };
    %outputPrefix    = 'sensorymotor';
    %fsROIAreas      = { 'precentral', 'postcentral'};

    % specify hull generation parameters
    numHullRefines                  = 1;        % number of times the hull is refined
    hullNormalSearchDepthFactor     = 5;        % how deep inward from the hull each face will search for the ROI
    hullNormalSearchOutwardFactor   = 8;        % how far the inner hull will search for a triangle of the outer hull


    
    %%
    % retrieve the root path and make sure dependencies can be found
    %
    addpath('../');
    gridSimRoot = gridSimRootPath;
    addpath([gridSimRoot, filesep, 'functions']);
    addpath(genpath([gridSimRoot, filesep, 'external']));
    
    

    %%
    %  build paths and read the files
    %
    
    % build the paths to the data
    bids_fsPath     = fullfile(bids_rootPath, 'derivatives', 'freesurfer', ['sub-' bids_sub]);
    bids_pialPath   = fullfile(bids_fsPath, 'surf', [hemi, '.pial.gii']);
    bids_annotPath  = fullfile(bids_fsPath, 'label', [hemi, '.aparc.Ext.annot']);

    % read the surface files (pial)
    gSurfPial = gifti(bids_pialPath);

    % read the annotations for the pial brain
    [~, annotVertexLabels, annotColortable] = read_annotation(bids_annotPath);



    %%
    %  isolate the areas of interest on the pial brain
    %

    % relabel the vertex annotation (color) labels to the ROI area indices
    roiVertexLabels = fsRelabelToAreas(fsROIAreas, annotColortable, annotVertexLabels);

    % fuse the labels to make them dichotome (0 not in ROI, 1 if it is)
    % 
    % beneficial for next step, as triangles will more easily be included in
    % cases where it's vertices border multiple areas
    roiVertexLabelsBin(~isnan(roiVertexLabels)) = 1;
    roiVertexLabelsBin(isnan(roiVertexLabels)) = 0;

    % get the indices of the vertices which are in the ROI
    roiVertexIDs = find(roiVertexLabelsBin == 1);

    % determine the ROI only triangles (whether whether a face consists of
    % three vertices that belong to the ROI vertices)
    roiFacesLabels = ismember(gSurfPial.faces, roiVertexIDs);
    roiFacesLabels = all(roiFacesLabels, 2);

    % determine the indices of the faces that should be extracted
    extractRoiFaceIndices = find(roiFacesLabels == 1);

    % recreate a gifti with only the faces and vertices of the ROI
    [vertexMatrix, facesMatrix, vertexConversion] = extract3DFaces(gSurfPial, extractRoiFaceIndices);
    gROISurfPial = gSurfPial;
    gROISurfPial.vertices = vertexMatrix;
    gROISurfPial.faces = facesMatrix;



    %%
    %  create a hull version of the pial brain and determine their inverted normals
    %

    % create a raw version of the pial
    rawSurfPial 		= [];
    rawSurfPial.vert 	= double(gSurfPial.vertices);
    rawSurfPial.tri 	= gSurfPial.faces;

    % compute the hull
    cortexcoarser 		= coarserModel(rawSurfPial, 10000);   % approach 2 -> nicer result than approach 1 in terms of triangles but smaller hull
    origin			 	= [0 20 40];
    smoothrad 			= 25;
    mult 				= 0.5;
    [cortexsmoothed] 	= smoothModel(cortexcoarser, smoothrad, origin, mult);
    hullCortex 			= hullModel(cortexsmoothed);

    % create a high resolution version of the hullcortex (increased number of triangles)
    highHullCortex = refinepatch(gifti(hullCortex));
    for iRef = 2:numHullRefines
        highHullCortex  = refinepatch(highHullCortex);
    end

    % calculate the centers, face normals of the high resolution hull triangles
    hullTr              = triangulation(double(highHullCortex.faces), double(highHullCortex.vertices));
    hullInvFaceNormals	= faceNormal(hullTr);
    hullFaceCenters     = incenter(hullTr);



    %%
    %  determine which hull triangles (and hull vertices) are over the ROI
    %

    % create a vector to store which high-res full triangles are over the ROI
    hullROIFaces        = zeros(size(highHullCortex.faces, 1), 1);

    % loop through each face in the high-res hull
    for iFace = 1:size(highHullCortex.faces, 1)

        % determine one point at the center of the face and another
        % at x times the normal distance
        p1 = hullFaceCenters(iFace, :);
        p2 = p1 + hullInvFaceNormals(iFace, :) * hullNormalSearchDepthFactor;

        % determine the search distance
        %distance = sqrt(sum((p1 - p2).^2, 2));
        distance = hullNormalSearchDepthFactor;

        % the distances of the pial vertices to the hull triangle's center
        vertexDistances = gROISurfPial.vertices;
        vertexDistances = sqrt(sum((p1 - vertexDistances).^2, 2));

        % determine which vertices are within the search distance
        withinVertices = vertexDistances <= distance;

        % if no vertex is within search distance, go to next triangle
        if nnz(withinVertices) == 0
           continue; 
        end

        % determine the indices of the vertices within search distance
        withinVertexIndices = find(withinVertices);

        % determine the faces that have at least one vertex within search distance
        withinFaces = ismember(gROISurfPial.faces, withinVertexIndices);
        withinFaces = any(withinFaces, 2);

        % determine the indices of the faces within search distance
        withinFaceIndices = find(withinFaces);

        % loop through the triangles within the search distance
        for k = 1:length(withinFaceIndices)

            % retrieve the vertices that make up the triangle
            tr = gROISurfPial.faces(withinFaceIndices(k), :);
            t = [   gROISurfPial.vertices(tr(1), :); ...
                    gROISurfPial.vertices(tr(2), :); ...
                    gROISurfPial.vertices(tr(3), :)]';

            % determine whether the line (inverted hull normal) intersects with this pial triangle
            [inside, ~] = triangle_contains_line_exp_3d (t, p1, p2);
            if inside
                % intersects

                % flag the hull triangle as over ROI
                hullROIFaces(iFace) = 1;

                % continue to the next
                break;

            end

        end     % end of pial triangle loop

    end     % end of high-res hull triangle loop

    % determine the indices of the faces in the hull that are over the ROI
    hullROIFaceIndices = find(hullROIFaces);

    % determine the indices of the vertices in the hull that are over the ROI
    hullROIVertexIndices = highHullCortex.faces(hullROIFaceIndices, :);
    hullROIVertexIndices = hullROIVertexIndices(:);
    hullROIVertexIndices = unique(hullROIVertexIndices);


    %
    % appearantly some faces in the middle are not included. however, all their surrounding
    % faces are, which causes a strange situation where a triangle is not
    % indexed as part of the ROI (hullROIFaceIndices/hullROIFaces) but all of it's vertices are
    % (hullROIVertexIndices/hullROIVertex) considered part of the ROI and when
    % the face is plotted looks like they are. Upon extraction (by faces) it is
    % then revealed to have holes
    %
    % to fix this, we include the faces of which all vertices are ROI
    hullROIFaces = ismember(highHullCortex.faces, hullROIVertexIndices);
    hullROIFaces = all(hullROIFaces, 2);
    hullROIFaceIndices = find(hullROIFaces);

    % create a seperate ROI gifti
    clear gROI;
    [gROI.vert, gROI.tri, ~] = extract3DFaces(highHullCortex, hullROIFaceIndices);
    gROI = gifti(gROI);



    %%
    %  project each vertex outward (along it's normal) to the outer hull (the hull created using approach 1)
    %

    % calculate the hull using approach 1
    outwardHullCortex = hullModel(rawSurfPial);         % approach 1 -> actually covers the brain (is used later to project the smaller hull out to)

    % calculate vertex normals of the ROI
    roiTr = triangulation(double(gROI.faces), double(gROI.vertices));
    roiVertexNormals = vertexNormal(roiTr);

    % variable to hold the new positions of the vertices
    % (initially holds the old positions)
    projectedVertices = gROI.vertices;

    % loop through each vertex in the ROI and project it
    for iVertex = 1:size(gROI.vertices, 1)

        % determine one point at the vertex position and another at x times the vertex it's normal distance
        p1 = gROI.vertices(iVertex, :);
        p2 = p1 - roiVertexNormals(iVertex, :) * hullNormalSearchOutwardFactor;

        % determine the search distance
        distance = hullNormalSearchOutwardFactor * 2;

        % the distances of the outer hull vertices to the current inner hull's vertex
        vertexDistances = outwardHullCortex.vert;
        vertexDistances = sqrt(sum((p1 - vertexDistances) .^ 2, 2));

        % determine which vertices are within the search distance
        withinVertices = vertexDistances <= distance;

        % check if no vertex is within search distance
        if nnz(withinVertices) == 0

            % message
            warning(['Could not find any vertices on the outer hull within the search distance of inner hull vertex ', num2str(iVertex), ', skipping']);

            % go to the next vertex
            continue; 

        end

        % determine the indices of the vertices within search distance
        withinVertexIndices = find(withinVertices);

        % determine the faces that have at least one vertex within search distance
        withinFaces = ismember(outwardHullCortex.tri, withinVertexIndices);
        withinFaces = any(withinFaces, 2);

        % determine the indices of the faces within search distance
        withinFaceIndices = find(withinFaces);

        % loop through the outer hull's triangles within the search distance
        interPoint = [];
        for k = 1:length(withinFaceIndices)

            % retrieve the vertices that make up the triangle
            tr = outwardHullCortex.tri(withinFaceIndices(k), :);
            t = [   outwardHullCortex.vert(tr(1), :); ...
                    outwardHullCortex.vert(tr(2), :); ...
                    outwardHullCortex.vert(tr(3), :)]';

            % determine whether the line (inner hull's vertex normal) intersects with this outer hull's triangle
            [inside, interPoint] = triangle_contains_line_exp_3d (t, p1, p2);

            % break the loop on intersection
            if inside, break; end

        end     % end of outer hull triangle loop

        if isempty(interPoint)
            % no intersection found

            % message
            warning(['Could not find an intersection point for inner hull vertex ', num2str(iVertex), ', skipping']);

            % go to the next vertex
            continue; 

        else

            % store the coordinates
            projectedVertices(iVertex, :) = interPoint;

        end

    end     % end of inner hull vertex loop


    % create a new outer ROI gifti
    clear gROIOuter;
    gROIOuter.vert = projectedVertices;
    gROIOuter.tri = gROI.faces;
    gROIOuter = gifti(gROIOuter);

    % 
    %{
    giiResultFilepath = [subjData.subjectBaseFolder, hemi, '_', outputPrefix, '_hull.gii'];
    save(gROIOuter, giiResultFilepath);
    %}
    
end