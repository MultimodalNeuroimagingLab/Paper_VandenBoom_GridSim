%   
%   Step 1 - generate random 3D sample-points on a given hull
%   [samplesFilename] = s1_generateSamples(bids_rootPath, bids_sub, hemi)
% 
%       bids_rootPath    = path to the BIDS root-directory where the data is located
%       bids_sub         = the subject to create a hull for
%       hemi             = the hemisphere that this step is applied on
%
%   Returns: 
%       samplesFilename  = The filename of the sample-set that was generated and stored
%
%   Copyright (C) 2019 Max van den Boom  (Multimodal Neuroimaging Lab, Mayo Clinic, Rochester MN)

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [samplesFilename] = s1_generateSamples(bids_rootPath, bids_sub, hemi)

	numSamples = 1000000;

    
    %%
    % retrieve the root path and make sure dependencies can be found
    %
    addpath('../');
    gridSimRoot = gridSimRootPath;
    addpath([gridSimRoot, filesep, 'functions']);
    addpath(genpath([gridSimRoot, filesep, 'external']));

    

    %%
    %  Build paths and read the files
    %
    
    % build the paths to the data
    bids_simPath     = fullfile(bids_rootPath, 'derivatives', [hemi, '_simulations'], ['sub-' bids_sub]);

	% read the hull
	gHull = gifti(fullfile(bids_simPath, [hemi, '_ext_hull.gii']));



	%%
	%  Generate samples
	%

	% create an empty sample-set
	SS = [];

	% generate random 3D sample points on a given 3D object
	[SS.samplePositions, SS.sampleTriangles] = genRandomSamplePoints3DObject(gHull, numSamples);
    
    
    
	%%
	%  Generate sample rotations
	%

	% generate a sample-list of grid rotations
	SS.sampleGridRotations = rand(size(SS.samplePositions, 1), 1) * 180;

    
    
    %%
    %  save the samples
    %
    
    % generate a unique samples filename
    samplesFilename = ['sampleSet_', datestr(now,'yyyymmdd_HHMMSS'), '.mat'];

    % save the samples
    outputFilename   = fullfile(bids_simPath, samplesFilename);
    save(outputFilename, 'SS');
    
end