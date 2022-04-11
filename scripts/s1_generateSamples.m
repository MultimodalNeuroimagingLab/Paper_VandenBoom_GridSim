%   
%   Step 1 - generate random 3D sample-points on a given hull
%   [SS, filenameSuffix] = s1_generateSamples(bids_rootPath)
% 
%       hullPath         = path to the hull gifti file, this is where the 3D sample point will be generated on top of
%
%
%   Returns: 
%       SS               = A structure that holds the sample-points
%       filenameSuffix   = A suggestion for a filename suffix (effectively the date/time in the format: yyyymmdd_HHMMSS)
%
%   Copyright (C) 2019 Max van den Boom  (Multimodal Neuroimaging Lab, Mayo Clinic, Rochester MN)

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [SS, filenameSuffix] = s1_generateSamples(hullPath)

	numSamples = 1000000;

    
    %%
    % retrieve the root path and make sure dependencies can be found
    %
    addpath('../');
    gridSimRoot = gridSimRootPath;
    addpath([gridSimRoot, filesep, 'functions']);
    addpath(genpath([gridSimRoot, filesep, 'external']));

    

    %%
    %  Read the hull
    %
	gHull = gifti(hullPath);



	%%
	%  Generate samples
	%

	% create an empty sample-set
	SS = [];

	% generate random 3D sample points on a given 3D object
	[SS.samplePositions, SS.sampleTriangles] = mx.three_dimensional.genRandomSamplePoints3DObject(gHull, numSamples);
    
    
    
	%%
	%  Generate sample rotations
	%

	% generate a sample-list of grid rotations
	SS.sampleGridRotations = rand(size(SS.samplePositions, 1), 1) * 180;
    
    
    %%
    % generate and return a file suffix
    filenameSuffix = [datestr(now,'yyyymmdd'), 'T', datestr(now,'HHMMSS')];
    
end