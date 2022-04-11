%
%   Step 3 - Select the 3D sample-points that have a searchlight classification score above a given threshold
%   [SS, filenameSuffix] = s3_sampleSelection(SS)
%
%       SS              = Structure that holds the sample-points and the searchlight classification results
%
%
%   Returns: 
%       SS              = The input structure with the only the sample-points remaining that were "higher" than the threshold
%       filenameSuffix  = Depending on the thresholding method a suggestion for a filename suffix
%
%
%   Copyright 2020, Max van den Boom (Multimodal Neuroimaging Lab, Mayo Clinic, Rochester MN)

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [SS, filenameSuffix] = s3_sampleSelection(SS)

    %focalOrNonFocal = 'foc';       % use the focal of non-focal scores
    focalOrNonFocal = 'nofoc'; 
    %threshMethod = 'static';       % theshold at a static value
    %threshStaticThreshold = 44;    % all samples with scores below this will be removed (if the method is set to static
    threshMethod = 'P95';           % threshold at the 95th percentile (all samples with scores below will be removed)


    %%
    %  Retrieve the root path and make sure dependencies can be found
    %
    addpath('../');
    gridSimRoot = gridSimRootPath;
    addpath([gridSimRoot, filesep, 'functions']);
    addpath(genpath([gridSimRoot, filesep, 'external']));



    %%
    %  Remove the empty (nan) samples
    %

    % check to use the focal or non-focal results
    if strcmpi(focalOrNonFocal, 'foc') == 1
        % focal

        remSamples = find(isnan(SS.sampleSearchFocalResult));

    elseif strcmpi(focalOrNonFocal, 'nofoc') == 1
        % non-focal

        remSamples = find(isnan(SS.sampleSearchNonFocalResult));

    else
        fprintf(2, 'Error: unknown set, choose between focal (''foc'') or non-focal (''nofoc'')\n');
        return;
    end

    % message
    disp (['Samples in set: ', num2str(size(SS.samplePositions, 1))]);
    disp (['Samples to be removed (nan values): ', num2str(length(remSamples))]);

    % remove the nan-samples
    SS = removeFromSampleSetByIndices(SS, remSamples);

    % message remaining
    disp (['Samples remaining: ', num2str(size(SS.samplePositions, 1))]);



    %%
    %  Remove the samples below the threshold
    %

    % determine the menthod of thresholding

    if strcmpi(threshMethod, 'static') == 1
        % static

        % use the given static threshold
        threshold = threshStaticThreshold;

        filenameSuffix = ['St', num2str(threshStaticThreshold)];
        
    elseif strcmpi(threshMethod, 'P95') == 1
        % 95th percentile

        % calculate cutoff
        if strcmpi(focalOrNonFocal, 'foc') == 1
            m = mean(SS.sampleSearchFocalResult);
            sd = std(SS.sampleSearchFocalResult);
            
            filenameSuffix = 'P95foc';
        
        else
            m = mean(SS.sampleSearchNonFocalResult);
            sd = std(SS.sampleSearchNonFocalResult);
            
            filenameSuffix = 'P95nofoc';
            
        end
        threshold = m + 1.64485 * sd;

    else
        fprintf(2, 'Error: unknown thresholding method was set\n');
        return;
    end

    % message
    disp (['Threshold is : ', num2str(threshold)]);

    % check to use the focal or non-focal results
    if strcmpi(focalOrNonFocal, 'foc') == 1
        % focal

        remSamples = find(SS.sampleSearchFocalResult < threshold);

    else
        % non-focal

        remSamples = find(SS.sampleSearchNonFocalResult < threshold);

    end


    % message
    disp (['Samples in set: ', num2str(size(SS.samplePositions, 1))]);
    disp (['Samples to be removed (by threshold): ', num2str(length(remSamples))]);

    % remove the nan-samples
    SS = removeFromSampleSetByIndices(SS, remSamples);

    % message remaining
    disp (['Samples remaining: ', num2str(size(SS.samplePositions, 1))]);

end

