%   Normalize and detrend a dataset
%
%   signalOut = cls_normalizeAndDetrendVolumeData(config, signal, tr) gives the normalized 
%   and detrended data according to the given configuration applied on the given data.
%
%   config = the configuration which holds the config.normAndDetrend
%            variable, the normAndDetrend variable holds the normalizing/detrend
%            actions (as a string cell array) that consequetively have to be applied
%            to the signal data. The following values can be given:
%               NormMatrixDemean        -> normalize the entire matrix around 0 (zero-center, std will stay the same)
%               NormMatrixPercSignal    -> normalize the entire matrix around 100 and scale to percentage
%                                          signal change from the overall mean. Also called 'Grand Mean Scaling'
%               NormMatrixStd           -> standardize the entire matrix
%               NormVolumesDemean       -> normalize each volume around 0 (zero-center, std per volume will stay the same)
%               NormVolumesPercSignal   -> normalize each volume to the mean of that volume and scale to percentage signal
%                                          change from volume mean. Also called 'Global scaling' or 'Proportional scaling'
%               NormVolumesStd          -> standardize each volume to the mean and std of that volume
%               NormVoxelsDemean        -> normalize each voxel around 0 (zero-center, std per voxel will stay the same)
%               NormVoxelsPercSignal    -> normalize each voxel (feature) to the mean of that voxel (feature) and scale to
%                                          the percentage signal change per voxel. (normalization scheme widely used, also in brainvoyager)
%               NormVoxelsStd           -> standardize each voxel (feature) to the mean and std of that voxel (feature)
%               MatlabDetrend           -> matlab linear detrending over volumes/time. Will return the matrix to the input grand mean
%               RTDetrend           -   -> Realtime software linear and nonlinear (quadratic, cubic) detrending over
%                                          volumes/time (original by P. Andersson). Will return the matrix to the input grand mean
%               SPMDetrend              -> SPM detrending, removes linear and nonlinear (quadratic, cubic) trends. Will keep the grand mean. Recommended
%               SPMDetrend2             -> SPM detrending, removes linear and nonlinear (quadratic, cubic) trends.
%                                          Will return the matrix to the input grand mean. Recommended method of detrending
%
%
%   signal = data matrix to normalize/detrend. The rows represent the 
%            voxels and the columns represent the volumes/timepoints
%
%       tr = the repetition time (in seconds) of the fmri sequence, only relevant when using
%            'SPMDetrend', otherwise leave empty
%
%
%
%   Example:
% 
%		t = 0:.3:150;
%		data = 2*sin(t) + t;
% 		data = [data;data] + rand(2, length(data)) * 20;
%
%		classConfig = [];
% 		classConfig.normAndDetrend = {'NormVoxelsPercSignal', 'SPMDetrend'};
%		detrendedSignal = cls_normalizeAndDetrendVolumeData(classConfig, data, 1.5);
%
%       plot(data'); figure; plot(detrendedSignal')
%
%
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [signalOut, constantsOut] = cls_normalizeAndDetrendVolumeData(config, signal, tr, constantsIn)
    if ~exist('constantsIn', 'var') || isempty(constantsIn)
        constantsIn = [];
    end
    if ~exist('tr', 'var') || isempty(tr)
        tr = [];
    end
	
    % create a copy of the dataset to apply the actions on
    signalOut = signal;
    constantsOut = constantsIn;

    % loop through the actions
    for i = 1:length(config.normAndDetrend)

        % pick the current action
        action = config.normAndDetrend{i};

        
        if strcmpi(action, 'NormMatrixDemean') == 1

            % normalize the entire matrix around 0 (zero-center, std will stay the same)
            signalMean = mean(signalOut(:));
            signalOut = signalOut - signalMean;

        elseif strcmpi(action, 'NormMatrixPercSignal') == 1

            % normalize the entire matrix around 100 and scale the matrix to percentage signal change from the overall mean
            % std will change. Also called 'Grand Mean Scaling'
            signalMean = mean(signalOut(:));
            signalOut = signalOut / signalMean * 100;

        elseif strcmpi(action, 'NormMatrixStd') == 1

            % standardize the entire matrix
            signalMean = mean(signalOut(:));
            signalStd = std(signalOut(:));
            signalOut = (signalOut - signalMean) / signalStd;

        elseif strcmpi(action, 'NormVolumesDemean') == 1

            % normalize each volume around 0 (zero-center, std per volume will stay the same)
            signalMeans = mean(signalOut, 1);
            signalOut = bsxfun(@minus, signalOut, signalMeans);

            % other way to do the same thing
            %signalMeans = mean(signalOut, 1);
            %signalMeans = repmat(signalMeans, size(signalOut,1), 1);
            %signalOut = signalOut - signalMeans;


        elseif strcmpi(action, 'NormVolumesPercSignal') == 1

            % normalize each volume to the mean of that volume and scale to signal change from volume mean
            % also called 'Global scaling' or 'Proportional scaling'
            signalMeans = mean(signalOut, 1);
            signalOut = bsxfun(@rdivide, signalOut, signalMeans);
            signalOut = signalOut * 100;

        elseif strcmpi(action, 'NormVolumesStd') == 1

            % standardize each volume (to the mean and std of that volume)
            signalMeans = mean(signalOut, 1);
            signalStds = std(signalOut);
            signalStds(signalStds == 0) = 1;
            signalOut = bsxfun(@minus, signalOut, signalMeans);
            signalOut = bsxfun(@rdivide, signalOut, signalStds);

            % other way to do the same thing
            %timeMean = mean(signalOut, 1);
            %timeMean = repmat(timeMean, size(signalOut,1), 1);
            %timeStd = std(signalOut);
            %timeStd = repmat(timeStd, size(signalOut,1), 1);
            %signalOut = (signalOut - timeMean) ./ timeStd;

        elseif strcmpi(action, 'NormVoxelsDemean') == 1

            % normalize each voxel around 0 (zero-center, std per voxel will stay the same)
            signalMeans = mean(signalOut, 2);
            signalOut = bsxfun(@minus, signalOut, signalMeans);

        elseif strcmpi(action, 'NormVoxelsPercSignal') == 1

            % normalize each voxel (feature) to the mean of that voxel (feature) and scale to the percentage signal change per voxel
            % (normalization scheme widely used, also in brainvoyager)
            signalMeans = mean(signalOut, 2);
            signalOut = bsxfun(@rdivide, signalOut, signalMeans);
            signalOut = signalOut * 100;
            
            %{
            if ~isempty(constantsOut)
                for iOther = 1:length(constantsOut)
                    constantsOut(iOther) = (constantsOut(iOther) / signalMeans(iOther)) * 100;
                end
            end
            %}

        elseif strcmpi(action, 'NormVoxelsStd') == 1

            % standardize each voxel (feature) to the mean and std of that voxel (feature)
            signalMeans = mean(signalOut, 2);
            signalStds = std(signalOut, [], 2);
            %signalStds(signalStds == 0) = 1;
            signalOut = bsxfun(@minus, signalOut', signalMeans')';
            signalOut = bsxfun(@rdivide, signalOut', signalStds')';

        elseif strcmpi(action, 'MatlabDetrend') == 1

            % store the grand mean before detrending
            signalMean = mean(signalOut(:));

            % matlab linear detrending (the matlab detrend function removes linear trends from each column, hence the transposing)
            % Note: this function normalizes the values so that the overall mean is 0, which is why we add the grand mean later
            % Note 2: The std will change as a result of the detrending but should still be valid, meaning that if a matrix with percentage
            % signal change goes in, the detrended output can still be considered in percentage signal change
            signalOut = detrend(signalOut')';

            % return back to the grand mean that it had before detrending
            signalOut = signalOut + signalMean;

        elseif strcmpi(action, 'RTDetrend') == 1 

            % store the grand mean before detrending
            signalMean = mean(signalOut(:));

            % detrending (realtime software - Andersson; removes linear and nonlinear trends from each row)
            % Note: adjusts the mean a bit from the input
            % Note 2: The std will change as a result of the detrending but should still be valid, meaning that if a matrix with percentage
            % signal change goes in, the detrended output can still be considered in percentage signal change
            % Note 3: lamda may vary, using 200 now
            signalOut = RT_detrend(signalOut, true, 200);      

            % return back to the grand mean that it had before detrending
            signalOut = signalOut - mean(signalOut(:)) + signalMean;


        elseif strcmpi(action, 'SPMDetrend') == 1

            % check if TR is given
            if isempty(tr) || ~isnumeric(tr)
               fprintf(2, 'Error: a valid TR (repetition time) is required when using SPMDetrend\n');
               return;
            end
            
            % spm detrending, removes linear and nonlinear (quadratic, cubic) trends
            % Note: this function keeps the mean at the input mean
            hpf = 128;
            K.RT = tr;
            K.row=ones(size(signalOut, 2));
            K.HParam=hpf;
            nK=spm_filter(K);
            fdata=zeros(size(signalOut));
            for i=1:size(fdata, 1) 
                fdata(i,:)=spm_filter(nK,signalOut(i,:)');
            end
            signalOut = fdata;
            clear fdata;

        elseif strcmpi(action, 'SPMDetrend2') == 1

            % store the grand mean before detrending
            signalMean = nanmean(signalOut(:));

            % spm detrending, removes linear and nonlinear (quadratic, cubic) trends
            % Note: this function normalizes the values so that the overall mean is 0
            signalOut = spm_detrend(signalOut', 3)';

            % return back to the grand mean that it had before detrending
            signalOut = signalOut + signalMean;

        end
        
    end
    
end
