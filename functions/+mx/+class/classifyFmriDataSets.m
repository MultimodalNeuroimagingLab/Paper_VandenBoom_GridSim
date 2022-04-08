%
%   Perform classification based on one or more fMRI datasets, either per set (within), between sets (between)
%   or concatenated sets (combine)
%
%   Usage: [results, success, signalSets] = classifyFmriDataSets(config, featureVoxels, sets, setMethod)
%
%   config                  = 
%   featureVoxels           = Indicates which voxels to use as features from the datasets, this 
%                             parameter can be used in several ways:
%                               [indices]    = one dimensional array that holds the linear indices of
%                                              the voxels in the input datasets to use as features
%                               [x-y-z]      = three-dimensional array (binary mask) to indicate which 
%                                              features to use for classification. The mask should have
%                                              the exact same dimension sizes as the first three dimensions
%                                              of the input datasets
%                               [cell array] = array specifying multiple classifications should
%                                              be performed using different feature-sets. Each cell
%                                              will result in classification using a specific set of 
%                                              features. The features to use for each cell (classification) 
%                                              are defined as the options above ([], [indices] or [x-y-z] mask)
%   sets                    = struct array of fmri datasets
%   sets(#).volumes         = The volume input data to pick the features from and classify on. The data 
%                             should hold 3D volume data in the first 3 dimensions, and time in the 4th
%   sets(#).conditions      = The experimental conditions per trial
%   sets(#).units           = The units the onsets and durations are in (can either be 'seconds' or 'scans')
%   sets(#).onsets          = A vector with the trial onsets
%   sets(#).durations       = A vector with the trial durations
%   sets(#).tr              = The repetition time used for the dataset
%   setMethod               = The dataset classification method, either:
%                               'within': classify per dataset, cross-validation
%                                         will be performed within each set
%                               'between': train on one dataset, classify on the
%                                          other. Exactly two input sets are expected
%                               'combine': concatinate all input sets and perform
%                                          cross-validation over the joint set
% 
%   Returns: 
%       results             = struct array with classification result(s). The first dimension (rows)
%                             represents the different feature-sets. The second dimension will represent
%                             the input datasets, the size of the second dimension depends on the classification
%                             method: 'within' will give one result per dataset (for each featureset); 'combined' and
%                             'between' will give one result per 2 or more datasets (for each featureset)
%                             Note: empty results structs can be returned on failure, use the success return
%                             variable to check
%       success             = vector which - for each selection of features (if featureVoxels is not cellarray
%                             then 1 flag; if featureVoxels is a cellarray then a flag for each cell) - indicates
%                             success (1) or failure (0) given that featureset
%       signalSets          = struct-array with the signal-sets and meta-info (onsets, durations, 
%                             conditions) that were used for classification. The first dimension (rows)
%                             represents the different feature-sets, the second dimension (columns) holds
%                             the input datasets.
%                             
%
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [results, success, signalSets] = classifyFmriDataSets(config, featureVoxels, sets, setMethod)

    splitConfig = [];
    splitConfig.numPoints = 1;
    splitConfig.threads = 1;
    splitConfig.silent = 'yes';
    
    % forward to split version of this function
    if nargout > 2
        [results, success, signalSets] = mx.class.classifyFmriDataSets_split(   config, ...
                                                                                featureVoxels, ...
                                                                                sets, ...
                                                                                setMethod, ...
                                                                                splitConfig);
    else
        [results, success] = mx.class.classifyFmriDataSets_split(   config, ...
                                                                    featureVoxels, ...
                                                                    sets, ...
                                                                    setMethod, ...
                                                                    splitConfig);
    end

end
