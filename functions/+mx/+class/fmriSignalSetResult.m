%   
%   Structure to store a fmri signal classification result in
%	h = cls_fmriSignalSetResult(varargin)
%
%
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function h = cls_fmriSignalSetResult(varargin)

    switch nargin
        case 0
            h   = struct(   'signalFormatting', [], ...
                            'epochedScanIndices', [], ...
                            'linearScanIndices', [], ...
                            'conditions', [], ...
                            'durations', [], ...
                            'numCorrect', [], ...
                            'numIncorrect', [], ...
                            'arrCorrect', [], ...
                            'correctPerCondition', [], ...
                            'accTotal', [], ...
                            'accPerCondition', [], ...
                            'durationResults', [], ...
                            'numFeatinput', [], ...             // number of features given in the call to classification
                            'numFeatStart', [], ...             // number of features at the start of classification (the 'minVoxelMeanForFeature' and 'minFeatureMean' options in config might have removed features)
                            'numFeatNo0OrNan', [], ...          // number of features after detrending and another round of removing 0 and nan features
                            'numFeatUsed', [] ...               // the actual number of features used in classification
                        );
            
        case 1
            h = cls_fmriSignalSetResult();
            
            % try to transfer the fields
            
            
        otherwise
            error('Unknown input');
            
    end
    
end
