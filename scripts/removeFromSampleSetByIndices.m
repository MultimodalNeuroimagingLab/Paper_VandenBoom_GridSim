%   
%   Remove samples (and their corresponding grids) from a sample-set struct
%   [SS] = removeFromSampleSetByIndices(SS, remSampleIndices)
% 
%       SS                  = the sample-set from which the samples need to be removed
%       remSampleIndices    = the indices of the samples that need to be removed
%
%   Returns: 
%       SS                  = the resulting sample-set from which the samples are removed
%
%   Copyright (C) 2020 Max van den Boom  (Multimodal Neuroimaging Lab, Mayo Clinic, Rochester MN)

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [SS] = removeFromSampleSetByIndices(SS, remSampleIndices)

    % remove the samples
    SS.samplePositions(remSampleIndices, :) = [];
    SS.sampleTriangles(remSampleIndices) = [];
    SS.sampleGridRotations(remSampleIndices) = [];
    if isfield(SS, 'sampleSearchFocalResult')
        SS.sampleSearchFocalResult(remSampleIndices) = [];
    end
    if isfield(SS, 'sampleSearchNonFocalResult')
        SS.sampleSearchNonFocalResult(remSampleIndices) = [];
    end
    if isfield(SS, 'sampleSearchNonFocalResultCounter')
        SS.sampleSearchNonFocalResultCounter(remSampleIndices) = [];
    end

    % remove the grids
    if isfield(SS, 'sampleGrids') || isfield(SS, 'sampleGridsProjected') || isfield(SS, 'sampleGridsTriangles')
        
        if ~isfield(SS, 'gridConfig')
            disp('Error: could not remove grids from sampleset, gridConfig not available');
            return;
        end
        if ~isfield(SS.gridConfig, 'elecTotal')
            disp('Error: could not remove grids from sampleset, the variable ''elecTotal'' is not available in the gridConfig');
            return;
        end
        
        % determine the electodes of the sample(grids) that need to be removed
        sampleGridElectrodes = ((remSampleIndices - 1) * SS.gridConfig.elecTotal) + 1;
        [~, Y] = meshgrid(sampleGridElectrodes, 0:SS.gridConfig.elecTotal - 1);
        if iscolumn(sampleGridElectrodes), sampleGridElectrodes = sampleGridElectrodes';    end
        sampleGridElectrodes = Y + sampleGridElectrodes;
        sampleGridElectrodes = sampleGridElectrodes(:);
        clear Y;
        
        % remove the electrodes of the sample(grids)
        if isfield(SS, 'sampleGrids')
            SS.sampleGrids(sampleGridElectrodes, :) = [];
        end
        if isfield(SS, 'sampleGridsProjected')
            SS.sampleGridsProjected(sampleGridElectrodes) = [];
        end
        if isfield(SS, 'sampleGridsTriangles')
            SS.sampleGridsTriangles(sampleGridElectrodes) = [];
        end
        
    end
    
    % remove the samples before transformation 
    if isfield(SS, 'preTrans_samplePositions')
        SS.preTrans_samplePositions(remSampleIndices, :) = [];
    end
    if isfield(SS, 'preTrans_sampleGrids')
        SS.preTrans_sampleGrids(sampleGridElectrodes, :) = [];
    end
    
end