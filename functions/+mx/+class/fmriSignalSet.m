%   
%   Structure to store fmri signal data in
%	h = fmriSignalSet(varargin)
%
%
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function h = fmriSignalSet(varargin)

    switch nargin
        case 0
            h   = struct(   'inputUnits', [], ...
                            'inputOnsets', [], ...
                            'inputDurations', [], ...
                            'inputConditions', [], ...
                            'inputTr', [], ...
                            'units', [], ...
                            'onsets', [], ...
                            'durations', [], ...
                            'conditions', [], ...
                            'onsetsScans', [], ...
                            'signalData', [], ...
                            'inputVolumesDims', [] ...
                        );
            

        case 1
            h = mx.class.fmriSignalSet();
            
            % try to transfer the fields
            
            
            
        otherwise
            error('Unknown input');
            
    end
    
end
