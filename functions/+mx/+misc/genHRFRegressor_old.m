%
%   Generate a regressor by convolving the HRF with the given design (onsets/durations). 
%	hrf = genHRFRegressor_old(TR, units, onsets, durations, numberOfVols)
%   
%       TR           		 = repetition time
%       units           	 = onset and duration units; either 'seconds' or 'scans'
%       onsets           	 = stimulus onsets
%       durations            = stimulus durations
%       numberOfVols         = total number of volumes/scans
%
%
%   Adapted from SPM8 routines
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function hrf = genHRFRegressor_old(TR, units, onsets, durations, numberOfVols)

    % if there are less durations than onsets
    if length(durations) < length(onsets)

        % message
        warning('Warning: less durations than onsets are given, replicating the first duration as the duration for each onset');
        
        % replicate the first duration as the duration for each onset
        durations = repmat(durations(1), 1, length(onsets)); 
        
    end
    
    % durations of 0 (scans or seconds) will be set to 1
    durations(durations == 0) = 1;
    
    % convert onsets and duration to scans (if in seconds)
    if strcmpi(units, 'seconds') == 1
        
        % convert scan onsets in seconds to scans
        onsets_Scans = round(onsets / TR);
        durations_Scans = round(durations / TR);
        
        % the durations could be rounded to 0 on a high TR, this is not possible, then assume 1 scan
        durations_Scans(durations_Scans < 1) = 1;
        
        % the first scan could be rounded to 0, there is no 0th
        % scan only the 1st scan, so correct this
        if onsets_Scans(1) == 0
           onsets_Scans(1) = 1; 
        end

    elseif strcmpi(units, 'scans') == 1
        
        % convert scan onsets in seconds to scans
        onsets_Scans = onsets;
        durations_Scans = durations;
        
    else
        
        % message
        fprintf(2, ['Unknown type of units: ', units, '\n']);
        return;
        
    end
    
    % create a binary regressor
    reg = zeros(numberOfVols, 1);
    for i=1:length(onsets_Scans)
        if onsets_Scans(i) < numberOfVols && (onsets_Scans(i) + durations_Scans(i)) < numberOfVols
            reg(onsets_Scans(i):onsets_Scans(i) - 1 + durations_Scans(i)) = 1;
        end
    end
    
        
    % calculate the HRF
    xBF.dt = TR;
    xBF.name = 'hrf';
    xBF.T = 16;
    bf = spm_get_bf(xBF);
    U.u = reg;
    U.name = {'act'};
    hrf = spm_Volterra(U, bf.bf);
    hrf = hrf';

end
