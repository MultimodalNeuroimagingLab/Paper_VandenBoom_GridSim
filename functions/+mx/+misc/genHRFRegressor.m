%
%   Generate a regressor by convolving the HRF with the given design (onsets/durations). 
%   hrf = genHRFRegressor(TR, units, onsets, durations, numVols)
%   
%       TR           		 = repetition time
%       units           	 = onset and duration units; either 'seconds' or 'scans'
%       onsets           	 = stimulus onsets
%       durations            = stimulus durations
%       numVols         	 = total number of volumes/scans
%
%
%   Adapted from SPM12 routines
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function hrf = genHRFRegressor(TR, units, onsets, durations, numVols)

    % if there are less durations than onsets
    if length(durations) < length(onsets)

        % message
        warning('Warning: less durations than onsets are given, replicating the first duration as the duration for each onset');
        
        % replicate the first duration as the duration for each onset
        durations = repmat(durations(1), 1, length(onsets)); 
        
    end
    
    if isrow(onsets)
        onsets = onsets';
    end
    
    if isrow(durations)
        durations = durations';
    end
    
    fMRI_T = 16;
    fMRI_T0 = 8;


    %
    dt     = TR / fMRI_T;

    % 
    if strcmpi(units, 'scans') == 1
        %TR = fMRI_T * dt;
    elseif strcmpi(units, 'seconds') == 1
        TR = 1;
    else
        error('Unknown unit "%s".', units);
    end
    
    %-Canonical hemodynamic response function
    [bf, p]      = spm_hrf(dt, [], fMRI_T);
    bf = spm_orth(bf);

    % 
    u     = onsets .^ 0;
    u     = spm_orth(u);
    if ~any(durations)
        u     = u / dt;
    end



    %-Create stimulus functions (32 bin offset)
    %======================================================================
    ton       = round(onsets * TR / dt) + 33;               % onsets
    tof       = round(durations * TR / dt) + ton + 1;          % offset
    sf        = sparse((numVols * fMRI_T + 128), size(u, 2));
    ton       = max(ton,1);
    tof       = max(tof,1);
    for j = 1:length(ton)
        if size(sf,1) > ton(j)
            sf(ton(j), :) = sf(ton(j), :) + u(j, :);
        end
        if size(sf,1) > tof(j)
            sf(tof(j), :) = sf(tof(j), :) - u(j, :);
        end
    end
    sf        = cumsum(sf);                             % integrate
    sf        = sf(1:(numVols * fMRI_T + 32), :);       % stimulus



    %-1st order terms
    %==========================================================================
    X     = [];
    Fc    = [];

    ind   = [];
    ip    = [];
    for k = 1:size(sf,2)
        for p = 1:size(bf,2)
            x = sf(:,k);
            d = 1:length(x);
            x = conv(full(x),bf(:,p));
            x = x(d);
            X = [X x];

            %-Indices and regressor names
            %------------------------------------------------------------------
            ind(end + 1)   = size(X,2);
            ip(end + 1)    = k;
        end
    end
    Fc(end + 1).i = ind;
    Fc(end).p     = ip;


    if ~isempty(X)
        X = X((0:(numVols - 1)) * fMRI_T + fMRI_T0 + 32,:);
    end

    %
    hrf = X';

end