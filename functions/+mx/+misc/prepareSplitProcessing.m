%	
%   Prepare for split processing, splitting up the input for speed
%   optimization (multi-threading) or memory maximization
%
%   [retVal, splitConfig, numSets, numPoints, ranges] = prepareSplitProcessing(totalPoints, splitConfig)
%
%   totalPoints           = the total amount of datapoints to be processed
%   splitConfig           = the configuration on splitting the input
%   splitConfig.numSets   = (optional) split input into x sets of datapoints [0 = no splitting into sets]
%   splitConfig.numPoints = (optional) splits the input in sets of x datapoints [0 = no splitting into sets]
%   splitConfig.threads   = (optional) number of threads to run the sets on (requires the
%                           set to be split using either numSets or numPoints); if left empty will use
%                           the number of cores as the number of threads
%   splitConfig.silent    = don't show messages (warnings and info), errors will still be shown
% 
%
%   Returns: 
%       retVal            = return value; 1 = succes, 0 = failure
%       splitConfig       = the updated split configuration
%       numSets           = the number of sets the datapoint have been split into
%       numPoints         = the number of datapoints in each set (the last set may
%                           be smaller, to accomodate the remainder of the total)
%       range             = the start and end indices of the subsets (referring to the index in the total set) 
%
%
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [retVal, splitConfig, numSets, numPoints, ranges] = prepareSplitProcessing(totalPoints, splitConfig)
    if exist('splitConfig', 'var') == 0
        splitConfig = [];
    end
    if ~isfield(splitConfig, 'numSets') || isempty(splitConfig.numSets) || splitConfig.numSets < 0
        splitConfig.numSets = 0;
    end
    if ~isfield(splitConfig, 'numPoints') || isempty(splitConfig.numPoints) || splitConfig.numPoints < 0
        splitConfig.numPoints = 0;
    end
    if ~isfield(splitConfig, 'threads') || isempty(splitConfig.threads) || splitConfig.threads < 1
        splitConfig.threads = feature('numcores');
    end
    if ~isfield(splitConfig, 'silent') || isempty(splitConfig.silent)
        splitConfig.silent = 'no';
    end
    
    % make sure only one of the two split methods is set
    if splitConfig.numSets > 0 && splitConfig.numPoints > 0

        % error message
        fprintf(2, 'Error: make sure to select only one of the two split methods. Either split into a\nnumber of sets (use numSets), or split based on number of datapoints per set (use numPoints)\n');
        
        % set failure and return
        retVal = 0;
        numSets = 0;
        numPoints = 0;
        ranges = [];
        return;    
        
    end
    
    % check if split by number of point or number of sets
    if splitConfig.numPoints > 0
        % split by number of points per set
        
        % calculate the number of sets and store the number of points per set
        numSets = ceil(totalPoints / splitConfig.numPoints);
        numPoints = splitConfig.numPoints;
        
    else
        % split by number of sets
        
        % calculate the number of point per set and store the number of sets
        numPoints = ceil(totalPoints / splitConfig.numSets);
        numSets = splitConfig.numSets;
        
    end
    
    % messages
    if strcmpi(splitConfig.silent, 'yes') == 0
        disp(['Processing ', num2str(totalPoints), ' datapoints, split into ', num2str(numSets), ' sets with a setsize of ', num2str(numPoints), ' datapoints']);
    end
    
    % check if there is just one single subset
    if numSets == 1
        
        % message
        if strcmpi(splitConfig.silent, 'yes') == 0
            warning('Split resulted in just one single subset, calling the split version of this function was unnecessary. Reverting to single thread processing');
        end
        
        % adjust the split config to require one thread
        splitConfig.threads = 1;
        
    end
    
    % set the point ranges for each set
    ranges = nan(numSets, 2);
    iCounter = 1;
    for i = 1:numSets
        if iCounter + numPoints <= totalPoints
            ranges(i, :) = [iCounter, iCounter + numPoints - 1];
        else
            ranges(i, :) = [iCounter, totalPoints];
        end
        iCounter = iCounter + numPoints;
    end
    
    % check if multiple threads are required
    if splitConfig.threads > 1
        % multiple threads

        % check if a pool with the required number of workers is already running
        pool = gcp('nocreate');
        if isempty(pool) || pool.NumWorkers ~= splitConfig.threads

            % message
            if strcmpi(splitConfig.silent, 'yes') == 0
                disp(['No pool found or pool does not have the required number of workers']);
            end
            
            % shutdown parallel pool if already running
            delete(gcp('nocreate'));

            % setup a pool of threads
            parpool('local', splitConfig.threads);

        end
        
    end
    
    % set success
    retVal = 1;
    
end
