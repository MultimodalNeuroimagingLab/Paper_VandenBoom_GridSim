%
%   Perform SVM classification with multiple (>2) classes
%	[group, distance, ties] = multisvmclassify_new(multisvm, features)
%
%   	multisvm    = 
%   	features    = 
%
%
%   Returns:
%       group  	    =
%       distance  	=
%       ties  	    =
%
%
%	Adapted from Mark Bruurmijn
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [group, distance, ties] = multisvmclassify_new(multisvm, features)
	
	% Get the type from the model. (The function multisvmtrain stores this in the structure.)
	type = multisvm{2, 1};
	
    % create the resulting classification array
	group = NaN( size(features,1), 1 );
	
    % 
	switch type
		case '1vsall'
            
            % loop through each row of features to classify
			for f = 1 : size(features,1)
                
				% calculate for each model what the sample's distance would be 
                % to the hyperplane (invert the distance).
                % Each model is one class versus the other classes
                distance = [];
                for m = 1:size(multisvm, 2)
                   [label, score] = predict(multisvm{1,m}, features(f,:));
                   distance(m) = -score(1, 1);
                end
                
                % it is unlikely/impossible that there will be a tie
                % with this type of classification
                ties = NaN;
                
                % pick the class with the highest distance value
                % the model/class where it is most likely to be that class
                % instead of the other classes in that model is picked
				[~,highest] = max( distance );
                
                % store the corresponding group label
				group(f) = setdiff(unique(multisvm{1, highest}.Y), 0);
                
			end
			
			
		case '1vs1'
            
            %
            [group, ~, similarFrequency] = mode( cell2mat(cellfun( @(m) predict(m, features), {multisvm{1,:}}, 'Un', false )), 2 );
            
            %
			distance = [];
            
            % determine the number of ties
			ties = nnz( cellfun(@length,similarFrequency) > 1 );
            
	end
	
end