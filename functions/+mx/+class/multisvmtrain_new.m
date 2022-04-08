%
%   Train a SVM for classification with multiple (>2) classes
%	multisvm = multisvmtrain_new(features, groups, type, varargin )
%
%   	features    = 
%   	groups      = 
% 	    type 		= 
%
%
%   Returns:
%       multisvm 	=
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
function multisvm = multisvmtrain_new(features, groups, type, varargin)
	
	% Check the type. This can be one-versus-all ('1vsall') or one-versus-one ('1vs1'). In a '1vsall' approach, there are as many
	% classifiers created as there are classes. Each classifier makes a distinction between one class and all the others. In the
	% '1vs1' approach, there are k over n classifiers created. Each classifier makes a distinction between one class and one other
	% class.
	if nargin < 3, type = '1vs1'; end
	
	
	% Extract the possible labels.
	labels = unique(groups);
	
	switch lower(type)
		case '1vsall'
            
			% Create as many classifiers as there are labels. This is for 1-vs-all classification.
			multisvm = cell(2,length(labels));


			for group = 1 : length(labels)
				% Train the classifier for the current class 'group' against all other points that do not belong to the 'group'. (This is
				% the 'one-versus-others' approach.
				multisvm{1, group} = fitcsvm( features, groups .* (groups == group), varargin{:} );
                
				
				% Add the multisvm type ('1vs1' or '1vsall'). This is used by the function multisvmclassify.
				multisvm{2, group} = type;
                
			end
			
			
		case '1vs1'
			% Get all possible combinations of binary classifiers using the presented groups.
			allCombinations = combnk( labels, 2 );


			% Empty classifier array.
			multisvm = cell(2,length(allCombinations));


			% Train the classifiers.
			for combination = 1 : length(allCombinations)
                
				% Select the two groups (0 and 1) for the current classifier to train on.
				groupsForThisCombination = NaN( size(groups) );
				groupsForThisCombination( groups == allCombinations(combination,1) ) = allCombinations(combination,1);
				groupsForThisCombination( groups == allCombinations(combination,2) ) = allCombinations(combination,2);

				% Train the classifier on the two selected groups.
				multisvm{1, combination} = fitcsvm( features, groupsForThisCombination, varargin{:} );

				% Add the multisvm type ('1vs1' or '1vsall'). This is used by the function multisvmclassify.
				multisvm{2, combination} = type;
                
                
			end
			
			
		otherwise
			error('The ''type'' should be ''1vs1'' or ''1vsall')
	end
end