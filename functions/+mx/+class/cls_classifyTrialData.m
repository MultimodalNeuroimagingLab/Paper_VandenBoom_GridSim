%
%   Classify a set of trials (with training and test data)
%	[predictedLabels, predictedCorrect, outTrainData, outTestData] = cls_classifyTrialData(config, trainData, trainConditions, testData, testConditions)
%
%   	config           = the configuration structure for the classification, with the following fields:
%   						  classificationMethod         		 = the method ['svm' or 'templatematching']
%   						  lateFeatureSelectionMethod   		 = additional feature selection/reduction ['none' or 'ANOVA']
%   						  classificationMethodTemplateMetric = only used with template matching, matching by correlation ['corr'] or sum of squares distance ['ssd']
%   	trainData       = matrix with the data to be used for training. The first
%           	          dimension(rows) should be the features and the second 
%               	      dimension (columns) should be the trials
% 	    trainConditions = vector with the conditions for the trianing trials
%   	testData        = matrix with the data to be classified (tested). The first
%               	      dimension(rows) should be the features and the second 
%                   	  dimension (columns) should be the trials
%   	testConditions  = vector with the conditions for the testing trials
%       	              (used to return the predicted correct vector)
%
%
%   Returns:
%       predictedLabels 	=
%		predictedCorrect 	= 
%		outTrainData		=
%		outTestData			=
%
%
%   Copyright 2019, Max van den Boom

%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function [predictedLabels, predictedCorrect, outTrainData, outTestData] = cls_classifyTrialData(config, trainData, trainConditions, testData, testConditions)

     % late feature selection based on trainingset
    if strcmpi(config.lateFeatureSelectionMethod, 'ANOVA') == 1

        % use ANOVA on all the features to determine the significant
        % features (voxels)
        featureIndices = mx.class.cls_selectFeaturesAnova_Epoch(config, trainData, trainConditions);

        if length(featureIndices) > 0
        
            % only use the features
            trainData = trainData(featureIndices, :);
            testData = testData(featureIndices, :);
    
        end
        
    end

    % select the classification method
    if strcmpi(config.classificationMethod, 'svm') == 1
        
        %{
        % recode condition (to ensure they are 0 based)
        recoded = unique(trainConditions);
        reTrainCondition = nan(1, length(trainConditions))
        reTestCondition = nan(1, length(testConditions))
        for r = 1:length(recoded)
            reTrainCondition(trainConditions == recoded(r)) = r - 1;
            reTestCondition(testConditions == recoded(r)) = r - 1;
        end
        %}
        
        % train the SVM model
        %SVMModel = mx.class.multisvmtrain(trainTrials', (trainTrialsConditions + 1)', '1vsall', 'autoscale', true, 'kernel_function', 'linear');
        %SVMModel = mx.class.multisvmtrain_new(trainData', (reTrainCondition + 1)', '1vsall', 'Standardize', true, 'KernelFunction', 'linear');      % updated to fitcsvm to work in Matlab 2018
        SVMModel = mx.class.multisvmtrain_new(trainData', (trainConditions + 1)', '1vsall', 'Standardize', true, 'KernelFunction', 'linear');      % updated to fitcsvm to work in Matlab 2018

        % test 
        %[predictedLabels, ~] = mx.class.multisvmclassify(SVMModel, testTrial');
        [predictedLabels, ~] = mx.class.multisvmclassify_new(SVMModel, testData');         % updated to fitcsvm to work in Matlab 2018   
        %predictedLabels = recoded(predictedLabels);
        predictedLabels = (predictedLabels - 1)';

        
        
        
    elseif strcmpi(config.classificationMethod, 'templatematching') == 1
        
        % create templates for all conditions and store in the results
        templates = nan(size(trainData, 1), max(trainConditions) + 1);
        for c = 0:max(trainConditions)
            templates(:, c + 1) = mean(trainData(:, (trainConditions == c)), 2);
        end
        
        % determine the matching method
        if strcmpi(config.classificationMethodTemplateMetric, 'corr') == 1
            % correlation
            
            % calculate the correlations between the templates and the test trials
            scores = corr(testData, templates);

            % make all correlations positieve
            scores = abs(scores);

            % take the highest correlation as the winner
            [~, predictedLabels] = max(scores');
            
        elseif strcmpi(config.classificationMethodTemplateMetric, 'ssd') == 1
            % sum squares distance (ssd)
            
            % calculate the SSD between the templates and the test trial
            scores = [];
            for cTrial = 1:size(testData, 2)
                for cTemplate = 1:size(templates, 2)
                   scores(cTrial, cTemplate) = mx.class.cls_ssd(testData(:, cTrial), templates(:, cTemplate));
                end                
            end
            
            % take the lowest ssd as the winner
            [~, predictedLabels] = min(scores);
            
        else

            % message
            fprintf(2, ['Unknown template classification metric: ', config.classificationMethodTemplateMetric, '\n']);
            return;

        end
        

        predictedLabels = predictedLabels - 1;
        
        
        % classify within run using template matching (leave on out)
        %[result] = mx.class.cls_classifyWithinUsingTemplates(signal, signalIndices, indicesConditions, indicesDurations, 'corr');
        
    else
        
        % message
        fprintf(2, ['Unknown classification method: ', config.classificationMethod, '\n']);
        return;
        
    end
    
    % test whether the predicted labels were correct to the given conditions
    predictedCorrect = predictedLabels == testConditions;

    % output the train and test data (might have changed after late feature selection)
    outTrainData = trainData;
    outTestData = testData;

end
