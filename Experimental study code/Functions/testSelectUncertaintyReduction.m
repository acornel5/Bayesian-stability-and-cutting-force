% Required named arguments:
%  - ChatterFreqUncertainty: What standard deviation is assigned to the chatter frequency likelihood

% Optional arguments:
%  - VariableIndex: Which variable to optimize for, specified by the index in the btheta array. Mutually exclusive with VariableName. If neither is set, the user is prompted to select one.
%  - VariableName: Which variable to optimize for, specified by the name.
%  - SampleCountToUse: How many samples should be used to calculate the uncertainty reduction. Cycle time increases with the square of this value. If not specified, all samples are used.


function [recommendedTest, extraMetrics] = testSelectUncertaintyReduction(iteration, sampledVals, testResults, acceptableTests, param, sampleParser, NameValueArgs, options)
    arguments
        iteration
        sampledVals
        testResults
        acceptableTests
        param
        sampleParser
        NameValueArgs.ChatterFreqUncertainty
        options.VariableIndex = ''
        options.VariableName = '';
        options.SampleCountToUse = 0;
        options.OverrideConstraints = false;
    end

    % Trim the list of sampled values to the correct length
    if options.SampleCountToUse == 0, options.SampleCountToUse = length(sampledVals); end
    sampleInds = randi(length(sampledVals), options.SampleCountToUse, 1);
    sampledVals = sampledVals(sampleInds);
	
	
    if ~xor(isempty(options.VariableIndex), isempty(options.VariableName))
        disp('The variable name for uncertainty test selection was not specified in the SimulationParameters file.')
        disp('The valid values are: ')
        for i = 1:length(sampleParser.fieldNames)
            disp(strcat(num2str(i), ": ", sampleParser.fieldNames(i)))
        end
        options.VariableIndex = input("What variable do you want to use? ");
        options.VariableName = sampleParser.fieldNames(options.VariableIndex);
    end
    if isempty(options.VariableIndex)
        options.VariableIndex = sampleParser.GetFieldIndex(options.VariableName);
    end
    if isempty(options.VariableName)
        options.VariableName = sampleParser.fieldNames{options.VariableIndex};
    end
	
    startingWeights = zeros(size(sampledVals));
    vars = zeros(size(sampledVals));
    for i = 1:length(sampledVals)
        startingWeights(i) = sampledVals(i).Weight;
        vars(i) = sampledVals(i).Values(options.VariableIndex);
    end

    chatterProbs = GetChatterProbabilities(sampledVals);
    % chatterProbs = repmat(chatterProbs, [1 1 1 length(param.FeedPerTooth) 1]);
    chatterLikelihoods = extractCalculatedStructField(sampledVals, 'Stability', [0, 0, 0, 0, 0], WeightedVals=false);
    % chatterLikelihoods = repmat(chatterLikelihoods, [1 1 1 length(param.FeedPerTooth) 1]);

    allChatterFrequencies = extractCalculatedStructField(sampledVals, 'ChatterFrequency', [0, 0, 0, 0, 0], WeightedVals=false);

    startingUncertainty = std(vars, startingWeights);


    % Initialize output    
    global taskList
    taskList.InitializeETA(length(chatterProbs(:)));
    expectedUncertainties = zeros(size(chatterProbs));
    for i = 1:length(expectedUncertainties(:))
        % Find the index of the test point and calculate the weights based
        % on likelihood of stability
        [inds(1), inds(2), inds(3), inds(4), inds(5)] = ind2sub(size(chatterProbs), i);
        stabilityInds = sub2ind(size(chatterLikelihoods), inds(1)*ones(size(sampledVals)), inds(2)*ones(size(sampledVals)), inds(3)*ones(size(sampledVals)), inds(4)*ones(size(sampledVals)), inds(5)*ones(size(sampledVals)), 1:length(sampledVals));
        updatedStableWeights = startingWeights .* squeeze(1 - chatterLikelihoods(stabilityInds));
        updatedUnstableWeights = startingWeights - updatedStableWeights;

         % If all samples predict that a cut will be stable, then you won't learn anything from the uncertainty 
        if all(updatedUnstableWeights == 0)
            expectedUncertainties(i) = startingUncertainty;
            continue
        end

        % If the chatter frequency learning is turned off
        if NameValueArgs.ChatterFreqUncertainty == realmax
            stableUncertainty = std(vars, updatedStableWeights);
            unstableUncertainty = std(vars, updatedUnstableWeights);        
            if isnan(stableUncertainty), stableUncertainty = startingUncertainty; end
            if isnan(unstableUncertainty), unstableUncertainty = startingUncertainty; end
            expectedUncertainties(i) = chatterProbs(i) * unstableUncertainty + (1-chatterProbs(i)) * stableUncertainty;
            continue
        end

        % Calculate the uncertainty reduction due to chatter frequency
        chatterFreqInds = min(inds, size(allChatterFrequencies, 1:5));
        chatterFreqInds = sub2ind(size(allChatterFrequencies), chatterFreqInds(1)*ones(size(sampledVals)), chatterFreqInds(2)*ones(size(sampledVals)), chatterFreqInds(3)*ones(size(sampledVals)), chatterFreqInds(4)*ones(size(sampledVals)), chatterFreqInds(5)*ones(size(sampledVals)), 1:length(sampledVals));
        predictedChatterFreqs = allChatterFrequencies(chatterFreqInds);
        weights = normpdf(predictedChatterFreqs, predictedChatterFreqs', NameValueArgs.ChatterFreqUncertainty) .* updatedUnstableWeights';
        chatterFreqUncertainties = weightedStd(vars', weights);
        unstableUncertainty = sum(chatterFreqUncertainties.*startingWeights)/sum(startingWeights);

        % Calculate the updated uncertainties
        stableUncertainty = std(vars, updatedStableWeights);
        if isnan(stableUncertainty), stableUncertainty = startingUncertainty; end
        if isnan(unstableUncertainty), unstableUncertainty = startingUncertainty; end
        expectedUncertainties(i) = chatterProbs(i) * unstableUncertainty + (1-chatterProbs(i)) * stableUncertainty;
    
        taskList.StepETA;
    end
    taskList.ClearETA;
	
    expectedReductions = ((expectedUncertainties - startingUncertainty)./startingUncertainty*100);

    disp(strcat("Max reduction: ", num2str(max(abs(expectedReductions), [], 'all'),"%.2f"), "%"))
    % Apply test constraints
	if ~options.OverrideConstraints
		selectionReductions = expectedReductions;
		selectionReductions(~acceptableTests) = NaN;
    else
        selectionReductions = expectedReductions;
    end

	minimizeDepth = input("Minimize axial depth? 0 for no, 1 for yes: ");
	if minimizeDepth == 1
		desiredLearningAmount = input("Desired reduction in uncertainty in %? ");
	end

	% Select test
	if minimizeDepth == 1 % Select the test with minimum axial depth which achieves the desired learning
		viableInds = selectionReductions <= -desiredLearningAmount;
		viableTestDepths = param.CutParamGrid.b;
		viableTestDepths(~viableInds) = NaN;
		[~, bestInd] = min(viableTestDepths(:));
	else % Select the test with the highest reduction in uncertainty
		[~, bestInd] = min(selectionReductions(:));
	end
	recommendedTest.RPM = param.CutParamGrid.rpm(bestInd);
	recommendedTest.b = param.CutParamGrid.b(bestInd);
	recommendedTest.a = param.CutParamGrid.a(bestInd);
	recommendedTest.fz = param.CutParamGrid.fz(bestInd);
	recommendedTest.CutDirection = param.CutParamGrid.cutDir(bestInd);
	recommendedTest.MRR = recommendedTest.RPM * recommendedTest.b * recommendedTest.a * recommendedTest.fz * param.FluteCount;
	recommendedTest.CutParamInd = bestInd;

        % Create the plot grid
    [~, aInd] = min(abs(param.RadialWidth - recommendedTest.a));
    [~, fzInd] = min(abs(param.FeedPerTooth - recommendedTest.fz));
    [~, cutDirInd] = min(abs(param.CutDirection - recommendedTest.CutDirection));
    rpmMatrix = squeeze(param.CutParamGrid.rpm(:,1,:,1,1));
    apMatrix = squeeze(param.CutParamGrid.b(:,1,:,1,1))*1000;
    expectedReductions = squeeze((expectedUncertainties(:,aInd,:,fzInd,cutDirInd) - startingUncertainty)./startingUncertainty*100);


    % Plot the uncertainty reduction
    figure(1)
    clf
    hold on
    surf(rpmMatrix, apMatrix, zeros(size(expectedReductions)), expectedReductions)
    contour(rpmMatrix, apMatrix, squeeze(acceptableTests(:,aInd,:,fzInd,cutDirInd)), [0.1 0.1], 'r')
    contour(rpmMatrix, apMatrix, squeeze(acceptableTests(:,aInd,:,fzInd,cutDirInd)), [0.9 0.9], 'b')
    view(2)

    plot3(recommendedTest.RPM, recommendedTest.b*1000, 200, "diamond", MarkerEdgeColor=[0.9290 0.6940 0.1250], MarkerFaceColor=[0.9290 0.6940 0.1250], MarkerSize=10)
	
        % Plot tested points
    for i = 1:length(testResults)
        if abs((testResults(i).a - recommendedTest.a)/recommendedTest.a)>0.01, continue, end
        if testResults(i).Broken == 1
            scatter3(testResults(i).RPM, testResults(i).b*1000, 2, 100, '*', 'MarkerEdgeColor', 'red', 'LineWidth', 2); % Chatter points
        elseif testResults(i).Stable == 1
            scatter3(testResults(i).RPM, testResults(i).b*1000, 2, 100, 'filled', 'MarkerEdgeColor', [100 200 0]/255, 'MarkerFaceColor', [100 200 0]/255);
        elseif testResults(i).Stable == 0
            scatter3(testResults(i).RPM, testResults(i).b*1000, 2, 100, 'x', 'MarkerEdgeColor', 'red', 'LineWidth', 2); % Chatter points
        elseif testResults(i).Stable == 0.5
            scatter3(testResults(i).RPM, testResults(i).b*1000, 2, 100, 'filled', 'MarkerEdgeColor', [255 200 0]/255, 'MarkerFaceColor', [255 200 0]/255);
        else
            error('Invalid stability classification')
        end
    end

    shading interp
    xlabel('Spindle speed (rpm)')
    ylabel('Axial depth (mm)')
    c = colorbar();
    set(get(c,'label'),'string', strcat("Expected \% improvement in ", options.VariableName, " uncertainty"));
    set(get(c,'label'),'interpreter','latex');
    c.TickLabelInterpreter = 'latex';
    colormap(flip(parula))
    title(strcat("a=", num2str(recommendedTest.a*1000,4), ...
        " mm, fz=", num2str(recommendedTest.fz*1000,4), ...
        " mm, ", inlineIfElse(recommendedTest.CutDirection==2, "climb", "conventional"), " milling"))


    extraMetrics{1} = strcat("Expected reduction in uncertainty: ", num2str((startingUncertainty - expectedUncertainties(bestInd))/startingUncertainty*100), "%");
end

% Calculate the weighted standard deviation for a matrix of weights with 
% different sets of weights. The stadnadr deviation will always be taken in the first (column) direction. The vars must be input as a column vector or matrix, and
% the weights must be a matrix with the same height
function stds =  weightedStd(vars, weights)
    weightSum = sum(weights, 1);
    weightedMean = sum(vars .* weights, 1) ./ weightSum;
    stds = (sum(weights .* (vars - weightedMean).^2, 1) ./ weightSum).^0.5;
end

function u = predictChatterUncertainty(predictedChatterFreqs, actualChatterFreq, chatterFreqStdDev, unstableWeights, variableValues)
    weights = normpdf(predictedChatterFreqs, actualChatterFreq, chatterFreqStdDev);
    u = std(variableValues, weights .* unstableWeights);
end


function u = expectedUncertainty(i, vars, chatterLikelihoods, chatterProbs, startingWeights, chatterFreqUncertainty, allChatterFrequencies, startingUncertainty)
    [inds(1), inds(2), inds(3), inds(4), inds(5)] = ind2sub(size(chatterProbs), i);
    stabilityInds = sub2ind(size(chatterLikelihoods), inds(1)*ones(size(vars)), inds(2)*ones(size(vars)), inds(3)*ones(size(vars)), inds(4)*ones(size(vars)), inds(5)*ones(size(vars)), 1:length(vars));
    updatedStableWeights = startingWeights .* squeeze(1 - chatterLikelihoods(stabilityInds));
    updatedUnstableWeights = 1 - updatedStableWeights;

     % If all samples predict that a cut will be stable, then you won't learn anything 
    if all(updatedUnstableWeights == 0)
        u = startingUncertainty;
        return
    end

    % If the chatter frequency learning is turned off
    if chatterFreqUncertainty == realmax
        stableUncertainty = std(vars, updatedStableWeights);
        unstableUncertainty = std(vars, updatedUnstableWeights);        
        if isnan(stableUncertainty), stableUncertainty = startingUncertainty; end
        if isnan(unstableUncertainty), unstableUncertainty = startingUncertainty; end
        u = chatterProbs(i) * unstableUncertainty + (1-chatterProbs(i)) * stableUncertainty;
        return
    end

    chatterFreqInds = min(inds, size(allChatterFrequencies, 1:5));
    chatterFreqInds = sub2ind(size(allChatterFrequencies), chatterFreqInds(1)*ones(size(vars)), chatterFreqInds(2)*ones(size(vars)), chatterFreqInds(3)*ones(size(vars)), chatterFreqInds(4)*ones(size(vars)), chatterFreqInds(5)*ones(size(vars)), 1:length(vars));
    predictedChatterFreqs = allChatterFrequencies(chatterFreqInds);
    weights = normpdf(predictedChatterFreqs, predictedChatterFreqs', chatterFreqUncertainty) .* updatedUnstableWeights';
    chatterFreqUncertainties = weightedStd(vars', weights);
    unstableUncertainty = sum(chatterFreqUncertainties.*startingWeights)/sum(startingWeights);


    stableUncertainty = std(vars, updatedStableWeights);
    % if isnan(stableUncertainty), stableUncertainty = startingUncertainty; end
    % if isnan(unstableUncertainty), unstableUncertainty = startingUncertainty; end
    u = chatterProbs(i) * unstableUncertainty + (1-chatterProbs(i)) * stableUncertainty;
end