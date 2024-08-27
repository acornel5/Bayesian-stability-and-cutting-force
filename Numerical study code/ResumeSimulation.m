% ResumeSimulation.m
% Takes the output of the GenerateInitialSamples file and starts running
% cutting updates. Every iteration, it will suggest cutting parameters to
% test and then wait for the response

global taskList;
taskList.Clear();
taskList.Add('Performing Bayesian updating');
taskList.Add('Iteration 0');

%%
for iteration = iteration:param.MaxIterations
    taskList.Update(['Iteration ' num2str(iteration) ' / ' num2str(param.MaxIterations)]);
    taskList.Add('Starting iteration');
    
    % The first run-through with iteration=0 does not add a test point
    if iteration ~= 0	
        % Select test
        taskList.Update('Selecting test point');
        recommendedTest = testSelector(iteration, sampledVals, testResults);
        
        % Classify test
        taskList.Update('Classifying the test point');
        testResults(end+1) = testClassifier(iteration, recommendedTest);

		%% Draw new samples
		taskList.Update('Sampling new distribution');
        sampledVals = updateSampleGenerator(sampledVals, testResults(end), testResults);
    end

    %% Calculate the physics models for all the samples, if necessary
    sampledVals = GenerateProbStabilityMap(sampledVals, physicsDelegate);

    %% Save figures and data
    taskList.Update('Saving figures');
    figureGeneration(iteration, sampledVals, nominalVals, testResults);

    % Save data
    taskList.Update('Saving current data');    
    iteration = iteration+1; % Increment the index so that restarting from the saved file will start properly 
    if param.SaveCalculatedData
        save(['Figures\\' num2str(iteration)]);
    else % Remove the CalculatedValues field to save space.
        save(['Figures\\' num2str(iteration)], '-regexp', ['^(?!', 'sampledVals','$).']);
        sampledValsTemp = sampledVals;
        sampledVals = rmfield(sampledVals, 'CalculatedValues');
        save(['Figures\\' num2str(iteration)], 'sampledVals', '-append');
        sampledVals = sampledValsTemp;
        clear sampledValsTemp
    end

    taskList.Remove;
end

save 'Figures\\Endstate'