% GenerateProbStabilityMap.m
% A. Cornelius
% Runs the physics algorithm if necessary

function [samples] = GenerateProbStabilityMap(samples, generator, param)
    global taskList

    % Check if the values have already been calculated
    if isfield(samples, 'CalculatedValues')
        return
    end

    % Start parpool
    gcp;

    % Calculate stability maps
    taskList.InitializeETA(length(samples(:)));
%     calculatedSamples(size(samples, 1), size(samples, 2), size(samples, 3)) = generator(samples(end, end, end).Values);
    parfor i = 1:length(samples(:))
        samples(i).CalculatedValues = generator(samples(i).Values);
        taskList.StepETA;
    end
    taskList.ClearETA;

    % Transfer values to the original sample array
%     [samples.CalculatedValues] = deal(calculatedSamples.CalculatedValues);
end