% GetSampledField.m
% 2022-6-22
% Cornelius

% Extracts a specific field from the sampled values, expanding from the
% weighted values to a full grid.

function [outputVals] = GetSampledField(sampledVals, fieldName)
    % Get field size
    d = size(sampledVals(1).CalculatedValues.(fieldName));
    dims = ones(length(d), 1);
    dims(end+1) = 1;

    % Create the output
    outputVals = [];
    for i = 1:length(sampledVals)
        if isempty(sampledVals(i).Weight), continue, end
        dims(end) = sampledVals(i).Weight;
        expandedMatrix = repmat(sampledVals(i).CalculatedValues.(fieldName), dims');
        outputVals = cat(length(d)+1, outputVals, expandedMatrix);
    end


    % Resize so that the first dimension is non-unity
    while size(outputVals, 1) == 1 && any(size(outputVals)~=1)
        currentSize = size(outputVals);
        newSize = [currentSize(2:end), 1];
        outputVals = reshape(outputVals, newSize);
    end
end