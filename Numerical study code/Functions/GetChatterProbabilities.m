function probs = GetChatterProbabilities(sampledVals)
    probs = zeros(size(sampledVals(1).CalculatedValues.Stability));
    weight = 0;
    for i = 1:length(sampledVals)
        weight = weight + sampledVals(i).Weight;
        probs = probs + (sampledVals(i).CalculatedValues.Stability * sampledVals(i).Weight);
    end

    probs = probs / weight;
%     stabilityMaps = GetSampledField(sampledVals, 'StabilityMap');

%     probs = sum(stabilityMaps, 4) / length(stabilityMaps);
end