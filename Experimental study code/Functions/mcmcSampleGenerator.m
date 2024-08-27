function [newSamples, diagnostics] = mcmcSampleGenerator(priorSamples, newTestResults, allTestResults, priorDist, pdfFormula, NameValueArgs, options)
    arguments
        priorSamples struct
        newTestResults struct
        allTestResults struct
        priorDist
        pdfFormula function_handle
        NameValueArgs.SampleCount double
        options.SamplesPerThread int32 = 5
        options.Beta double = 0
        options.MinimumSampleCount = 100
    end

    retainedSamples = rejectionSample(priorSamples, @(sample) pdfFormula(sample, [], newTestResults), options.MinimumSampleCount);

    mapProbabilityFormula = @(sample) pdfFormula(sample, priorDist, allTestResults);
    [newSamples, chainInds, acceptRate] = ...
        mcmc_Parallel( ...
        retainedSamples, ...
        mapProbabilityFormula, ...
        NameValueArgs.SampleCount);

    diagnostics = [];
end


function retainedSamples = rejectionSample(priorSamples, pdfFunc, minimumSampleCount)
    for i = 1:length(priorSamples)
        priorSamples(i) = pdfFunc(priorSamples(i));
        disp(strcat(num2str(i), "/", num2str(length(priorSamples))))
    end

    % Scale the probabilities to always accept some minimum number of
    % samples
    for i = 1:length(priorSamples)
        allWeights(i) = priorSamples(i).LogProbability;
    end

    % Find nth highest probability
    sortedWeights = sort(allWeights, 2, "descend");
    nthWeight = sortedWeights(minimumSampleCount);
    allWeights = allWeights - nthWeight;
    for i = 1:length(priorSamples)
        priorSamples(i).LogProbability = allWeights(i);
    end

    % Filter old samples
    for i = 1:length(priorSamples)
        prob = exp(priorSamples(i).LogProbability);
        pass = 0;
        for j = 1:priorSamples(i).Weight
            if rand < prob
                pass = pass + 1;
            end
        end
        priorSamples(i).Weight = pass;
    end
    retainedSamples = priorSamples([priorSamples.Weight] ~= 0);

    if length(retainedSamples) < 2
        error('ERROR: Number of remaining samples is insufficient.')
    end

end

function [newSamples, diagnostics] = mcmc(priorSamples, pdfFunc, options)

end