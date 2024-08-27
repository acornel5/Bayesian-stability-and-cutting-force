% acceptanceFormulaChatterFreq.m

% Calculates the probability for a specific sample. The sample must be
% input as an array with the correct size for the prior


% This returns an anonymous function which calculates the probability for a
% given sample. The input/output for that sample is:

% Input: a struct with the following fields
%  - .Values: the set of randomly sampled values to calculate the
%    probability for
%  - .CalculatedValues: an OPTIONAL input with all the pre-calculated
%    physics model information calculated by the delegate for this sample

% Output:
%  - .Values: the set of randomly sampled values for this sample
%  - .LogProbability: the log of the PDF for this sample. Note that this is
%    not necissarily normalized.
%  - .Diagnostics: Extraneous probabilities for debugging purposes
%     - .LogPriorProb: 
%     - .LogStableProbs:
%     - .LogChatterFreqProbs: 
%  - .CalculatedValues: the calculated physics model info for this sample
%  - Any other fields present on the original input will also be passed to
%    the output (e.g. Weight for the MCMC chain).

function sample = acceptanceFormulaLikelihood(sample, prior, testedPoints, physicsAlgorithm, param, namedArgs, options)
    arguments
        sample
        prior
        testedPoints
        physicsAlgorithm
        param
        namedArgs.StabilityUncertainty
        namedArgs.chatterFreqStdDev
        namedArgs.breakageMean
        namedArgs.breakageStdDev
        options.LearnFromUnbrokenCuts logical = false
        options.PowerStdDev = realmax;
        options.LearnPowerFromUnstableCuts = true;
    end

    if strcmp(namedArgs.StabilityUncertainty, 'inParam')
        namedArgs.StabilityUncertainty = param.StabilityUncertainty;
    end

    % Check if it conforms to empirical tests
    priorProb = 0;
    stableProbs = zeros(size(testedPoints));
    chatterFreqProbs = zeros(size(testedPoints));
    brokenProbs = zeros(size(testedPoints));
    PowerProbs = zeros(size(testedPoints));

    if ~isstruct(sample)
        sample = struct('Values', sample);
    end
    % Check if the provided sample has already been calculated
    if ~isfield(sample, 'CalculatedValues')
        sample.CalculatedValues = physicsAlgorithm(sample.Values);
    end
    
    % Calculate the probability on the initial prior distribution
    % If the probability on the prior is 0 (e.g. on a uniform
    % distribution), the calculation can be aborted here.
    if ~isempty(prior)
        priorProb = log(prior.PDF(sample.Values));
        if priorProb == -Inf
			sample.LogProbability = -Inf;
			sample.Diagnostics.LogPriorProb = -Inf;
			sample.Diagnostics.LogStableProbs = 0;
			sample.Diagnostics.LogChatterFreqProbs = 0;
            return
        end
    end

    % If the stability likelihoods aren't precalculated, then they'll be
    % calculated from the 1d stability vector
    if ~strcmp(namedArgs.StabilityUncertainty, 'precalculated')
        rpmVec1D = linspace(min(param.Rpms), max(param.Rpms), size(sample.CalculatedValues.Stability1D, 1));
    end

    for i = 1:length(testedPoints)
        cutParamInd = testedPoints(i).CutParamInd;

        % Stability likelihood
        if strcmp(namedArgs.StabilityUncertainty, 'precalculated')
            stabilityMap = repmat(sample.CalculatedValues.Stability, [1 1 1 length(param.FeedPerTooth) 1]);
            stabilityLikelihood = stabilityMap(cutParamInd);
            stableProbs(i) = log(inlineIfElse(testedPoints(i).Stable, 1-stabilityLikelihood, stabilityLikelihood));
        else
            [~, rpmInd] = min(abs(rpmVec1D - testedPoints(i).Rpm));
            [~, aInd, ~, ~, cutDirInd] = ind2sub(size(sample.CalculatedValues.Stability),cutParamInd);
            blim = sample.CalculatedValues.Stability1D(rpmInd, aInd, 1, 1, cutDirInd);
            stabilityLikelihood = 1 - normcdf(testedPoints(i).b, blim, namedArgs.StabilityUncertainty);
            stableProbs(i) = log(inlineIfElse(testedPoints(i).Stable, 1-stabilityLikelihood, stabilityLikelihood));
        end

        % Chatter frequency likelihood
        if testedPoints(i).ChatterFrequency ~= 0
            predictedChatterFrequency = getPartialField(sample.CalculatedValues.ChatterFrequency, cutParamInd, size(param.CutParamGrid.a));
            actualChatterFrequency = testedPoints(i).ChatterFrequency;
            chatterFreqProbs(i) = log(normpdf(predictedChatterFrequency, actualChatterFrequency, namedArgs.chatterFreqStdDev));
        else % Chatter is not observed
            chatterFreqProbs(i) = 0;
        end

        % Power
        if (testedPoints(i).Stable || options.LearnPowerFromUnstableCuts) && testedPoints(i).Power ~= 0 && options.PowerStdDev ~= realmax
            if options.PowerStdDev < 0 % Negative values use the uncertainty as a percentage of the prediction
                powerUncertainty = sample.CalculatedValues.SpindlePower(cutParamInd) * (-options.PowerStdDev);
            else
                powerUncertainty = options.PowerStdDev;
            end
            PowerProbs(i) = log(normpdf(testedPoints(i).Power, sample.CalculatedValues.SpindlePower(cutParamInd), powerUncertainty));
        else
        end

        if options.LearnFromUnbrokenCuts
            if ~testedPoints(i).Broken
                brokenProbs(i) = log(1-sample.CalculatedValues.StableBreakageLikelihood(cutParamInd));
            end
        end
    end
	
    sample.LogProbability = priorProb + sum(stableProbs) + sum(chatterFreqProbs) + sum(brokenProbs) + sum(PowerProbs);
    sample.Diagnostics.PriorProb = priorProb;
	sample.Diagnostics.StableProbs = stableProbs;
	sample.Diagnostics.ChatterFreqProbs = chatterFreqProbs;
    sample.Diagnostics.BrokenProbs = brokenProbs;
    sample.Diagnostics.PowerProbs = PowerProbs;
end