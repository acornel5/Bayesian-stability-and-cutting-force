% mcmcWithChatterFrequency.m
% A. Cornelius
% 2021-8-31
% This is a flexible function for calculating stability maps via Metropolis
% Hastings. It uses a list of standard deviations and a starting point to
% generate candidate samples. Those samples are then passed to an input 
% delegate function to generate the stability map. That delegate allows
% this function to be generalized over different simulations.
% 
% This method uses the parallel MCMC framework described in "A general
% construction for parallelizing Metropolis−Hastings algorithms", applied
% to Haario's adaptive Metropolis algorithm.


function [rSamples, rChain, acceptanceRate] = mcmc_Parallel( ...
	startingSamples, ...
	probabilityAlgorithm, ...
    SampleCount, ...
    options)
    arguments
        startingSamples struct
        probabilityAlgorithm
        SampleCount int32
        options.SamplesPerThread = 5
        options.Beta = 0
    end

    % Initialize the matfile
    if isfile('Figures\\MCMC_Diagnostics.mat'), delete('Figures\\MCMC_Diagnostics.mat'), end
    m = matfile('Figures\\MCMC_Diagnostics.mat', 'Writable', true);
    
	% Update task list
    global taskList;

	% Define the different functions used to calculate the probability of the sample being accepted
	d = length(startingSamples(1).Values);
    fullSampleMatrix = [];
    for i = 1:length(startingSamples)
        samples = repmat(startingSamples(i).Values, startingSamples(i).Weight, 1);
        fullSampleMatrix = [fullSampleMatrix; samples];
    end
    covMatrix = cov(fullSampleMatrix);
		
	% allSamples is initialized to the maximum expected length to minimize time spent copying variables if the length is exceeded.
    rChain = 1;
    
	% Initialize parallel calculation
	threadCount = str2double(getenv('NUMBER_OF_PROCESSORS'));
	N = threadCount * options.SamplesPerThread;
	    
    % Initialize the output array
    rSamples(1) = probabilityAlgorithm(startingSamples(1).Values);
    rSamples(1).Weight = 0;
    lastChainInd = 1;
		
    % Start parpool
    gcp;
    
    taskList.InitializeCustomETA(ETA_FunctionMCMC(SampleCount));
    
    % Save starting memory
    [userMem, systemMem] = memory;
    m.Memory(1,1) = struct('UserMem', userMem, 'SystemMem', systemMem);
    m.MemoryVars(1,1) = struct('Variables', whos);

    m.iterationTotalTime(1,1) = 0;
    m.iterationTimeBreakdown(1,1) = struct('Times', 0);
    while length(rSamples) < SampleCount
        iterationTic = tic;
		iterationTimeBreakdown = toc(iterationTic);
		% Generate N samples in parallel
        currentValues = rSamples(lastChainInd).Values;
		candidates = (1-options.Beta) .* mvnrnd(currentValues, 2.38^2 * covMatrix / d, N) + options.Beta .* mvnrnd(currentValues, 0.01*eye(d)./d.*currentValues, N);
		iterationTimeBreakdown(end+1) = toc(iterationTic);

        parfor i = 1:N
			candSamples(i) = probabilityAlgorithm(candidates(i,:));
        end
		iterationTimeBreakdown(end+1) = toc(iterationTic);

        % Calculate and normalize acceptance probabilities. The probability
        % of staying at the last point is the LAST element in the list
        candAcceptProbs = [candSamples.LogProbability];
		candAcceptProbs(end+1) = rSamples(lastChainInd).LogProbability; 
		candAcceptProbs = candAcceptProbs - max(candAcceptProbs); % Normalize the exponential acceptance probabilities so that they are centered about 0
        candAcceptProbs = exp(candAcceptProbs); % Since these are normalized, we can shift back to non-exponential space with minimal risk of exceeding double capacity
        iterationTimeBreakdown(end+1) = toc(iterationTic);

		% Draw samples
		inds = randsample(N+1, N, true, candAcceptProbs);
        newChainInds = inds;
        iterationTimeBreakdown(end+1) = toc(iterationTic);
        
        % Add weight to the most recent sample
        rSamples(lastChainInd).Weight = rSamples(lastChainInd).Weight + nnz(inds == N+1);
        newChainInds(newChainInds == N+1) = lastChainInd;
        iterationTimeBreakdown(end+1) = toc(iterationTic);
        
        % Add new elements to the chain
        for i = 1:N
            if ~any(inds == i), continue, end
            tempSample = candSamples(i);
            tempSample.Weight = nnz(inds == i);
%             candSamples(i).Weight = nnz(inds == i);
            
%             if candSamples(i).Weight ~= 0
            samples = repmat(tempSample.Values, tempSample.Weight, 1);
            fullSampleMatrix = [fullSampleMatrix; samples];
            rSamples(end+1) = tempSample;
													  
            newChainInds(newChainInds == i) = length(rSamples);
%             end
        end
        iterationTimeBreakdown(end+1) = toc(iterationTic);
        
        rChain = [rChain; newChainInds];
        lastChainInd = newChainInds(end);
        iterationTimeBreakdown(end+1) = toc(iterationTic);

        % Update covariance matrix
        covMatrix = cov(fullSampleMatrix);
        iterationTimeBreakdown(end+1) = toc(iterationTic);

        % Update progress bar
        taskList.StepETA([length(rSamples), length(rChain)]);
		
        % Clear candidates
%         candSamples = rmfield(candSamples, 'Weight');
%         clear candSamples

        % Save diagnostics
        [userMem, systemMem] = memory;
        m.Memory(end+1,1) = struct('UserMem', userMem, 'SystemMem', systemMem);
        % m.MemoryVars(end+1,:) = struct('Variables', whos);
        m.FullSampleMatrix = fullSampleMatrix;
        m.LastChainInd = lastChainInd;
        m.rChain = rChain;
        m.lastChainInd = lastChainInd;
        m.iterationTotalTime(end+1,1) = toc(iterationTic);
        iterationTimeBreakdown(end+1) = toc(iterationTic);
        m.iterationTimeBreakdown(end+1,1) = struct('Times', iterationTimeBreakdown);
    end

    acceptanceRate = length(rSamples) / length(rChain);
    taskList.Remove;

    m.Completed = true;
end

function S = CalcCovMatrix(startingSamples, currentSamples)
    fullMatrix = [];
    for i = 1:length(startingSamples)
        samples = repmat(startingSamples(i).Values, startingSamples(i).Weight, 1);
        fullMatrix = [fullMatrix; samples];
    end
    for i = 1:length(currentSamples)
        samples = repmat(currentSamples(i).Values, currentSamples(i).Weight, 1);
        fullMatrix = [fullMatrix; samples];
    end
    S = cov(fullMatrix);
end