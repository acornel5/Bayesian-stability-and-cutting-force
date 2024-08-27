% Calculates the stability map and chatter frequency using the zero-order
% approximation frequency domain stability algorithm. Rather than returning
% a binary stability map, this function returns a grid of stability
% likelihoods based on the Gaussian CDF stability likelihood function.

% Required previous calculations:
%  - FrequenciesFRF: the frequencies that the FRF are sampled at, in Hz
%  - FRFxx: the X-axis FRF, in m/N
%  - FRFyy: the Y-axis FRF, only used if AssumeSymmetric is not set to true

% Required named arguments:
%  - StabilityUncertainty: standard deviation for the Gaussian CDF
%    likelihood calculation. If this is set to 'inParam', then it will look
%    for the 'param.StabilityUncertainty' field instead.

% Possible outputs:
%  - Stability: the stability likelihoods for the grid
%  - ChatterFrequency: the chatter frequencies for every spindle speed
%  - Stability1D: The blim values for spindle speeds, calculated at the
%    resolution RpmIncrement1D

% Options:
% - RpmIncrement1D: Specifies the higher rpm-resolution used for
%   calculating the stability map and saving to the Stability1D output.
%   Defaults to 1 rpm.
% - AssumeSymmetric: Assume that the X and Y dynamics are symmetrical
%   (i.e., FRFxx = FRFyy). Defaults to false.
% - ReduceSymmetricFRF: Approximate the FRF using a smaller number of 
%   points to reduce calculation time. Only works if AssumeSymmetric is on.

function rStruct = physicsFunZOA_StabilityLikelihood(inputVals, parser, previousCalcs, param, namedArgs, options)
    arguments
        inputVals double
        parser
        previousCalcs struct
        param struct
        namedArgs.StabilityUncertainty;
        options.RpmIncrement1D double = 1 % The ZOA will run calculations at a higher rpm-resolution than the stability grid and can optionally save 
        options.AssumeSymmetric logical = false %
        options.ReduceSymmetricFRF logical = true % Approximate the FRF using a smaller number of points to reduce calculation time. Only works for 
    end

    if namedArgs.StabilityUncertainty == 'inParam'
        namedArgs.StabilityUncertainty = param.StabilityUncertainty;
    end

    cuttingForce = parser.forceKsBeta(inputVals);
    cuttingForce(2) = rad2deg(cuttingForce(2));

    f = previousCalcs.FrequenciesFRF;
    w = f * 2 * pi;
    if options.AssumeSymmetric
        frfxx = previousCalcs.FRFxx;
        if options.ReduceSymmetricFRF, [w, frfxx] = dynamicResampleFRF(w, frfxx); end
        frfyy = frfxx;
    else
        frfxx = previousCalcs.FRFxx;
        frfyy = previousCalcs.FRFyy;
    end
    
    % Preallocate the outputs
    lengthOfOneDVec = length(min(param.Rpms):options.RpmIncrement1D:max(param.Rpms));
    stability1D = zeros(lengthOfOneDVec, length(param.RadialWidth), 1, 1, length(param.CutDirection));
    stability = single(zeros(length(param.Rpms), length(param.RadialWidth), length(param.Aps), 1, length(param.CutDirection)));
    chatterFreqs = zeros(length(param.Rpms), length(param.RadialWidth), 1, 1, length(param.CutDirection));

    % The cutting parameter grid is in the order: (rpm, radialWidth,
    % axialDepth, feedPerFlute, cutDirection
    for cutDirInd = 1:length(param.CutDirection)    
        for radialWidthInd = 1:length(param.RadialWidth)
            [stability1D(:,radialWidthInd,1,1,cutDirInd), chatterFreqs1D] = AltintasSLD(...
                w', frfxx, frfyy, ... % Dynamics
                cuttingForce, ... % Force model
                param.Rpms, param.Aps, param.RadialWidth(radialWidthInd), param.CutDirection(cutDirInd), ... % Cutting parameters
                param, ... % Tool info etc
                RpmIncrement1D = options.RpmIncrement1D); % Options

            % Convert from the 1D stability map to a 2d stable/unstable map
            stability(:,radialWidthInd,:,1,cutDirInd) = zeros(length(param.Rpms), length(param.Aps));
            rpms1D = min(param.Rpms):options.RpmIncrement1D:max(param.Rpms);
            outRpms = param.Rpms;
            outAps  = param.Aps;

            for outRpmInd = 1:length(param.Rpms)
                [~, inRpmInd] = min(abs(rpms1D - outRpms(outRpmInd)));
                apLim = stability1D(inRpmInd, radialWidthInd, 1, 1, cutDirInd);
                stability(outRpmInd, radialWidthInd, :, 1, cutDirInd) = normcdf(outAps, apLim, namedArgs.StabilityUncertainty);
			    chatterFreqs(outRpmInd,radialWidthInd,1,1,cutDirInd) = chatterFreqs1D(inRpmInd);
            end
        end
    end


    rStruct.Stability = stability;
    rStruct.ChatterFrequency = chatterFreqs;
    rStruct.Stability1D = stability1D;
end

