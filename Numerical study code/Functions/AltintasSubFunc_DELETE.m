function rStruct = AltintasSubFunc(inputVals, parser, previousCalcs, param, options)
    arguments
        inputVals double
        parser
        previousCalcs struct
        param struct
        options.RpmIncrement1D double = 1
        options.AssumeSymmetric logical = false
        options.ReduceSymmetricFRF logical = true
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
            rpms1D = min(param.Rpms):param.RpmIncrement1D:max(param.Rpms);
            outRpms = param.Rpms;
            outAps  = param.Aps;

            for outRpmInd = 1:length(param.Rpms)
                [~, inRpmInd] = min(abs(rpms1D - outRpms(outRpmInd)));
                apLim = stability1D(inRpmInd, radialWidthInd, 1, 1, cutDirInd);
                stableInds = outAps > apLim;
                stability(outRpmInd, radialWidthInd, stableInds, 1, cutDirInd) = 1;
			    chatterFreqs(outRpmInd,radialWidthInd,1,1,cutDirInd) = chatterFreqs1D(inRpmInd);
            end
        end
    end


    rStruct.Stability = stability;
    rStruct.ChatterFrequency = chatterFreqs;
    rStruct.Stability1D = stability1D;
end

