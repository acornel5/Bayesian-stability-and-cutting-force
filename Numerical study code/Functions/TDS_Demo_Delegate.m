% AltintasDelegate1DoF.m
% A. Cornelius
% 2021-7-27
% Maps a symmetric N DoF dynamic system and ks/beta pair to calculate a
% stability map and chatter frequencies. The dynamics inputs should be in
% the form of fn, k, and zeta values.

% Example input:
% [ks, beta, cr, fn1, k1, zeta1, fn2, k2, zeta2, ...]

function [calculatedSample] = TDS_Demo_Delegate(inputVals, param)
    % Resphape the input vals into the appropriate value arrays
    cuttingForce = inputVals(1:3);
    ks = cuttingForce(1);
    beta = cuttingForce(2);
    cr = cuttingForce(3);
    tdsMat = material(ks, beta);
    tdsMat.Cr = cr;

    tdsDynamics = param.DynamicModel;
    
%     ;inputVals(4:end);
%     dynamics = transpose(reshape(dynamics, 3, []));
%     tdsDynamics = dynamicModel.FromFnKZeta(dynamics, 'Symmetric');
   
    rpms = linspace(param.RpmMin, param.RpmMax, param.RpmSteps);
    aps  = linspace(param.ApMin,  param.ApMax,  param.ApSteps);


    tdsTool = param.tdsTool; %simpleEndmill(param.ToolDiameter, param.FluteCount, param.HelixAngle);

    calculatedSample.StabilityMap = logical(zeros(param.RpmSteps, param.ApSteps));
    for rpmInd = 1:param.RpmSteps
        for apInd = 1:param.ApSteps
%             disp(strcat(num2str((rpmInd-1) * param.ApSteps + apInd), ' / ', num2str(param.RpmSteps * param.ApSteps)))
            results = fastTimeDomainSimulation(...
                rpms(rpmInd), ...
                aps(apInd), ...
                param.RadialWidth, ...
                param.FeedDirection, ...
                param.FeedPerTooth, ...
                tdsTool, ...
                tdsMat, ...
                tdsDynamics, ...
                param);
    
            [calculatedSample.StabilityMap(rpmInd, apInd), ~] = classifier(results.xDisplacement, results.time, rpms(rpmInd), param.FluteCount, param.StabilityCutoff);
        end
    end
    
% 	calculatedSample.ChatterFrequencies = zeros(size(calculatedSample.StabilityMap));
end

function [stable, mmetric] = classifier(signal, times, rpm, flutes, cutoff)
    % Determine how many steps there are per flute pass
    flutePassingFreq = rpm / 60 * flutes;
    timeStep = times(2)-times(1);
    timePerFlutePass = 1 / flutePassingFreq;
    stepsPerFlutePass = timePerFlutePass / timeStep;
    
    % Extract once per rev indices
    oncePerRevIndices = round(1:stepsPerFlutePass:length(signal));
    oncePerRevSamples = signal(oncePerRevIndices);
    
    % Classify as stable
    mmetric = mean(abs(diff(oncePerRevSamples)));
    stable = ~(mmetric < cutoff); % Inverted because 0 is stable
end