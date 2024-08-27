% Returns a grid of logical values representing whether a test is
% acceptable or not according to the constraints

function acceptableTests = testConstraints(samples, param, options)
    arguments
        samples
        param
        options.ApRange = [0 realmax]
        options.ChatterProbability = [0 1];
        options.BreakageProbability = [0 1];
        options.SpindleBendingMomentProbability = [0 1];
        options.SpindleBendingMomentLimit = 375;
        options.SpindleTorqueProbability = [0 1];
    end

    acceptableTests = true(size(param.CutParamGrid.rpm));

    % Axial depth limit
    acceptableTests(~InRange(param.CutParamGrid.b, options.ApRange(1), options.ApRange(2))) = false;

    % Chatter probability
    if isfield(samples(1).CalculatedValues, 'Stability')
        chatterProb = CalculateAverageLikelihood(samples, 'Stability');
        chatterProb = repmat(chatterProb, [1 1 1 length(param.FeedPerTooth)]);
        acceptableTests(~InRange(chatterProb, options.ChatterProbability(1), options.ChatterProbability(2))) = false;
    end
    % Stable breakage likelihood
    if isfield(samples(1).CalculatedValues, 'StableBreakageLikelihood')
        breakageProbability = CalculateAverageLikelihood(samples, 'StableBreakageLikelihood');
        acceptableTests(~InRange(breakageProbability, options.BreakageProbability(1), options.BreakageProbability(2))) = false;
    end
    % Bending moment
    if isfield(samples(1).CalculatedValues, 'BendingMomentStable')
        weightSum = 0;
        P = zeros(size(samples(1).CalculatedValues.Stability));
        for i = 1:length(samples)
            P = P + samples(i).Weight .* (1-samples(i).CalculatedValues.Stability) .* samples(i).CalculatedValues.BendingMomentStable > options.SpindleBendingMomentLimit;
            P = P + samples(i).Weight .* samples(i).CalculatedValues.Stability .* samples(i).CalculatedValues.BendingMomentUnstable > options.SpindleBendingMomentLimit;
            weightSum = weightSum + samples(i).Weight;
        end
        P = P/weightSum;
        acceptableTests(~InRange(P, options.SpindleBendingMomentProbability(1), options.SpindleBendingMomentProbability(2))) = false;
    end

    % Spindle torque
    if isfield(samples(1).CalculatedValues, 'MeanTorque')
        torqueLimit = getTorqueLimit(param.CutParamGrid.rpm);
        weightSum = 0;
        P = zeros(size(samples(1).CalculatedValues.Stability));
        for i = 1:length(samples)
            P = P + samples(i).Weight * (samples(i).CalculatedValues.MeanTorque > torqueLimit);
            weightSum = weightSum + samples(i).Weight;
        end
        P = P/weightSum;
        acceptableTests(~InRange(P, options.SpindleTorqueProbability(1), options.SpindleTorqueProbability(2))) = false;
    end
end


        
function P = CalculateAverageLikelihood(samples, fieldName)
    weightSum = 0;
    P = zeros(size(samples(1).CalculatedValues.(fieldName)));

    for i = 1:length(samples)
        P = P + samples(i).Weight * samples(i).CalculatedValues.(fieldName);
        weightSum = weightSum + samples(i).Weight;
    end

    P = P / weightSum;
end

function P = CalculateToleranceLikelihood(samples, fieldName, maxVal)
    weightSum = 0;
    P = zeros(size(samples(1).CalculatedValues.(fieldName)));

    for i = 1:length(samples)
        P = P + samples(i).Weight * samples(i).CalculatedValues.(fieldName)<maxVal;
        weightSum = weightSum + samples(i).Weight;
    end

    P = P / weightSum;
end

% Check if a nubmer is in a specified range, always inclusive of both ends
function vals = InRange(input, min, max)
    vals = input >= min & input <= max;
end


function torque = getTorqueLimit(rpms)
    torque = min(47*ones(size(rpms)), 60*130e3/2/pi ./ rpms);
end