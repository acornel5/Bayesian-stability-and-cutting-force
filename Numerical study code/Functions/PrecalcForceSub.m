function rStruct = PrecalcForceSub(precalcForces, sample, param, options)
arguments
    precalcForces
    sample double
    param struct
    options.UnstableBendingMomentFactor double = 4
end

    % Get force coefficients from the sample
    % Right now it only considers ktc/knc
    cuttingForceKsBeta = sample(1:2);
    forceCoefs = [...
        cuttingForceKsBeta(1)*sin(cuttingForceKsBeta(2)), ... % ktc
        cuttingForceKsBeta(1)*cos(cuttingForceKsBeta(2)), ... % knc
        0, ... % kte
        0]; % kne

    % Force is assumed to be independent of spindle speed
    aGrid = param.CutParamGrid.a(1,:,:,:,:);
    bGrid = param.CutParamGrid.b(1,:,:,:,:);
    fzGrid = param.CutParamGrid.fz(1,:,:,:,:);

    % Initialize outputs
    stableBendingMoment = zeros(size(aGrid)); % Peak stable bending moment
    unstableBendingMoment = zeros(size(aGrid)); % Peak unstable bending moment
    meanTorque = zeros(size(aGrid)); % Mean cutting torque

    % This won't properly handle multiple feed directions
    for i = 1:numel(aGrid)        
        [stableBendingMoment(i), unstableBendingMoment(i), meanTorque(i)] = forceEstimateSub(...
            precalcForces{i}, ... % Precalculated forces
            forceCoefs, ... % Cutting force coefficients
            bGrid(i), ... % Axial depth
            fzGrid(i), ... % Feed per flute
            param.ToolGageLength, ...  % Tool gage length
            options.UnstableBendingMomentFactor);
    end

    rStruct.BendingMomentStable = stableBendingMoment;
    rStruct.BendingMomentUnstable = unstableBendingMoment;
    rStruct.MeanTorque = meanTorque;
end

function [bendingMomentStable, bendingMomentUnstable, meanTorque] = forceEstimateSub(precalcForces, forceCoefs, b, fz, gageLength, unstableBendingMomentFactor)
    ktc = forceCoefs(1);
    knc = forceCoefs(2);
    kte = forceCoefs(3);
    kne = forceCoefs(4);

    % Calculate instantaneous bending moments
    forceX = ...
        precalcForces.ktcX * ktc * fz + ...
        precalcForces.kncX * knc * fz + ...
        precalcForces.kteX * kte + ...
        precalcForces.kneX * kne;
    forceY = ...
        precalcForces.ktcY * ktc * fz + ...
        precalcForces.kncY * knc * fz + ...
        precalcForces.kteY * kte + ...
        precalcForces.kneY * kne;
    forceR = sqrt(forceX.^2 + forceY.^2);
    bendingMoments = forceR * (gageLength - b/2);

    % Output peak bending moments
    bendingMomentStable = max(bendingMoments);
    bendingMomentUnstable = bendingMomentStable * unstableBendingMomentFactor;

    % Calculate mean torque
    torques = ...
        precalcForces.ktcTorque * ktc * fz + ...
        precalcForces.kncTorque * knc * fz + ...
        precalcForces.kteTorque * kte + ...
        precalcForces.kneTorque * kne;
    meanTorque = mean(torques);
end