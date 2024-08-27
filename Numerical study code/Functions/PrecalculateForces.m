% PrecalculateForces.m
% Calculates the mechanistic force model for unit values for ktc, knc, kte,
% kne. This allows the cutting forces for each sample to be calculated
% rapidly without needing to calculate the entire force model.
% The output is a struct with the following fields:
%  - dt: Time step for the simulation, at 1 rpm
%  - timesPerRpm: Time vector for the simulation, at 1 rpm
%  - freqsPerRpm: Frequency vector for the simulation, at 1 rpm
%  - fluteCuttingWall: Booleans for the time vector indicating at what
%    times there is a flute in contact with the wall
% The other fields are generally self-explanatory

function precalcVals = PrecalculateForces(param, options)
arguments
    param struct
    options.StepsPerRev = 2000
end

    % global taskList
    % taskList.Add('Precalculating forces');

    % Force is assumed to be independent of spindle speed
    aGrid = param.CutParamGrid.a(1,:,:,:,:);
    bGrid = param.CutParamGrid.b(1,:,:,:,:);
    cutDirGrid = param.CutParamGrid.cutDir(1,:,:,:,:);
    outputSize = size(aGrid);
    outputSize(1) = 1; % Compress along spindle speed direction

    % Initialize output grid
    precalcVals = cell(outputSize);
    % precalcVals{outputSize(1), outputSize(2), outputSize(3), outputSize(4), outputSize(5)} = 0;
    
    % Run precalculation    
    % taskList.InitializeETA(numel(precalcVals));
    for i = 1:numel(precalcVals)
        precalcVals{i} = precalcSub(...
            aGrid(i), ... % Radial width
            bGrid(i), ... % Axial depth
            1, ... % Feed per flute
            cutDirGrid(i), ... % Cut direction, currently not a variable
            param, ...
            options.StepsPerRev); % General parameters
        % taskList.StepETA;
    end
    % taskList.ClearETA;
    % taskList.Remove;
end

% Calculates the unit force values for a single set of cutting parameters
% forceFields is a struct cell array defining the cutting force for
% different values of a and b, as configured in param. Force is
% assumed to be independent of spindle speed. They should be accessed
% as forceFields(rpmInd, aInd, bInd, cutDirInd). Each struct has the following
% fields:
%     'dt'
%     'timesPerRpm'
%     'freqsPerRpm'
%     'kncX'
%     'kncY'
%     'kncTorque'
%     'kncXFFT'
%     'kncYFFT'
%     'kncTorqueFFT'
%     'ktcX'
%     'ktcY'
%     'ktcTorque'
%     'ktcXFFT'
%     'ktcYFFT'
%     'ktcTorqueFFT'
%     'kteX'
%     'kteY'
%     'kteTorque'
%     'kteXFFT'
%     'kteYFFT'
%     'kteTorqueFFT'
%     'kneX'
%     'kneY'
%     'kneTorque'
%     'kneXFFT'
%     'kneYFFt'
%     'kneTorqueFFT'
%     'fluteCuttingWall'
function precalcVals = precalcSub(a, b, fz, cutDir, param, stepsPerRev)
    arguments
        a double
        b double
        fz double
        cutDir int32
        param struct
        stepsPerRev int32
    end

    % Ktc
    [ktcX, ktcY, ktcZ, ktcTorque, origTimesPerRpm, fluteCuttingWall] = ForceSimSLE(1, b, a, fz, cutDir, [1, 0, 0, 0, 0, 0], param, StepsPerRev=stepsPerRev);
    ktcXFFT = trimmedFFT(ktcX);
    ktcYFFT = trimmedFFT(ktcY);
    ktcZFFT = trimmedFFT(ktcZ);
    ktcTorqueFFT = trimmedFFT(ktcTorque);

    % Knc
    [kncX, kncY, kncZ, kncTorque, ~, ~] = ForceSimSLE(1, b, a, fz, cutDir, [0, 1, 0, 0, 0, 0], param, StepsPerRev=stepsPerRev);
    kncXFFT = trimmedFFT(kncX);
    kncYFFT = trimmedFFT(kncY);
    kncZFFT = trimmedFFT(kncZ);
    kncTorqueFFT = trimmedFFT(kncTorque);

    % Kac
    [kacX, kacY, kacZ, kacTorque, ~, ~] = ForceSimSLE(1, b, a, fz, cutDir, [0, 0, 1, 0, 0, 0], param, StepsPerRev=stepsPerRev);
    kacXFFT = trimmedFFT(kacX);
    kacYFFT = trimmedFFT(kacY);
    kacZFFT = trimmedFFT(kacZ);
    kacTorqueFFT = trimmedFFT(kacTorque);

    % Kte
    [kteX, kteY, kteZ, kteTorque, ~, ~] = ForceSimSLE(1, b, a, fz, cutDir, [0, 0, 0, 1, 0, 0], param, StepsPerRev=stepsPerRev);
    kteXFFT = trimmedFFT(kteX);
    kteYFFT = trimmedFFT(kteY);
    kteZFFT = trimmedFFT(kteZ);
    kteTorqueFFT = trimmedFFT(kteTorque);

    % Kne
    [kneX, kneY, kneZ, kneTorque, ~, ~] = ForceSimSLE(1, b, a, fz, cutDir, [0, 0, 0, 0, 1, 0], param, StepsPerRev=stepsPerRev);
    kneXFFT = trimmedFFT(kneX);
    kneYFFT = trimmedFFT(kneY);
    kneZFFT = trimmedFFT(kneZ);
    kneTorqueFFT = trimmedFFT(kneTorque);

    % Kae
    [kaeX, kaeY, kaeZ, kaeTorque, ~, ~] = ForceSimSLE(1, b, a, fz, cutDir, [0, 0, 0, 0, 0, 1], param, StepsPerRev=stepsPerRev);
    kaeXFFT = trimmedFFT(kaeX);
    kaeYFFT = trimmedFFT(kaeY);
    kaeZFFT = trimmedFFT(kaeZ);
    kaeTorqueFFT = trimmedFFT(kaeTorque);

    % Frequency/time vectors
    dt = origTimesPerRpm(2) - origTimesPerRpm(1);
    freqsPerRpm = [0:1/dt/length(kncY):(1-1/(2*length(kncY)))*1/dt]';
    freqsPerRpm = freqsPerRpm(1:length(freqsPerRpm)/2+1);
    
    timesPerRpm = origTimesPerRpm(1:length(origTimesPerRpm)/2+1);
    timesPerRpm = timesPerRpm * 2;
    fluteCuttingWall = interp1(origTimesPerRpm, fluteCuttingWall, timesPerRpm, 'nearest', 'extrap');

    % Assemble output matrix
    precalcVals = struct(...
        'dt', dt, ...
        'timesPerRpm', timesPerRpm, ...
        'freqsPerRpm', freqsPerRpm, ...
        'ktcX', ktcX, ...
        'ktcY', ktcY, ...
        'ktcZ', ktcZ, ...
        'ktcTorque', ktcTorque, ...
        'ktcXFFT', ktcXFFT, ...
        'ktcYFFT', ktcYFFT, ...
        'ktcZFFT', ktcZFFT, ...
        'ktcTorqueFFT', ktcTorqueFFT, ...
        'kncX', kncX, ...
        'kncY', kncY, ...
        'kncZ', kncZ, ...
        'kncTorque', kncTorque, ...
        'kncXFFT', kncXFFT, ...
        'kncYFFT', kncYFFT, ...
        'kncZFFT', kncZFFT, ...
        'kncTorqueFFT', kncTorqueFFT, ...
        'kacX', kacX, ...
        'kacY', kacY, ...
        'kacZ', kacZ, ...
        'kacTorque', kacTorque, ...
        'kacXFFT', kacXFFT, ...
        'kacYFFT', kacYFFT, ...
        'kacZFFT', kacZFFT, ...
        'kacTorqueFFT', kacTorqueFFT, ...
        'kteX', kteX, ...
        'kteY', kteY, ...
        'kteZ', kteZ, ...
        'kteTorque', kteTorque, ...
        'kteXFFT', kteXFFT, ...
        'kteYFFT', kteYFFT, ...
        'kteZFFT', kteZFFT, ...
        'kteTorqueFFT', kteTorqueFFT, ...
        'kneX', kneX, ...
        'kneY', kneY, ...
        'kneZ', kneZ, ...
        'kneTorque', kneTorque, ...
        'kneXFFT', kneXFFT, ...
        'kneYFFT', kneYFFT, ...
        'kneZFFT', kneZFFT, ...
        'kneTorqueFFT', kneTorqueFFT, ...
        'kaeX', kaeX, ...
        'kaeY', kaeY, ...
        'kaeZ', kaeZ, ...
        'kaeTorque', kaeTorque, ...
        'kaeXFFT', kaeXFFT, ...
        'kaeYFFT', kaeYFFT, ...
        'kaeZFFT', kaeZFFT, ...
        'kaeTorqueFFT', kaeTorqueFFT, ...
        'fluteCuttingWall', fluteCuttingWall);
end

function rVal = trimmedFFT(signal)
    N = length(signal);
    mu = mean(signal);
    signal = signal - mu;
    rVal = fft(signal);%/(N/2);
    rVal = rVal(1:length(signal)/2+1);
end


function [FFT, frequencies] = signalFFT(signal, samplingRate)
N = length(signal);
mu = mean(signal);
signal = signal - mu;
rawFFT = fft(signal) / (N/2); % Calculate FFT and correct amplitude
warning('off')
FFT = rawFFT(1:N/2+1);
warning('on')
FFT(1) = mu; % Replace DC value with mean

frequencies = [0:samplingRate/N:(1-1/(2*N))*samplingRate]';
warning('off')
frequencies = frequencies(1:N/2+1);                 % frequency, Hz
warning('on')
end