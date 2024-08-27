% manualTestClassifier.m
% Cornelius A.

% Prompts the user to classify the test as stable/unstable using whatever
% metric they prefer.

function test = tdsTestClassifier( ...
    iteration, ...
    test, ...
    param, ...
    dynamics, ...
    forceModel)

    tdsDynamics = dynamicModel.FromFnKZeta(dynamics, 'Symmetric');
    % tdsDynamics = dynamicModel.GetRigidModel();
    tdsMat = material(forceModel(1), forceModel(2));
    tdsMat.Kte = 10000;

    tdsTool = simpleEndmill(param.ToolDiameter, param.FluteCount, param.HelixAngle);

    rpm = test.RPM;
    ap = test.b;
    ae = test.a;
    fz = test.fz;
    cutDirection = test.CutDirection;

    results = fastTimeDomainSimulation(...
        rpm, ...
        ap, ...
        ae, ...
        cutDirection, ...
        fz, ...
        tdsTool, ...
        tdsMat, ...
        tdsDynamics, ...
        param);

    [test.Stable, mMetric, oncePerRevTimes, oncePerRevSamples] = classifier(results.xDisplacement, results.time, rpm, param.FluteCount, param.StabilityCutoff);
    [chatterFreq, freqs, stableFFT, chatterFFT, chatterFreqMags] = IdentifyChatterFreq(results, 'xDisplacement', rpm, ap, test.Stable, param, 5000, 1);
    % Plot the frequency-domain data
    f = figure(1);
    clf
    hold on
    plot(freqs, abs(chatterFFT)*1e6, 'r')
    plot(freqs, abs(stableFFT)*1e6, 'b')
    % stem(freqs, abs(stableFFT)*1e6, 'b', MarkerEdgeColor='none', MarkerFaceColor='none')
    if ~test.Stable
        plot(chatterFreq, chatterFreqMags*1e6, 'or')
    end
    xlabel('Frequency (Hz)')
    ylabel('Magnitude ($\mu$m)')
    xlim([0 3000])
    legend(["Chatter data" "Flute passing freqs" strcat("Chatter freq: ", num2str(round(chatterFreq)), " Hz")], Location="northeast")
    folderName = strcat('Figures\\TestClassificationFreqDomain');
    if not(isfolder(folderName)), mkdir(folderName), end
    fileName = [num2str(iteration)];

    print('-dpng', strcat(folderName, '\\', fileName), '-r600');
    close(f);

    if ~test.Stable
        test.ChatterFrequency = chatterFreq;
    else
        test.ChatterFrequency = 0;
    end
    torque = mean(results.torque);
    test.Power = torque * rpm * 2 * pi/ (60);

    % Plot the time-domain data
    f = figure(1);
    clf
    hold on
    plot(results.time-results.time(1), results.xDisplacement*1e6);
    plot(oncePerRevTimes-results.time(1), oncePerRevSamples*1e6, 'or', MarkerSize=4);
    title(strcat("n=", num2str(round(rpm)), " rpm, b=", num2str(ap*1000,4), " mm, fz=", num2str(fz*1000,4), " mm"))
    xlabel('Time (s)')
    ylabel('X displacement ($\mu$m)')
    legend(["X displacement", strcat("m=", num2str(mMetric*1e6,3), " $\mu$m")])
    folderName = strcat('Figures\\TestClassificationTimeDomain');
    if not(isfolder(folderName)), mkdir(folderName), end
    fileName = [num2str(iteration)];

    print('-dpng', strcat(folderName, '\\', fileName), '-r600');
    close(f)
end

function [stable, mmetric, oncePerRevTimes, oncePerRevSamples] = classifier(signal, times, rpm, flutes, cutoff)
    % Determine how many steps there are per flute pass
    flutePassingFreq = rpm / 60 * flutes;
    timeStep = times(2)-times(1);
    timePerFlutePass = 1 / flutePassingFreq;
    stepsPerFlutePass = timePerFlutePass / timeStep;
    
    % Extract once per rev indices
    oncePerRevIndices = round(1:stepsPerFlutePass:length(signal));
    oncePerRevSamples = signal(oncePerRevIndices);
    oncePerRevTimes = times(oncePerRevIndices);
    
    % Classify as stable
    mmetric = mean(abs(diff(oncePerRevSamples)));
    stable = (mmetric < cutoff); % Inverted because 0 is stable
end

function [chatterFrequencies, freqs, stableFFT, chatterFFT, chatterFreqMags] = IdentifyChatterFreq(tdsResults, field, rpm, ap, stable, param, maxFreq, peakCount)
    tolerance = 0.01;

    signal = tdsResults.(field);
    flutePassingFreq = rpm / 60; % Includes runout freqs %rpm * param.FluteCount / 60;
    
    freq = 1/(tdsResults.time(2)-tdsResults.time(1));

    [fft, fftFreqs] = signalFFT(signal, tdsResults.sampleFrequency);
    filteredFFT = combFilterFFT(fft, fftFreqs, flutePassingFreq, 0.01);
    
    freqs = fftFreqs;
    chatterFFT = filteredFFT;
    stableFFT = fft - chatterFFT;

    chatterFrequencies = zeros(peakCount, 1);
    chatterFreqMags = zeros(peakCount, 1);
    if ~stable
        for i = 1:peakCount
            [~, chatterIndex] = max(abs(filteredFFT));
            chatterFrequencies(i) = fftFreqs(chatterIndex);
            chatterFreqMags(i) = abs(filteredFFT(chatterIndex));
            filteredFFT(round((1-tolerance)*chatterIndex):round((1+tolerance)*chatterIndex)) = 0;
        end
    end
end