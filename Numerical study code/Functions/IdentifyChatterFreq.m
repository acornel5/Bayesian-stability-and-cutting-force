function [chatterFrequencies, fig] = IdentifyChatterFreq(tdsResults, field, rpm, ap, stable, param, maxFreq, peakCount)
    tolerance = 0.01;

    signal = tdsResults.(field);
    flutePassingFreq = rpm * param.FluteCount / 60;
    
    [fft, fftFreqs] = signalFFT(signal, tdsResults.sampleFrequency);
    filteredFFT = combFilterFFT(fft, fftFreqs, flutePassingFreq, 0.01);
    
    
    
    title([num2str(rpm) ' rpm, ' num2str(ap*1000) ' mm'])
    % fig = plotFFTFigure(fft, filteredFFT, fftFreqs, rpm, ap, flutePassingFreq, stable);
    % fig.CurrentAxes.XLim = [0 maxFreq];

    chatterFrequencies = zeros(peakCount, 1);
    if ~stable
        for i = 1:peakCount
            [~, chatterIndex] = max(abs(filteredFFT));
            chatterFrequencies(i) = fftFreqs(chatterIndex);
            filteredFFT(round((1-tolerance)*chatterIndex):round((1+tolerance)*chatterIndex)) = 0;
        end
    end
end