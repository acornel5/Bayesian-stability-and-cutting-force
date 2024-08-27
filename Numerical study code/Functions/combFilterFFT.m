% combFilterFFT.m
% A. Cornelius
% 2020-7-7
% This function sets the amplitude for an FFT to 0 near harmonics of some
% specified frequency. The width of the filter region is set as a
% percentage of the input frequency

function combFilteredData = combFilterFFT(fft, frequencies, flutePassingFreq, tolerance)
    % Calculate the tolerance index offset
    resolution = frequencies(2)-frequencies(1);
    toleranceOffset = round(flutePassingFreq / resolution * tolerance);

    % Determine which indexes should be filtered
    indexes = [];
    nextFrequency = 0;
    while nextFrequency < max(frequencies)
        [~, baseIndex] = min(abs(frequencies - nextFrequency));
        indexes = [indexes, (baseIndex-toleranceOffset):(baseIndex+toleranceOffset)];
        nextFrequency = nextFrequency + flutePassingFreq;
    end
    
    % Get rid of any indexes which are out of the range
    indexes(indexes < 1) = [];
    indexes(indexes > length(frequencies)) = [];
    
    % Set values to zero and returns the filtered FFT
    fft(indexes) = 0;
    combFilteredData = fft;
end