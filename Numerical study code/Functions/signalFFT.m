% signalFFT.m
% A. Cornelius
% 2020-7-6
% This function calculates the fourier transform frequency content of an 
% input signal

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