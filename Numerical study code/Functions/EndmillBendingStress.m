% SLE.m
% Cornelius A.
% 2022-5-2

% Calculates the maximum bending stress in a stable cut based on the tool 
% deflection, using the frequency domain method to estimate deflection and
% Euler-Bernoulli beam theory to estimate the stres

function rStruct = EndmillBendingStress(inputVals, previousCalcs, parser, param, precalcData, NameValueArgs, options)
    arguments
        inputVals double
        previousCalcs struct
        parser
        param struct
        precalcData
        NameValueArgs.BeamDiameter
        NameValueArgs.BeamLength
        NameValueArgs.UltimateTensileStrengthMean
        options.StructuralDamping double = 0.046
        options.ElasticModulus double = 600e9
        options.Density = 7800;
        options.StressConcentrationFactor double = 2.135
        options.UltimateTensileStrengthStdDev = 0
        options.UseFrequencyDomainApproximation = true
    end

    % Extract the ultimate tensile strength from the sample vector if
    % needed
    if isstring(NameValueArgs.UltimateTensileStrengthMean)
        NameValueArgs.UltimateTensileStrengthMean = parser.returnField(inputVals, NameValueArgs.UltimateTensileStrengthMean);
    end
    if isstring(options.UltimateTensileStrengthStdDev)
        options.UltimateTensileStrengthStdDev = parser.returnField(inputVals, options.UltimateTensileStrengthStdDev);
    end

    % Calculate the root stress for a unit tool-tip displacement
    lambda = 1.87510407;
    sigma = 0.734095514;
    rawMaxDisplacement = -(cosh(lambda) - cos(lambda) - sigma*(sinh(lambda) - sin(lambda)));
    rawMaxMoment = lambda^2/NameValueArgs.BeamLength^2 * 2;
    I = pi/4 * (NameValueArgs.BeamDiameter/2)^4;
    unitMaxMoment = rawMaxMoment / max(abs(rawMaxDisplacement)) * options.ElasticModulus * I;
    unitMaxStress = unitMaxMoment / I * NameValueArgs.BeamDiameter / 2 * options.StressConcentrationFactor;

    f = previousCalcs.FrequenciesFRF;
    frfxx = previousCalcs.FRFxx;
    frfyy = previousCalcs.FRFyy;
    f = reshape(f, numel(f), 1);
    frfxx = reshape(frfxx, numel(frfxx), 1);
    frfyy = reshape(frfyy, numel(frfyy), 1);


    % Resphape the input vals into the appropriate value arrays
    cuttingForce = parser.forceKtcKnc(inputVals);
    ktc = cuttingForce(1);
    knc = cuttingForce(2);
    kte = 0;
    kne = 0;


    bendingStress = NaN(size(param.CutParamGrid.a));


    for rpmInd = 1:length(param.Rpms)
        % Get the FRF for the relevant frequencies
        rpm = param.Rpms(rpmInd);
        freqs = precalcData{1,1,1}.freqsPerRpm * rpm;                
        frfxxInterp = reshape(qinterp1(f, frfxx, freqs, 0), size(precalcData{1,1,1}.ktcXFFT));
        frfyyInterp = reshape(qinterp1(f, frfyy, freqs, 0), size(precalcData{1,1,1}.ktcXFFT));
        frfxxInterp(isnan(frfxxInterp)) = 0;
        frfyyInterp(isnan(frfyyInterp)) = 0;

        for aInd = 1:length(param.RadialWidth)
            for bInd = 1:length(param.Aps)
                a = param.RadialWidth(aInd);
                b = param.Aps(bInd);
                fz = param.FeedPerTooth;
                cutDir = param.CutDirection;
        
                % Save time by omitting 0 depths
                if b==0
                    bendingStress(rpmInd, aInd, bInd) = 0;
                    continue
                end
                
                % Frequency-domain forces
                ktcFz = ktc * fz;
                kncFz = knc * fz;
                precalcValues = precalcData{1, aInd, bInd};
                F_x = precalcValues.ktcXFFT .* ktcFz + precalcValues.kncXFFT .* kncFz;% + precalcValues.kteXFFT * kte + precalcValues.kneXFFT * kne;
                F_y = precalcValues.ktcYFFT .* ktcFz + precalcValues.kncYFFT .* kncFz;% + precalcValues.kteYFFT * kte + precalcValues.kneYFFT * kne;

                % Frequency-domain response
                X_f = F_x .* frfxxInterp;
                Y_f = F_y .* frfyyInterp;

                % Calculate the maximum deflection. The frequency domain
                % approximation assumes that there is only one dominant
                % vibrational frequency: this is generally close to
                % accurate for the tools this is applied for, but does have
                % some approximation error.
                if options.UseFrequencyDomainApproximation
                    maxDeflection = max(sqrt(abs(X_f).^2 + abs(Y_f).^2)) / (length(F_x)-1);
                else
                    x_t = real(ifft(X_f));
                    y_t = real(ifft(Y_f));
                    maxDeflection = max(sqrt(x_t.^2 + y_t.^2));
                end

                maxStress = maxDeflection * unitMaxStress;
                bendingStress(rpmInd, aInd, bInd) = maxStress;
            end
        end
    end    

    rStruct.StableBendingStress = bendingStress;
    rStruct.StableBreakageLikelihood = normcdf(bendingStress, NameValueArgs.UltimateTensileStrengthMean, options.UltimateTensileStrengthStdDev);
end



function Yi = qinterp1(x,Y,xi,methodflag)
% Performs fast interpolation compared to interp1
%
% qinterp1 provides a speedup over interp1 but requires an evenly spaced
% x array.  As x and y increase in length, the run-time for interp1 increases
% linearly, but the run-time for
% qinterp1 stays constant.  For small-length x, y, and xi, qinterp1 runs about
% 6x faster than interp1.
%
%
% Usage:
%   yi = qinterp1(x,Y,xi)  - Same usage as interp1
%   yi = qinterp1(x,Y,xi,flag)
%           flag = 0       - Nearest-neighbor
%           flag = 1       - Linear (default)
%
% Example:
%   x = [-5:0.01:5];   y = exp(-x.^2/2);
%   xi = [-4.23:0.3:4.56];
%   yi = qinterp1(x,y,xi,1);
%
% Usage restrictions
%    x must be monotonically and evenly increasing
%    e.g.,  x=-36:0.02:123;
%
%    Y may be up to two-dimensional
%
% Using with non-evenly spaced arrays:
%   Frequently the user will wish to make interpolations "on the fly" from
%   a fixed pair of library (i.e., x and y) vectors.  In this case, the
%   user can generate an equally-spaced set of library data by calling
%   interp1 once, and then storing this library data in a MAT-file or
%   equivalent.  Because the speed of qinterp1 is independent of the length
%   of the library vectors, the author recommends over-sampling this
%   generated set untill memory considerations start limitting program speed.
%
%   If the user wishes to use two or more spacings (i.e., a closely-spaced
%   library in the region of fine features, and a loosely-spaced library in
%   the region of coarse features), just create multiple libraries, record
%   the switching points, and send the search data to different qinterp1
%   calls depending on its value.
%
%   Example:
%       x1 = [-5:0.01:5];   x2 = [-40:1:-5 5:1:40];
%       y1 = exp(-x1.^2/3); y2 = exp(-x2.^2/3);
%       xi = [-30:0.3:30];
%       in = xi < 5 & xi > -5;
%       yi(in) = qinterp1(x1,y1,xi(in));
%       yi(~in) = qinterp1(x2,y2,xi(~in));
% Author: N. Brahms
% Copyright 2006
% Forces vectors to be columns
x = x(:); xi = xi(:);
sx = size(x); sY = size(Y);
if sx(1)~=sY(1)
    if sx(1)==sY(2)
        Y = Y';
    else
        error('x and Y must have the same number of rows');
    end
end
if nargin>=4
    method=methodflag;
else
    method = 1;    % choose nearest-lower-neighbor, linear, etc.
                   % uses integer over string for speed
end
% Gets the x spacing
ndx = 1/(x(2)-x(1)); % one over to perform divide only once
xi = xi - x(1);      % subtract minimum of x
% Fills Yi with NaNs
s = size(Y);
if length(s)>2
    error('Y may only be one- or two-dimensional');
end
Yi = NaN*ones(length(xi),s(2));
switch method
    case 0 %nearest-neighbor method
        rxi = round(xi*ndx)+1;        % indices of nearest-neighbors
        flag = rxi<1 | rxi>length(x) | isnan(xi);
                                      % finds indices out of bounds
        nflag = ~flag;                % finds indices in bounds
        Yi(nflag,:) = Y(rxi(nflag),:);
    case 1 %linear interpolation method
        fxi = floor(xi*ndx)+1;          % indices of nearest-lower-neighbors
        flag = fxi<1 | fxi>length(x)-1 | isnan(xi);
                                        % finds indices out of bounds
        nflag = ~flag;                  % finds indices in bounds
        Yi(nflag,:) = (fxi(nflag)-xi(nflag)*ndx).*Y(fxi(nflag),:)+...
            (1-fxi(nflag)+xi(nflag)*ndx).*Y(fxi(nflag)+1,:);
end
end