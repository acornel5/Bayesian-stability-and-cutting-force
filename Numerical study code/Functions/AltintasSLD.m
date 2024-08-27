%  AltintasSLD.m
%  A. Cornelius
%  2021-7-1
%  Calculates the stability map and chatter frequencies for a given process, cutting force
%  coefficients, and X/Y frequency response functions.

function [SLD, chatterFreqs] = AltintasSLD(w, FRFxx, FRFyy, cuttingForce, rpms, aps, radialWidth, cutDirection, tool, options)
arguments
    w double
    FRFxx
    FRFyy
    cuttingForce double
    rpms double
    aps double
    radialWidth double
    cutDirection double
    tool struct
    options.RpmIncrement1D double = 1
end

% maxMagnitude = max([abs(FRFxx), abs(FRFyy)], [], 'all');
% acceptedIndices = (abs(FRFxx) > 0.01*maxMagnitude) | (abs(FRFyy) > 0.01*maxMagnitude);
% w = w(acceptedIndices);
% FRFxx = FRFxx(acceptedIndices);
% FRFyy = FRFyy(acceptedIndices);

% Define stability map limits
omega_vector = min(rpms):options.RpmIncrement1D:max(rpms);       % spindle speed, rpm
blimmax = max(aps)*100;                    % maximum axial depth of cut, m
maxLobeNumber = ceil(w(end)/2/pi*60/(min(rpms)*tool.FluteCount))+1;
minLobeNumber = max(1, floor(w(1)/2/pi*60/(max(rpms)*tool.FluteCount))-1);

% Define tool
d = tool.ToolDiameter;                        % tool diameter, m
Nt = tool.FluteCount;                             % number of teeth

% Define cut
rd = radialWidth;                            % radial depth of cut, m
ud = cutDirection;                             % 1 = up milling, 2 = down milling
if ud == 1                          % up milling
    phistart = 0;                   % starting angle, deg
    phiexit = acos((d/2 - rd)/(d/2))*180/pi; % exit angle, deg
else                                % down milling
    phistart = 180 - acos((d/2 - rd)/(d/2))*180/pi;
    phiexit = 180;
end
phis = phistart*pi/180;             % rad
phie = phiexit*pi/180;              % rad

% Load cutting force model from input
Ks = cuttingForce(1);
beta = cuttingForce(2);
Kn = 1/tan(beta*pi/180);
Kt = Ks/sqrt(1 + Kn^2);             % tangential force coefficient, N/m^2

% Orient frequency response function
alphaxx = 0.5*((cos(2*phie)-2*Kn*phie+Kn*sin(2*phie))-(cos(2*phis)-2*Kn*phis+Kn*sin(2*phis)));
alphaxy = 0.5*((-sin(2*phie)-2*phie+Kn*cos(2*phie))-(-sin(2*phis)-2*phis+Kn*cos(2*phis)));
alphayx = 0.5*((-sin(2*phie)+2*phie+Kn*cos(2*phie))-(-sin(2*phis)+2*phis+Kn*cos(2*phis)));
alphayy = 0.5*((-cos(2*phie)-2*Kn*phie-Kn*sin(2*phie))-(-cos(2*phis)-2*Kn*phis-Kn*sin(2*phis)));
FRF_or = ...
    [alphaxx*reshape(FRFxx, [1 1 length(w)]) alphaxy*reshape(FRFyy, [1 1 length(w)]); ...
    alphayx*reshape(FRFxx, [1 1 length(w)]) alphayy*reshape(FRFyy, [1 1 length(w)])];    % m/N

%Calculate eigenvalues
D = eig2(FRF_or);
lambda2 = D(2,:)';
lambda1 = D(1,:)';

% Calculate blim vectors
blim1 = (2*pi/Nt/Kt)./((real(lambda1)).^2 + (imag(lambda1)).^2) .* (real(lambda1) .* (1 + (imag(lambda1)./real(lambda1)).^2));  % m
blim2 = (2*pi/Nt/Kt)./((real(lambda2)).^2 + (imag(lambda2)).^2) .* (real(lambda2) .* (1 + (imag(lambda2)./real(lambda2)).^2));
[index1] = find(blim1 > 0);
blim1 = blim1(index1);
w1 = w(index1);
psi1 = atan2(imag(lambda1), real(lambda1));
psi1 = psi1(index1);
[index2] = find(blim2 > 0);
blim2 = blim2(index2);
w2 = w(index2);
psi2 = atan2(imag(lambda2), real(lambda2));
psi2 = psi2(index2);

% Initialize stability map 
blimvec = blimmax * ones(size(omega_vector));

% Stability map from lambda1
index = find(blim1 < blimmax);
blimMap1 = [];
chatterFreqMap1 = [];
if ~isempty(index)
    blim1 = blim1(index);
    w1 = w1(index);
    psi1 = psi1(index);
    epsilon1 = pi - 2*psi1;
    
	blimMap1 = realmax*ones(length(omega_vector), maxLobeNumber - minLobeNumber);
	chatterFreqMap1 = zeros(length(omega_vector), maxLobeNumber - minLobeNumber);

    for N = minLobeNumber:(maxLobeNumber - 1)
		[blimMap1(:, N), chatterFreqMap1(:, N)] = altintasMappingSub(Nt, w1, N, epsilon1, blim1, omega_vector, options.RpmIncrement1D);
    end
end

% Stability map from lambda2
index = find(blim2 < blimmax);
blimMap2 = [];
chatterFreqMap2 = [];
if ~isempty(index)
    blim2 = blim2(index);
    w2 = w2(index);
    psi2 = psi2(index);
    epsilon2 = pi - 2*psi2;

	blimMap2 = realmax*ones(length(omega_vector), maxLobeNumber - minLobeNumber);
	chatterFreqMap2 = zeros(length(omega_vector), maxLobeNumber - minLobeNumber);

    for N = minLobeNumber:(maxLobeNumber - 1)
		[blimMap2(:, N), chatterFreqMap2(:, N)] = altintasMappingSub(Nt, w2, N, epsilon2, blim2, omega_vector, options.RpmIncrement1D);
    end
    
end

% No chatter detected at any frequency less than depth of cut
if isempty(blimMap1) && isempty(blimMap2)
    SLD = omega_vector * 0;
    chatterFreqs = omega_vector * 0;
    return
end

% Identify the limiting values and chatter freqs at each rpm 
joinedBlim = [blimMap1, blimMap2];
joinedChatterMap = [chatterFreqMap1, chatterFreqMap2];
[SLD, colInds] = min(joinedBlim, [], 2);
rowInds = [1:length(colInds)]';
chatterFreqs = joinedChatterMap(sub2ind(size(joinedChatterMap), rowInds, colInds)) / (2*pi);


% 2j

SLD = fillmissing(SLD, 'linear');
end

function [blimOut, chatterFreqs] = altintasMappingSub(fluteCount, w_N, N, epsilon, blimStart, rpmsVec, rpmIncrement)
    w_N = reshape(w_N, size(epsilon));
    omega_N = round((60/fluteCount)*w_N./(epsilon + 2*(N-1)*pi) / rpmIncrement) * rpmIncrement;     % rpm
    blimOut = realmax*ones(size(rpmsVec));
    chatterFreqs = zeros(size(rpmsVec));

    minRPM = min(rpmsVec);
    maxRPM = max(rpmsVec);

    if (max(omega_N) > minRPM) && (min(omega_N) < maxRPM)
        for i = 2:length(omega_N)
            % Interpolate the omega vector
            if (omega_N(i) > maxRPM && omega_N(i-1) > maxRPM) || (omega_N(i) < minRPM && omega_N(i-1) < minRPM), continue, end

            initialRange = [omega_N(i-1) omega_N(i)];
            
            if omega_N(i) <= maxRPM && omega_N(i-1) <= maxRPM && omega_N(i) >= minRPM && omega_N(i-1) >= minRPM % Spindle speeds completely in range
                interpOmega = initialRange(1):sign(initialRange(2)-initialRange(1))*rpmIncrement:initialRange(2);
                interpBlim = blimStart(i-1):(blimStart(i)-blimStart(i-1))/(length(interpOmega)-1):blimStart(i);
                interpChatterFreq = w_N(i-1):(w_N(i)-w_N(i-1))/(length(interpOmega)-1):w_N(i);
            else % Spindle speeds partially in range, partial blim interpolation
                trimmedRange = min(initialRange, [maxRPM maxRPM]);
                trimmedRange = max(trimmedRange, [minRPM minRPM]);                
                startFrac = (trimmedRange(1) - initialRange(1)) / (initialRange(2) - initialRange(1));
                endFrac   = (trimmedRange(2) - initialRange(1)) / (initialRange(2) - initialRange(1));
                interpOmega = trimmedRange(1):sign(trimmedRange(2)-trimmedRange(1))*rpmIncrement:trimmedRange(2);
                interpLength = length(interpOmega);
                startBlim = (1-startFrac) * blimStart(i-1) + startFrac * blimStart(i);
                endBlim   = (1-endFrac) * blimStart(i-1) + endFrac * blimStart(i);
                startW_N = (1-startFrac) * w_N(i-1) + startFrac * w_N(i);
                endW_N   = (1-endFrac) * w_N(i-1) + endFrac * w_N(i);
                interpBlim = startBlim:(endBlim-startBlim)/(interpLength-1):endBlim;
                interpChatterFreq = startW_N:(endW_N-startW_N)/(interpLength-1):endW_N;
            end
                        
            % Map these blim and chatter frequency to the matrix
            rpmIndices = (interpOmega - minRPM) / rpmIncrement + 1;
            for j = 1:length(rpmIndices)
                if interpBlim(j) < blimOut(rpmIndices(j)) 
                    blimOut(rpmIndices(j)) = interpBlim(j);
                    chatterFreqs(rpmIndices(j)) = interpChatterFreq(j);
                end
            end
        end
    end
end

function D = eig2(A)
% function D = eig2(A)
%
% Compute in one shot the eigen-values of multiples (2 x 2) matrices
%
% INPUT:
%   A: (2 x 2 x n) array
% OUTPUT:
%   D: (2 x n). EIG2 returns in D(:,k) three eigen-values of A(:,:,k)
%
% See also: ParabolaRoots, eig3, eig
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% History:
%     Original 27-May-2010

if size(A,1) ~= 2 || size(A,2) ~= 2
    error('A must be [3x3xn] array');
end

A = reshape(A, 4, []).';

P3 = 1;
% Trace
P2 = -(A(:,1)+A(:,4));

% Determinant
P1 = A(:,1).*A(:,4) - A(:,2).*A(:,3);

% Find the roots of characteristic polynomials
D = ParabolaRoots(P3, P2, P1).';

end % eig2

function roots = ParabolaRoots(varargin)
% function roots = ParabolaRoots(P)
%
% Find roots of second order polynomials
%
% INPUT
%   P: (n x 2) array, each row corresponds to coefficients of each
%   polynomial, P(:,1)*x^2 + P(:,2)*x + P(:,3)
% OUTPUT
%   roots: (n x 2) array, each row correspond to the roots of P
%
% To adjust the parameter below which the the discriminant is considerered
% as nil, use
%   ParabolaRoots(P, tol)
% Adjusting tol is useful to avoid the real roots become complex due to
% numerical accuracy. The default TOL is 0
%
% See also: roots, CardanRoots, eig2
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% History:
%     Original 27-May-2010

% Adjustable parameter
tol = 0;

if nargin<3
    P = varargin{1};
    a = P(:,1);
    b = P(:,2);
    c = P(:,3);
    if nargin>=2
        tol = varargin{2};
    end
else
    [a b c] = deal(varargin{1:3});
    if nargin>=4
        tol = varargin{2};
    end
end

if ~isequal(a,1)
    b = b./a;
    c = c./a;
end

b = 0.5*b;

delta = b.^2 - c;
delta(abs(delta)<tol) = 0;
sqrtdelta = sqrt(delta);

roots = [sqrtdelta -sqrtdelta];
roots = bsxfun(@minus, roots, b);

end