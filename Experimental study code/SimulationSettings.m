% SimulationSettings.m
% This file contains all configuration options for running a simulation
% set. If any modifications are required which can't be done through this
% file, the appropriate files should be modified so that the full
% simulation options can be contained in this one file.

% Label info
testName = '12.7mm endmill, steel testing';

%% Configure global testing parameters. These can be overriden in specific functions if necessary
% Tool
param.ToolDiameter = 12e-3; % Tool diameter, meters
param.FluteCount = 4; % Flute count
param.HelixAngle = 0; % Flute helix angle, radians

% Cutting parameters
param.RadialWidth = param.ToolDiameter * [0.5]; % Radial widths of cut, meters
assert(all(param.ToolDiameter >= param.RadialWidth), 'Radial width of cut must be smaller than tool diameter')
param.CutDirection = [2]; % Feed directions. 1 = conventional milling, 2 = climb milling
param.FeedPerTooth = [0.03e-3 0.06e-3]; % Feed rate per flute, meters
param.Rpms = linspace(4000, 7000, 200); % List of spindle speed, rpm
param.Aps = linspace(0e-3, 1.5e-3, 100);
param.CutParamGrid = cuttingParameterGrid(... % Defines the grid of parameters using full-factorial combinations
    param.Rpms, ... % Spindle speeds
    param.RadialWidth, ... % Radial width
    param.Aps, ... % Axial depth
    param.FeedPerTooth, ... % Feed per flute
    param.CutDirection); % Cut direction

% General calculation settings
param.MaxIterations = 50;
param.SaveCalculatedData = false; % Remove calculated data to save file space
param.StabilityUncertainty = 0.05e-3;

% Plot settings
param.PlotNominals = true; % Plot nominal stability map and values on the diagrams.

%% Define the prior and true values
% [ks (N/m^2), beta (rad), kte, fn1 (Hz), k1 (N/m), zeta1, fn2, k2, zeta2]
forcePriorMean = [1800e6	deg2rad(68) 50e3];
forcePriorStdDevs = forcePriorMean .* [0.2 0.1 0.25]; % Define prior stdDevs as a percentage of the mean
forcePriorCov = diag(forcePriorStdDevs.^2);
forcePrior = Distribution_MVN_Truncated(forcePriorMean, forcePriorCov, [0 0 0], [realmax realmax realmax]);
fn1Prior = Distribution_Uniform(1000, 2000);
dist3 = Distribution_MVN_Truncated([5e6, 0.03],diag([1e6, 0.01]).^2, [0 0], [realmax realmax]);
dist4 = Distribution_Uniform(1000, 2000);
dist5 = Distribution_MVN_Truncated([5e6, 0.03],diag([1e6, 0.01]).^2, [0 0], [realmax realmax]);

penaltyFunc = @(btheta) exp(-1/2*(max(btheta(7)-btheta(4),0)).^2/100^2); % Insert an extra function like this
priorDist = Distribution_Mixed({forcePrior fn1Prior dist3 dist4 dist5}, PenaltyFunction=penaltyFunc);
%priorDist = Distribution_MVN_Truncated(forcePriorMean, forcePriorCov, zeros(size(forcePriorMean)), realmax*ones(size(forcePriorMean)));
% param.plotLim = [forcePriorMean - 2*forcePriorStdDevs; forcePriorMean + 2*forcePriorStdDevs]';
clear forcePriorMean forcePriorCov

param.sampleVariableNames = ["$k_s$" "$\beta$" "$k_{te}$" "$f_{n1}$" "$k_1$" "$\zeta_1$" "$f_{n2}$" "$k_2$" "$\zeta_2$"];
parser = sampleParser('KsBeta', [1 2], 'fnKZeta', [4 9], ["ks" "beta" "kte" "fn1" "k1" "zeta1" "fn2" "k2" "zeta2"]);

% Nominal values for comparison. Only plotted if PlotNominals is turned on.
trueStabilityMapVals = [2.16E+09	deg2rad(68)	5.00E+04 1243 2.92e6 0.0619 1747 1.2e7 0.0288];

%% Define functions
% Physics algorithms
physicsFuncs = {...
    @(sample) physicsFuncFrfFromSymmetricModes(sample, parser, FrequencyRange='auto'), ...
    @(sample, previousCalcs) physicsFunZOA_StabilityLikelihood(sample, parser, previousCalcs, param, RpmIncrement1D = 25, StabilityUncertainty='inParam', AssumeSymmetric=true, ReduceSymmetricFRF=true), ...
    @(sample) physicsMeanCuttingTorque(sample, parser, param)};
physicsDelegate = @(sample) GenericPhysicsDelegate(sample, physicsFuncs, ["Stability" "ChatterFrequency" "SpindlePower"]);
clear physicsFuncs

% Algorithm to select the next test
testSelectors = {
    @(iteration, sampledVals, testResults, acceptableTests) testSelectUncertaintyReduction2( ...
        iteration, ...
        sampledVals, ...
        testResults, ...
        acceptableTests, ...
        param, ...
        parser, ...
        ChatterFreqUncertainty = 20, SampleCountToUse=500, OverrideConstraints = true, SpindlePowerUncertainty=-0.25), ...
    @(iteration, sampledVals, testResults, acceptableTests) testSelectMaximizeMRR( ...
        iteration, ...
        sampledVals, ...
        testResults, ...
        acceptableTests, ...
        param, ...
        ChatterRiskTolerance = 0.5)};
selectorNames = ["Uncertainty reduction" "Chatter risk"];
constraints = @(samples) testConstraints(samples, param, ChatterProbability=[0 0.5]);

testSelector = @(iteration, sampledVals, testResults) ...
    testSelectSwitching(iteration, sampledVals, testResults, testSelectors, selectorNames, param, constraints);

% Algorithm to label tests as stable or unstable
testClassifier = @(iteration, test) manualTestClassifier(iteration, test, param);

% PDF formula. The prior is not injected here because it's not used while calculating updated probabilities for the old points
pdfFormula = @(sample, prior, testedPoints) acceptanceFormulaLikelihood(sample, prior, testedPoints, ...
    physicsDelegate, ... % Algorithm to generate the full stability maps
    param, ... % General parameters
    StabilityUncertainty='precalculated', ... % The stability uncertainty is calculated directly in the ZOA function to save time
    chatterFreqStdDev = 20, ... % Standard deviation for chatter frequency agreement, Hz
    PowerStdDev=-0.25, LearnPowerFromUnstableCuts=true); % Power learning uncertainty. Negative numbers indicate a fraction of the nominal prediction, positive values are fixed uncertainties in W. Set to realmax to remove power learning.

% Sample generators
sampleCount = 4000;
initialSampleGenerator = @() initialSampleGenerator(priorDist, pdfFormula, SampleCount = sampleCount, Parallelize = true);
updateSampleGenerator = @(priorSamples, newTestResults, oldTestResults) mcmcSampleGenerator(priorSamples, newTestResults, oldTestResults, priorDist, pdfFormula, SampleCount = sampleCount, Beta=0.02);
clear sampleCount

% Figure plotting
plotFuncs = {@saveProbabilityFigureBinaryStability, @PlotMatrixCorrelations, @chatterFrequencyPlot};
plotFuncNames = ["Probability of stability" "Variable plotmatrix" "Chatter frequencies"];
figureGeneration = @(iteration, sampledVals, nominalVals, testResults) PlotFiguresSub(plotFuncs, plotFuncNames, iteration, sampledVals, nominalVals, testResults, param);


% The cutting parameter grid is in the order: (rpm, radialWidth,
% axialDepth, feedPerFlute, cutDirection
function grid = cuttingParameterGrid(rpms, as, bs, fzs, cutDirs)
    baseGrid = zeros(length(rpms), length(as), length(bs), length(fzs), length(cutDirs));
    grid.rpm = baseGrid;
    grid.a = baseGrid;
    grid.b = baseGrid;
    grid.fz = baseGrid;
    grid.cutDir = baseGrid;

    for rpmInd = 1:length(rpms)
        for aInd = 1:length(as)
            for bInd = 1:length(bs)
                for fzInd = 1:length(fzs)
                    for cutDirInd = 1:length(cutDirs)
                        grid.rpm(rpmInd, aInd, bInd, fzInd, cutDirInd) = rpms(rpmInd);
                        grid.a(rpmInd, aInd, bInd, fzInd, cutDirInd) = as(aInd);
                        grid.b(rpmInd, aInd, bInd, fzInd, cutDirInd) = bs(bInd);
                        grid.fz(rpmInd, aInd, bInd, fzInd, cutDirInd) = fzs(fzInd);
                        grid.cutDir(rpmInd, aInd, bInd, fzInd, cutDirInd) = cutDirs(cutDirInd);
                    end
                end
            end
        end
    end
end