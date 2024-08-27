% SimulationSettings.m
% This file contains all configuration options for running a simulation
% set. If any modifications are required which can't be done through this
% file, the appropriate files should be modified so that the full
% simulation options can be contained in this one file.

% Label info
testName = '12.7mm endmill, numerical case study';

%% Configure global testing parameters. These can be overriden in specific functions if necessary
% Tool
param.ToolDiameter = 12.7e-3; % Tool diameter, meters
param.FluteCount = 3; % Flute count
param.HelixAngle = 0; % Flute helix angle, radians

% Cutting parameters
param.RadialWidth = param.ToolDiameter * [0.5]; % Radial widths of cut, meters
assert(all(param.ToolDiameter >= param.RadialWidth), 'Radial width of cut must be smaller than tool diameter')
param.CutDirection = [2]; % Feed directions. 1 = conventional milling, 2 = climb milling
param.FeedPerTooth = [0.05e-3 0.1e-3]; % Feed rate per flute, meters
param.Rpms = linspace(4000, 12000, 200); % List of spindle speed, rpm
param.Aps = linspace(0e-3, 8e-3, 100);
param.CutParamGrid = cuttingParameterGrid(... % Defines the grid of parameters using full-factorial combinations
    param.Rpms, ... % Spindle speeds
    param.RadialWidth, ... % Radial width
    param.Aps, ... % Axial depth
    param.FeedPerTooth, ... % Feed per flute
    param.CutDirection); % Cut direction

% General calculation settings
param.MaxIterations = 50;
param.SaveCalculatedData = false; % Remove calculated data to save file space
param.StabilityUncertainty = 0.01e-3;
param.StabilityCutoff = 1e-6;
param.ChatterFreqUncertainty = 20;
param.SpindlePowerUncertainty = 20;

% Plot settings
param.PlotNominals = true; % Plot nominal stability map and values on the diagrams.

%% Define the prior and true values
% [ks (N/m^2), beta (rad), kte, fn1 (Hz), k1 (N/m), zeta1, fn2, k2, zeta2]
priorMean = [600e6 deg2rad(68) 50e3 1500 10e6 0.03];
priorCov = diag((priorMean.*[0.2 0.1 0.2 0.1 0.2 0.3]).^2);
priorDist = Distribution_MVN_Truncated(priorMean, priorCov, [0 0 0 0 0 0], [realmax realmax realmax realmax realmax realmax]);

clear priorMean priorCov

param.sampleVariableNames = ["K_s (N/mm^2)" "\beta (deg)" "k_{te} (N/mm)" "f_n (Hz)" "k (MN/m)" "\zeta (%)"];
parser = sampleParser('KsBeta', [1 2], 'fnKZeta', [4 6], ["ks" "beta" "kte" "fn" "k" "zeta"]);

% Nominal values for comparison. Only plotted if PlotNominals is turned on.
trueStabilityMapVals = [700e6 deg2rad(68) 10e3 1200 8e6 0.02];

%% Define functions
% Physics algorithms
physicsFuncs = {...
    @(sample) physicsFuncFrfFromSymmetricModes(sample, parser, FrequencyRange='auto'), ...
    @(sample, previousCalcs) physicsFunZOA_StabilityLikelihood(sample, parser, previousCalcs, param, RpmIncrement1D = 25, StabilityUncertainty='inParam', AssumeSymmetric=true, ReduceSymmetricFRF=true), ...
    @(sample) physicsMeanCuttingTorque(sample, parser, param)};
physicsDelegate = @(sample) GenericPhysicsDelegate(sample, physicsFuncs, ["Stability" "ChatterFrequency" "SpindlePower" "FRFxx" "FrequenciesFRF"]);
clear physicsFuncs

% Algorithm to select the next test
testSelectors = {
    @(iteration, sampledVals, testResults, acceptableTests) testSelectMaximizeMRR( ...
        iteration, ...
        sampledVals, ...
        testResults, ...
        acceptableTests, ...
        param, ...
        ChatterRiskTolerance = 0.5), @(iteration, sampledVals, testResults, acceptableTests) testSelectUncertaintyReduction2(iteration, sampledVals, testResults, acceptableTests, param, parser, ChatterFreqUncertainty='inParam', SpindlePowerUncertainty='inParam', SampleCountToUse=1000)};
selectorNames = ["Chatter risk" "Uncertainty reduction"];
constraints = @(samples) testConstraints(samples, param, ChatterProbability=[0 0.5]);

testSelector = @(iteration, sampledVals, testResults) ...
    testSelectSwitching(iteration, sampledVals, testResults, testSelectors, selectorNames, param, constraints);

% Algorithm to label tests as stable or unstable
% testClassifier = @(iteration, test) manualTestClassifier(iteration, test, param);

testClassifier = @(iteration, test) tdsTestClassifier(iteration, ...
    test, ...
    param, ...
    trueStabilityMapVals(4:6), ...
    trueStabilityMapVals(1:2));


% PDF formula. The prior is not injected here because it's not used while calculating updated probabilities for the old points
pdfFormula = @(sample, prior, testedPoints) acceptanceFormulaLikelihood(sample, prior, testedPoints, ...
    physicsDelegate, ... % Algorithm to generate the full stability maps
    param, ... % General parameters
    StabilityUncertainty='precalculated', ... % The stability uncertainty is calculated directly in the ZOA function to save time
    chatterFreqStdDev = 10, ... % Standard deviation for chatter frequency agreement, Hz
    LearnFromUnbrokenCuts = false, PowerStdDev='inParam', LearnPowerFromUnstableCuts=false); % Whether the PDF should learn from cuts that don't break the tool

% Sample generators
sampleCount = 4000;
initialSampleGenerator = @() initialSampleGenerator(priorDist, pdfFormula, SampleCount = sampleCount, Parallelize = true);
updateSampleGenerator = @(priorSamples, newTestResults, oldTestResults) mcmcSampleGenerator(priorSamples, newTestResults, oldTestResults, priorDist, pdfFormula, SampleCount = sampleCount, Beta=0.02);
clear sampleCount

% Figure plotting
plotFuncs = {@saveProbabilityFigureBinaryStability, @PlotMatrixCorrelations, @chatterFrequencyPlot, @PlotPowerPredictions};
plotFuncNames = ["Probability of stability" "Variable plotmatrix" "Chatter frequencies" "Power estimates"];
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