% StartSimulation.m
% This file generates the initial samples for the cutting test updates and
% initializes other relevant variables.

%% Clear all data to start with blank slate
clear
close all

%% Add the function source code folder to the path
addpath("Functions");

%% Create output directory
if not(isfolder('Figures'))
    mkdir('Figures')
end

%% Define the time tracker
global taskList;
taskList = TaskList();

%% Load the simulation parameters
SimulationSettings();
taskList.Add(testName);

%% Generate starting stability maps
% Generating the nominal stability map
nominalVals = pdfFormula(trueStabilityMapVals, priorDist, []);
taskList.Add("Generating initial samples");
sampledVals = initialSampleGenerator();

%% Initialize tested point array
taskList.Update('Save initial data');
testResults = struct('RPM', {}, 'b', {}, 'a', {}, 'fz', {}, 'MRR', {}, 'CutDirection', {}, 'Stable', {}, 'ChatterFrequency', {}, 'CutParamInd', {}, 'Power', {});
iteration = 0;
chatterFrequency = 0;
% save('Figures\\StartingData')

taskList.Remove;
taskList.Remove;

ResumeSimulation();