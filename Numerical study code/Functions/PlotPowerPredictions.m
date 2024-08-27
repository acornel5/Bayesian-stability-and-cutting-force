function PlotPowerPredictions(iteration, samples, nominalVals, testResults, param)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'defaultLegendInterpreter','latex');
    set(groot, 'defaultTextInterpreter', 'latex');

    if length(testResults) == 0
        return
    end

    f = figure(1);
    clf
    hold on

    % Plot actual results
    for i = 1:length(testResults)
        plot(i, testResults(i).Power, 'ob')
        powerPredictions = [];
        for j = 1:length(samples)
            powerPredictions = [powerPredictions; repmat(samples(j).CalculatedValues.SpindlePower(testResults(i).CutParamInd), samples(j).Weight, 1)];
        end
        meanPrediction = mean(powerPredictions);
        lowerBound = meanPrediction - prctile(powerPredictions, 2.5);
        upperBound = prctile(powerPredictions, 97.5) - meanPrediction;
        errorbar(i,meanPrediction,lowerBound, upperBound,"vertical",'r')
    end

    folderName = 'Figures\\PowerPredictions';
    if not(isfolder(folderName)), mkdir(folderName), end
    fileName = [num2str(iteration)];
    print(f, '-dpng', strcat(folderName, '\\', fileName), '-r600');

    close(f)

end