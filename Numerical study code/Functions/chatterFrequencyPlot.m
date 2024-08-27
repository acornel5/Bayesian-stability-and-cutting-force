function chatterFrequencyPlot(iteration, sampledVals, nominalVals, testResults, param)
    requiredFields = ["ChatterFrequency"];
    CheckForFields(nominalVals.CalculatedValues, requiredFields, 'chatterFrequencyPlot')

    set(0, 'DefaultAxesFontSize', 14, ...
      'DefaultAxesTitleFontWeight', 'normal', ...
      'DefaultAxesTitleFontSizeMultiplier', 1) ;
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'defaultLegendInterpreter','latex');
    set(groot, 'defaultTextInterpreter', 'latex');

    rpms = param.Rpms;

    allChatterFrequencies = GetSampledField(sampledVals, 'ChatterFrequency');
    allTrueChatterFreqs = nominalVals.CalculatedValues.ChatterFrequency;

    % Define bin limits
    if any(allTrueChatterFreqs ~= 0)
        binEdges = linspace(min(allTrueChatterFreqs, [], 'all')*0.8, max(allTrueChatterFreqs, [], 'all')*1.2, 100);
    else
        binEdges = linspace(min(allChatterFrequencies, [], 'all')*0.8, max(allChatterFrequencies, [], 'all')*1.2, 100);
    end

    % Create figure for each radial width
    for aeInd = 1:length(param.RadialWidth)
        
        cutDirInd = 1; % HACK

        chatterFrequenciesSlice = extractCalculatedStructField(sampledVals, 'ChatterFrequency', [0, aeInd, 1, 1, cutDirInd]);
        chatterFrequencies = reshape(chatterFrequenciesSlice, [size(chatterFrequenciesSlice,1), size(chatterFrequenciesSlice,6)]);
        % chatterFrequencies = squeeze(allChatterFrequencies(:,aeInd,1,1,:));
        % trueChatterFreqs = squeeze(allTrueChatterFreqs(:,aeInd,1,1,:));
        trueChatterFreqs = extractCalculatedStructField(nominalVals, 'ChatterFrequency', [0, aeInd, 1, 1, cutDirInd]);

        counts = zeros(length(rpms), length(binEdges)-1);
        for i = 1:length(rpms)
            counts(i, :) = histcounts(chatterFrequencies(i, :), binEdges);
        end
        counts = counts ./ size(chatterFrequencies,2);
        
        % Calculate the X and Y matrices
        binCenters = binEdges(1:end-1) + (binEdges(2) - binEdges(1))/2;
        [freqMatrix, rpmMatrix] = meshgrid(binCenters, rpms);
    
        % Create the figure
        fig = figure('visible', 'off');
        surf(rpmMatrix, freqMatrix, zeros(size(rpmMatrix)), counts);
        hold on
        
        % Plot nominal chatter freqs
        b = plot3(rpmMatrix(:,1), trueChatterFreqs, ones(size(trueChatterFreqs)), 'b');
        
        % Plot actual chatter freqs
        for i = 1:length(testResults)
            if ((testResults(i).a - param.RadialWidth(aeInd)) > param.RadialWidth(aeInd)) > 0.01, continue, end
            if testResults(i).ChatterFrequency == 0, continue, end
            if length(testResults(i).ChatterFrequency) == 1
                scatter3(testResults(i).RPM, testResults(i).ChatterFrequency, 1, 'o', 'MarkerEdgeColor', 'blue')
            else
                rpmvec = testResults(i).RPM*ones(size(testResults(i).ChatterFrequency));
                sizes = linspace(30, 5, length(testResults(i).ChatterFrequency));
                scatter3(rpmvec, testResults(i).ChatterFrequency, ones(size(sizes)), sizes, 'o', 'MarkerEdgeColor', 'blue')
            end
        end
        
        view(2)
        shading interp
        colormap(flipud(gray))
        c=colorbar;
        % c=colorbar('Ticks',[0,0.2,0.4,0.6,0.8, 1],'TickLabels',{'0\%', '20\%', '40\%', '60\%', '80\%', '100\%'});
        set(get(c,'label'),'string','Density of predicted chatter frequencies (1/Hz)');
        set(get(c,'label'),'interpreter','latex');
        c.TickLabelInterpreter = 'latex';
    
        xlim([min(rpmMatrix, [], 'all'), max(rpmMatrix, [], 'all')])
        ylim([min(binEdges), max(binEdges)])
        xlabel('Spindle speed (rpm)')
        ylabel('Chatter frequency (Hz)')
		
		legend([b], ["Nominal dominant chatter freq."], Location="southwest")
    

        folderName = strcat('Figures\\ChatterFrequencies\\RadialWidth', num2str(param.RadialWidth(aeInd)*1000), 'mm');
        if not(isfolder(folderName)), mkdir(folderName), end
        fileName = [num2str(iteration)];
        print(fig, '-dpng', strcat(folderName, '\\', fileName), '-r250');
        close(fig)

    end
end


function arr = smashArray(arr)
    newSize = size(arr);
    newSize(newSize==1) = [];
    if length(newSize) == 1, newSize = [newSize 1]; end
    arr = reshape(arr, newSize);
end