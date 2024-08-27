function saveProbabilityFigureBinaryStability(iteration, sampledVals, nominalVals, testedPoints, param)
    set(0, 'DefaultAxesFontSize', 14, ...
      'DefaultAxesTitleFontWeight', 'normal', ...
      'DefaultAxesTitleFontSizeMultiplier', 1) ;
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'defaultLegendInterpreter','latex');
    set(groot, 'defaultTextInterpreter', 'latex');

    requiredFields = ["Stability"];
    CheckForFields(nominalVals.CalculatedValues, requiredFields, 'saveProbabilityFigureBinaryStability')

    % rpmMatrix = GetRpmMatrix(param);
    rpmMatrix = squeeze(param.CutParamGrid.rpm(:,1,:,1,1));
    % apMatrix = GetApMatrix(param)*1000; % Converted to mm
    apMatrix = squeeze(param.CutParamGrid.b(:,1,:,1,1))*1000;
    allChatterProbabilities = GetChatterProbabilities(sampledVals);

    %% Doesn't support different cut directions!
    for aeInd = 1:length(param.RadialWidth)
        cutDirInd = 1; % HACK
        chatterProbabilities = squeeze(allChatterProbabilities(:,aeInd,:));
	    % Draw image
        f = figure('visible', 'off');
        clf(f)
        surf(rpmMatrix, apMatrix, zeros(size(chatterProbabilities)), chatterProbabilities)
        hold on
        view(2)
        shading interp
        if ~isempty(nominalVals)
            if exist('nominalVals.CalculatedValues.StabilityMap1D', 'var')
                rpmVec = param.RpmMin:param.RpmIncrement1D:param.RpmMax;
                plot3(rpmVec, nominalVals.CalculatedValues.StabilityMap1D(:,aeInd,1,:,cutDirInd)*1000, ones(size(rpmVec)), 'color', 'blue', 'LineWidth', 1);
            else
                contour(rpmMatrix, apMatrix, squeeze(nominalVals.CalculatedValues.Stability(:,aeInd,:,1,cutDirInd))*1000,[500 500], 'blue', 'LineWidth', 1);
            end

            % contour(rpmMatrix, apMatrix, nominalVals.CalculatedValues.StableBendingStress/1e6, [340 340])
        end
    
        %% Plot tested points
        % Extract the stable and chatter history points
        for i = 1:length(testedPoints)
            if abs((testedPoints(i).a - param.RadialWidth(aeInd))/param.RadialWidth(aeInd))>0.01, continue, end
            if testedPoints(i).Broken == 1
                scatter3(testedPoints(i).RPM, testedPoints(i).b*1000, 2, 100, '*', 'MarkerEdgeColor', 'red', 'LineWidth', 2); % Chatter points
            elseif testedPoints(i).Stable == 1
                scatter3(testedPoints(i).RPM, testedPoints(i).b*1000, 2, 100, 'filled', 'MarkerEdgeColor', [100 200 0]/255, 'MarkerFaceColor', [100 200 0]/255);
            elseif testedPoints(i).Stable == 0
                scatter3(testedPoints(i).RPM, testedPoints(i).b*1000, 2, 100, 'x', 'MarkerEdgeColor', 'red', 'LineWidth', 2); % Chatter points
            elseif testedPoints(i).Stable == 0.5
                scatter3(testedPoints(i).RPM, testedPoints(i).b*1000, 2, 100, 'filled', 'MarkerEdgeColor', [255 200 0]/255, 'MarkerFaceColor', [255 200 0]/255);
            else
                error('Invalid stability classification')
            end
        end
    
        %% Legend
        stableForLegend  = scatter3(NaN, NaN, NaN, 100, 'filled', 'MarkerEdgeColor', [100 200 0]/255, 'MarkerFaceColor', [100 200 0]/255);
        unstableForLegend = scatter3(NaN, NaN, NaN, 100, 'x', 'MarkerEdgeColor', 'red', 'LineWidth', 2);
        b = plot([-2 -1], [-2, -1], 'blue', 'LineWidth', 1);
        legend([stableForLegend, unstableForLegend, b], ["Stable" "Unstable" "Tap test stability map"], 'Location', 'northwest')
    
        %% Set axis limits, labels, etc
        title(strcat("Radial width = ", num2str(param.RadialWidth(aeInd)*1000), " mm"))
        xlim([min(rpmMatrix, [], 'all') max(rpmMatrix, [], 'all')])
        ylim([0 max(apMatrix, [], 'all')])
        xlabel('Spindle speed (rpm)')
        ylabel('Axial depth of cut (mm)')
        colormap(flipud(gray))
        c=colorbar('Ticks',[0,0.2,0.4,0.6,0.8, 1],'TickLabels',{'100\%', '80\%', '60\%', '40\%', '20\%', '0\%'});
        set(get(c,'label'),'string','Expected probability of stability');
	    set(get(c,'label'),'interpreter','latex');
	    c.TickLabelInterpreter = 'latex';
	    caxis([0 1])
        hold off;
        folderName = strcat('Figures\\ChatterProbability\\RadialWidth', num2str(param.RadialWidth(aeInd)*1000), 'mm');
        if not(isfolder(folderName)), mkdir(folderName), end
        fileName = [num2str(iteration)];
    
        print('-dpng', strcat(folderName, '\\', fileName), '-r250');
        close(f)
    end
end