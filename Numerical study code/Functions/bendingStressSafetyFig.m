function bendingStressSafetyFig(iteration, sampledVals, nominalVals, testedPoints, param)
    set(0, 'DefaultAxesFontSize', 14, ...
      'DefaultAxesTitleFontWeight', 'normal', ...
      'DefaultAxesTitleFontSizeMultiplier', 1) ;
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'defaultLegendInterpreter','latex');
    set(groot, 'defaultTextInterpreter', 'latex');

    requiredFields = ["StableBreakageLikelihood"];
    CheckForFields(nominalVals.CalculatedValues, requiredFields, 'bendingMomentSafetyFig')

    rpmMatrix = squeeze(param.CutParamGrid.rpm(:,1,:,1,1));
    apMatrix = squeeze(param.CutParamGrid.b(:,1,:,1,1))*1000;

    safetyProbs = zeros(size(sampledVals(1).CalculatedValues.StableBreakageLikelihood));
    cnt = 0;
    for i = 1:length(sampledVals)
        safetyProbs = safetyProbs + sampledVals(i).CalculatedValues.StableBreakageLikelihood * sampledVals(i).Weight;
        cnt = cnt + sampledVals(i).Weight;
    end
    safetyProbs = 1 - safetyProbs / cnt; % Flip because the likelihoods are for breaking not safety!

    %% Doesn't support different cut directions!
    for aeInd = 1:length(param.RadialWidth)
        cutDirInd = 1; % HACK
	    % Draw image
        f = figure('visible', 'off');
        surf(rpmMatrix, apMatrix, zeros(size(apMatrix)), 1-squeeze(safetyProbs(:,aeInd,:)))
        hold on
        view(2)
        shading interp
        if ~isempty(nominalVals)
            contour(rpmMatrix, apMatrix, squeeze(nominalVals.CalculatedValues.StableBreakageLikelihood(:,aeInd,:)), [0.5 0.5], 'b')
        end
    
        %% Plot tested points
        % Extract the stable and chatter history points
        for i = 1:length(testedPoints)
            if ((testedPoints(i).a - param.RadialWidth(aeInd))/param.RadialWidth(aeInd))>0.01, continue, end
            if testedPoints(i).Stable == 1
                scatter3(testedPoints(i).RPM, testedPoints(i).b*1000, 2, 100, 'filled', 'MarkerEdgeColor', [100 200 0]/255, 'MarkerFaceColor', [100 200 0]/255);
            elseif testedPoints(i).Stable == 0
                scatter3(testedPoints(i).RPM, testedPoints(i).b*1000, 2, 100, 'x', 'MarkerEdgeColor', 'red', 'LineWidth', 3); % Chatter points
            elseif testedPoints(i).Stable == 0.5
                scatter3(testedPoints(i).RPM, testedPoints(i).b*1000, 2, 100, 'filled', 'MarkerEdgeColor', [255 200 0]/255, 'MarkerFaceColor', [255 200 0]/255);
            else
                error('Invalid stability classification')
            end
        end
    
        %% Legend
        stableForLegend  = scatter3(NaN, NaN, NaN, 100, 'filled', 'MarkerEdgeColor', [100 200 0]/255, 'MarkerFaceColor', [100 200 0]/255);
        unstableForLegend = scatter3(NaN, NaN, NaN, 100, 'x', 'MarkerEdgeColor', 'red', 'LineWidth', 2);
        legend([stableForLegend, unstableForLegend], ["Stable" "Unstable"], 'Location', 'northwest')
    
        %% Set axis limits, labels, etc
        title(strcat("Radial width = ", num2str(param.RadialWidth(aeInd)*1000), " mm"))
        xlim([min(rpmMatrix, [], 'all') max(rpmMatrix, [], 'all')])
        ylim([0 max(apMatrix, [], 'all')])
        xlabel('Spindle speed (rpm)')
        ylabel('Axial depth of cut (mm)')
        colormap(flipud(gray))
        c=colorbar('Ticks',[0,0.2,0.4,0.6,0.8, 1],'TickLabels',{'100\%', '80\%', '60\%', '40\%', '20\%', '0\%'});
        set(get(c,'label'),'string','Probability of safe bending stress');
	    set(get(c,'label'),'interpreter','latex');
	    c.TickLabelInterpreter = 'latex';
	    caxis([0 1])
        hold off;
        folderName = strcat('Figures\\StableBendingStress\\RadialWidth', num2str(param.RadialWidth(aeInd)*1000), 'mm');
        if not(isfolder(folderName)), mkdir(folderName), end
        fileName = [num2str(iteration)];
    
        print('-dpng', strcat(folderName, '\\', fileName), '-r250');
        close(f)
    end
end