function [recommendedTest, extraMetrics] = testSelectMaximizeMRR(iteration, sampledVals, testResults, acceptableTests, param, NameValueArgs)
    arguments
        iteration
        sampledVals
        testResults
        acceptableTests
        param
        NameValueArgs.ChatterRiskTolerance
    end
    
    chatterProbs = GetChatterProbabilities(sampledVals);
    mrrMatrix    = param.CutParamGrid.a .* param.CutParamGrid.b .* param.CutParamGrid.rpm .* param.CutParamGrid.fz .* param.FluteCount;

    mrrMatrix(chatterProbs > NameValueArgs.ChatterRiskTolerance) = 0;

    mrrMatrix(~acceptableTests) = 0;

    [~, ind] = max(mrrMatrix(:));

    testRpm = param.CutParamGrid.rpm(ind);
    testa   = param.CutParamGrid.a(ind);
    testb   = param.CutParamGrid.b(ind);
    testfz  = param.CutParamGrid.fz(ind);
    testDir = param.CutParamGrid.cutDir(ind);
    testMRR = mrrMatrix(ind);

    recommendedTest = struct("RPM", testRpm, "a", testa, "b", testb, "fz", testfz, "CutDirection", testDir, "MRR", testMRR, "CutParamInd", ind);

    bestPreviousMRR = 0;
    for i = 1:length(testResults)
        if testResults(i).Stable && (testResults(i).MRR > bestPreviousMRR), bestPreviousMRR = testResults(i).MRR; end
    end
    extraMetrics = strcat(num2str(((testMRR/bestPreviousMRR)-1)*100,'%.2f'), "% MRR improvement");

    [~,aeInd,~,fzInd,cutDirInd] = ind2sub(size(param.CutParamGrid.rpm), ind);
    rpmGrid = squeeze(param.CutParamGrid.rpm(:,1,:,1,1));
    apGrid  = squeeze(param.CutParamGrid.b(:,1,:,1,1))*1000;
    figure(1)
    clf
    hold on
    chatterProbs = GetChatterProbabilities(sampledVals);
    chatterProbs = squeeze(chatterProbs(:,aeInd,:,1,cutDirInd));
    surface(rpmGrid, apGrid, zeros(size(rpmGrid)), chatterProbs)
    contour(rpmGrid, apGrid, squeeze(acceptableTests(:,aeInd,:,fzInd,cutDirInd)), [0.1 0.1], 'r')
    contour(rpmGrid, apGrid, squeeze(acceptableTests(:,aeInd,:,fzInd,cutDirInd)), [0.9 0.9], 'b')        
	
    plot3(recommendedTest.RPM, recommendedTest.b*1000, 200, "diamond", MarkerEdgeColor=[0.9290 0.6940 0.1250], MarkerFaceColor=[0.9290 0.6940 0.1250], MarkerSize=10)
	
        % Plot tested points
    for i = 1:length(testResults)
        if abs((testResults(i).a - recommendedTest.a)/recommendedTest.a)>0.01, continue, end
        if testResults(i).Broken == 1
            scatter3(testResults(i).RPM, testResults(i).b*1000, 2, 100, '*', 'MarkerEdgeColor', 'red', 'LineWidth', 2); % Chatter points
        elseif testResults(i).Stable == 1
            scatter3(testResults(i).RPM, testResults(i).b*1000, 2, 100, 'filled', 'MarkerEdgeColor', [100 200 0]/255, 'MarkerFaceColor', [100 200 0]/255);
        elseif testResults(i).Stable == 0
            scatter3(testResults(i).RPM, testResults(i).b*1000, 2, 100, 'x', 'MarkerEdgeColor', 'red', 'LineWidth', 2); % Chatter points
        elseif testResults(i).Stable == 0.5
            scatter3(testResults(i).RPM, testResults(i).b*1000, 2, 100, 'filled', 'MarkerEdgeColor', [255 200 0]/255, 'MarkerFaceColor', [255 200 0]/255);
        else
            error('Invalid stability classification')
        end
    end

    shading interp
    title(strcat("Radial width = ", num2str(param.RadialWidth(aeInd)*1000), " mm"))
    xlim([min(rpmGrid, [], 'all') max(rpmGrid, [], 'all')])
    ylim([0 max(apGrid, [], 'all')])
    xlabel('Spindle speed (rpm)')
    ylabel('Axial depth of cut (mm)')
    colormap(flipud(gray))
    c=colorbar('Ticks',[0,0.2,0.4,0.6,0.8, 1],'TickLabels',{'100\%', '80\%', '60\%', '40\%', '20\%', '0\%'});
    set(get(c,'label'),'string','Expected probability of stability');
    set(get(c,'label'),'interpreter','latex');
    c.TickLabelInterpreter = 'latex';
    caxis([0 1])
    hold off;
    folderName = 'Figures\\TestSelection';
    if not(isfolder(folderName)), mkdir(folderName), end
    fileName = [num2str(iteration)];
    print('-dpng', strcat(folderName, '\\', fileName), '-r300');
end
