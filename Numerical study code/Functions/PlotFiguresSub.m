function PlotFiguresSub(plotFuncs, funcNames, iteration, sampledVals, nominalVals, testResults, param)
    global taskList;
    taskList.Add('Outputting figures');

    for i = 1:length(plotFuncs)
        taskList.Update(funcNames(i));
        plotFuncs{i}(iteration, sampledVals, nominalVals, testResults, param);
    end

    taskList.Remove;
end