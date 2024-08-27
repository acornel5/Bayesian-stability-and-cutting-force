function recommendedTest = testSelectSwitching(iteration, sampledVals, testResults, functionList, functionNames, param, constraintFunc, options)
    arguments
        iteration
        sampledVals
        testResults
        functionList
        functionNames
        param
        constraintFunc
        options.AskForConfirmation = true
    end


    if exist('constraintFunc', 'var')
        acceptableTests = constraintFunc(sampledVals);
    else
        acceptableTests = true(size(param.CutParamGrid.rpm));
    end


    pass = false;
    while ~pass
        % Select the test
        if length(functionNames) ~= 1
            disp('The following test selectors are available: ')
            for i = 1:length(functionNames)
                disp(strcat(num2str(i), ": ", functionNames(i)))
            end
            disp('Or enter 0 to input a manual test.')
            index = input('Which do you want to use? ');
        else
            index = 1;
        end

        if index ~= 0
            [recommendedTest, extraMetrics] = functionList{index}(iteration, sampledVals, testResults, acceptableTests);
        else
            recommendedTest = manualTestSelect(param);
            extraMetrics = {};
        end

        % Confirm the test
        disp("Run at test at:")
        disp(strcat("Spindle speed: ", num2str(recommendedTest.RPM), " rpm"))
        disp(strcat("Axial depth: ", num2str(recommendedTest.b*1000), " mm, ", num2str(recommendedTest.b*1000/25.4), " in"))
        disp(strcat("Radial depth: ", num2str(recommendedTest.a*1000), " mm, ", num2str(recommendedTest.a*1000/25.4), " in"))
        disp(strcat("Feed per flute: ", num2str(recommendedTest.fz*1000), " mm, ", num2str(recommendedTest.fz*1000/25.4), " in"))
        disp(strcat("Feedrate: ", num2str(recommendedTest.RPM*recommendedTest.fz*param.FluteCount*1000), " mm/min, ", num2str(recommendedTest.RPM*recommendedTest.fz*param.FluteCount*1000/25.4), " in/min"))
        if recommendedTest.CutDirection == 1
            disp("Cut direction: conventional")
        else
            disp("Cut direction: climb")
        end

        disp(strcat("MRR: ", num2str(recommendedTest.MRR*1e6), " cm^3/min, ", num2str(recommendedTest.MRR*1e6/2.54^3), " in^3/min"))
        for i = 1:length(extraMetrics)
            disp(extraMetrics{i})
        end
        disp("")

        if options.AskForConfirmation
            accept = input("Enter 1 to accept, 0 to override this test: ");
            if accept ~= 0 && accept ~= 1
                continue
            end
            if accept == 1
                pass = true;
                continue
            end
        else
            pass = true;
        end
    end
	
	close all % Close any diagnostic figures that were opened
end


function test = manualTestSelect(param)
    % Get the entry
    if length(param.Rpms) > 1
        test.RPM = input("Enter spindle speed in rpm: ");
    else
        test.RPM = param.Rpms;
    end
    if length(param.Aps) > 1
        test.b   = input("Enter axial depth in mm: ")/1000;
    else
        test.b = param.Aps;
    end
    if length(param.RadialWidth) > 1
    test.a = input("Enter radial depth in mm: ")/1000;
    else
        test.a = param.RadialWidth;
    end
    if length(param.FeedPerTooth) > 1
    test.fz = input("Enter feed per flute in mm: ")/1000;
    else
        test.fz = param.FeedPerTooth;
    end
    if length(param.CutDirection) > 1
        test.CutDirection = input("Enter 1 for conventional or 2 for climb: ");
    else
        test.CutDirection = param.CutDirection;
    end
    test.MRR = test.RPM * test.b * test.a * test.fz * param.FluteCount;

    % Find the closest test index
    [~, rpmInd] = min(abs(test.RPM-param.Rpms));
    [~, aInd] = min(abs(test.a-param.RadialWidth));
    [~, bInd] = min(abs(test.b - param.Aps));
    [~, fzInd] = min(abs(test.fz - param.FeedPerTooth));
    [~, cutDirInd] = min(abs(test.CutDirection-param.CutDirection));

    test.CutParamInd = sub2ind(size(param.CutParamGrid.rpm), rpmInd, aInd, bInd, fzInd, cutDirInd);
end