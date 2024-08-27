function [samples, diagnostics] = initialSampleGenerator(priorDist, probFunc, NameValueArgs, options)
arguments
    priorDist
    probFunc
    NameValueArgs.SampleCount int32
    options.Parallelize logical = true
end

    global taskList

    timer = tic;
    taskList.InitializeETA(NameValueArgs.SampleCount);
    if options.Parallelize
        parfor i = 1:NameValueArgs.SampleCount
            samples(i) = probFunc(priorDist.Sample(), priorDist, []);
            taskList.StepETA;
        end
    else
        for i = 1:NameValueArgs.SampleCount
            samples(i) = probFunc(priorDist.Sample(), priorDist, []);
            taskList.StepETA;
        end
    end
    taskList.ClearETA;

    [samples(:).Weight] = deal(1); % All initial samples have a weight of 1
    diagnostics.calcTime = toc(timer);
end