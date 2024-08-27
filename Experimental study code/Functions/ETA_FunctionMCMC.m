% A simple progress tracker for MCMC calculations
% Calling this function returns a queue. Sending a value to this queue via
% send(queue, val) will set the counter based on the input values and
% update the display. The input should be a 2-element tuple:
%  - The first element should be the number of accepted samples
%  - The second element should be the total number of calculated samples

function [queue] = ETA_FunctionMCMC(numberToGenerate)
    % Start timer
    timer = tic;
    tracker = 0;
    numberToGenerate = double(numberToGenerate);
    
    global taskList;
    
    function func(reportedProgress)
        vals = reportedProgress{1};
        tracker = vals(1);
        totalCalculated = vals(2);
        elapsedTime = toc(timer);
        eta = seconds(elapsedTime / tracker * (numberToGenerate - tracker));
        eta.Format = 'hh:mm:ss';
        elapsedTime = seconds(elapsedTime);
        elapsedTime.Format = 'hh:mm:ss';
        fracComplete = tracker / numberToGenerate;
        
        taskList.DisplayTasks;
        
        % Output progress text
        if tracker ~= numberToGenerate
            disp([num2str(tracker) ' / ' num2str(totalCalculated) ' / ' num2str(numberToGenerate)])
            disp(['Acceptance rate: ' num2str(100*tracker / totalCalculated, 2) '%'])
            disp(['Elapsed: ' char(elapsedTime) '     ETA: ' char(eta)])
        else
            disp(['COMPLETE. Elapsed time: ' char(elapsedTime)])
        end
        
        % Output progressbar
        progressBarLength = 50;
        filledCount = round(progressBarLength * fracComplete);
        unfilledCount = progressBarLength - filledCount;
        progressBar = ['[' repmat('*', 1, filledCount) repmat('-', 1, unfilledCount) ']'];
        disp(progressBar)
    end

    queue = parallel.pool.DataQueue;

    afterEach(queue, @func);
end