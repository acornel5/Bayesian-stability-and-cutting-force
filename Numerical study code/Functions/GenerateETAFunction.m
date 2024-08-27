% A simple progress tracker for single or multi-threaded calculations
% Calling this function returns a queue. Sending a value to this queue via
% send(queue, val) will increment the counter by one and update the
% displayed progress. Currently, it does not matter what value you send.

function [queue, timer] = GenerateETAFunction(numberToGenerate)
    % Start timer
    timer = tic;
    tracker = 0;
    numberToGenerate = double(numberToGenerate);
    
    global taskList;
    
    function func(num)
        tracker = tracker + 1;
        % Reset the timer on the first iteration to help mitigate timing
        % startup error
        if tracker == 1
            timer = tic;
        end
        elapsedTime = toc(timer);
        eta = seconds(elapsedTime / tracker * (numberToGenerate - tracker));
        eta.Format = 'hh:mm:ss';
        elapsedTime = seconds(elapsedTime);
        elapsedTime.Format = 'hh:mm:ss';
        fracComplete = tracker / numberToGenerate;
        
        taskList.DisplayTasks;
        
        % Output progress text
        if tracker ~= numberToGenerate
            disp([num2str(tracker) ' / ' num2str(numberToGenerate)])
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