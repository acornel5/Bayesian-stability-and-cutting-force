% dynamicResampleFRF.m
% 2021-10-7
% A. Cornelius

% This function resamples an input FRF to reduce the number of points
% while maintaining accuracy. The function checks how much the accuracy
% would be degraded if any given point was omitted, and if accuracy is
% maintained the point is removed.


function [x, y] = dynamicResampleFRF(x,y)
    % Constants
    maxSingleIterationClear = 10;
    maxDeviation = 0.01;
    
    
    clearInds = false(size(x));
    startInd = 1;
	while startInd <= length(x)-2
        endInd = startInd + 2;
 		while endInd <= min(length(x), startInd + 10)
			startX = x(startInd);
			startY = y(startInd);
			endX = x(endInd);
			endY = y(endInd); 
			nominal = y(endInd-1);
			interpolated = linInterp(x(startInd+1), startX, startY, endX, endY);
            error = interpolated - nominal;
            if (abs(real(error) / real(nominal)) >= maxDeviation) || (abs(imag(error) / imag(nominal)) >= maxDeviation)
				break
            end
            clearInds(endInd-1) = true;
            endInd = endInd+1;
        end
        startInd = endInd;
    end
    x(clearInds) = [];
    y(clearInds) = [];
    if sum(clearInds) > 10
        [x,y] = dynamicResampleFRF(x,y);
    end
end


function y = linInterp(x, x1, y1, x2, y2)
	frac = (x-x1)/(x2-x1);
	y = y1*(1-frac) + y2*frac;
end