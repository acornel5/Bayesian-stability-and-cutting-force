classdef Distribution_Uniform
    %DISTRIBUTION_MVN   A fixed uniform distribution
    % Contains various properties for sampling or evaluating from the
    % target probability distribution
    
    properties
        LowerLimit
        UpperLimit
    end
    
    methods
        function obj = Distribution_Uniform(lowerLimit, upperLimit)
            if ~(isvector(lowerLimit) && isvector(upperLimit) && length(lowerLimit) == length(upperLimit))
                error('Limits must be vectors of the same size')
            end
            if any(lowerLimit > upperLimit)
                error('Upper limit elements must be uniformly greater than lower limit')
            end
            obj.LowerLimit = lowerLimit;
            obj.UpperLimit = upperLimit;
        end
        
        function prob = PDF(obj, theta)
            if all(theta >= obj.LowerLimit) && all(theta <= obj.UpperLimit)
                % Trim out any variables where the upper and lower value
                % are the same
                varWidths = obj.UpperLimit - obj.LowerLimit;
                varWidths(varWidths == 0) = [];
                
                prob = 1 / prod(varWidths);
            else
                prob = 0;
            end
        end
		
		function prob = logPDF(obj, theta)
			if all(theta >= obj.LowerLimit) && all(theta <= obj.UpperLimit)
                % Trim out any variables where the upper and lower value
                % are the same
                varWidths = obj.UpperLimit - obj.LowerLimit;
                varWidths(varWidths == 0) = [];
                
                prob = -log(prod(varWidths));
            else
                prob = -inf;
            end
		end
        
        function theta = Sample(obj)
            theta = rand(size(obj.LowerLimit)) .* (obj.UpperLimit - obj.LowerLimit) + obj.LowerLimit;
        end
        
        function d = Dimension(obj)
            d = length(obj.LowerLimit);
        end
    end
end