classdef Distribution_Mixed
    %DISTRIBUTION_MVN   A distribution composed of mixed sub-distributions
    % Contains various properties for sampling or evaluating from the
    % target probability distribution
    
    properties
        SubDistributions
        dimension
		PenaltyFunction
    end
    
    methods
        function obj = Distribution_Mixed(subDistributions, options)
			arguments
				subDistributions
				options.PenaltyFunction = [];
			end
            obj.SubDistributions = subDistributions;
			obj.PenaltyFunction = options.PenaltyFunction;
            
            obj.dimension = 0;
            for i = 1:length(obj.SubDistributions)
                obj.dimension = obj.Dimension + obj.SubDistributions{i}.Dimension;
            end
        end
        
        function prob = PDF(obj, theta)
            prob = 1;
            counter = 1;
            for i = 1:length(obj.SubDistributions)
                range = counter:counter+obj.SubDistributions{i}.Dimension-1;
                prob = prob * obj.SubDistributions{i}.PDF(theta(range));
                counter = range(end)+1;
            end
			
			if ~isempty(obj.PenaltyFunction)
				prob = prob .* obj.PenaltyFunction(theta);
			end
        end
		
		function prob = logPDF(obj, theta)
            prob = 0;
            counter = 1;
            for i = 1:length(obj.SubDistributions)
                range = counter:counter+obj.SubDistributions{i}.Dimension-1;
                prob = prob + obj.SubDistributions{i}.logPDF(theta(range));
                counter = range(end)+1;
            end
			
			if ~isempty(obj.PenaltyFunction)
				prob = prob + obj.PenaltyFunction(theta);
			end
        end

        
        function theta = Sample(obj)
            theta = zeros(1, obj.Dimension);
            
            counter = 1;
            for i = 1:length(obj.SubDistributions)
                range = counter:counter+obj.SubDistributions{i}.Dimension-1;
                theta(range) = obj.SubDistributions{i}.Sample();
                counter = range(end)+1;
            end

            % Rejection sample with the penalty function
            penalty = obj.PenaltyFunction(theta);
            a = rand();
            if a > penalty, theta = obj.Sample(); end
        end
        
        function d = Dimension(obj)
            d = obj.dimension;
        end

        
    end
end