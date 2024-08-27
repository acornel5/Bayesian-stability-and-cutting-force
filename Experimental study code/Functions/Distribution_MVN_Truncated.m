classdef Distribution_MVN_Truncated
    %DISTRIBUTION_MVN_Truncated   A fixed truncated multivariate Gaussian distribution
    % Contains various properties for sampling or evaluating from the
    % target probability distribution
    
    properties
        Mean
        Covariance
		LowerLim
		UpperLim
        ZeroCov
    end
    
    methods
        function obj = Distribution_MVN_Truncated(mean, covariance, lowerLim, upperLim)
            if ~ismatrix(covariance) || size(covariance, 1) ~= size(covariance, 2)
                error('Covariance matrix must be square')
            end
            if ... 
				~isvector(mean) || ...
				length(mean) ~= size(covariance, 1) || ...
				length(mean) ~= length(lowerLim) || ...
				length(mean) ~= length(upperLim)
				
                error('Mean and limit vectors must be the same size as the covariance matrix')
            end
			if any(upperLim <= lowerLim)
				error('Upper limit must be universally greater than lower limit')
			end
            
			obj.Mean = mean;
            obj.Covariance = covariance;
			obj.LowerLim = lowerLim;
			obj.UpperLim = upperLim;
            
            % Check if the covariance matrix is all zeros, indicating a
            % fixed value
            obj.ZeroCov = all(obj.Covariance == 0);
        end
        
        function prob = PDF(obj, theta)
            if ~obj.ZeroCov
			    if any(theta < obj.LowerLim) || any(theta > obj.UpperLim)
				    prob = 0; % Probability is out of the truncated limit
			    else
				    prob = mvnpdf(theta, obj.Mean, obj.Covariance);
                end
            else
                if all(theta==obj.Mean)
                    prob = 1;
                else
                    prob = 0;
                end
            end
        end
		
		function prob = logPDF(obj, theta)
            if ~obj.ZeroCov
			    if any(theta < obj.LowerLim) || any(theta > obj.UpperLim)
				    prob = -inf; % Probability is out of the truncated limit
			    else
					prob = logmvnpdf(theta, obj.Mean, obj.Covariance);
                end
            else
                if all(theta==obj.Mean)
                    prob = 0;
                else
                    prob = -inf;
                end
            end
		end
        
        function theta = Sample(obj)
            if all(obj.ZeroCov) % If the entire covariance matrix is zero, then just return the mean value
                theta = obj.Mean;
                return
            end
            if any(obj.ZeroCov)
                error('There is a zero entry in the covariance matrix, which will lead to erratic sampling results.')
            end
            theta = mvnrnd(obj.Mean, obj.Covariance, 1);
	        if obj.PDF(theta) == 0 % If outside the truncation then sample again
			    theta = obj.Sample();
            end
        end
        
        function d = Dimension(obj)
            d = length(obj.Mean);
        end
    end
end

function [logp] = logmvnpdf(x,mu,Sigma)
	% outputs log likelihood array for observations x  where x_n ~ N(mu,Sigma)
	% x is NxD, mu is 1xD, Sigma is DxD
    warning off
    [N,D] = size(x);
	const = -0.5 * D * log(2*pi);
	xc = bsxfun(@minus,x,mu);
	term1 = -0.5 * sum((xc / Sigma) .* xc, 2); % N x 1
	term2 = const - 0.5 * logdet(Sigma);    % scalar
	logp = term1' + term2;
    warning on
end

function y = logdet(A)
	U = chol(A);
	y = 2*sum(log(diag(U)));
end
