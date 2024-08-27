classdef dynamicModel
    % Defines the X and Y modal parameters of a tool or workpiece
    
    properties
        XDynamicMatrix % Each row should be a mode, with stiffness, mass, and damping in the columns respectively
        YDynamicMatrix
    end
    
	methods(Static)
		function obj = GetRigidModel()
			obj = dynamicModel('Rigid', 'Rigid');
        end
        
        % Get a new s
        function obj = FromFnKZeta(xMatrix, yMatrix)
            if strcmp(xMatrix, 'Rigid')
                xMatrix = [realmax realmax realmax];
            else
                xMatrix = reshape(xMatrix, [], 3);
                fnx = xMatrix(:,1);
                wnx = fnx * (2*pi);
                kx = xMatrix(:,2);
                zetax = xMatrix(:,3);
                mx = kx ./ wnx.^2;
                cx = zetax .* 2 .* sqrt(kx.*mx);
                obj.XDynamicMatrix = reshape(xMatrix, [], 3);
                xMatrix = [kx mx cx];
            end
            
            if strcmp(yMatrix, 'Rigid')
                yMatrix = [realmax realmax realmax];
            elseif strcmp(yMatrix, 'Symmetric')
                yMatrix = xMatrix;
            else
                yMatrix = reshape(yMatrix, [], 3);
                fny = yMatrix(:,1);
                wny = fny * (2*pi);
                ky = yMatrix(:,2);
                zetay = yMatrix(:,3);
                my = ky ./ wny.^2;
                cy = zetay .* 2 .* sqrt(ky.*my);
                obj.XDynamicMatrix = reshape(yMatrix, [], 3);
                yMatrix = [ky my cy];
            end
            
            obj = dynamicModel(xMatrix, yMatrix);
        end
	end
	
    methods
        function obj = dynamicModel(xMatrix, yMatrix)
            if strcmp(xMatrix, 'Rigid')
                obj.XDynamicMatrix = [realmax realmax realmax];
            else
                obj.XDynamicMatrix = reshape(xMatrix, [], 3);
            end
            if strcmp(yMatrix, 'Rigid')
                obj.YDynamicMatrix = [realmax realmax realmax];
            elseif strcmp(yMatrix, 'Symmetric')
                obj.YDynamicMatrix = obj.XDynamicMatrix;
            else
                obj.YDynamicMatrix = reshape(yMatrix, [], 3);
            end
            %DYNAMICMODEL Construct an instance of this class
            %   Detailed explanation goes here
        end
       
        
        function [k, m, c, w, n] = XParameters(obj)
            % Extract the X axis vibrations in one output
            k = obj.XDynamicMatrix(:, 1);
            m = obj.XDynamicMatrix(:, 2);
            c = obj.XDynamicMatrix(:, 3);
            w = (k ./ m) .^ 0.5;
            n = length(k);
        end
        
        function [k, m, c, w, n] = YParameters(obj)
            % Extract the X axis vibrations in one output
            k = obj.YDynamicMatrix(:, 1);
            m = obj.YDynamicMatrix(:, 2);
            c = obj.YDynamicMatrix(:, 3);
            w = (k ./ m) .^ 0.5;
            n = length(k);
        end

    end
end

