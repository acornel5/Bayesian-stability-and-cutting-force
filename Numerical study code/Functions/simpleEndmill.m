classdef simpleEndmill
    % A simple endmill with constant flute spacing and helix angle, and no runout
    
    properties
        Diameter % Tool diameter in meters
        FluteCount % Number of flutes
        HelixAngle % Helix angle
    end
    
    methods
        function obj = simpleEndmill(diameter, fluteCount, helixAngle)
            % Constructor
            obj.Diameter = diameter;
            obj.FluteCount = fluteCount;
            obj.HelixAngle = helixAngle;
        end
        
        function helixVector = HelixVector(obj)
            helixVector = ones(obj.FluteCount, 1) * obj.HelixAngle; 
        end
        
        function fluteAngles = FluteAngleVector(obj)
            initialVec = (linspace(0, 360, 1+obj.FluteCount));
            fluteAngles = initialVec(1:end-1);
        end
    end
end

