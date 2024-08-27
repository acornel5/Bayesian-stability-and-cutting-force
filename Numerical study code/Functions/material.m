classdef material
    %MATERIAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Ks          % specific cutting force, N/m^2
        Beta        % force angle, deg
        Ktc         % tangential cutting coefficient, N/m^2 
        Krc         % radial cutting coefficient, N/m^2 
        Kte = 0     % tangential edge constant, N/m
        Kre = 0     % radial edge constant, N/m
        Ct  = 0     % tangential process damping coefficient, N s/m
        Cr  = 0     % radial process damping coefficient, N s/m
    end
    
    methods
        function obj = material(ks, beta)
            %MATERIAL Construct an instance of this class
            obj.Ks = ks;
            obj.Beta = beta;
            obj.Ktc = ks*sin(beta);
            obj.Krc = ks*cos(beta);
        end
        
        function [ktc, krc, kte, kre, ct, cr] = MaterialMatrix(obj)
            % Outputs all the parameters in a matrix for easy 1 line
            % assignment
            ktc = obj.Ktc;
            krc = obj.Krc;
            kte = obj.Kte;
            kre = obj.Kre;
            cr  = obj.Cr;
            ct  = obj.Ct;
        end
    end
end

