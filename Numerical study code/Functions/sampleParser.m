classdef sampleParser
    % The sample parser assists with extracting specific fields from the
    % theta vector, since elements of the vector may be in different
    % positions depending on how exactly the vector is constructed. For
    % example, a function can request a named value which is returned
    % regardless of its position in the vector, or the parser can convert
    % between different forms of the force model.
    
    properties
        forceModelForm
        dynamicsModelForm
        forceModelRange
        dynamicsModelRange
        fieldNames
    end
    
    methods
        %% Constructor
        function obj = sampleParser(forceModel, forceModelRange, dynamicsModel, dynamicsModelRange, fieldNames)
            obj.forceModelForm = forceModel;
            obj.forceModelRange = forceModelRange;
            obj.dynamicsModelForm = dynamicsModel;
            obj.dynamicsModelRange = dynamicsModelRange;
            obj.fieldNames = fieldNames;
        end

        %% Force model parsing
        % Returns the 2 element force model [Ks beta]
        function vals = forceKsBeta(obj, sample)
            sixElementModel = obj.force6Element(sample);
            vals = [sqrt(sum(sixElementModel([1 2]).^2)), atan(sixElementModel(1)/sixElementModel(2))];
        end

        % Returns the 2 element force model [ktc knc]
        function vals = forceKtcKnc(obj, sample)
            sixElementModel = obj.force6Element(sample);
            vals = sixElementModel([1, 2]);
        end

        % Returns the 4 element force model [ktc knc kte kne]
        function vals = force4Element(obj, sample)
            sixElementModel = obj.force6Element(sample);
            vals = sixElementModel([1:2, 4:5]);
        end

        % Returns the 6 element force model [ktc knc kac kte kne kae]
        function vals = force6Element(obj, sample)
            switch obj.forceModelForm
                case 'KsBeta'
                    modelVals = sample(obj.forceModelRange(1):obj.forceModelRange(2));
                    vals = [modelVals(1)*sin(modelVals(2))  modelVals(1)*cos(modelVals(2)) 0 0 0 0];
                case 'KtcKnc'
                    modelVals = sample(obj.forceModelRange(1):obj.forceModelRange(2));
                    vals = [modelVals(1) 0 0 modelVals(2) 0 0];
                case '4Element'
                    modelVals = sample(obj.forceModelRange(1):obj.forceModelRange(2));
                    vals = [modelVals(1:2) 0 modelVals(3:4) 0];
                case '6Element'
                    vals = sample(obj.forceModelRange(1):obj.forceModelRange(2));
                otherwise, error('Invalid force model type')
            end
        end

        %% Dynamics model parsing
        % Return the dynamics as a matrix of modes. Each row contains one
        % natural frequency, stiffness, and damping ratio value
        function vals = dynamicsFnKZeta(obj, sample)
            switch obj.dynamicsModelForm
                case 'fnKZeta'
                    vals = sample(obj.dynamicsModelRange(1):obj.dynamicsModelRange(2));
                    assert(mod(length(vals),3)==0)
                    vals = reshape(vals, [3, numel(vals)/3])';
                case 'kmc'
                    tempVals = sample(obj.dynamicsModelRange(1):obj.dynamicsModelRange(2));
                    assert(mod(length(tempVals),3)==0)
                    tempVals = reshape(tempVals, [3, numel(tempVals)/3])';
                    k = tempVals(:,1);
                    m = tempVals(:,2);
                    c = tempVals(:,3);

                    wn = sqrt(k/m);
                    fn = wn/(2*pi);
                    ccrit = 2*sqrt(k*m);
                    zeta = c/ccrit;
                    vals = [fn k zeta];
                otherwise, error('Invalid dynamics model type')
            end
        end

        % Return the dynamics as a matrix of modes. Each row contains one
        % stiffness, mass, and viscous damping coefficient
        function vals = dynamicsKMC(obj, sample)
            switch obj.dynamicsModelForm
                case 'fnKZeta'
                    tempVals = sample(obj.dynamicsModelRange(1):obj.dynamicsModelRange(2));
                    assert(mod(length(tempVals),3)==0)
                    tempVals = reshape(tempVals, [3, numel(tempVals)/3])';
                    fn = tempVals(:,1);
                    k = tempVals(:,2);
                    zeta = tempVals(:,3);

                    wn = fn * 2 * pi;
                    m = k ./ wn.^2;
                    ccrit = 2*sqrt(k.*m);
                    c = ccrit .* zeta;
                    vals = [k m c];
                case 'KMC'
                    vals = sample(obj.dynamicsModelRange(1):obj.dynamicsModelRange(2));
                    assert(mod(length(vals),3)==0)
                    vals = reshape(vals, [3, numel(vals)/3])';
                otherwise, error('Invalid dynamics model type')
            end
        end

        function vals = rcsaCouplingParam(obj, sample)
            switch obj.dynamicsModelForm
                case 'RCSA_KtfCtf'
                    vals = [0 0 sample(obj.dynamicsModelRange(1)) 0 0 sample(obj.dynamicsModelRange(2))];
                otherwise, error('Invalid dynamics model type')
            end
        end

        % Extract a single named field from the sample 
        function val = returnField(obj, sample, fieldName)
            indices = ismember(obj.fieldNames, fieldName);

            % Error catching
            if sum(indices) > 1, error(strcat("Invalid value. The named value '", fieldName, "' appears multiple times in the theta vector.")), end
            if ~any(indices), error(strcat("Invalid value. The named value '", fieldName, "' is not defined in the theta vector.")), end
            
            val = sample(indices);
        end

        % Find the index for a specific field name
        function ind = GetFieldIndex(obj, fieldName)
            indices = ismember(obj.fieldNames, fieldName);
            % Error catching
            if sum(indices) > 1, error(strcat("Invalid value. The named value '", fieldName, "' appears multiple times in the theta vector.")), end
            if ~any(indices), error(strcat("Invalid value. The named value '", fieldName, "' is not defined in the theta vector.")), end
            ind = find(indices);
        end

    end
end