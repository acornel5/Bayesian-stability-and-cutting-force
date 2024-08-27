% GenericPhysicsDelegate.m
% 2023-8-30
% 
% This function calls a series of other functions on an input sample

function output = GenericPhysicsDelegate(sample, functions, desiredFields)
    output = struct('Values', sample);

    % Invoke all of the sub-functions
    for i = 1:length(functions)
        % Invoke the physics function. If the function only takes one
        % input, then it will only get the sample values as an input. If it
        % takes two inputs, it will also get the previously calculated
        % values (e.g. the stability map)
        switch nargin(functions{i})
            case 1, subOutput = functions{i}(sample);
            case 2, subOutput = functions{i}(sample, output);
            otherwise, error(strcat("Physics function requested an invalid number of inputs"))
        end
        
        % Append the returned fields to the output
        fields = fieldnames(subOutput);
        for j = 1:length(fields)
            output.(fields{j}) = subOutput.(fields{j});
        end
    end

    % Get rid of all non-desired fields
    fields = fieldnames(output);
    for i = 1:length(fields)
        if ~any(strcmp(fields{i}, desiredFields))
            output = rmfield(output, fields{i});
        end
    end
end