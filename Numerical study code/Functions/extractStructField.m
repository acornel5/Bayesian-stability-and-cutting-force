% extractStructField.m
% 2021-11-18
% A. Cornelius

% Takes a list of structs and extracts a specific field from all of them,
% converting it into an array. The field must have the same dimensions on
% all of them.

function vars = extractStructField(arr, fieldName)
    sampleVals = eval(['arr(1).', fieldName]);
    vars = zeros(length(arr), length(sampleVals));
    
    for i = 1:length(arr)
        vars(i,:) = eval(['arr(', num2str(i), ').', fieldName]);
    end
end

