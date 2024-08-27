function CheckForFields(S, requiredFields, callingFuncName)
    throwError = false;
    missingFields = [];
    for i = 1:length(requiredFields)
        if ~isfield(S, requiredFields(i))
            throwError = true;
            missingFields(end+1) = requiredFields(i);
        end
    end

    if throwError
        disp(strcat('The following calculated fields are missing for function ', callingFuncName, ':'))
        for i = 1:length(missingFields)
            disp(missingFields(i))
        end
    end
end

