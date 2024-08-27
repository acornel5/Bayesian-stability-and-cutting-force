% getPartialField.m
% 2023-11-2

% Some fields may not be calculated for every set of cutting parameters in
% the full grid: for example, chatter frequency is currently only
% saved for one axial depth. This takes a nominal full-grid index and
% returns the appropriate value from the flattened grid.

function rVal = getPartialField(array, nominalInd, nominalGridSize)
    [ind(1), ind(2), ind(3), ind(4), ind(5)] = ind2sub(nominalGridSize, nominalInd);

    partialGridSize = size(array);
    while length(partialGridSize) < 5
        partialGridSize = [partialGridSize 1];
    end
    ind(partialGridSize==1) = 1;

    rVal = array(sub2ind(partialGridSize, ind(1), ind(2), ind(3), ind(4), ind(5)));
end
