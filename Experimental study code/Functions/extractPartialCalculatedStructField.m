% Extracts an appropriate slice of a squeezed set of dimensions, even when
% the target grid is not of the nominal grid size.
% Nominal inds are the 1d index

function rArr = extractPartialCalculatedStructField(samples, fieldName, nominalInd, nominalSize, opts)
    arguments
        samples
        fieldName string
        nominalInd
        nominalSize
        opts.WeightedVals logical = true
    end

    % Determine the true field size, including trailing singleton
    % dimensions
    fieldSize = [1 1 1 1 1];
    rawFieldSize = size(samples(1).CalculatedValues.(fieldName));
    fieldSize(1:length(rawFieldSize)) = rawFieldSize;

    [nominalInds(1), nominalInds(2), nominalInds(3), nominalInds(4), nominalInds(5)] = ...
        ind2sub(nominalSize, nominalInd);
    for i = 1:length(nominalInds)
        nominalInds(i) = inlineIfElse(fieldSize(i)~=1, nominalInds(i), 1);
    end

    rpmInd = inlineIfElse(nominalInds(1)~=0, nominalInds(1), 1:fieldSize(1));
    aeInd = inlineIfElse(nominalInds(2)~=0, nominalInds(2), 1:fieldSize(2));
    bInd = inlineIfElse(nominalInds(3)~=0, nominalInds(3), 1:fieldSize(3));
    fzInd = inlineIfElse(nominalInds(4)~=0, nominalInds(4), 1:fieldSize(4));
    cutDirInd = inlineIfElse(nominalInds(5)~=0, nominalInds(5), 1:fieldSize(5));

    % Calculate the total size of the array
    weights = zeros(size(samples));
    for i = 1:length(samples)
        if opts.WeightedVals
            if isfield(samples, 'Weight')
                weights(i) = samples(i).Weight;
            else
                weights(i) = 1;
            end
        else
            weights(i) = 1;
        end
    end
    rowCount = sum(weights);

    rArr = zeros([length(rpmInd), length(aeInd), length(bInd), length(fzInd), length(cutDirInd), rowCount]);
    currentInd = 1;
    for i = 1:length(samples)
        newVals = samples(i).CalculatedValues.(fieldName)(rpmInd, aeInd, bInd, fzInd, cutDirInd);
        rArr(:,:,:,:,:,currentInd:currentInd+weights(i)-1) = repmat(newVals, [1 1 1 1 1 weights(i)]);
        currentInd = currentInd + weights(i);
        % rArr = cat(6, rArr, repmat(newVals, [1 1 1 1 1 weight]));
    end

    function val = inlineIfElse(condition, trueVal, falseVal)
        if condition
            val = trueVal;
        else
            val = falseVal;
        end
    end
end
