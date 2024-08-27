function rStruct = physicsFuncFrfFromSymmetricModes(inputVals, parser, options)
    arguments
        inputVals double
        parser
        options.FrequencyRange = 'auto'
    end

    % Set the frequency range
    if strcmp(options.FrequencyRange, 'auto') % Automatically span the natural frequency range
        modesFnKZeta = parser.dynamicsFnKZeta(inputVals);
        minFn = min(modesFnKZeta(:,1));
        maxFn = max(modesFnKZeta(:,1));
        f = minFn/2:0.5:maxFn*2;
    else % Manually specified natural frequencies
        f = options.FrequencyRange;
    end

    modesKMC = parser.dynamicsKMC(inputVals);

    w = f*2*pi;
    frfxx = zeros(size(f));
    for i = 1:size(modesKMC,1)
        frfxx = frfxx + 1 ./ (-w.^2 * modesKMC(i,2) + 1i * w * modesKMC(i,3) + modesKMC(i,1));
    end
    frfyy = frfxx;

    rStruct.FRFxx = frfxx;
    rStruct.FRFyy = frfyy;
    rStruct.FrequenciesFRF = f;
end