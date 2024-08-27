function output = physicsMeanCuttingTorque(inputVals, parser, param, options)
    arguments
        inputVals
        parser
        param
        options.EdgeCoefficientFactor = 0
    end
    
    % Get the cutting force
    ktcKnc = parser.forceKtcKnc(inputVals);
    ktc = ktcKnc(1);
    % kte = options.EdgeCoefficientFactor * ktc;

    try
        kte = parser.returnField(inputVals, 'kte');
    catch
        kte = 0;
    end
    
    engagement = acos((param.ToolDiameter/2 - param.CutParamGrid.a)/(param.ToolDiameter/2)); % exit angle, rad
    % engagementFrac = engagement / (2*pi);



    % output.SpindlePower = param.CutParamGrid.a .* param.CutParamGrid.fz .* param.CutParamGrid.b .* param.CutParamGrid.rpm .* param.FluteCount .* ktc/60;
    cuttingPower = param.CutParamGrid.a .* param.CutParamGrid.fz .* param.CutParamGrid.b .* param.CutParamGrid.rpm .* param.FluteCount .* ktc/60;
    edgePower = param.ToolDiameter * param.CutParamGrid.rpm / 120.* param.CutParamGrid.b .* kte .* engagement .* param.FluteCount;
    output.SpindlePower = cuttingPower + edgePower;
    output.SpindleTorque = output.SpindlePower ./ param.CutParamGrid.rpm / 2 / pi * 60;

    %output.SpindleTorque = single(param.FluteCount .* param.CutParamGrid.b .* param.CutParamGrid.fz .* ktc .* (param.ToolDiameter .* acos(1-2*param.CutParamGrid.a/param.ToolDiameter))/(4*pi));
    %output.SpindlePower = param.CutParamGrid.rpm.* output.SpindleTorque * 2 * pi / 60;
end