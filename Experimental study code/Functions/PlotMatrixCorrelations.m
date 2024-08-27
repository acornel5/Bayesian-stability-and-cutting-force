function PlotMatrixCorrelations(iteration, samples, nominalVals, testResults, param)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'defaultLegendInterpreter','latex');
    set(groot, 'defaultTextInterpreter', 'latex');


    % TODO: Setting these limits manually is a nasty hack
    if isfield(param, 'plotLim')
        limits = param.plotLim;
        tickLabels = [
            limits(:,1), 0.5*(limits(:,1)+limits(:,2)), limits(:,2)];
    end

	
    data = AssembleSampledVals(samples);
    N = size(data,2);

    % Scale stiffness. Nasty hack, need to rework this
	% Scale data to +-100% of measured value
% 	for i = 1:N
% 		data(:,i) = (data(:,i) - trueVals(i)) ./ trueVals(i);
% 	end
    
    % Plot data
    fig = figure('visible', 'off');
    clf
    [S,Ax,BigAx,H,HAx] = plotmatrix(data);

    trueVals = [1650 68 50 1254 3.33e6 4.34 1745 1.26e7 2.53]

    % Set marker size
    for i=1:N
        for j=1:N
            if i==j
                continue;
            end
            S(i,j).MarkerSize = 1;
        end
    end

    % Plot true values and correlations, set axis limits
    if ~isempty(nominalVals.Values)
        % Scale stiffness. Nasty hack, need to rework this
%         trueVals(1:3) = trueVals(1:3)*1e-6;
        
        % Diagonal
        for i = 1:N
            hold(HAx(i), 'on');
%             plot(HAx(i),[0 0], HAx(i).YLim, 'Color', 'red')
            if exist('limits','var') == 1
                HAx(i).XLim = [limits(i,1) limits(i,2)];
            end
        end

        % Non-diagonal
        corrVals = corr(data);
        for i = 1:N
            for j = 1:N
                if i==j || j>=i % Skip diagonals
                    continue;
                end
                
                % True values
                hold(Ax(i,j), 'on')
%                 plot(Ax(i,j), 0, 0, 'o', 'Color', 'red')

                % Plot orthogonal fit
%                 [~, orthFitFunc] = OrthogonalLinearFit(data(:, [j i]));
%                 fplot(Ax(i,j), orthFitFunc, limits(j,:), 'color', 'blue')
                
                % Axis limits
                if exist('limits','var') == 1
                    Ax(i,j).XLim = [limits(j,1) limits(j,2)];
                    Ax(i,j).YLim = [limits(i,1) limits(i,2)];
                end
            end
        end
    end

    
    for i = 1:N
%         Ax(i,1).YTick = [-0.5 0 0.5];
%         Ax(i,1).YTickLabels = ["-50\%" "0\%" "50\%"];
        Ax(i,1).FontSize = 10;
        
%         Ax(N,i).XTick = [-0.5 0 0.5];
%         Ax(N,i).XTickLabels = ["-50\%" "0\%" "50\%"];
        Ax(N,i).FontSize = 10;
    end

    
for i = 1:N
    Ax(i,i).XTickLabel = {};
    Ax(i,i).YTickLabel = {};
end

    Ax(1,1).YTick = [];
    Ax(1,1).YTickLabels = {};
    Ax(1,1).XTick = [];
    HAx(1).YTick = [];
    HAx(1).YTickLabels = {};
    % HAx(N).XTick = tickLabels(N,:);
%     HAx(N).XTick = [-0.5 0 0.5];
%     HAx(N).XTickLabels = ["-50\%" "0\%" "50\%"];        
    HAx(N).FontSize = 10;

%     Ax(N,N).XTickLabels = [" " " " " "];


    % Clear values above the diagonal
    for i = 1:N
        for j = i+1:N
            delete(Ax(i,j))
        end
    end


    % Set axis labels
    for i = 1:N
%         Ax(i,1).YLabel.Position(1) = -100;
%         HAx(i).YLabel.Position(1) = -100;
        Ax(i,1).YLabel.String = param.sampleVariableNames{i};
        % Ax(i,1).YTick = tickLabels(i,:);
        
        if i~=N
            Ax(N,i).XLabel.String =  param.sampleVariableNames{i};
            % Ax(N,i).XTick = tickLabels(i,:);
        end
    end

    HAx(N).XLabel.String = param.sampleVariableNames{N};
    
    % For some reason the first YLabel needs to be 50pt further right
    Ax(1,1).YLabel.Position(1) = Ax(1,1).YLabel.Position(1) - 2;
%     Ax(1,1).YLabel.String = '$k_{tc}$';
    % Create legend
    hold(BigAx, 'on');
    plot(BigAx, NaN, NaN, '-o', 'color', 'red', 'MarkerEdgeColor', 'red')
    plot(BigAx, NaN, NaN, 'color', 'blue')
    plot(BigAx, NaN, NaN, 'color', [255 130 0]/255)
%     legend(BigAx, ["Measured values" "Orthogonal best fit"], 'location', 'northeast')

    % Save figure
    folderName = 'Figures\\Histogram';
    if not(isfolder(folderName)), mkdir(folderName), end
    fileName = [num2str(iteration)];
    print(fig, '-dpng', strcat(folderName, '\\', fileName), '-r250');

    close(fig);
end

function matrix = AssembleSampledVals(samples)
    matrix = zeros(sum([samples.Weight]), length(samples(1).Values));
    sampleCounter = 1;
    for i = 1:length(samples)
        if isempty(samples(i).Weight), continue, end

        matrix(sampleCounter:sampleCounter+samples(i).Weight-1,:) = repmat(samples(i).Values, samples(i).Weight, 1);
        sampleCounter = sampleCounter+samples(i).Weight;
    end
end

% function matrix = AssembleSampledVals(samples)
%     for i = 1:length(samples)
%         matrix(i,:) = samples(i).Values;
%     end
% end