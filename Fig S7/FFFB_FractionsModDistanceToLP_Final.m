%% Fig S7: fraction of modulated cells vs distance from looped cells

saveDir_FigS7 = fullfile(Params.MainDir, 'plot', 'Fig S7');
if ~exist(saveDir_FigS7, 'dir')
    mkdir(saveDir_FigS7)
end

%% Settings
bin_size    = 10;                 % width of distance bins (um)
nBins       = 30;
numShuffles = 1000;
% rng(1);                        % uncomment for a reproducible shaded CI

FractionModBinAll = cell(1, 2);
CI_low_all        = cell(1, 2);
CI_high_all       = cell(1, 2);

for F = 1:2

    %% Collect min distance to a looped cell + modulation flag for every non-looped cell
    allMinDist = [];   % column: min distance to any looped cell
    allIsMod   = [];   % column: cell is significantly excited or inhibited

    for ani = 1:Params.n_ani{F}
        for pos = 1:Params.n_pos{F}(ani)

            if isempty(data{F}{ani,pos}), continue; end
            stat = CellStat{F}{ani,pos};
            if ~any(stat.Loop), continue; end

            % Keep responsive or looped cells only
            keep   = stat.responsiveAD | stat.Loop;
            xy     = stat.Coordyx_um(keep, :);
            isLoop = logical(stat.Loop(keep));
            isLoop = isLoop(:);
            isMod  = logical(stat.SigExcAD | stat.SigInhAD);
            isMod  = isMod(keep);
            isMod  = isMod(:);

            % Distance of each non-looped cell to each looped cell
            loopXY  = xy(isLoop, :);
            otherXY = xy(~isLoop, :);
            dx = otherXY(:,1) - loopXY(:,1)';
            dy = otherXY(:,2) - loopXY(:,2)';
            minDist = min(sqrt(dx.^2 + dy.^2), [], 2);

            allMinDist = [allMinDist; minDist];
            allIsMod   = [allIsMod; isMod(~isLoop)];
        end
    end

    %% Fraction of modulated cells per distance bin
    % First bin is 0-20 um, then bin_size-wide bins, extending to cover ALL cells
    binEdges = [0, 20 : bin_size : max(allMinDist) + bin_size];
    nAll     = numel(binEdges) - 1;
    binIdx   = discretize(allMinDist, binEdges);
    isMod    = double(allIsMod);

    % NaN for empty bins; keep the first nBins for the figure
    fracAll = accumarray(binIdx, isMod, [nAll 1], @mean, NaN)';
    FractionModBinAll{F} = fracAll(1:nBins);

    %% Null distribution: shuffle modulation labels across bins
    % (pool includes all cells, as in the original script)
    shuffFrac = zeros(numShuffles, nBins);
    for s = 1:numShuffles
        shuffled = isMod(randperm(numel(isMod)));
        fracAll  = accumarray(binIdx, shuffled, [nAll 1], @mean, NaN)';
        shuffFrac(s,:) = fracAll(1:nBins);
    end

    % 95% confidence interval of the shuffled fractions
    CI = prctile(shuffFrac, [2.5 97.5]);
    CI_low_all{F}  = CI(1,:);
    CI_high_all{F} = CI(2,:);
end


  %% Plot for figure  

figure;
for F=1:2
    prange=[1:19];
    subplot(1,2,3-F)
    xvect=[bin_size*2:bin_size:bin_size*size(FractionModBinAll{F},2)-1];
    
    % Plot the shaded confidence interval
    s1= fill([xvect(prange) fliplr(xvect(prange))], ...
        [CI_low_all{F}(1,prange) fliplr(CI_high_all{F}(1,prange))], ...
        [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'none');  hold on
    
    p1=plot(xvect(prange), FractionModBinAll{F}(prange),'-o',  'LineWidth', 1, ...
        'MarkerSize', 2, 'MarkerEdgeColor', Params.colors{F}, 'MarkerFaceColor', Params.colors{F},'Color', Params.colors{F});
    
    legend([p1,s1],{'Real data','95% CI from shuffled bins'},'location','southwest'); legend boxoff
if F==2
    ylabel('Fraction of Modulated cells');
    end
    xlabel('Minimum distance to a looped cell (um)');

xlim([10 205])
    title([Params.ExpName{F}])
    box off
  set(gca,'FontSize',4.6,'LineWidth',0.6,'TickDir','out','TickLength',[0.01, 0.01]);
% despline(gca);
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';  
    
end
set(gcf,'Units','centimeters','Position',[10 10 9 3.5]);

PlotTitle='Fraction of modulated cells vs Distance from looped cells vs shuffled bins  - All Cells' ;
basePath = fullfile(saveDir_FigS7,[PlotTitle]);
saveFigAll(gcf, basePath);