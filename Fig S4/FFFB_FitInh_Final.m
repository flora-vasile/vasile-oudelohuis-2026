saveDir_FigS4 = fullfile(Params.MainDir, 'plot', 'Fig 4');
if ~exist(saveDir_FigS4,'dir')
    mkdir(saveDir_FigS4)
end

minRespRange = 0.04;

for F = 1:2

    %% Build per-session stimulus-response arrays
    for ani = 1:Params.n_ani{F}
        for pos = 1:Params.n_pos{F}(ani)
            if isempty(data{F}{ani,pos})
                continue
            end

            raw = data{F}{ani,pos}.Raw;              % [cell x time x stim x dir x laser x rep]

            % Visual-only stimuli (conditions 2-9) and photostim-only condition (1,1)
            visTrials{F}{ani,pos}   = raw(:,:,2:9,:,:,:);
            photoTrials{F}{ani,pos} = squeeze(raw(:,:,1,1,:,:));

            % Collapse stim-type x direction into a single "condition" dimension
            [nCells, nTime, nStim, nDir, nLaser, nRep] = size(visTrials{F}{ani,pos});
            visTrialsCat{F}{ani,pos} = reshape(visTrials{F}{ani,pos}, ...
                nCells, nTime, nStim*nDir, nLaser, nRep);

            % Prepend the photostim-only condition as condition 1
            allTrialsCat{F}{ani,pos}(:,:,1,:,:)   = photoTrials{F}{ani,pos};
            allTrialsCat{F}{ani,pos}(:,:,2:17,:,:) = visTrialsCat{F}{ani,pos};

            % Average over the stimulus-on window
            condResp{F}{ani,pos}    = squeeze(nanmean(allTrialsCat{F}{ani,pos}(:,Params.StimONResponse,:,:,:), 2));
            visCondResp{F}{ani,pos} = squeeze(nanmean(visTrialsCat{F}{ani,pos}(:,Params.StimONResponse,:,:,:), 2));
            condMean{F}{ani,pos}    = squeeze(nanmean(condResp{F}{ani,pos}, 4));

            stimResp{F}{ani,pos} = squeeze(nanmean(visTrials{F}{ani,pos}(:,Params.StimONResponse,:,:,:,:), 2));

            %% Find directions with a significant Laser-ON response
            visTrialsCatOn{F}{ani,pos} = squeeze(visTrialsCat{F}{ani,pos}(:,:,:,2,:));

            base   = squeeze(mean(visTrialsCatOn{F}{ani,pos}(:,1:15,:,:), 2));
            stimOn = squeeze(mean(visTrialsCatOn{F}{ani,pos}(:,Params.StimONResponse,:,:), 2));
            baseMean   = squeeze(mean(base,3));
            stimOnMean = squeeze(mean(stimOn,3));

            dirResponsive{F}{ani,pos} = false(size(visTrialsCatOn{F}{ani,pos},1), size(visTrialsCatOn{F}{ani,pos},3));

                        alphaVal=0.01;
            for ROI = 1:size(visTrialsCatOn{F}{ani,pos},1)
                for Dir = 1:size(visTrialsCatOn{F}{ani,pos},3)
                    h = ttest(squeeze(base(ROI,Dir,:)), squeeze(stimOn(ROI,Dir,:)), alphaVal);
                    if h & stimOnMean(ROI,Dir) > baseMean(ROI,Dir)
                        dirResponsive{F}{ani,pos}(ROI,Dir) = h;
                    end
                end
            end

            %% Find inhibited cells using stimuli responses that are non-zero
            visCondRespOn{F}{ani,pos} = NaN(size(visCondResp{F}{ani,pos}));

            for ROI = 1:size(visCondRespOn{F}{ani,pos},1)
                for Dir = 1:size(visCondRespOn{F}{ani,pos},2)
                    if dirResponsive{F}{ani,pos}(ROI,Dir)
                        visCondRespOn{F}{ani,pos}(ROI,Dir,:,:) = visCondResp{F}{ani,pos}(ROI,Dir,:,:);
                    end
                end
            end

            inhibitedCell{F}{ani,pos} = false(size(visTrialsCatOn{F}{ani,pos},1), 1);

            permuted{F}{ani,pos}    = permute(visCondRespOn{F}{ani,pos}, [1 3 2 4]);
            reshapedResp{F}{ani,pos} = reshape(permuted{F}{ani,pos}, size(permuted{F}{ani,pos},1), size(permuted{F}{ani,pos},2), []);

            for ROI = 1:size(reshapedResp{F}{ani,pos},1)
                [h, p] = ttest2(reshapedResp{F}{ani,pos}(ROI,1,:), reshapedResp{F}{ani,pos}(ROI,2,:), 'Alpha', alphaVal);
                if ~isnan(nanmean(reshapedResp{F}{ani,pos}(ROI,2,:),3)) && h && ...
                        nanmean(reshapedResp{F}{ani,pos}(ROI,2,:),3) < nanmean(reshapedResp{F}{ani,pos}(ROI,1,:),3)
                    inhibitedCell{F}{ani,pos}(ROI) = true;
                end
            end

            % Mask condition-mean responses to keep only ON-responsive directions
            % (condition 1, the photostim-only condition, is always kept)
            condMeanOn{F}{ani,pos} = NaN(size(condMean{F}{ani,pos}));
            condMeanOn{F}{ani,pos}(:,1,:) = condMean{F}{ani,pos}(:,1,:);

            for ROI = 1:size(condMeanOn{F}{ani,pos},1)
                for Dir = 2:size(condMeanOn{F}{ani,pos},2)
                    if dirResponsive{F}{ani,pos}(ROI,Dir-1)
                        condMeanOn{F}{ani,pos}(ROI,Dir,:) = condMean{F}{ani,pos}(ROI,Dir,:);
                    end
                end
            end

            cellAniID{F}{ani,pos}   = ani * ones(numel(inhibitedCell{F}{ani,pos}), 1);   % animal label per cell
cellPosID{F}{ani,pos}  = pos * ones(numel(inhibitedCell{F}{ani,pos}), 1);   % session (pos) label per cell
            
        end
    end

    %% Pool sessions/animals together
    condMeanAll{F} = cat(1, condMeanOn{F}{:});
    inhibited{F}   = logical(cat(1, inhibitedCell{F}{:}));
    
        
cellAniAll{F} = cat(1, cellAniID{F}{:});
cellPosAll{F} = cat(1, cellPosID{F}{:});
cellAniinh{F} = cellAniAll{F}(inhibited{F});
cellPosinh{F} = cellPosAll{F}(inhibited{F});

    inhibitedResp{F} = condMeanAll{F}(inhibited{F},:,:);

    stimRespAll{F} = cat(1, stimResp{F}{:});
    stimRespInh{F} = stimRespAll{F}(inhibited{F},:,:,:,:);

    % Flag cells whose Laser-OFF response range inheeds threshold
    rangeOK{F} = range(condMeanAll{F}(:,:,1), 2) > minRespRange;


    %% Linear regression Inhited cells

    nInhibited  = size(inhibitedResp{F}, 1);
    nConditions = size(inhibitedResp{F}, 2);

    % Preallocate with NaN so skipped neurons are marked
    intercept{F} = nan(nInhibited,1);
    slope{F} = nan(nInhibited,1);

    for n = 1:nInhibited
        x = squeeze(inhibitedResp{F}(n,:,1));  % Laser OFF
        y = squeeze(inhibitedResp{F}(n,:,2));  % Laser ON

        % Skip if fewer than 3 non-NaN values in y
        if sum(~isnan(y)) < 3
            continue
        end


        if range(x) < minRespRange | max(x) < minRespRange
            continue  % Skip neuron
        end

        validIdx = ~isnan(y);

        % Fit
        p = polyfit(x(validIdx), y(validIdx), 1);
        yfit = polyval(p, x(validIdx));

        slope{F}(n) = p(1);
        intercept{F}(n) = p(2);
    end

    valid{F} = ~isnan(slope{F}); %&  slope{F}>=1;

end

%% Number of mice and sessions contributing valid inhibited cells
for F = 1:2
    aniValid  = cellAniinh{F}(valid{F});
    posValid  = cellPosinh{F}(valid{F});

    nMice     = numel(unique(aniValid));
    nSess     = size(unique([aniValid posValid], 'rows'), 1);   % unique (ani,pos) pairs

    nMiceinh  = numel(unique(cellAniinh{F}));
    nSessinh  = size(unique([cellAniinh{F} cellPosinh{F}], 'rows'), 1);

    fprintf(['%s: %d mice, %d sessions with inhibited cells passing range threshold ' ...
             '(%d mice, %d sessions with any inhibited cell); n = %d cells\n'], ...
        Params.ExpName{F}, nMice, nSess, nMiceinh, nSessinh, sum(valid{F}));
end


%% Plot cumulative distribution

% Slope statistics
if kstest(slope{1}(valid{1}))|kstest(slope{2}(valid{2}))
    pSlope = ranksum(slope{1}(valid{1}),slope{2}(valid{2}));
    testSlope='Wilcoxon';
else
    [~,pSlope]=ttest2(slope{1}(valid{1}),slope{2}(valid{2}));
    testSlope='ttest';
end

% Intercept statistics
if kstest(intercept{1}(valid{1}))|kstest(intercept{2}(valid{2}))
    pIntercept = ranksum(intercept{1}(valid{1}),intercept{2}(valid{2}));
    testIntercept='Wilcoxon';
else
    [~,pIntercept]=ttest2(intercept{1}(valid{1}),intercept{2}(valid{2}));
    testIntercept='ttest';
end

for F=1:2
    [intCdfY{F}, intCdfX{F}] = ecdf(intercept{F}(valid{F}));
    
    % cap at 0.3
    % Find where x <= 0.3
    idx = intCdfX{F} <= 0.2;
    
    % Keep only up to 0.3
    intCdfX{F} = intCdfX{F}(idx);
    intCdfY{F} = intCdfY{F}(idx);
    
    % Append final point at x = 0.3, y = 1
    intCdfX{F}(end+1) = 0.2;
    intCdfY{F}(end+1) = 1;
    
    [slopeCdfY{F}, slopeCdfX{F}] = ecdf(slope{F}(valid{F}));
end

%% Example fit
F = 1;
n = 153;

x = squeeze(inhibitedResp{F}(n,:,1));
y = squeeze(inhibitedResp{F}(n,:,2));

validIdx = ~isnan(y);

p = polyfit(x(validIdx), y(validIdx), 1);
yfit = polyval(p, x(validIdx));

figure;
scatter(x(validIdx), y(validIdx),3, 'filled','MarkerFaceColor',[0.4 0.4 0.4]); hold on;
plot(x(validIdx), yfit, 'r-', 'LineWidth', 1,'color',Params.colors{F});
title(sprintf('Neuron %d', n));
% Add slope and offset text at bottom center
txt = sprintf('Slope = %.2f\nOffset = %.2f', p(1), p(2));
text(0.7, 0.06, txt, ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 3, ...
    'Color', 'k');

xlabel('Laser OFF'); ylabel('Laser ON');
xLimits = xlim; yLimits = ylim;
xlim([min(xLimits(1),yLimits(1)) max(xLimits(2),yLimits(2))]); ylim([min(xLimits(1),yLimits(1)) max(xLimits(2),yLimits(2))]);
plot([0 max(xLimits(2),yLimits(2))], [0 max(xLimits(2),yLimits(2))],'color',[0.5 0.5 0.5],'LineStyle', ':');
axis equal
xlim([0 0.15]); ylim([0 0.15])
set(gca,'FontSize',4.6,'LineWidth',0.6,'TickDir','out','TickLength',[0.02, 0.01]);
% despline(gca,0.1);
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';

set(gcf,'Units','centimeters','Position',[10 10 2.5 2]);

basePath = fullfile(saveDir_FigS4, ['Example Inhited cell ' num2str(n) ' ' Params.ExpName{F}]);
saveFigAll(gcf, basePath);

%% Example tuning curve

figure;

% Manual positions: [left bottom width height]
ax1 = subplot('Position',[0.10 0.55 0.25 0.35]); % top left
ax2 = subplot('Position',[0.10 0.10 0.25 0.35]); % bottom left
ax3 = subplot('Position',[0.52 0.10 0.35 0.80]); % right (full height)

angles = 0:45:315;
sfs = [1 2];
axesList = [ax1 ax2];

for i = 1:2
    
    axes(axesList(i)); hold on
    sf = sfs(i);
    
    % ---- Extract ----
    resp_off = squeeze(stimRespInh{F}(n,:,sf,1,:));
    resp_on  = squeeze(stimRespInh{F}(n,:,sf,2,:));
    
    % ---- Mean & SEM ----
    off_mean = mean(resp_off,2);
    off_sem  = std(resp_off,0,2) ./ sqrt(size(resp_off,2));
    on_mean  = mean(resp_on,2);
    on_sem   = std(resp_on,0,2) ./ sqrt(size(resp_on,2));
    
    off_mean = off_mean(:)';
    off_sem  = off_sem(:)';
    on_mean  = on_mean(:)';
    on_sem   = on_sem(:)';
    
    % ---- Plot ----
    errorbar(angles, off_mean, off_sem, 'o-', ...
        'Color',[0.7 0.7 0.7], ...
        'MarkerFaceColor',[0.7 0.7 0.7], ...
        'MarkerEdgeColor','none', ...
        'LineWidth',0.8, ...
        'MarkerSize',2, ...
        'CapSize',0);
    
    errorbar(angles, on_mean, on_sem, 'o-', ...
        'Color',[0.2 0.2 0.2], ...
        'MarkerFaceColor',[0.2 0.2 0.2], ...
        'MarkerEdgeColor','none', ...
        'LineWidth',0.8, ...
        'MarkerSize',2, ...
        'CapSize',0);
    
    xlim([-15 330])
    ylim([0 0.45])
    ylabel('Response')
    
    if i == 2
        xlabel('Orientation (°)')
        xticks(angles(1:2:end))
        %                   ylim([0 0.45])
        
    else
        xticks([])
    end
    text(0.03, 0.95, sprintf('SFTF = %d', sf), ...
        'Units','normalized', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'FontSize',4, ...
        'Color','k');
    
    box off
    set(gca,'FontSize',4.6,'LineWidth',0.6,'TickDir','out')
end



%  subplot(1,2,2); hold on
axes(ax3); hold on
scatter(x(validIdx), y(validIdx),3, 'filled','MarkerFaceColor',[0.4 0.4 0.4]); hold on;
plot(x(validIdx), yfit, 'r-', 'LineWidth', 1,'color',Params.colors{F});
% Add slope and offset text at bottom center
txt = sprintf('Slope = %.2f\nOffset = %.2f', p(1), p(2));
text(0.7, 0.06, txt, ...
    'Units', 'normalized', ...      % use axis-relative positioning
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 3, ...
    'Color', 'k');

xlabel('Laser OFF'); ylabel('Laser ON');
xLimits = xlim; yLimits = ylim;
xlim([min(xLimits(1),yLimits(1)) max(xLimits(2),yLimits(2))]); ylim([min(xLimits(1),yLimits(1)) max(xLimits(2),yLimits(2))]);
plot([0 max(xLimits(2),yLimits(2))], [0 max(xLimits(2),yLimits(2))],'color',[0.5 0.5 0.5],'LineStyle', ':');
axis equal
xlim([0 0.45]); ylim([0 0.45])
set(gca,'FontSize',4.6,'LineWidth',0.6,'TickDir','out','TickLength',[0.02, 0.01]);
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';

set(gcf,'Units','centimeters','Position',[10 10 4.5 2]);

PlotTitle=['Example Inhited cell Tuning ' num2str(n) ' ' Params.ExpName{F}];
basePath = fullfile(saveDir_FigS4,[PlotTitle]);
saveFigAll(gcf, basePath);

fprintf('--- Slope ---\n')
fprintf('FF slope: mean = %.4f, std = %.4f\n', nanmean(slope{2}(valid{2})), nanstd(slope{2}(valid{2})))
fprintf('FB slope: mean = %.4f, std = %.4f\n', nanmean(slope{1}(valid{1})), nanstd(slope{1}(valid{1})))
fprintf('FF vs FB slope p = %.4g\n', pSlope)

fprintf('\n--- Intercept ---\n')
fprintf('FF intercept: mean = %.4f, std = %.4f\n', nanmean(intercept{2}(valid{2})), nanstd(intercept{2}(valid{2})))
fprintf('FB intercept: mean = %.4f, std = %.4f\n', nanmean(intercept{1}(valid{1})), nanstd(intercept{1}(valid{1})))
fprintf('FF vs FB intercept p = %.4g\n', pIntercept)


figure;
subplot(1,2,1); hold on
line([1 1], [0 1], 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1);
h1=plot(slopeCdfX{1}, slopeCdfY{1}, 'LineWidth', 1, 'Color', Params.colors{1});
xMean = nanmean(slope{1}(valid{1}));
plot(xMean, 1,'v','MarkerSize', 6,'MarkerFaceColor', Params.colors{1},'MarkerEdgeColor', 'none');
h2=plot(slopeCdfX{2}, slopeCdfY{2}, 'LineWidth', 1, 'Color', Params.colors{2});
xMean = nanmean(slope{2}(valid{2}));
plot(xMean, 1,'v','MarkerSize', 6,'MarkerFaceColor', Params.colors{2},'MarkerEdgeColor', 'none');
% title(['Slope ' testSlope ' p = ' num2str(pSlope)]);   box off
box off
xlabel('Multiplicative Slope','Fontsize', 12);ylabel('Cumulative frequency')
ylim([0 1]); xlim([min([slopeCdfX{1}; slopeCdfX{2}]) 3])
set(gca,'FontSize',4.6,'LineWidth',0.6,'TickDir','out','TickLength',[0.01, 0.01]);
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';

subplot(1,2,2); hold on
line([0 0], [0 1], 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1);
h1=plot(intCdfX{1}, intCdfY{1}, 'LineWidth', 1, 'Color', Params.colors{1});
xMean = nanmean(intercept{1}(valid{1}));
plot(xMean, 1,'v','MarkerSize', 6,'MarkerFaceColor', Params.colors{1},'MarkerEdgeColor', 'none');
h2=plot(intCdfX{2}, intCdfY{2}, 'LineWidth', 1, 'Color', Params.colors{2});
xMean = nanmean(intercept{2}(valid{2}));
plot(xMean, 1,'v','MarkerSize', 6,'MarkerFaceColor', Params.colors{2},'MarkerEdgeColor', 'none');
legend([h1, h2], {['FB n= ' num2str(sum(valid{1}))], ['FF n= ' num2str(sum(valid{2}))]}); legend boxoff

box off
xlabel('Additive Offset','Fontsize', 12);ylabel('Cumulative frequency')
ylim([0 1]); xlim([min([intCdfX{1}; intCdfX{2}]) 0.2])
set(gca,'FontSize',4.6,'LineWidth',0.6,'TickDir','out','TickLength',[0.01, 0.01]);
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';
set(gcf,'Units','centimeters','Position',[10 10 9 3.5]);

PlotTitle='_Cumulative distributions Laser OFF vs Laser ON Inh ';
basePath = fullfile(saveDir_FigS4,[PlotTitle]);
saveFigAll(gcf, basePath);