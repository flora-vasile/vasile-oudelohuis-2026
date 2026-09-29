function plotDecoderComparison_mixedANOVA_fig(OFF_perf, ON_perf, varargin)
% Compare decoding performance for FF (F=2) and FB (F=1) on a single plot
% with 4 columns: FF OFF, FF ON, FB OFF, FB ON (paired lines)
% and perform a 2-way mixed ANOVA (Laser repeated, Group independent).

%% --- Parse inputs ---
p = inputParser;
addParameter(p, 'plotColors', [1 0.4 0.4; 0.2 0.6 1], @isnumeric); % FF, FB
addParameter(p, 'groupLabels', {'FF','FB'}, @iscell); % FF first
addParameter(p, 'ExpName', '', @ischar);
addParameter(p, 'saveDir', '', @ischar);
parse(p, varargin{:});
opt = p.Results;

sem = @(x) nanstd(x) / sqrt(sum(~isnan(x)));

%% --- Reorder data: FF first, FB second ---
F_order = [2 1]; % FF first
FF_OFF = OFF_perf{2}(:); FF_ON = ON_perf{2}(:);
FB_OFF = OFF_perf{1}(:); FB_ON = ON_perf{1}(:);

%% --- Normalize Laser ON to Laser OFF ---
FF_norm = FF_ON ./ FF_OFF; % Laser ON as % of OFF
FB_norm = FB_ON ./ FB_OFF;

% For plotting: OFF columns = 100, ON columns = normalized values
FF_OFF_plot = ones(size(FF_OFF)); % OFF baseline = 100%
FB_OFF_plot = ones(size(FB_OFF));


%% --- Prepare table for mixed ANOVA ---
% Get number of subjects in each group
nSubj_FF = numel(FF_OFF);
nSubj_FB = numel(FB_OFF);

% --- Prepare response vector (Y) ---
Y = [FF_OFF; FF_ON; FB_OFF; FB_ON];

% --- Subject IDs (unique across groups) ---
subject = [(1:nSubj_FF)'; (1:nSubj_FF)'; ...         % FF
           (1:nSubj_FB)'+nSubj_FF; (1:nSubj_FB)'+nSubj_FF]; % FB

% --- Laser factor (within-subject) ---
laser = [zeros(nSubj_FF,1); ones(nSubj_FF,1); ...   % FF OFF/ON
         zeros(nSubj_FB,1); ones(nSubj_FB,1)];      % FB OFF/ON

% --- Group factor (between-subject) ---
group = [ones(nSubj_FF*2,1); 2*ones(nSubj_FB*2,1)];

% --- Build table ---
tbl_lme = table(Y, laser, group, subject, ...
    'VariableNames', {'Accuracy','Laser','Group','Subject'});

% --- Fit linear mixed-effects model ---
lme = fitlme(tbl_lme, 'Accuracy ~ Laser*Group + (1|Subject)');

% --- Get ANOVA ---
anovaTbl = anova(lme);

laser_p = anovaTbl.pValue(2); % main effect of Laser (within)
group_p = anovaTbl.pValue(3); % main effect of Group (between)
int_p   = anovaTbl.pValue(4); % Laser × Group interaction


%% t-test
% FF
idxFF = ~isnan(FF_OFF) & ~isnan(FF_ON);
FF_off = FF_OFF(idxFF);
FF_on  = FF_ON(idxFF);

% Normality tests
isNonNormalFF = kstest(zscore(FF_off)) || kstest(zscore(FF_on));

if isNonNormalFF
    pFF = signrank(FF_off, FF_on);
    testNameFF = 'signrank';
else
    [~, pFF] = ttest(FF_off, FF_on);
    testNameFF = 'ttest';
end

% FB
idxFB = ~isnan(FB_OFF) & ~isnan(FB_ON);
FB_off = FB_OFF(idxFB);
FB_on  = FB_ON(idxFB);

isNonNormalFB = kstest(zscore(FB_off)) || kstest(zscore(FB_on));

if isNonNormalFB
    pFB = signrank(FB_off, FB_on);
    testNameFB = 'signrank';
else
    [~, pFB] = ttest(FB_off, FB_on);
    testNameFB = 'ttest';
end


%% --- Plotting ---
fig = figure('Position',[100 100 700 600]); hold on;

% Define x positions for the 4 columns
xPos = [0.5 1.5 2.5 3.5]; % 1,2 = FF OFF/ON; 4,5 = FB OFF/ON
xOffset = [-0.15 0.15]; % mean line offsets

allData = {FF_OFF_plot, FF_norm, FB_OFF_plot, FB_norm};

% Plot paired lines
for iGroup = 1:2
    OFF_all = allData{(iGroup-1)*2+1};
    ON_all  = allData{(iGroup-1)*2+2};
    xOFF = xPos((iGroup-1)*2+1);
    xON  = xPos((iGroup-1)*2+2);
    
    if iGroup==1
      nSubj=nSubj_FF;
    else
      nSubj=nSubj_FB;

    end
  
    for iSubj = 1:nSubj
        if ~isnan(OFF_all(iSubj)) && ~isnan(ON_all(iSubj))
            plot([xOFF xON], [OFF_all(iSubj) ON_all(iSubj)], '-', 'Color', [opt.plotColors(iGroup,:) 0.5], 'LineWidth', 0.5);
        end
    end
    
    % Plot mean ± SEM
    means = [nanmean(OFF_all), nanmean(ON_all)];
    sems  = [sem(OFF_all), sem(ON_all)];
    
% Plot mean ± SEM with connected means
xPair = xPos((iGroup-1)*2+1 : (iGroup-1)*2+2);

% Plot mean line connecting OFF and ON
plot(xPair, means, '-', ...
    'Color', opt.plotColors(iGroup,:), ...
    'LineWidth', 1);

% Plot mean markers
plot(xPair, means, 'o', ...
    'Color', opt.plotColors(iGroup,:), ...
    'MarkerFaceColor', opt.plotColors(iGroup,:), ...
    'MarkerSize', 3);

% Plot SEM bars
for j = 1:2
    line([xPair(j) xPair(j)], ...
        [means(j)-sems(j), means(j)+sems(j)], ...
        'Color', opt.plotColors(iGroup,:), ...
        'LineWidth', 1);
end

end

% Formatting
xlim([0 4]); %ylim([0 100]);
set(gca, ...
    'XTick', [0.5 1.5 2.5 3.5], ...
    'XTickLabel', {'Laser OFF','Laser ON','Laser OFF','Laser ON'});
ylabel('Relative Decoding Performance');
box off;

yMax = 1.6;
xFF = mean(xPos(1:2)); % FF OFF–ON
xFB = mean(xPos(3:4)); % FB OFF–ON

% FF significance
% plot(xPos(1:2), [yMax yMax], 'k', 'LineWidth', 1.5)
text(xFF, yMax*1.02, p2star(pFF), ...
    'HorizontalAlignment', 'center', 'FontSize', 6)

% FB significance
% plot(xPos(3:4), [yMax yMax], 'k', 'LineWidth', 1.5)
text(xFB, yMax, p2star(pFB), ...
    'HorizontalAlignment', 'center', 'FontSize', 10)

% Interaction mixed anova
text((xFB+xFF)/2, yMax+0.2, ['Interaction: ' p2star(int_p)], ...
    'HorizontalAlignment', 'center', 'FontSize', 6)

% plot(xlim,[12.5 12.5],'k:'); text(0.2,15,'Chance');

% Figure title with ANOVA results
% suptitle(sprintf('Decoder Performance: %s\nLaser p=%.3f | Pathway p=%.3f | Interaction p=%.3f', ...
%     opt.ExpName, laser_p, group_p, int_p));

% despline(gca,0.5);
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';
set(gca,'FontSize',5,'LineWidth',0.6,'TickDir','out','TickLength',[0.01, 0.01]);
xtickangle(20)
set(gcf,'Units','centimeters','Position',[10 10 4.2 3.5]);

%% --- Save figure if requested ---
if ~isempty(opt.saveDir)
    filename = sprintf('DecoderComparison_MixedANOVA_%s.png', opt.ExpName);
    %saveas(fig, fullfile(opt.saveDir, filename));
    saveas(gcf,[opt.saveDir filesep '_DecoderComparison_MixedANOVA ' ],'fig');
        saveas(gcf,[opt.saveDir filesep '_DecoderComparison_MixedANOVA ' ],'svg');
    fprintf('Figure saved to %s\n', fullfile(opt.saveDir, filename));
end

end
