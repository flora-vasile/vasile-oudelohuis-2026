saveDir_Fig3 = [Params.MainDir filesep 'plot\Fig 3'];
if ~exist(saveDir_Fig3,'dir')
    mkdir(saveDir_Fig3)
end

%% dprime
for F=1:2
    for ani=1:Params.n_ani{F}
        for pos=1:Params.n_pos{F}(ani)
            if ~isempty(data{F}{ani,pos})
                %Coordinates
                ROICoordyxResp{F}{ani,pos}=CellStat{F}{ani,pos}.Coordyx_um(CellStat{F}{ani,pos}.responsiveAD,:);
                
                AllDirRespAllRepSes{F}{ani,pos}=squeeze(nanmean(data{F}{ani,pos}.StimAv(:,StimONResponse,:,:),2));
                AllDirRespSes{F}{ani,pos}=squeeze(nanmean(nanmean(data{F}{ani,pos}.StimAv(:,StimONResponse,:,:),2),4));
                
                dprime{F}{ani,pos}=(AllDirRespSes{F}{ani,pos}(CellStat{F}{ani,pos}.responsiveAD,2)...
                    -AllDirRespSes{F}{ani,pos}(CellStat{F}{ani,pos}.responsiveAD,1))./...
                    sqrt(1/2*(nanvar(AllDirRespAllRepSes{F}{ani,pos}(CellStat{F}{ani,pos}.responsiveAD,2,:),0,3)...
                    +nanvar(AllDirRespAllRepSes{F}{ani,pos}(CellStat{F}{ani,pos}.responsiveAD,1,:),0,3)));
                
            end
        end
    end
end


%% Plot Map
% Create scatter plot with color mapping
figure;
F=2; ani=18; pos=3;
subplot(2,1,1)
scatter(ROICoordyxResp{F}{ani,pos}(:,2), ROICoordyxResp{F}{ani,pos}(:,1), 5, dprime{F}{ani,pos}, 'filled'); hold on

% Set colormap
Orange_to_blue_colormap
colormap(cmap)

% Set color limits
maxAbsValue = max(max(abs(dprime{F}{ani,pos})));
caxis([-0.4, 0.4]) %caxis([-maxAbsValue, maxAbsValue])

% Add a colorbar to show the color scale
colorbar;
% Add axis labels and title
set(gca, 'Fontsize', 6); set(gca,'linewidth',1);

xlim([min(ROICoordyxResp{F}{ani,pos}(:,2))-10, max(ROICoordyxResp{F}{ani,pos}(:,2))+10])
ylim([min(ROICoordyxResp{F}{ani,pos}(:,1))-10, max(ROICoordyxResp{F}{ani,pos}(:,1))+10])

set(gca, 'Color', [0.5 0.5 0.5]);
set(gcf, 'Color', [1 1 1], 'InvertHardcopy', 'off');

ax = gca;

ax.XColor = 'none';   % hides x-axis line + ticks
ax.YColor = 'none';   % hides y-axis line + ticks
ax.TickLength = [0 0]; % removes tick marks

hold on

% Length of scale bar (in same units as your coordinates)
scaleLength = 50;  % 50 µm

% Get axis limits
xl = xlim;
yl = ylim;

% Position the bar (bottom-left corner with some padding)
xStart = xl(1) + 50;        % 10 units from left edge
yStart = yl(1) + 50;        % 10 units from bottom edge

% Draw scale bar
plot([xStart xStart + scaleLength], ...
    [yStart yStart], ...
    'k', 'LineWidth', 0.5)

% Optional: add text label
text(xStart + scaleLength/2, ...
    yStart + 5, ...
    '50\mum', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 4, ...
    'Color', 'k')


F=1; ani=9; pos=1;
subplot(2,1,2)
scatter(ROICoordyxResp{F}{ani,pos}(:,2), ROICoordyxResp{F}{ani,pos}(:,1), 5, dprime{F}{ani,pos}, 'filled'); hold on

% Set colormap
Orange_to_blue_colormap
colormap(cmap)

% Set color limits
maxAbsValue = max(max(abs(dprime{F}{ani,pos})));
caxis([-maxAbsValue, maxAbsValue])
% Add a colorbar to show the color scale
colorbar;
% Add axis labels and title
set(gca, 'Fontsize', 6); set(gca,'linewidth',1);
xlabel('X position (um)');
ylabel('Y position (um)');
xlim([min(ROICoordyxResp{F}{ani,pos}(:,2))-10, max(ROICoordyxResp{F}{ani,pos}(:,2))+10])
ylim([min(ROICoordyxResp{F}{ani,pos}(:,1))-10, max(ROICoordyxResp{F}{ani,pos}(:,1))+10])
set(gca, 'Color', [0.5 0.5 0.5]);
set(gcf, 'Color', [1 1 1], 'InvertHardcopy', 'off');

ax = gca;

ax.XColor = 'none';   % hides x-axis line + ticks
ax.YColor = 'none';   % hides y-axis line + ticks
ax.TickLength = [0 0]; % removes tick marks

set(gcf,'Units','centimeters','Position',[10 10 4 5.5])

PlotTitle='_dprime maps';
basePath = fullfile(saveDir_Fig3,[tv PlotTitle]);
saveFigAll(gcf, basePath);


