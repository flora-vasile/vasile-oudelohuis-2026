
saveDir_Fig3 = [Params.MainDir filesep 'plot\Fig 3'];
if ~exist(saveDir_Fig3,'dir')
    mkdir(saveDir_Fig3)
end

%% Plot raw traces
figure; hold on
F=1; ani=11 ; pos=3; SF=1; ROI=49; % FB
% F=2; ani=7 ; pos=3; SF=1; ROI=1; % FF
datatmp=data{F}{ani,pos}.Raw;
datatmp(:,:,1,2,:,:)=datatmp(:,:,1,1,:,:);
meanRepdataAll = squeeze(nanmean(datatmp(:,1:60,:,:,:,:),6));
semRepdataAll = squeeze(nanstd(datatmp(:,1:60,:,:,:,:),[],6)./sqrt(size(data{F}{ani,pos}.Raw,6)));

plotOrder=[5 4 1 2 3 6 9 8 7];
x = [1:60];
a =  meanRepdataAll(:,:,:,:,1);
b = semRepdataAll(:,:,:,:,1);
c =  meanRepdataAll(:,:,:,:,2);
d = semRepdataAll(:,:,:,:,2);

for Dir = 1:size(meanRepdataAll,3)
    tightSubplot(3,3,plotOrder(Dir),0.01, 0.01); hold on
        p1 = plot(movmean(meanRepdataAll(ROI,:,Dir,SF,1),3), 'Color', [0.1255 1 0.7804], 'linewidth', 0.8); hold on
        p2 = plot(movmean(meanRepdataAll(ROI,:,Dir,SF,2),3), 'Color', [0 0.3804 0.2863], 'linewidth', 0.8);
        h1 = fill([x flip(x)],[movmean(a(ROI,:,Dir,SF)-b(ROI,:,Dir,SF),3) flip(movmean(a(ROI,:,Dir,SF)+b(ROI,:,Dir,SF),3))],[0.1255 1 0.7804],'linestyle','none','facealpha',.3);
        h2 = fill([x flip(x)],[movmean(c(ROI,:,Dir,SF)-d(ROI,:,Dir,SF),3) flip(movmean(c(ROI,:,Dir,SF)+d(ROI,:,Dir,SF),3))], [0 0.3804 0.2863],'linestyle','none','facealpha',.3);
ylim([-0.02, 0.4])
    y1=ylim;
    h=fill([15,30,30,15],[y1(1),y1(1),y1(2),y1(2)], [0.9020 0.9020 0.9020],'linestyle','none');
    
    xlim([7 45])
    
    
    ax = gca;
    
    ax.XColor = 'none';   % hides x-axis line + ticks
    ax.YColor = 'none';   % hides y-axis line + ticks
    ax.TickLength = [0 0]; % removes tick marks
    
    if plotOrder(Dir) == 9
        
        % Length of scale bar (in same units as your coordinates)
        scaleLength_x = 7.5;
        scaleLength_y = 0.1;
        
        % Get axis limits
        xl = xlim;
        yl = ylim;
        
        % Position the bar (bottom-left corner with some padding)
        xStart = 37.5;       
        yStart = -0.03;        
        
        % Draw scale bar
        plot([xStart xStart + scaleLength_x], ...
            [yStart yStart], ...
            'k', 'LineWidth', 0.5)
        
        % Optional: add text label
        text(xStart + scaleLength_x/2, ...
            yStart + 5, ...
            '0.5s', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 4, ...
            'Color', 'k')
        
        % Position the bar (bottom-left corner with some padding)
        xStart = xStart+scaleLength_x;        % 10 units from left edge
        
        % Draw scale bar
        plot([xStart xStart], ...
            [yStart yStart + scaleLength_y], ...
            'k', 'LineWidth', 0.5)
        
        % Optional: add text label
        text(xStart + scaleLength_y/2, ...
            yStart + 5, ...
            '0.1 Events/F0', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 4, ...
            'Color', 'k')
    end
    
    
    set(gca, 'Fontsize', 6)
    set(gca,'linewidth',1)
end
set(gcf,'Units','centimeters','Position',[10 10 4 4])
basePath = fullfile(saveDir_Fig3,[tv '_ Responsive cells - ' Params.ExpName{F}]);
saveFigAll(gcf, basePath);


