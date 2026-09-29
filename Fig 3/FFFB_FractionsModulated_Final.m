saveDir_Fig3 = [Params.MainDir filesep 'plot\Fig 3'];
if ~exist(saveDir_Fig3,'dir')
    mkdir(saveDir_Fig3)
end
colors = [
    0.8 0.8 0.8;
    1 0.6 0.4;
    0.6 0.8 1];

n_ani=Params.n_ani;
n_pos=Params.n_pos;


figure;
for F=1:2
    
    
    responsive = {};
    excited = {};
    inhibited = {};
    
    for ani=1:n_ani{F}
        for pos=1:n_pos{F}(ani)
            if ~isempty(data{F}{ani,pos})
                responsive{end+1} = CellStat{F}{ani,pos}.responsiveAD;
                excited{end+1} = CellStat{F}{ani,pos}.SigExcAD;
                inhibited{end+1} = CellStat{F}{ani,pos}.SigInhAD;
                
            end
        end
    end
    
    responsive = vertcat(responsive{:});
    excited = vertcat(excited{:});
    inhibited = vertcat(inhibited{:});
    
    
    nROIs=sum(responsive);
    FracExc=(sum(excited)./nROIs)*100;
    FracInh=(sum(inhibited)./nROIs)*100;
    FracNoChange=((nROIs-sum(excited)-sum(inhibited))./nROIs)*100;
    
    
    subplot(2,1,3-F); hold on
    h=pie([FracNoChange FracExc FracInh]);
    
    % lgd=legend(labels, 'Location', 'southwest');
    axis off
    legend off
    % Loop through wedges and set color
    colorIdx = 1;
    for i = 1:length(h)
        h(i).LineWidth=1;
        if isgraphics(h(i), 'patch')
            h(i).FaceColor = colors(colorIdx, :);
            colorIdx = colorIdx + 1;
        end
        if isgraphics(h(i), 'patch') && isgraphics(h(i+1), 'text')
            % Get the centroid of the slice
            sliceCenter = mean(h(i).Vertices, 1);
            
            % Move the text to the centroid (inside the slice)
            h(i+1).Position = sliceCenter;
            h(i+1).FontSize = 6;
            h(i+1).FontWeight = 'bold';
            h(i+1).Color = 'black';
        end
    end
    
    
    title([Params.ExpName{F} ' - n = ' num2str(sum(responsive))])
    set(gca, 'Fontsize',6);
    
end

set(gcf,'Units', 'centimeters','Position',[20 10 3 7])
PlotTitle='_Fraction modulated cells pie chart';
basePath = fullfile(saveDir_Fig3,[PlotTitle]);
saveFigAll(gcf, basePath);


