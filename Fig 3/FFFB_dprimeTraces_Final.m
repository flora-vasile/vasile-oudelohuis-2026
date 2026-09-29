saveDir_Fig3 = [Params.MainDir filesep 'plot\Fig 3'];
if ~exist(saveDir_Fig3,'dir')
    mkdir(saveDir_Fig3)
end

% concatenate data for each animal
for F=1:2

dataStimAv{F} = cell(size(data{F},1),1);
responsive{F} = cell(size(data{F},1),1);

for ani = 1:size(data{F},1)
    tmp = {};
    tmpresp = {};
    for pos = 1:size(data{F},2)
        if ~isempty(data{F}{ani,pos})
            tmp{end+1} = data{F}{ani,pos}.StimAv;
            tmpresp{end+1} = CellStat{F}{ani,pos}.responsiveAD;
        end
    end
    
    if ~isempty(tmp)
        dataStimAv{F}{ani} = vertcat(tmp{:});
        responsive{F}{ani} = vertcat(tmpresp{:});
    end
end
end

% Calculate d´
for F=1:2
    for ani=1:n_ani{F}
        if ~isempty(dataStimAv{F}{ani})
meanROI{F}{ani}=mean(dataStimAv{F}{ani},4);

dprimeTraceAni{F}{ani}=(meanROI{F}{ani}(responsive{F}{ani},:,2)-meanROI{F}{ani}(responsive{F}{ani},:,1))./...
sqrt(1/2*(var(dataStimAv{F}{ani}(responsive{F}{ani},:,2,:),0,4)...
+var(dataStimAv{F}{ani}(responsive{F}{ani},:,1,:),0,4)));


for i=1:size(dprimeTraceAni{F}{ani},1)
PeakOffAni{F}{ani}(i)=max(dprimeTraceAni{F}{ani}(i,:));
DINormtoPeakTraceAni{F}{ani}(i,:) = dprimeTraceAni{F}{ani}(i,:)./PeakOffAni{F}{ani}(i);
end


meanRespDINormtoPeakTraceAni{F}{ani}=mean(dprimeTraceAni{F}{ani}(:,Params.StimONResponse),2);

[B,I]=sort(meanRespDINormtoPeakTraceAni{F}{ani}, 'descend');

for i=1: size(dprimeTraceAni{F}{ani},1)
    dprimeTraceAni2{F}{ani}(i,:) = dprimeTraceAni{F}{ani}(I(i),:);
end   

        end
    end
end



figure;
subplot(1,2,1)
HeatMapPlot(dprimeTraceAni2{2}{14}(:,1:60),'ROI #','Time (s)',[' '],[-0.4 0.4])
Orange_to_blue_colormap
colormap(cmap)
a=colorbar;
set(gca, 'Fontsize', 12,'LineWidth', 2,'TickDir','out','TickLength',[0.008, 0.01])
box off

subplot(1,2,2)
HeatMapPlot(dprimeTraceAni2{1}{11}(:,1:60),' ','Time (s)',[' '],[-0.4 0.4])
colormap(cmap)
a=colorbar;
a.Label.String = ' ';
set(gca, 'Fontsize', 12,'LineWidth', 2,'TickDir','out','TickLength',[0.008, 0.01])
set(gcf,'Units','centimeters','Position',[10 10 15 10])
box off

PlotTitle='_Raster plot d';
basePath = fullfile(saveDir_Fig3,[PlotTitle]);
saveFigAll(gcf, basePath);
